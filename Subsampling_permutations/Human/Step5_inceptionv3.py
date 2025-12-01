#!/usr/bin/env python3

#although not used for the study, this is a useful code for running inception v3 using gpu on python rather than cpu on R

import os, sys, logging
from glob import glob
import h5py
import numpy as np
import tensorflow as tf

# Quieter TF logs
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
logging.getLogger("tensorflow").setLevel(logging.ERROR)

# ---------------- config ----------------
try:
    perm = int(os.environ["SLURM_ARRAY_TASK_ID"])
except Exception:
    sys.exit("Error: SLURM_ARRAY_TASK_ID not set to a valid integer")

IMG_BASE_DIR = "pngLs/ml_pngLs/H"
OUT_DIR      = "features_h5/H"
LOG_DIR      = "logs"
BATCH_SIZE   = 32
IMG_SIZE     = (299, 299)
NUM_SPLITS   = 4
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

print(f"→ Starting extraction for perm {perm}")
print("TensorFlow version:", tf.__version__)
gpus = tf.config.list_physical_devices("GPU")
if gpus:
    print("GPUs detected:", gpus)
    for g in gpus:
        try:
            tf.config.experimental.set_memory_growth(g, True)
        except Exception as e:
            print("  Warning: could not set memory growth:", e)
else:
    print("No GPU detected — will run on CPU")

# Load model
model = tf.keras.applications.InceptionV3(weights="imagenet", include_top=False)
AUTOTUNE = tf.data.AUTOTUNE

# ---------- utilities ----------
def can_decode_png(path: str) -> bool:
    """Return True if file exists, non-empty, and decodes as 3‑ch PNG."""
    try:
        if not os.path.isfile(path):
            return False
        if os.path.getsize(path) < 8:  # PNG header is 8 bytes
            return False
        # quick decode test in eager mode
        b = tf.io.read_file(path)
        _ = tf.image.decode_png(b, channels=3)
        return True
    except Exception:
        return False

def process_images(valid_paths):
    """Build dataset from validated paths and return features array."""
    ds = tf.data.Dataset.from_tensor_slices(valid_paths)
    ds = ds.map(lambda p: tf.image.decode_png(tf.io.read_file(p), channels=3),
                num_parallel_calls=AUTOTUNE)
    ds = ds.map(lambda img: tf.image.resize(img, IMG_SIZE),
                num_parallel_calls=AUTOTUNE)
    ds = ds.map(tf.keras.applications.inception_v3.preprocess_input,
                num_parallel_calls=AUTOTUNE)
    # if you want a last-resort guard (may silently drop a bad example):
    # ds = ds.apply(tf.data.experimental.ignore_errors())
    ds = ds.batch(BATCH_SIZE).prefetch(AUTOTUNE)

    parts = []
    for batch in ds:
        conv = model(batch, training=False)                 # [bs, 8,8,2048]
        flat = tf.reshape(conv, (conv.shape[0], -1)).numpy()
        parts.append(flat)
    return np.vstack(parts)

# ---------- gather, validate, split ----------
img_dir = os.path.join(IMG_BASE_DIR, f"perm_{perm:02d}")
all_imgs = sorted(glob(os.path.join(img_dir, "*.png")))
if not all_imgs:
    print(f"[WARN] No PNGs found in {img_dir}, exiting.")
    sys.exit(0)

# Validate files once up-front
valid_imgs = []
bad_imgs   = []
for p in all_imgs:
    if can_decode_png(p):
        valid_imgs.append(p)
    else:
        bad_imgs.append(p)

print(f"Found {len(all_imgs)} PNGs; valid={len(valid_imgs)}, bad={len(bad_imgs)}")
if bad_imgs:
    bad_log = os.path.join(LOG_DIR, f"bad_pngs_perm_{perm:02d}.txt")
    with open(bad_log, "w") as fh:
        fh.write("\n".join(bad_imgs))
    print(f"[WARN] Wrote list of bad PNGs to {bad_log}")

if not valid_imgs:
    print("[ERROR] All images failed validation.")
    sys.exit(1)

# Balanced splits (no empty subsets)
chunks = np.array_split(np.array(valid_imgs, dtype=object), NUM_SPLITS)

# ---------- run inference ----------
all_feats, all_cells = [], []
for idx, subset in enumerate(chunks, start=1):
    subset = subset.tolist()
    if len(subset) == 0:
        print(f"  Split {idx}/{NUM_SPLITS} is empty — skipping")
        continue
    print(f"  Perm {perm}: split {idx}/{NUM_SPLITS}, {len(subset)} images")
    feats = process_images(subset)
    names = [os.path.splitext(os.path.basename(p))[0] for p in subset]
    all_feats.append(feats)
    all_cells.extend(names)

features_all = np.vstack(all_feats).astype("float32")
cells_all    = np.array(all_cells, dtype=object)

# ---------- save ----------
out_file = os.path.join(OUT_DIR, f"H_perm{perm:02d}.h5")
with h5py.File(out_file, "w") as f:
    f.create_dataset("features", data=features_all, compression="gzip")
    dt = h5py.string_dtype(encoding="utf-8")
    f.create_dataset("cells", data=cells_all, dtype=dt)
print(f"✅ Saved features for perm {perm} → {out_file}")

