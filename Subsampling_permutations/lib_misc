# Originally written by Lucia Gomez

# LIBRARY FOR MISCELANEOUS HUMOUS FUNCTIONS

library(Seurat) ; library(SeuratObject) ; library(bmrm) ; library(parallel) ; library(scales) ; library(diptest) ; library(fitdistrplus)


# gene selection for landscapes, expressed in all groups (datasets) an more than 10 cells overall
gene_selector <- function(seuratobj,grouping,threshold_grouping,threshold_ncells){
  
  df <- data.frame(cells=colnames(seuratobj),grouping=grouping)
  
  crossgroup_genes <- local({
    crossgroup_genes <- list()
    for(i in levels(as.factor(df$grouping))){
      
      sumsvect <- rowSums(seuratobj@assays$RNA@counts[,df$cells[df$grouping==i]])
      crossgroup_genes[[i]] <- names(sumsvect)[which(rowSums(seuratobj@assays$RNA@counts[,df$cells[df$grouping==i]]) > threshold_grouping)] 
    }
    crossgroup_genes <- Reduce(intersect,crossgroup_genes) 
  })
  
  sumsvect <- rowSums(seuratobj@assays$RNA@counts)
  expr_genes10 <- names(sumsvect)[which(rowSums(seuratobj@assays$RNA@counts) > threshold_ncells)]
  
  selected_genes <- intersect(crossgroup_genes,expr_genes10)
  
  return(selected_genes)
  
}


# dataset integration parwise (if reference_dataset is not provided) or with reference ----
dataset_integration <- function(list_seuratobjs,nfeats,FindIntDims,FindInt_kScore,IntKweight,reference_dataset=NULL){
  
  if(is.null(reference_dataset)){
    
    for (i in 1:length(list_seuratobjs)) {
      list_seuratobjs[[i]] <- NormalizeData(list_seuratobjs[[i]], verbose = FALSE)
      list_seuratobjs[[i]] <- FindVariableFeatures(list_seuratobjs[[i]], selection.method = "vst",nfeatures=nfeats, verbose = FALSE)
    }
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = list_seuratobjs)
    mer.anchors <- FindIntegrationAnchors(object.list = list_seuratobjs, anchor.features=features,dims=FindIntDims,k.score = FindInt_kScore)
    mer.integrated <- IntegrateData(anchorset = mer.anchors, features = features,k.weight = IntKweight)
    return(mer.integrated)
    
  }else{
    
  for (i in 1:length(list_seuratobjs)) {
    list_seuratobjs[[i]] <- NormalizeData(list_seuratobjs[[i]], verbose = FALSE)
    list_seuratobjs[[i]] <- FindVariableFeatures(list_seuratobjs[[i]], selection.method = "vst",nfeatures=nfeats, verbose = FALSE)
  }
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = list_seuratobjs)
  reference_data <- which(names(full.list) == reference_dataset)
  mer.anchors <- FindIntegrationAnchors(object.list = list_seuratobjs, reference=reference_data,anchor.features=features,dims=FindIntDims,k.score = FindInt_kScore)
  mer.integrated.ref <- IntegrateData(anchorset = mer.anchors, features = features,k.weight = IntKweight)
  return(mer.integrated.ref)
  }
}


# integrate datasets into an already integrated object ----
double_dataset_integration <- function(list_seuratobjs,nfeats,FindIntDims,FindInt_kScore,IntKweight,reference_dataset,features){
  reference_data <- which(names(full.list) == reference_dataset)
  mer.anchors <- FindIntegrationAnchors(object.list = list_seuratobjs, reference=reference_data,anchor.features=features,dims=FindIntDims,k.score = FindInt_kScore)
  mer.integrated.ref <- IntegrateData(anchorset = mer.anchors, features = features,k.weight = IntKweight)
  return(mer.integrated.ref)
}



# cell selection depending on grouping ----
cell.selector <- function(seuratobject,cellnames,grouping,n){set.seed(1234);
  df <- data.frame(cellnames=cellnames,celltypes=grouping) ; library(dplyr) ; df_f <- df %>% group_by(celltypes) %>% sample_n(size = n)
  S_f <- subset(seuratobject,cells=df_f$cellnames)}



# cost matrix function ----
costM <- function(X,y,...) {
  C <- costMatrix(y)
  C <- C / tabulate(y)[col(C)] / tabulate(y)
  C <- C/sum(C)
}



# ordinal SVM full model training function ----
training_ordi_full <- function(datamatrix,target,cost,lambda_full,epsilon_full,maxiter){
  set.seed(1234)
  # train full model
  full <- nrbm(ordinalRegressionLoss(datamatrix,target,C=cost),LAMBDA=lambda_full,EPSILON_TOL=epsilon_full,MAX_ITER=maxiter)
}



# ordinal SVM reduced model training-testing CV function ----
custom_red_and_pred <- function(fullmodel,xtrain,target,xtest,ngenesselect,
                                lambda_red,epsilon_red,maxiter,nfolds){
  set.seed(1234)
  print("extract feat weights from full model and select top X depending on the dataset for testing")
  w <- t(attr(fullmodel,"gradient")) 
  # check which genes are common between w and xtest and, based on this, select top X model genes (most negative or positive weights)
  print(paste0("select top ",ngenesselect," model genes"))
  common <- intersect(colnames(xtest),rownames(w))
  w_f <- as.matrix(w[common,])
  
  # feat selection on weights
  w_f <- w_f[order(w_f),] ; w_f <- w_f[c(1:ngenesselect,(length(w_f)-(ngenesselect-1)):length(w_f))] ; print(length(w_f))
  # subset xtrain and xtest by feats selected
  xtrain_red <- xtrain[,names(w_f)] ; print(dim(xtrain_red))
  xtest_red <- xtest[,names(w_f)] ; print(dim(xtest_red))
  
  costM <- function(X,y,...) {
    C <- costMatrix(y)
    C <- C / tabulate(y)[col(C)] / tabulate(y)
    C <- C/sum(C)
  }
  cost=costM(X=xtrain_red,y=target) ; print(length(target)) ; print(dim(cost))
  
  # train reduced model without crossval (for future predictions)
  redmodel <- nrbm(ordinalRegressionLoss(xtrain_red,target,C=cost),LAMBDA = lambda_red,EPSILON_TOL =epsilon_red ,MAX_ITER = maxiter)
  
  # run reduced model with crossval and predict
  print("train reduced model with CV and predict")
  folds <- balanced.cv.fold(target,nfolds)
  pred <- simplify2array(mclapply(levels(folds),mc.cores=16,function(f) {
    w <- nrbm(ordinalRegressionLoss(xtrain_red[folds!=f,],target[folds!=f],C=cost),LAMBDA = lambda_red,EPSILON_TOL =epsilon_red ,MAX_ITER = maxiter)
    Y <- predict(w,xtest_red) ; Y[folds!=f] <- NA ; Y
  }))
  # prediction result aggregated across CVs
  pred <- rowSums(pred,na.rm=TRUE)
  # scale prediction from 0 to 1
  #pred <- rescale(pred,to=c(0,1))
  
  # return reduced model (without fold) and prediction vector
  return(list(redmodel=redmodel,pred=pred))
  
}


ordi_activate <- function(ordiscore){
  
  # detect if not unimodal, if so, enter the loop
  if(diptest::dip.test(ordiscore,B=100)$p.value < 0.05){
    
      print("is not unimodal, check if is uniform")
      
      if( ks.test( ordiscore, "punif")$p.value > 0.05 ) {
        print("is uniform, nothing to do")
        return(ordiscore)
        
      }else{
        print("not uniform, smooth tanh*1")
        ordiscore <- tanh((ordiscore - median(ordiscore))*1)
      }
      
  }else{
    print("is unimodal, strong tahn")
    sd <- sd(ordiscore) ; print(sd)
    if( sd<0.15 ){
      print("sd<0.15, strong tanh*5")
      ordiscore <- tanh((ordiscore - median(ordiscore))*5)  
    }else if(sd<0.175){
      print("sd<0.175, strong tanh*2.5")
      ordiscore <- tanh((ordiscore - median(ordiscore))*2.5)
    }else{
      print("sd>0.175, strong tanh*2")
      ordiscore <- tanh((ordiscore - median(ordiscore))*2)
    }
  }
  
}
    

ordi_normalize <- function(ordiscore,ordigroups=NULL,ordigroups_other=NULL,limits1=NULL,limits2=NULL,limits3=NULL,applyWinsor=TRUE){
  
  ordiscore <- rescale(ordi_activate(ordiscore),to=c(0,1))
  
  if(is.null(limits1)){
    print("limits not set, letting default and applying tanh activation")
  }else{
    print("applying limits")
    for (i in levels(as.factor(ordigroups))){
      print(i)
      if(i==1){
        ordiscore[ordigroups==i] <- rescale(ordi_activate(ordiscore[ordigroups==i]),to=limits1)
      }else if(i==2){
        ordiscore[ordigroups==i] <- rescale(ordi_activate(ordiscore[ordigroups==i]),to=limits2)
      }else{
        ordiscore[ordigroups==i] <- rescale(ordi_activate(ordiscore[ordigroups==i]),to=limits3)
      }
    }
  }
  
  if(is.null(ordigroups_other)){
    print("2D correction not needed")
  }else{
    print("applying 2D correction")
    for (i in levels(as.factor(ordigroups_other))){ordiscore[ordigroups_other==i] <- rescale(ordi_activate(ordiscore[ordigroups_other==i]),to=c(0,1))}
  }

  if(applyWinsor==TRUE){
    print("applying winsor transformation")
    ordiscore <- Winsorize(ordiscore,val = quantile(ordiscore, probs = c(0.01, 0.99)))
  }else{
    print("not applying winsor transformation, letting default")
    return(ordiscore)
  }
  
  return(ordiscore)
}
    














