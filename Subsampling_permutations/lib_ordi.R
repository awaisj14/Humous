

# LIBRARY FOR ORDINALS
################################################################################################################################################################

library(bmrm) ; library(Seurat) ; library(scales) ; library(parallel)

# CELL SELECTION FUNCTION
#####################################
cell.selector <- function(seuratobject,cellnames,grouping,n){set.seed(1234);
  df <- data.frame(cellnames=cellnames,celltypes=grouping) ; library(dplyr) ; df_f <- df %>% group_by(celltypes) %>% sample_n(size = n)
  S_f <- subset(seuratobject,cells=df_f$cellnames)}
#####################################


# COST MATRIX FUNCTION
#####################################
costM <- function(X,y,...) {
  C <- costMatrix(y)
  C <- C / tabulate(y)[col(C)] / tabulate(y)
  C <- C/sum(C)
}
#####################################


# TRAINING FUNCTION - OLD ONE (humous v2)
#####################################

# training_ordi_red_old <- function(datamatrix,target,cost,lambda_full,epsilon_full,maxiter,
#                           ngenesselect,ncrossvals,lambda_red,epsilon_red){
#   set.seed(1234)
#   # train full model and extract feature weights
#   full <- nrbm(ordinalRegressionLoss(datamatrix,target,C=cost),LAMBDA=lambda_full,EPSILON_TOL=epsilon_full,MAX_ITER=maxiter)
#   w <- t(attr(full,"gradient")) 
#   
#   # feat selection and run reduced model with crossval
#   w <- w[order(w),] ; w <- w[c(1:ngenesselect,(length(w)-(ngenesselect-1)):length(w))] ; print(head(w)) ; print(tail(w))
#   datamatrix_red <- datamatrix[,names(w)]
#   
#   folds <- balanced.cv.fold(target,num.cv=ncrossvals) ; table(folds,target)
#   red_list <- list()
#   for (f in levels(folds)){
#     # seed will be a random integer between 1000 and 2000
#     set.seed(round(runif(1, 1000, 2000)))
#     red_list[[f]] <- nrbm(ordinalRegressionLoss(datamatrix_red[folds!=f,],target[folds!=f],C=cost),LAMBDA=lambda_red,EPSILON_TOL=epsilon_red,MAX_ITER=maxiter)
#   }
# 
#   w_red <- local({
#     w_red <- red_list$`1`
#     bindedweights <- rlist::list.cbind(red_list) # this generates a 2D array with genes in rows, folds in columns and weights as content
#     # each column of bindedweights is a vector of model weights
#     ordi.meanw <- rowMeans(bindedweights)
#     w_red[1:length(colnames(datamatrix_red))] <- ordi.meanw
#     w_red
#   })
#   attr(w_red,"genenames") <- colnames(datamatrix_red)
#   
#   w_red
# }
#####################################



# TRAINING FUNCTION - FULL MODEL
#####################################
training_ordi_full <- function(datamatrix,target,cost,lambda_full,epsilon_full,maxiter){
  set.seed(1234)
  # train full model
  full <- nrbm(ordinalRegressionLoss(datamatrix,target,C=cost),LAMBDA=lambda_full,EPSILON_TOL=epsilon_full,MAX_ITER=maxiter)
}
#####################################



# CUSTOM RED-TRAINING AND PREDICTION
#####################################
custom_red_and_pred <- function(fullmodel,xtrain,target,xtest,ngenesselect,
                                lambda_red,epsilon_red,maxiter,nfolds){
  set.seed(1234)
  # extract feat weights from full model and select top X depending on the dataset for testing
  w <- t(attr(fullmodel,"gradient")) 
  # check which genes are common between w and xtest and, based on this, select top X model genes (most negative or positive weights)
  common <- intersect(colnames(xtest),rownames(w))
  w_f <- as.matrix(w[common,])

  # feat selection on weights
  w_f <- w_f[order(w_f),] ; w_f <- w_f[c(1:ngenesselect,(length(w_f)-(ngenesselect-1)):length(w_f))] ; print(length(w_f))
  # subset xtrain and xtest by feats selected
  xtrain_red <- xtrain[,names(w_f)] ; print(dim(xtrain_red))
  xtest_red <- xtest[,names(w_f)] ; print(dim(xtest_red))
  
  # another option for adding cols with 0s in case genes in test are not enough
  # x2 <- cbind(t(as.matrix(FerretInt@assays$integrated@scale.data))[,common], t(matrix(0, length(attr(red_AGE_diff1[[2]],"genenames"))-length(common), length(colnames(FerretInt) ))  ))
  
  # # run reduced model within training species (known target) and find prediction medians and sds to correct prediction in testing data (unknown targets)
  # pred <- nrbm(ordinalRegressionLoss(xtrain_red,target,C=cost),LAMBDA = 0.1,EPSILON_TOL =1e-7 ,MAX_ITER = 1000)
  # df <- data.frame(target=age_sel_diff1$y_age,pred=predict(pred,xtrain_red) )
  # df <- df %>% dplyr::group_by(target) %>% summarise_at(vars(pred), list(min=min,max=max,median=median,sd=sd))
  
  costM <- function(X,y,...) {
    C <- costMatrix(y)
    C <- C / tabulate(y)[col(C)] / tabulate(y)
    C <- C/sum(C)
  }
  cost=costM(X=xtrain_red,y=target) ; print(length(target)) ; print(dim(cost))
  
  # run reduced model with crossval and predict
  folds <- balanced.cv.fold(target,nfolds)
  pred <- simplify2array(mclapply(levels(folds),mc.cores=16,function(f) {
    w <- nrbm(ordinalRegressionLoss(xtrain_red[folds!=f,],target[folds!=f],C=cost),LAMBDA = lambda_red,EPSILON_TOL =epsilon_red ,MAX_ITER = maxiter)
    Y <- predict(w,xtest_red) ; Y[folds!=f] <- NA ; Y
  }))
  # return prediction result aggregated across CVs
  pred <- rowSums(pred,na.rm=TRUE)
  
  # correct prediction weights by medians and sds on training data
  
  
}

#####################################
# age_sel_diff1$ordi_age_diff1 <- custom_red_and_pred(fullmodel=full_AGE_diff1,xtrain=x,target=age_sel_diff1$y_age,xtest=x,ngenesselect=50,
#                                                     lambda_red=0.05,epsilon_red=1e-7,maxiter=1000,nfolds=5)
# 
# folds <- balanced.cv.fold(target,nfolds)
# pred <- simplify2array(mclapply(levels(folds),mc.cores=16,function(f) {
#   w <- nrbm(ordinalRegressionLoss(xtrain_red[folds!=f,],target[folds!=f],C=cost),LAMBDA = lambda_red,EPSILON_TOL =epsilon_red ,MAX_ITER = maxiter)
#   Y <- predict(w,xtest_red) ; Y[folds!=f] <- NA ; Y
# }))

# 
# hinge <- function(Xtrain,logi,Xtest){
#   folds= balanced.cv.fold(logi,num.cv=2)
#   hinge_precs <- mclapply(levels(folds),mc.cores=1,function(f) {
#     w <- nrbm(hingeLoss(x=Xtrain[folds!=f,],y=logi[folds!=f]),LAMBDA=0.001,MAX_ITER=5000)
#     Y <- predict(w,Xtest)
#     Y <- attr(Y,"decision.value") #; Y[folds!=f] <- NA 
#   })
# }
# 
# res <- list()
# for (i in levels(as.factor(target))){
#   res[[i]] <- hinge(Xtrain=xtrain_red,
#                     logi=ifelse(age_sel_diff1$y_age==i,TRUE,FALSE),
#                     Xtest=xtest[,rownames(w_f)])
# }
# 
# 
# 
# res <- list()
# for (i in levels(as.factor(target))){
#   res[[i]] <- hinge(Xtrain=xtrain_red,
#                     logi=ifelse(age_sel_diff1$y_age==i,TRUE,FALSE),
#                     Xtest=xtrain_red)
# }
# 
# 
# a <- as.data.frame(res$`12`)
# table(a$w.1>0)


