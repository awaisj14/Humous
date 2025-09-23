

# LIBRARY FOR LANDSCAPES 

library(IRanges)
library(ggplot2)
library(DelayedArray)
library(SummarizedExperiment)
library(rasterVis)
library(ggplot2)
library(ggmap)
library(scales)


rowmean <- function(M,f,na.rm=TRUE) {
  f <- as.factor(f)
  M <- rowsum(M,f)
  M[levels(f),,drop=TRUE] / as.vector(table(f))
}

#' return a lowres integer array of nearest neighbors - grid
knn_array_lowres <- function(x,y,k,gx = seq(min(x),max(x),length.out=50L),gy = seq(min(y),max(y),length.out=50L)) { # grid of size gx*gy averaging k cells (neighbors)
  g <- matrix(NA,length(gy),length(gx))
  knn <- FNN::get.knnx(cbind(x,y),cbind(gx[col(g)],gy[row(g)]),k=k)$nn.index
  knn <- array(knn,c(dim(g),k))
  attr(knn,"x") <- gx
  attr(knn,"y") <- gy
  knn
}



#' return a mediumres integer array of nearest neighbors - grid
knn_array_medres <- function(x,y,k,gx = seq(min(x),max(x),length.out=250L),gy = seq(min(y),max(y),length.out=250L)) { # grid of size gx*gy averaging k cells (neighbors)
  g <- matrix(NA,length(gy),length(gx))
  knn <- FNN::get.knnx(cbind(x,y),cbind(gx[col(g)],gy[row(g)]),k=k)$nn.index
  knn <- array(knn,c(dim(g),k))
  attr(knn,"x") <- gx
  attr(knn,"y") <- gy
  knn
}



#' return a highres integer array of nearest neighbors - grid
knn_array_highres <- function(x,y,k,gx = seq(min(x),max(x),length.out=500L),gy = seq(min(y),max(y),length.out=500L)) { # grid of size gx*gy averaging k cells (neighbors)
  g <- matrix(NA,length(gy),length(gx))
  knn <- FNN::get.knnx(cbind(x,y),cbind(gx[col(g)],gy[row(g)]),k=k)$nn.index
  knn <- array(knn,c(dim(g),k))
  attr(knn,"x") <- gx
  attr(knn,"y") <- gy
  knn
}


#' Split rows into batches of 100 genes
#' Compute average expression of each gene in m according to index structure found in i
#'@param m an expression matrix (row=gene,col=cell)
#'@param i an array of column indexes, to be reduced over the last dimension
knn_rowMeans <- function(m,i) {
  B <- local({
    batches <- split(seq(nrow(m)),seq(nrow(m)) %/% 100L)
    B <- lapply(batches,function(b){
      print(paste("batch",b,"out of", (length(batches)*100) )) # progress report
      M <- m[b,as.vector(i)]
      M <- array(M,c(length(b),dim(i)))
      DelayedArray(rowMeans(M,dims=length(dim(i))))
    })
    do.call(rbind,B)
  })       
}


# plot the average batch of gene
gg_landscape <- function(B,x=seq(dim(B)[3]),y=seq(dim(B)[2])) {
  D <- reshape2::melt(as.array(B))
  D$x <- x[D$Var3]
  D$y <- y[D$Var2]
  ggplot(D) + facet_wrap(~Var1) +
    geom_tile(aes(x=x,y=y,fill=value)) +
    scale_fill_gradientn(colors = c("#35978f", "#80cdc1", "#c7eae5", "#f5f5f5","#f6e8c3", "#dfc27d", "#bf812d", "#8c510a","#543005", "#330000")) +
    #scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11,"RdYlBu"))) +
    coord_equal() +
    #guides(fill="none") +
    xlab("") + ylab("") + theme(panel.grid = element_blank())
}

gg_landscape_image <- function(m,x=seq(dim(m)[2]),y=seq(dim(m)[1])) {
  D <- reshape2::melt(m)
  D$x <- x[D$Var2]
  D$y <- y[D$Var1]
  ggplot(D) + 
    geom_tile(aes(x=x,y=y,fill=value)) +
    scale_fill_gradientn(colors = c("#35978f", "#80cdc1", "#c7eae5", "#f5f5f5","#f6e8c3", "#dfc27d", "#bf812d", "#8c510a","#543005", "#330000")) +
    #scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11,"RdYlBu"))) +
    coord_equal() +
    theme_void() + theme(legend.position = "none")
}


gg_landscape_blues <- function(B,x=seq(dim(B)[3]),y=seq(dim(B)[2])) {
  D <- reshape2::melt(as.array(B))
  D$x <- x[D$Var3]
  D$y <- y[D$Var2]
  ggplot(D) + facet_wrap(~Var1) +
    geom_tile(aes(x=x,y=y,fill=value)) +
    #scale_fill_gradientn(colors = c("#35978f", "#80cdc1", "#c7eae5", "#f5f5f5","#f6e8c3", "#dfc27d", "#bf812d", "#8c510a","#543005", "#330000")) +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9,"Blues")) +
    coord_equal() +
    #guides(fill="none") +
    xlab("") + ylab("") + theme(panel.grid = element_blank())
}

# plot landscapes for one gene on multispecies
multispecies_L <- function(gene){
  h <- as.array(landscapes_MR_RNA$H_M[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster() + scale_y_continuous(breaks = c(-0.95,-0.7,-0.4,-0.15),labels = c("aRG","oRG/bRG/IPC","iN","mN"))  + scale_x_continuous(breaks = c(0,1),labels = c("10pcw","24pcw"))
  m <- as.array(landscapes_MR_RNA$M_M[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster() + scale_y_continuous(breaks = c(-0.95,-0.7,-0.4,-0.15),labels = c("aRG","IPC","iN","mN")) + scale_x_continuous(breaks = c(0,1),labels = c("E12","E17"))
  o <- as.array(landscapes_MR_RNA$O_M[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() %>% base_ggraster() + scale_y_continuous(breaks = c(-0.95,-0.7,-0.4,-0.15),labels = c("aRG","oRG/bRG/IPC","iN","mN")) + scale_x_continuous(breaks = c(0,1),labels = c("4wiv","37w"))
  p <- gridExtra::grid.arrange(h,m,o,nrow=1)
  p
}

# plot substraction maps for one gene on multispecies
multispecies_SM <- function(gene){
  a <- levelplot(flip(
    raster::overlay( (as.array(landscapes_MR_RNA$H_M[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,
                     ( as.array(landscapes_MR_RNA$M_M[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,fun=function(x,y){ (scale(x)-scale(y)) })
    ,2),at=seq(-5,5,0.1),xlab="Developmental stage (time)",ylab="Differentiation stage (cell type)",scales=list(
      x=list(at=c(0,1),labels=c("min","max")),y=list(at=c(0,1),labels= c("RG","N"))),
    par.settings = rasterTheme(region = c("#F7D452","#F7D452","white","white","#B797AE","#B797AE")))
  b <- levelplot(flip(
    raster::overlay( (as.array(landscapes_MR_RNA$H_M[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,
                     ( as.array(landscapes_MR_RNA$O_M[gene,,]) %>% scales::rescale(to=c(0,1)) %>% raster::raster() ) ,fun=function(x,y){ (scale(x)-scale(y)) })
    ,2),at=seq(-5,5,0.1),xlab="Developmental stage (time)",ylab="Differentiation stage (cell type)",scales=list(
      x=list(at=c(0,1),labels=c("min","max")),y=list(at=c(0,1),labels= c("RG","N"))),
    par.settings = rasterTheme(region = c("#EA7D94","#EA7D94","white","white","#B797AE","#B797AE")))
  p <- gridExtra::grid.arrange(a,b,nrow=1)
  p
}

# generate a random dataset
random_SE <- function(ngene=5000,ncell=1000) {
  x <- runif(ncell)
  y <- runif(ncell)
  m <- matrix(runif(ngene*ncell),ngene,ncell)
  SummarizedExperiment(m,colData=DataFrame(x=x,y=y))
}


# transform landscapes object to matrix
landscape_to_matrix <- function(landscape){
  M <- as.array(landscape)
  dim(M) <- c(dim(M)[1],prod(dim(M)[-1]))
  M <- t(scale(t(M)))
  M[is.na(M)] <- 0
  M
}


# function to plot one raster landscape as gplot object
base_ggraster <- function(obj){
  gplot(obj) + geom_tile(aes(x,-y,fill = -value)) +
    coord_equal() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9,"Blues"))) + theme_bw() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),
      axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"))),
      axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"))),
      axis.text.y = element_text(angle=90,vjust=0,hjust=0.5),
      axis.ticks.y = element_blank()
    ) + 
    labs(fill='') + 
    xlab("Developmental stage (time)") + 
    ylab("Differentiation stage (cell type)") 
}



# function to generate landscapes for website
website_ggraster <- function(obj,title){
  gplot(obj) + geom_tile(aes(x,-y,fill = value)) +
    coord_equal() + scale_fill_gradientn(colors = (RColorBrewer::brewer.pal(9,"Blues"))) + theme_bw(base_size = 30) +
    theme(
      legend.position="none", #axis.title=element_text(size=30 ),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),
      axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"))), 
      axis.text.x=element_blank(), axis.text.y=element_blank(),
      axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"))),
      axis.ticks.y = element_blank(),axis.ticks.x = element_blank()) +
    labs(fill='') + xlab("Age (time)") + ylab("Differentiation (cell type)") 
}


# function to generate landscapes for deeplearning
ML_ggraster <- function(obj){
  par(mar=c(0,0,0,0))
  gplot(obj) + geom_tile(aes(x,-y,fill = -value)) +
    coord_equal() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9,"Blues"))) + theme_nothing() + labs(x = NULL, y = NULL) +scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) 
}



# function to calculate correlation between two landscapes using boostraping and thresholding
    # returns a list containing correlation coeficient, correlation sd and n pixels considered after thresholding (area occupied by the gene in the landscape)
L_TH_BoosCorr <- function(gene1,gene2,threshold,n_boosts,boostsize,L1,L2){
  tryCatch({
    df <- data.frame(L1=as.vector(L1[gene1,,])%>%rescale(to=c(0,1)),L2=as.vector(L2[gene2,,])%>%rescale(to=c(0,1)),stringsAsFactors=FALSE)
    df$touse <- ifelse(df$L1>threshold | df$L2>threshold,"use","notuse") ; df <- df[df$touse=="use",]  # keep just pixels to be used 
    set.seed(1234) ; boos_idx <- replicate(n_boosts,sample(nrow(df),boostsize,replace=TRUE)) # corrs are influenced by sample size so we set fixed size and do 10 boostraps
    boostcorrs <- vector(mode="integer",length=ncol(boos_idx))
    for (i in 1:ncol(boos_idx)){
      df_f <- df[boos_idx[,i],]
      boostcorrs[i] <- cor(df_f$L1,df_f$L2)
    } 
    corr <- median(boostcorrs) 
    return(list(corr=corr,corrsd=sd(boostcorrs),dimTH_L=nrow(df)))
    #p <- ggplot(df,aes(L1,L2)) + geom_point() + geom_smooth(method="lm")
    #lmcoef <- coef(lm(df$L1~df$L2))[2] # slope of linear fit
    print(corr) ; print(dim(df)) #; print(summary(boostcorrs)) #; print(as.numeric(lmcoef)) ; print(p)
  },error = function(e) {
    message("something went wrong")
    message(e)
  } )
}


