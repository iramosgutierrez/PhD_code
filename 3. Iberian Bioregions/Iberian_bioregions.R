library(ape)
library(phyloregion)
library(terra)
library(rgdal)
library(stringr)

#Find optimal K####
resolutions <- c("20_t", "20_e")

for(resolution in resolutions){
  
  res <- stringr::word(resolution, 1, sep="_")
  group <-stringr::word(resolution, 2, sep="_")
  
  
  if(group=="e"){matrix <- read.csv(paste0("3. Iberian Bioregions/lists/Matrix_endemic_",res,"km.csv"), check.names = F)}
  if(group=="t"){matrix <- read.csv(paste0("3. Iberian Bioregions/lists/Matrix_total_",res,"km.csv"), check.names = F)}
  
  
  comm.grids<- as.vector(matrix[,1])
  comm.species<- as.vector(names(matrix[-1]))
  commvector<-vector(mode = "numeric")
  for(i in 2:ncol(matrix)){
    vect<-as.numeric( matrix[,i])
    commvector<- c(commvector, vect)}
  dense_comm <- matrix(commvector , length(comm.grids),length(comm.species),
                       dimnames=list(comm.grids,comm.species))
  
  
  x <-phyloregion::dense2sparse( dense_comm)
  pbc<- phyloregion::beta_diss(x,  index.family = "s")
  
  optimum <- phyloregion::optimal_phyloregion(pbc[[1]], k = 50)
  k <- optimum$optimal$k

  print(paste0("Resolution ", resolution, "(taxo). Optimum k=", k))
  
  
  reglist <- list()
  for( t in 1:100){
    tree <- read.tree(paste0("3. Iberian Bioregions/phylogenies/rtree_",t,"_ed.tre"))
    tree<- phyloregion::phylobuilder(species = comm.species, tree = tree, extract = T)
    x <-phyloregion::dense2sparse( dense_comm)
    subphy <- phyloregion::match_phylo_comm(phy=tree, comm=x)$phy
    submat <- phyloregion::match_phylo_comm(tree, x)$com
    pbc <- phyloregion::phylobeta(submat, subphy, index.family = "s")
    reglist[[t]]<- pbc[[1]]
  }
  
  
  mean.distable <- reglist[[1]] 
  dists <- rowMeans(sapply(reglist, '['))
  mean.distable[1:length(mean.distable)] <- dists[1:length(mean.distable)]
  opt <- optimal_phyloregion(mean.distable, k=50)
  k.n <- opt$optimal$k
  
  print(paste0("Resolution ", resolution, "(phylo). Optimum k=", k.n))
  
}


#Obtain regions####
resolutions <- c("20_e", "20_t")#"50_e", "50_t", "20_e", "20_t"

for(resolution in resolutions){
  res <- stringr::word(resolution, 1, sep="_")
  group <-stringr::word(resolution, 2, sep="_")
  
  if(resolution=="50_e") {k.tax<-8 ; k.phy<-13}
  if(resolution=="50_t") {k.tax<-13; k.phy<-14}
  if(resolution=="20_e") {k.tax<-8 ; k.phy<-10}
  if(resolution=="20_t") {k.tax<-4 ; k.phy<-11}
  if(resolution=="10_e") {k.tax<-40; k.phy<-11}
  # if(resolution=="10_t") {k.tax<-  ; k.phy<-  }
  
  grid <- readOGR(paste0("3. Iberian Bioregions/grids/UTM",res,"x",res,".shp"))
  names(grid@data)<-"grids"
  
  
  if(group=="e"){matrix <- read.csv(paste0("3. Iberian Bioregions/lists/Matrix_endemic_",res,"km.csv"), check.names = F)}
  if(group=="t"){matrix <- read.csv(paste0("3. Iberian Bioregions/lists/Matrix_total_",res,"km.csv"), check.names = F)}
  
  
  comm.grids<- as.vector(matrix[,1])
  comm.species<- as.vector(names(matrix[-1]))
  commvector<-vector(mode = "numeric")
  for(i in 2:ncol(matrix)){
    vect<-as.numeric( matrix[,i])
    commvector<- c(commvector, vect)}
  dense_comm <- matrix(commvector , length(comm.grids),length(comm.species),
                       dimnames=list(comm.grids,comm.species))
  
  
  x <-phyloregion::dense2sparse( dense_comm)
  pbc<- phyloregion::beta_diss(x,  index.family = "s")
  
  y <- phyloregion::phyloregion(x = pbc[[1]], shp = grid, method="average", k = k.tax)
  # save(y, file=paste0("results/regions_", resolution, "_taxo.Rda"))
  unis <- as.numeric(colSums(table(y$membership)))
  cexuni <- unis; cexuni[unis<6]<- 2;   cexuni[unis>5]<- 3.5
  pchuni <- unis; pchuni[unis<6]<- 1.2; pchuni[unis>5]<- 2
  
  
  png(paste0("3. Iberian Bioregions/results/",res,"_",group,"_taxo.png"), width = 210, 
      height = 210, units = "mm", res=600)
  par(mar=c(2,0.75,0.75,0.75), mfrow=c(1,1))
  layout(mat = matrix(c(1,1,1,1,1,1, 1,1,1,1,1,1, 
                        1,1,1,1,1,1, 1,1,1,1,1,1, 
                        3,3,1,1,2,2, 3,3,1,1,2,2),
                      nrow = 6, ncol = 6))
  
  phyloregion::plot.phyloregion(y, cex=cexuni[y$shp@data$cluster]/2)
  
  par(mar=c(2,0,0,2))
  cols <- y$region.df[!(duplicated(y$region.df$cluster)),]
  dendro <- as.phylo(hclust(y$region.dist))
  plot.phylo(dendro,show.tip.label = F, direction = "downwards")
  tiplabels(text="",frame = "none", pch=16, cex=cexuni*1.25, col= cols$COLOURS)
  tiplabels(frame = "none", cex=cexuni/2.5)
  
  
  par(mar=c(6,6,3,3))
  phyloregion::plot_NMDS(y, cex=pchuni*2,xlab="", ylab="", xaxt="n", yaxt="n")
  axis(3);axis(4)
  phyloregion::text_NMDS(y, cex=pchuni*0.8, font=3)
  
  dev.off()
  
  
  reglist <- list()
  for( t in 1:100){
    tree <- read.tree(paste0("3. Iberian Bioregions/phylogenies/rtree_",t,"_ed.tre"))
    tree<- phyloregion::phylobuilder(species = comm.species, tree = tree, extract = T)
    x <-phyloregion::dense2sparse( dense_comm)
    subphy <- phyloregion::match_phylo_comm(phy=tree, comm=x)$phy
    submat <- phyloregion::match_phylo_comm(tree, x)$com
    pbc <- phyloregion::phylobeta(submat, subphy, index.family = "s")
    reglist[[t]]<- pbc[[1]]
  }
  mean.distable <- reglist[[1]] 
  dists <- rowMeans(sapply(reglist, '['))
  mean.distable[1:length(mean.distable)] <- dists[1:length(mean.distable)]
  
  k.n <- k.phy
  y <- phyloregion(mean.distable, k=k.n,  shp=grid)
  # save(y, file=paste0("results/regions_", resolution, "_phylo.Rda"))
  unis <- as.numeric(colSums(table(y$membership)))
  cexuni <- unis; cexuni[unis<6]<- 2;   cexuni[unis>5]<- 3.5
  pchuni <- unis; pchuni[unis<6]<- 1.2; pchuni[unis>5]<- 2
  
  png(paste0("results/", resolution, "_phylo.png"), 
      height = 210,width = 210, units="mm", res=600 )
  par(mar=c(2,0.75,0.75,0.75), mfrow=c(1,1))
  
  layout(mat = matrix(c(1,1,1,1,1,1, 1,1,1,1,1,1, 
                        1,1,1,1,1,1, 1,1,1,1,1,1, 
                        3,3,1,1,2,2, 3,3,1,1,2,2),
                      nrow = 6, ncol = 6))
  
  
  phyloregion::plot.phyloregion(y, cex=cexuni[y$shp@data$cluster]/2)
  
  par(mar=c(2,0,0,2))
  cols <- y$region.df[!(duplicated(y$region.df$cluster)),]
  dendro <- as.phylo(hclust(y$region.dist))
  plot.phylo(dendro,show.tip.label = F, direction = "downwards")
  tiplabels(text="",frame = "none", pch=16, cex=cexuni*1.25, col= cols$COLOURS)
  tiplabels(frame = "none", cex=cexuni/2.5)
  
  
  par(mar=c(6,6,3,3))
  phyloregion::plot_NMDS(y, cex=pchuni*2,xlab="", ylab="", xaxt="n", yaxt="n")
  axis(3);axis(4)
  phyloregion::text_NMDS(y, cex=pchuni*0.8, font=3)
  
  
  dev.off()
}

#Comparison with variables####

# see variable_comparisons.R

load("rasters/lda_20_e_taxo.Rda")
copyldatable <- function(ld){
  
  print(paste0("AUC: ", round(ld$mean.AUC, digits = 3)))
  restable <- data.frame(matrix(nrow=7, ncol=2, dimnames = list(NULL, c("var", "vals"))))
  
  table <- ld$summary.table
  rownames(table) <- table$var
  
  # importances <- abs(table$mean)
  # import.cum<- cumsum(importances)
  # import.sum <- sum(importances)
  # thres <- which(import.cum>import.sum*0.7)[1]
  # if(thres==0){thres<-1}
  # table$var[1:thres] <- paste0(table$var[1:thres], "*")
  
  
  
  if(!any(c("isolation","isolation*")%in% table$var)){
    isotable <- data.frame("var"="isolation", "mean"=NA, "sd"=NA)
    rownames(isotable)<-"isolation"
    table <- rbind(table, isotable)
  }
  
  # table <- table[c("isolation","alt.heter", "ph.mn", "bio1.mn", "bio4.mn", "bio12.mn", "bio15.mn"),]
  
  restable$var <- table$var
  
  restable$vals <- paste0(round(table$mean, digits = 3 ), "\U00B1", round(table$sd, digits = 3 ))
  
  
  write.table(restable, "clipboard", row.names = F, col.names = F, sep="\t", quote=F)
  
}
copyldatable(ld1)
copyldatable(ld2)
copyldatable(ld3)
copyldatable(ld4)
rm(list=ls())

load("rasters/lda_20_e_phylo.Rda")
copyldatable <- function(ld){
  
  print(paste0("AUC: ", round(ld$mean.AUC, digits = 3)))
  restable <- data.frame(matrix(nrow=7, ncol=2, dimnames = list(NULL, c("var", "vals"))))
  
  table <- ld$summary.table
  rownames(table) <- table$var
  
  # importances <- abs(table$mean)
  # import.cum<- cumsum(importances)
  # import.sum <- sum(importances)
  # thres <- which(import.cum>import.sum*0.7)[1]
  # if(thres==0){thres<-1}
  # table$var[1:thres] <- paste0(table$var[1:thres], "*")
  
  
  
  if(!any(c("isolation","isolation*")%in% table$var)){
    isotable <- data.frame("var"="isolation", "mean"=NA, "sd"=NA)
    rownames(isotable)<-"isolation"
    table <- rbind(table, isotable)
  }
  
  # table <- table[c("isolation","alt.heter", "ph.mn", "bio1.mn", "bio4.mn", "bio12.mn", "bio15.mn"),]
  
  restable$var <- table$var
  
  restable$vals <- paste0(round(table$mean, digits = 3 ), "\U00B1", round(table$sd, digits = 3 ))
  
  
  write.table(restable, "clipboard", row.names = F, col.names = F, sep="\t", quote=F)
  
}
copyldatable(ld1)
copyldatable(ld2)
copyldatable(ld3)
copyldatable(ld4)
copyldatable(ld5)
rm(list=ls())

load("rasters/lda_20_t_taxo.Rda")
copyldatable <- function(ld){
  
  print(paste0("AUC: ", round(ld$mean.AUC, digits = 3)))
  restable <- data.frame(matrix(nrow=7, ncol=2, dimnames = list(NULL, c("var", "vals"))))
  
  table <- ld$summary.table
  rownames(table) <- table$var
  
  # importances <- abs(table$mean)
  # import.cum<- cumsum(importances)
  # import.sum <- sum(importances)
  # thres <- which(import.cum>import.sum*0.7)[1]
  # if(thres==0){thres<-1}
  # table$var[1:thres] <- paste0(table$var[1:thres], "*")
  
  
  
  if(!any(c("isolation","isolation*")%in% table$var)){
    isotable <- data.frame("var"="isolation", "mean"=NA, "sd"=NA)
    rownames(isotable)<-"isolation"
    table <- rbind(table, isotable)
  }
  
  # table <- table[c("isolation","alt.heter", "ph.mn", "bio1.mn", "bio4.mn", "bio12.mn", "bio15.mn"),]
  
  restable$var <- table$var
  
  restable$vals <- paste0(round(table$mean, digits = 3 ), "\U00B1", round(table$sd, digits = 3 ))
  
  
  write.table(restable, "clipboard", row.names = F, col.names = F, sep="\t", quote=F)
  
}
copyldatable(ld1)
copyldatable(ld2)
copyldatable(ld3)
rm(list=ls())

load("rasters/lda_20_t_phylo.Rda")
copyldatable <- function(ld){
  
  print(paste0("AUC: ", round(ld$mean.AUC, digits = 3)))
  restable <- data.frame(matrix(nrow=7, ncol=2, dimnames = list(NULL, c("var", "vals"))))
  
  table <- ld$summary.table
  rownames(table) <- table$var
  
  # importances <- abs(table$mean)
  # import.cum<- cumsum(importances)
  # import.sum <- sum(importances)
  # thres <- which(import.cum>import.sum*0.7)[1]
  # if(thres==0){thres<-1}
  # table$var[1:thres] <- paste0(table$var[1:thres], "*")
  
  
  
  if(!any(c("isolation","isolation*")%in% table$var)){
    isotable <- data.frame("var"="isolation", "mean"=NA, "sd"=NA)
    rownames(isotable)<-"isolation"
    table <- rbind(table, isotable)
  }
  
  # table <- table[c("isolation","alt.heter", "ph.mn", "bio1.mn", "bio4.mn", "bio12.mn", "bio15.mn"),]
  
  restable$var <- table$var
  
  restable$vals <- paste0(round(table$mean, digits = 3 ), "\U00B1", round(table$sd, digits = 3 ))
  
  
  write.table(restable, "clipboard", row.names = F, col.names = F, sep="\t", quote=F)
  
}
copyldatable(ld1)
copyldatable(ld2)
copyldatable(ld3)


#.