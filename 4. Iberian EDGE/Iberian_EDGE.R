devtools::install_github("iramosgutierrez/RIRG")
library(RIRG)
library(openxlsx)
library(rgdal)
library(ape)
library(picante)
library(sf)
library(ggplot2)

#

#0. own functions####
plotcontinuous <- function(shp, col, pal="maxent", breaks=10, round=1, ...){
  library(RColorBrewer)
  library(rgdal)
  if(is.null(pal)){pal <- NA}
  if(!(pal %in% hcl.pals())){cols <- colorRampPalette(colors = c("purple", "blue", "green","yellow", "red"))}
  if(pal %in% hcl.pals()) {cols <- colorRampPalette(hcl.colors(n=10, palette = pal))}
  if(pal=="maxent"){cols <- colorRampPalette(colors = c("#2b83ba", "#abdda4", "#ffffbf", "#fdae61", "#d7191c"))}
  # Rank variable for colour assignment
  shp <- shp[!is.na(shp@data[,col]),]
  shp$order = findInterval(shp@data[,col], sort(shp@data[,col]))
  # Make plot
  cols.plot <- cols(nrow(shp))[shp$order]
  cols.plot[shp@data[,col]==0] <- "#DBDBDB"
  plot(shp, pch=19, col=cols.plot, ...)
  # Add a simple legend
  legend("topright", col=cols(breaks), pch=15,
         legend=round(seq(range(shp@data[,col])[1], range(shp@data[,col])[2],
                          length.out=breaks),round), cex=0.85 , bty='n' )
}
#3. Plotting PD and SR      ####

#tree preparation: see cluster analysis tree_preparation.R

#Phylogenetic indexes calculations: see cluster analyses Phylo.R 
res <- "10"

afliber.orig <- read.csv(paste0("4. Iberian EDGE/files/lists/afliber_matrix_", res,"x",res,".csv"))
colnames(afliber.orig) <- gsub("\\.", "-", colnames(afliber.orig))


for(i in 1:100){load(paste0("4. Iberian EDGE/forcluster/output/",res,"/phylo_",res,"_", i, ".rda"))}

if(res=="10"){
  pd.list <-   NULL; for(i in 1:100){pd.list <- c(pd.list, get(paste0("pd.10.list_",i )))
  rm(list = paste0("pd.10.list_",i ))}
  
  mpd.list <-  NULL; for(i in 1:100){mpd.list <- c(mpd.list, get(paste0("mpd.10.list_",i )))
  rm(list = paste0("mpd.10.list_",i ))}
  
  mntd.list <- NULL; for(i in 1:100){mntd.list <- c(mntd.list, get(paste0("mntd.10.list_",i )))
  rm(list = paste0("mntd.10.list_",i ))}
}


#obs PD
{
  pdo.obs <- sapply(pd.list, '[', 2)
  pdo.pval<- sapply(pd.list, '[', 7)
  
  pdo.obs.df <- as.data.frame(sapply(pdo.obs, '['))
  colnames(pdo.obs.df) <- paste0("pdo.obs_t", 1:length(pdo.obs.df))
  pdo.obs.df$mean <- rowMeans(pdo.obs.df)
  
  pdo.pval.df <- as.data.frame(sapply(pdo.pval, '['))
  colnames(pdo.pval.df) <- paste0("pdo.pval_t", 1:length(pdo.pval.df))
  pdo.pval.df$mean <- rowMeans(pdo.pval.df)
}

#ses PD
{
  pd.obs <- sapply(pd.list, '[', 6)
  pd.pval<- sapply(pd.list, '[', 7)
  
  pd.obs.df <- as.data.frame(sapply(pd.obs, '['))
  colnames(pd.obs.df) <- paste0("pd.obs_t", 1:length(pd.obs.df))
  pd.obs.df$mean <- rowMeans(pd.obs.df)
  
  pd.pval.df <- as.data.frame(sapply(pd.pval, '['))
  colnames(pd.pval.df) <- paste0("pd.pval_t", 1:length(pd.pval.df))
  pd.pval.df$mean <- rowMeans(pd.pval.df)
}

#MPD
{
  mpd.obs <- sapply(mpd.list, '[', 6)
  mpd.pval<- sapply(mpd.list, '[', 7)
  
  mpd.obs.df <- as.data.frame(sapply(mpd.obs, '['))
  colnames(mpd.obs.df) <- paste0("mpd.obs_t", 1:length(mpd.obs.df))
  mpd.obs.df$mean <- rowMeans(mpd.obs.df)
  
  mpd.pval.df <- as.data.frame(sapply(mpd.pval, '['))
  colnames(mpd.pval.df) <- paste0("mpd.pval_t", 1:length(mpd.pval.df))
  mpd.pval.df$mean <- rowMeans(mpd.pval.df)
}

#MNTD
{
  mntd.obs <- sapply(mntd.list, '[', 6)
  mntd.pval<- sapply(mntd.list, '[', 7)
  
  mntd.obs.df <- as.data.frame(sapply(mntd.obs, '['))
  colnames(mntd.obs.df) <- paste0("mntd.obs_t", 1:length(mntd.obs.df))
  mntd.obs.df$mean <- rowMeans(mntd.obs.df)
  
  mntd.pval.df <- as.data.frame(sapply(mntd.pval, '['))
  colnames(mntd.pval.df) <- paste0("mntd.pval_t", 1:length(mntd.pval.df))
  mntd.pval.df$mean <- rowMeans(mntd.pval.df)
}

phylo <- data.frame("id"=(afliber.orig[,1]), 
                    "pdo.obs.mn"= pdo.obs.df$mean, "pdo.pval.mn"= pdo.pval.df$mean, 
                    "pd.obs.mn"=  pd.obs.df$mean,  "pd.pval.mn"=  pd.pval.df$mean,  
                    "mpd.obs.mn"= mpd.obs.df$mean, "mpd.pval.mn"= mpd.pval.df$mean,
                    "mntd.obs.mn"=mntd.obs.df$mean,"mntd.pval.mn"=mntd.pval.df$mean  )
nacols <- which(apply(apply(phylo[,2:7], 2, is.na), 1, any))
if(length(nacols)>0){phylo <- phylo[-nacols,]}

grid <- rgdal::readOGR(paste0("files/grids/UTM", res,"x",res, ".shp"))
grid <- RIRG::joinAttributeTable(grid, phylo, names(grid)[1], "id")
plotcontinuous <- function(shp, col, pal=NULL, breaks=10, round=1, ...){
  library(RColorBrewer)
  library(rgdal)
  if(is.null(pal)){pal <- NA}
  if(!(pal %in% hcl.pals())){cols <- colorRampPalette(colors = c("purple", "blue", "green","yellow", "red"))}
  if(pal %in% hcl.pals()) {cols <- colorRampPalette(hcl.colors(n=10, palette = pal))}
  if(pal=="maxent"){cols <- colorRampPalette(colors = c("#2b83ba", "#abdda4", "#ffffbf", "#fdae61", "#d7191c"))}
  # Rank variable for colour assignment
  shp <- shp[!is.na(shp@data[,col]),]
  shp$order = findInterval(shp@data[,col], sort(shp@data[,col]))
  # Make plot
  plot(shp, pch=19, col=cols(nrow(shp))[shp$order], ...)
  # Add a simple legend
  legend("topright", col=cols(breaks), pch=15,
         legend=round(seq(range(shp@data[,col])[1], range(shp@data[,col])[2],
                          length.out=breaks),round), cex=0.85 , bty='n' )
}

{
  png(paste0("4. Iberian EDGE/results/1. phylo/PDobs_", res, ".png"), height = 150, width = 150, units = "mm",
      res=450)
  par(mar=c(.1, .1, .1, .1))
  plotcontinuous(grid, col="pdo.obs.mn",    pal="maxent", border="transparent" )#,main="PD OBSERVED",cex=0.75 
  dev.off()
}#png PD obs

{
  png(paste0("4. Iberian EDGE/results/1. phylo/PD_", res, ".png"), height = 150, width = 150, units = "mm",
      res=450)
  par(mar=c(.1, .1, .1, .1))
  plotcontinuous(grid, col="pd.obs.mn",    pal="maxent", border="transparent" )#,main="PD OBSERVED",cex=0.75 
  cent.pd.l <- rgeos::gCentroid(grid[grid@data$pd.pval.mn<0.025,], byid = T)
  cent.pd.h <- rgeos::gCentroid(grid[grid@data$pd.pval.mn>0.975,], byid = T)
  plot(cent.pd.l, add=T, pch=4, col="#00000088", cex=log10(as.numeric(res))-0.4)
  plot(cent.pd.h, add=T, pch=16 , col="#00000088", cex=log10(as.numeric(res))-0.4)
  dev.off()
}#png PD ses



#4. EDGE Results ####
res <- "10"

afliber.orig <- read.csv(paste0("4. Iberian EDGE/files/lists/afliber_matrix_", res,"x",res,".csv"))
colnames(afliber.orig) <- gsub("\\.", "-", colnames(afliber.orig))
grid.iberia <- rgdal::readOGR(paste0("4. Iberian EDGE/files/grids/UTM", res,"x",res, ".shp"))


edge <- read.csv("4. Iberian EDGE/files/lists/EDGE_list_angio.csv")
edge <- edge[edge$Species%in%colnames(afliber.orig),]

cols <- which(colnames(afliber.orig) %in% edge$Species)
afliber.filt <- afliber.orig[,c(1,cols)]
cols2 <- rowSums(afliber.filt[,-1])
cols2 <- which(cols2>0)
afliber.filt <- afliber.filt[c(1,cols2),]

grid.edge <- grid.iberia
grid.edge@data$richness <- as.numeric(0)
grid.edge@data$edge.val <- as.numeric(0)
grid.edge@data$ratio <- as.numeric(0)
for( i in 1:nrow(afliber.filt)){
  grid <- afliber.filt[i,1]
  afliber.filt.i <- afliber.filt[i,]
  cols3 <- which(afliber.filt.i!="0")
  cols3 <- cols3[cols3!=1]
  spp <- colnames(afliber.filt.i[cols3])
  
  edge.i <- edge[edge$Species%in%spp,]
  edge.sum <- sum(edge.i$EDGE.med)
  grid.edge@data$richness[grid.edge@data[1]==grid] <- nrow(edge.i)
  grid.edge@data$edge.val[grid.edge@data[1]==grid] <- edge.sum
  grid.edge@data$ratio[grid.edge@data[1]==grid] <- edge.sum/nrow(edge.i)
}


{
  png(paste0("results/2. edge/EDGEsum_", res, ".png"), height = 150, width = 150, units = "mm",
      res=450)
  par(mar=c(.1, .1, .1, .1))
  plotcontinuous(grid.edge, col="edge.val",  pal="maxent" , border="transparent" )
  dev.off()
}#png EDGE values sum

{
  png(paste0("results/2. edge/EDGErich_", res, ".png"), height = 150, width = 150, units = "mm",
      res=450)
  par(mar=c(.1, .1, .1, .1))
  plotcontinuous(grid.edge, col="richness",  pal="maxent" , border="transparent" )
  dev.off()
}#png EDGE richness

res <- "10"
afliber.orig <- read.csv(paste0("4. Iberian EDGE/files/lists/afliber_matrix_", res,"x",res,".csv"))
colnames(afliber.orig) <- gsub("\\.", "-", colnames(afliber.orig))
matrix <- afliber.orig[,-1]
rownames(matrix)<- afliber.orig[,1]
colnames(matrix) <- colnames(afliber.orig[,-1])

edge <- read.csv("4. Iberian EDGE/files/lists/EDGE_total_results_summary.csv") 
edge <- edge[edge$Species%in%colnames(afliber.orig),]

edge <- edge[, c("Species","EDGE.med")]
colnames(edge)<- c("taxon", "value")

ses.edge <- function(matrix, edge, niter=999){
  
  nonedgespp <- colnames(matrix)[!(colnames(matrix)%in% edge$taxon)]
  nonedgespp.n <- which(colnames(matrix) %in% nonedgespp)
  if(length(nonedgespp.n)!=0){matrix <- matrix[,-nonedgespp.n]}
  
  
  
  results <- data.frame("id"=1:length(rownames(matrix)),
                        "grid"=NA_character_,
                        "obs.edge"=NA_real_,
                        "ses.edge"=NA_real_,
                        "sd.edge"=NA_real_,
                        "position"=NA_real_,
                        "edge.obs.z"=NA_real_,
                        "pval"=NA_real_,
                        "runs"=niter )
  for( i in 1:length(rownames(matrix))){
    RIRG::progressbar(i, length(rownames(matrix)),1,"secs")
    grid <- rownames(matrix)[i]
    results[i, "grid"]<- grid
    
    
    matrix.i <- matrix[i,]
    
    spp.i <- colnames(matrix.i[,matrix.i!=0])
    
    obs.edge <- sum(edge[edge$taxon%in%spp.i, "value"])
    results[i, "obs.edge"] <- obs.edge
    
    edge.iter <- as.numeric(NULL)
    for(j in 1:niter){
      spp.iter.j <- sample(colnames(matrix), size = length(spp.i), replace = F)
      edge.iter[j] <- sum(edge[edge$taxon%in%spp.iter.j, "value"])
    }
    
    results[i, "ses.edge"] <- mean(edge.iter)
    results[i, "sd.edge"] <- sd(edge.iter)
    edge.z <- (obs.edge-mean(edge.iter))/sd(edge.iter)
    results[i, "edge.obs.z"]<- edge.z
    
    ord.pos <- which(order(c(obs.edge, edge.iter))==1)
    
    results[i,"position"] <- ord.pos
    results[i,"pval"] <- ord.pos/(niter+1)
  } 
  results[results$obs.edge==0, "edge.obs.z"] <- 0
  return(results)   
}    
sesedge <- ses.edge(matrix = matrix, edge = edge )

grid.iberia <- rgdal::readOGR(paste0("4. Iberian EDGE/files/grids/UTM", res,"x",res, ".shp"))
grid <- RIRG::joinAttributeTable(grid.iberia, sesedge, names(grid.iberia)[1], "grid")
{
  png(paste0("4. Iberian EDGE/results/2. edge/EDGE.ses_", res, "(all).png"), height = 150, width = 150, units = "mm",
      res=450)
  par(mar=c(.1, .1, .1, .1))
  plotcontinuous(grid, col="edge.obs.z",  pal="maxent" , border="transparent" )
  cent.pd.h <- rgeos::gCentroid(grid[grid@data$pval>0.95 & grid@data$ses.edge!=0,], byid = T)
  plot(cent.pd.h, add=T, pch=16 , col="#00000088", cex=log(as.numeric(res), base=20)-0.2)
  dev.off()
}#png EDGE ratio

#5. EDGE zone prioritization ####

#SEE FOR-CLUSTER PRIORITIZATION

grid.iberia <- sf::st_read(paste0("4. Iberian EDGE/files/grids/UTM10x10.shp"))
colnames(grid.iberia)[1] <-"MGRS"
ibpen.cont <- sf::st_read("4. Iberian EDGE/files/grids/PeninsulaIberica.shp")
PN <- sf::st_read("4. Iberian EDGE/files/ENP.shp")

nzones <- as.numeric(21)#in EDGE complete 2perc: 80% = 10 steps; 85%= 15 steps; 90%=21steps
reglist <- list()
for( i in 1:nzones){
  reglist[[i]] <-  sf::st_read(paste0("E:/UNI/4. DOCTORADO/5. EDGE/forcluster/4. prioritization/output/9b. Newdata2/output/regions/region_st",i,".dbf"))
}
selected <- c("30SVG70","30TUN57","30STF70","31TDG38","31SDE80","30TTK56","30SYJ22","30SWG19","29TNG65",
              "29SQA29","30TUK86","30SXG45","29TPH74","30TWK67","30TUK15","30TWN84","30SVH07","30SUF48",
              "31SFE02","30TXL18","29SNC11")

colvect <- c(RColorBrewer::brewer.pal(8, "Set1"), RColorBrewer::brewer.pal(nzones-7, "Dark2"))
colvect <- c(colvect, RColorBrewer::brewer.pal(12, "Set3"))
colvect<-colvect[-11]

png("4. Iberian EDGE/results/priorityareas.png", height = 210, width = 210, units="mm", res=600)
par(bg="aliceblue")
plot(grid.iberia[1], border="transparent", col="transparent", 
     main="\nEDGE zones of\nIberian angiosperm flora", reset=F)
plot(ibpen.cont[1], col="#fafafa", add=T)
plot(PN[1]    ,border="transparent", col="#b0b0b0ff", add=T)
plot(PNPort[1],border="transparent", col="#b0b0b0ff", add=T)

for( i in nzones:1){
  plot(reglist[[i]][1], add=T, col=paste0(colvect[i], "AA"), border=colvect[i])
  selgrid <- grid.iberia[grid.iberia$MGRS==selected[i],]
  centr <- sf::st_centroid(selgrid)
  plot(centr, pch=16, add=T, col="black", cex = 0.5)
}
for( i in nzones:1){
  cent <- sf::st_centroid(reglist[[i]])
  text(cent$geometry[[1]][1], cent$geometry[[1]][2]+30000, as.character(i),
       col="black", cex = 0.85, font=2)
}

plot(ibpen.cont[1], col="transparent", add=T)

dev.off()




