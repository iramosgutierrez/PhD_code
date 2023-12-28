
library(terra)
library(phytools)

grid <- vect("4. Iberian EDGE/grid_20x20_withvalues.shp")
plotregions <- function(y){
  par(mfrow=c(1,2))
  phyloregion::plot.phyloregion(y,)
  
  par(mar=c(2,0,0,2))
  cols <- y$region.df[!(duplicated(y$region.df$cluster)),]
  dendro <- as.phylo(hclust(y$region.dist))
  plot.phylo(dendro,show.tip.label = F, direction = "downwards", )
  tiplabels(text="",frame = "none", pch=16, cex=5, col= cols$COLOURS)
  tiplabels(frame = "none")
}
compareAICs <- function(table){
  mod.iso <- glm(sep~isolation, data=table , family = "binomial")
  # mod.alm <- glm(sep~alt.mn   , data=table , family = "binomial")
  mod.alh <- glm(sep~alt.heter, data=table , family = "binomial")
  mod.ph  <- glm(sep~ph.mn    , data=table , family = "binomial")
  # mod.con <- glm(sep~cont     , data=table , family = "binomial")
  mod.b01 <- glm(sep~bio1.mn  , data=table , family = "binomial")
  mod.b04 <- glm(sep~bio4.mn  , data=table , family = "binomial")
  mod.b12 <- glm(sep~bio12.mn , data=table , family = "binomial")
  mod.b15 <- glm(sep~bio15.mn , data=table , family = "binomial")
  
  tab <- AIC(mod.iso, mod.alh , mod.ph  ,#mod.alm , mod.con, 
             mod.b01 , mod.b04 , mod.b12 , mod.b15)
  print(tab[order(tab$AIC),])
}
compare.ld <- function(table){
  
  ret <- list()
  
  ind <- sample(2, nrow(table),
                replace = TRUE,
                prob = c(0.7, 0.3))
  training <- table[ind==1,-c(1,which(colnames(table)=="cluster"))]
  testing <-  table[ind==2,-c(1,which(colnames(table)=="cluster"))]
  
  linear <- MASS::lda(sep~., training, groping=sep)
  
  pred <- predict(linear, testing)$class
  tab1 <- table(Predicted = pred, Actual = testing$sep)
  
  coefs <- as.data.frame(coefficients(linear))
  coefs$var <- rownames(coefs)
  coefs <- coefs[order(abs(coefs$LD1), decreasing = T),]
  
  ret$coefs <- coefs[,2:1]

  ret$accur <- sum(diag(tab1))/sum(tab1)
  ret$mod   <- linear
  
  ret$training <- training
  ret$testing  <- testing
  
  
  ret$pred.train <- predict(linear, training)
  ret$pred.test  <- predict(linear, testing)
  
  
  return(ret)
  
}
ld.stats <- function(table, n=100){
  ldstats <- list()
  ldstats["results"]<-list()
  ldstats["AUC"]<-list()
  for(i in 1:n){
    comps <- compare.ld(table)
    ldstats[["results"]][[i]] <- comps
    names(ldstats[["results"]])[i]<-i
    
    pred <- predict(comps$mod, comps$training)
    suppressMessages(
      auc <-pROC::roc(response=as.factor(comps$training$sep), 
                      predictor=as.numeric(pred$class))
    )
   
    ldstats[["AUC"]][[i]] <- as.numeric(auc$auc)
    names(ldstats[["AUC"]])[i]<-i
    
  }

  
  mean.accur <- mean(unlist(lapply(ldstats[["results"]], "[", 2)))
  sd.accur <- sd(unlist(lapply(ldstats[["results"]], "[", 2)))
  mean.auc <- mean(unlist(lapply(ldstats[["AUC"]], "[", 1)))
  
  list.table <- as.data.frame(lapply(ldstats[["results"]], "[", 1)[[1]])
  list.table$iter <- 1
  for (i in 2:n){
    tab <- as.data.frame(lapply(ldstats[["results"]], "[", 1)[[i]])
    tab$iter <- i
    list.table <- rbind(list.table, tab)
  }
  
  summ.table <- data.frame("var"=unique(list.table$coefs.var), "mean"=NA, "sd"=NA)
  for(var in unique(list.table$coefs.var)){
    summ.table[summ.table$var==var, "mean"] <- mean(list.table$coefs.LD1[list.table$coefs.var==var])
    summ.table[summ.table$var==var, "sd"] <- sd(list.table$coefs.LD1[list.table$coefs.var==var])
  } 
  
  ldstats$mean.accur <- mean.accur
  ldstats$sd.accur <- sd.accur
  ldstats$mean.AUC <- mean.auc
  
  
  importances <- abs(summ.table$mean)
  import.cum<- cumsum(importances)
  import.sum <- sum(importances)
  thres <- which(import.cum>import.sum*0.7)[1]
  if(thres==0){thres<-1}
  
  
  ldstats$summary.table <- summ.table
  
  # message(paste0("Accuracy: ",round(ldstats$mean.accur, digits = 4), 
  #              " +- ", round(ldstats$sd.accur, digits = 4), "\n"))
  message(paste0("Mean AUC: ",round(ldstats$mean.AUC, digits = 4), "\n"))
  
  print(ldstats$summary.table)
  
  message(paste0("Important variables (",thres,"):\n",
                 paste0(summ.table$var[1:thres], collapse = ", "),
          "\n"))
  return(ldstats)
  
  
}
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
importantvars <- function(ld){
  
  summ.table <- ld$summary.table
  
  importances <- abs(summ.table$mean)
  import.cum<- cumsum(importances)
  import.sum <- sum(importances)
  thres <- which(import.cum>import.sum*0.7)[1]
  if(thres==0){thres<-1}
  
   message(paste0("Important variables (",thres,"):\n",
                 paste0(summ.table$var[1:thres], collapse = ", "),
                 "\n"))
}
normalize <- function(df){
  for(i in 1:ncol(df))
  x <- df[,i]
  x.norm <- x
  x.norm[complete.cases(x)] <- (x[complete.cases(x)] -min(x[complete.cases(x)] ))/
    (max(x[complete.cases(x)]) - min(x[complete.cases(x)] ))
  df[,i] <- x.norm
  
  return(df)
  
  }


resolutions <- c("20_e_taxo", "20_e_phylo","20_t_taxo", "20_t_phylo")

# 20_e_Taxo####
resolution <- resolutions[1]

res <- stringr::word(resolution, 1, sep="_")
group <-stringr::word(resolution, 2, sep="_")
dis <- stringr::word(resolution, 3, sep="_")


regions <- load(paste0("4. Iberian EDGE/results/regions_",resolution, ".Rda"))
plotregions(y)


table <- merge(values(grid), y$membership, by.x=colnames(values(grid))[1], by.y="grids")
table <- table[,-which(colnames(table)%in%c("alt.mn","cont"))]
table[,-c(1,ncol(table))] <- scale(table[,-c(1,ncol(table))] )#own function normalize()could also be used
#Baleares vs rest
table1 <- table
table1$sep <- NA
table1$sep[table$cluster==8] <- 0
table1$sep[table$cluster!=8] <- 1

compareAICs(table1) #Isolation 
ld1 <- ld.stats(table1)
copyldatable(ld1)



#2-5-7 vs 1-3-4-6
table2 <- table
table2 <- table[table$cluster!=8,]
table2$sep <- NA

table2$sep[table2$cluster%in%c(1,3,4,6)] <- 0
table2$sep[table2$cluster%in%c(7,2,5)  ] <- 1
compareAICs(table2) #Bio12, annual precipitation
ld2 <- ld.stats(table2[,-which(colnames(table2)=="isolation")])
copyldatable(ld2)

# 1-3 vs 4-6
table3 <- table
table3 <- table3[table3$cluster%in% c(1,3,4,6),]
table3$sep[table3$cluster%in%c(1,3)] <- 0
table3$sep[table3$cluster%in%c(4,6)  ] <- 1
compareAICs(table3) #pH
ld3 <- ld.stats(table3[,-which(colnames(table3)=="isolation")])
copyldatable(ld3)

# 7 vs 2-5
table4 <- table
table4 <- table4[table4$cluster%in% c(7,2,5),]
table4$sep[table4$cluster%in%c(7)] <- 0
table4$sep[table4$cluster%in%c(2,5)  ] <- 1
compareAICs(table4) #bio4, temperature seasonality, small differences
ld4 <- ld.stats(table4[,-which(colnames(table4)=="isolation")])
copyldatable(ld4)

# save.image(paste0("4. Iberian EDGE/lda_",resolution, ".Rda"))
rm(list=ls()[-(which(ls()%in%c("grid", "plotregions", "compareAICs", 
                               "resolutions","compare.ld","ld.stats",
                               "copyldatable")))])


# 20_e_Phylo####
resolution <- resolutions[2]

res <- stringr::word(resolution, 1, sep="_")
group <-stringr::word(resolution, 2, sep="_")
dis <- stringr::word(resolution, 3, sep="_")


regions <- load(paste0("4. Iberian EDGE/results/regions_",resolution, ".Rda"))
plotregions(y)


table <- merge(values(grid), y$membership, by.x=colnames(values(grid))[1], by.y="grids")
table <- table[,-which(colnames(table)%in%c("alt.mn","cont"))]
table[,-c(1,ncol(table))] <- scale(table[,-c(1,ncol(table))] )
#Gimnésicas vs rest
table1 <- table
table1 <- table1[!(table1$cluster%in%c(7,8)),]
table1$sep <- NA
table1$sep[table1$cluster==10] <- 0
table1$sep[table1$cluster!=10] <- 1

compareAICs(table1) #Isolation
ld1 <- ld.stats(table1)
copyldatable(ld1)

#1-4-2-5 vs 9-3-6
table2 <- table
table2 <- table2[!(table2$cluster%in%c(7,8, 10)),]
table2$sep <- NA
table2$sep[table2$cluster%in%c(1,2,4,5)] <- 0
table2$sep[table2$cluster%in%c(9,3,6)] <- 1

compareAICs(table2) #bio1/bio15;small differences
ld2 <- ld.stats(table2)
copyldatable(ld2)

#1-4 vs 2-5
table3 <- table
table3 <- table3[table3$cluster%in%c(1,4,2,5),]
table3$sep <- NA
table3$sep[table3$cluster%in%c(1,4)] <- 0
table3$sep[table3$cluster%in%c(2,5)] <- 1

compareAICs(table3) #small differences bio15/pH
ld3 <- ld.stats(table3)
copyldatable(ld3)

#1 vs 4
table4 <- table
table4 <- table4[table4$cluster%in%c(1,4),]
table4$sep <- NA
table4$sep[table4$cluster%in%c(1)] <- 0
table4$sep[table4$cluster%in%c(4)] <- 1

compareAICs(table4) 
ld4 <- ld.stats(table4[,-(which(colnames(table)=="isolation"))])
copyldatable(ld4)

#1 vs 4
table5 <- table
table5 <- table5[table5$cluster%in%c(2,5),]
table5$sep <- NA_integer_
table5$sep[table5$cluster%in%c(2)] <- 0
table5$sep[table5$cluster%in%c(5)] <- 1

compareAICs(table5) 
ld5 <- ld.stats(table5)
copyldatable(ld5)

# save.image(paste0("4. Iberian EDGE/lda_",resolution, ".Rda"))
rm(list=ls()[-(which(ls()%in%c("grid", "plotregions", "compareAICs", 
                               "resolutions","compare.ld","ld.stats",
                               "copyldatable")))])
# 20_t_taxo####
resolution <- resolutions[3]

res <- stringr::word(resolution, 1, sep="_")
group <-stringr::word(resolution, 2, sep="_")
dis <- stringr::word(resolution, 3, sep="_")


regions <- load(paste0("4. Iberian EDGE/results/regions_",resolution, ".Rda"))
plotregions(y)


table <- merge(values(grid), y$membership, by.x=colnames(values(grid))[1], by.y="grids")
table <- table[,-which(colnames(table)%in%c("alt.mn","cont"))]
table[,-c(1,ncol(table))] <- scale(table[,-c(1,ncol(table))] )

#4 vs 1-2-3
table1 <- table
table1$sep <- NA
table1$sep[table$cluster==4] <- 0
table1$sep[table$cluster!=4] <- 1

compareAICs(table1) #Bio1 (T mn) 
ld1 <- ld.stats(table1)
copyldatable(ld1)

#2 vs 1-3
table2 <- table
table2 <- table2[table2$cluster%in%c(1,2,3),]
table2$sep <- NA
table2$sep[table2$cluster==2] <- 0
table2$sep[table2$cluster!=2] <- 1

compareAICs(table2) #pH, near bio12 (anunual precip) 
ld2 <- ld.stats(table2[,-(which(colnames(table)=="isolation"))])
copyldatable(ld2)

#1 vs 3
table3 <- table
table3 <- table3[table3$cluster%in%c(1,3),]
table3$sep <- NA
table3$sep[table3$cluster==1] <- 0
table3$sep[table3$cluster==3] <- 1

compareAICs(table3) #Bio15:Precipitation seasonality
ld3 <- ld.stats(table3[,-(which(colnames(table)=="isolation"))])
copyldatable(ld3)

# save.image(paste0("4. Iberian EDGE/lda_",resolution, ".Rda"))
rm(list=ls()[-(which(ls()%in%c("grid", "plotregions", "compareAICs", 
                               "resolutions","compare.ld","ld.stats",
                               "copyldatable")))])
# 20_t_Phylo####
resolution <- resolutions[4]

res <- stringr::word(resolution, 1, sep="_")
group <-stringr::word(resolution, 2, sep="_")
dis <- stringr::word(resolution, 3, sep="_")


regions <- load(paste0("4. Iberian EDGE/results/regions_",resolution, ".Rda"))
plotregions(y)


table <- merge(values(grid), y$membership, by.x=colnames(values(grid))[1], by.y="grids")
table <- table[,-which(colnames(table)%in%c("alt.mn","cont"))]
table[,-c(1,ncol(table))] <- scale(table[,-c(1,ncol(table))] )


#3-8-11-6-5-9 vs 10-7-4-1-2
table1 <- table
table1$sep <- NA
table1$sep[table1$cluster%in%c(3,8,11,6,5,9)] <- 0
table1$sep[table1$cluster%in%c(10,7,4,1,2)] <- 1

compareAICs(table1) #Bio1 
ld1 <- ld.stats(table1)
copyldatable(ld1)

#3-8 vs 11-6-5-9
table2 <- table
table2 <- table2[table2$cluster%in%c(3,8,11,6,5,9),] 
table2$sep <- NA
table2$sep[table2$cluster%in%c(3,8)] <- 0
table2$sep[table2$cluster%in%c(11,6,5,9)] <- 1

compareAICs(table2) #pH
ld2 <- ld.stats(table2[,-(which(colnames(table)=="isolation"))])
copyldatable(ld2)


#1 vs 2
table3 <- table
table3 <- table3[table3$cluster%in%c(1,2),] 
table3$sep <- NA
table3$sep[table3$cluster==1] <- 0
table3$sep[table3$cluster==2] <- 1

compareAICs(table3) #bio12, near pH
ld3 <- ld.stats(table3[,-(which(colnames(table)=="isolation"))])
copyldatable(ld3)

# save.image(paste0("4. Iberian EDGE/lda_",resolution, ".Rda"))
rm(list=ls()[-(which(ls()%in%c("grid", "plotregions", "compareAICs", 
                               "resolutions","compare.ld","ld.stats",
                               "copyldatable")))])

