library(ape)
library(picante)

pd.10.list   <- list()
mpd.10.list  <- list()
mntd.10.list <- list()

i <- Sys.getenv("i")

afliber.orig <- read.csv("files/lists/afliber_matrix_10x10.csv")
colnames(afliber.orig) <- gsub("\\.", "-", colnames(afliber.orig))



tree <- read.tree(paste0("files/trees/tree_", i, ".tre"))

tips <- colnames(afliber.orig)[colnames(afliber.orig)%in%tree$tip.label]
afliber <- afliber.orig[,c(1,which(colnames(afliber.orig)%in%tips))]


rownames(afliber)<- afliber[,1]
afliber <- afliber[,-1]
afliber <- as.matrix(afliber)
cophtree <- cophenetic(tree)


pd <- ses.pd(afliber,tree , null.model = "taxa.labels", include.root = F, runs = 99)     
pd$id <- rownames(pd)


mpd <- ses.mpd(afliber,cophtree , null.model = "taxa.labels", runs = 99)  
mpd$id <- rownames(mpd)


mntd <- picante::ses.mntd(afliber,cophtree , null.model = "taxa.labels", runs = 99)
mntd$id <- rownames(mntd)

pd.10.list  [[1]]<- pd
mpd.10.list [[1]]<- mpd
mntd.10.list[[1]]<- mntd
print(i)



assign(paste0("pd.10.list_",  i), pd.10.list   )
assign(paste0("mpd.10.list_", i), mpd.10.list  )
assign(paste0("mntd.10.list_",i), mntd.10.list )

save(list = c(paste0("pd.10.list_",  i),
              paste0("mpd.10.list_", i),
              paste0("mntd.10.list_",i)),
     
     file = paste0("output/10/phylo_10_", i,".rda"))


