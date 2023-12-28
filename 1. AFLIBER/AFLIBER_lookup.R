#INICIO####

library(rgbif)
library(rgdal)
library(maptools)
library(sf)
library(plyr)
library(raster)

PIB <- readOGR("1. AFLIBER/files/PROVINCIAS_PIB+BAL.shp")
proj4string(PIB) <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
PIBtr <- spTransform(PIB, "+proj=longlat +ellps=GRS80 +no_defs")

malla <-  readOGR("1. AFLIBER/files/Iberian_Peninsula_10x10_Grid.shp")
mallatr <- spTransform(malla, "+proj=longlat +ellps=GRS80 +no_defs")
proj4string(mallatr)

#OBTENCI?N DE INFO####
especiesTOTAL<- readxl::read_xlsx("1. AFLIBER/files/LISTADO_FINAL.xlsx")
especies <- subset(especiesTOTAL, especiesTOTAL$Endemic!="Alóctona")
especies$FINAL_NAME <- as.character(especies$FINAL_NAME)
num.esp<- nrow(especies)
especies$Esp.Subesp <- sapply(strsplit(especies$FINAL_NAME, " "), length)



#DESCARGA DE DATOS DE GBIF####
for (i in 1:nrow(especies)) {
  esp <- toString(especies$FINAL_NAME[i])
  if(especies$Esp.Subesp[i]==2){
    keys <- name_suggest(q = esp , rank = "species")
  }else{
    keys <- name_suggest(q = esp , rank = "subspecies")}
  keys <- as.data.frame(keys[1])
  l <- length(keys$data.key)
  if(l==1){
    busq1 <- occ_search(taxonKey = keys$data.key[1] , limit= 250000, decimalLatitude="35,45", decimalLongitude="-10,5")
    resultados <- busq1$data
  }else{
    for (j in 1:l){
      busq.j <- occ_search(taxonKey = keys$data.key[j], limit= 250000, decimalLatitude="35,45", decimalLongitude="-10,5")
      unionESPT <- busq.j$data
      if (j==1){
        resultados <- unionESPT
      }else{
        df1 <- as.data.frame(resultados)
        df2 <- as.data.frame(unionESPT)
        resultados <- rbind.fill (df1, df2)
      }  
    }
  }
  #guardar datos de gbif como csv
  if( (resultados[1]=="no data found, try a different search")&&
      (resultados[2]=="no data found, try a different search"))
  {
    if(exists("especieserror"))
    {
      a<- data.frame(paste0("ERROR_",esp))
      especieserror <- rbind.fill (especieserror, a)}
    else
    {
      especieserror<- data.frame(paste0("ERROR_",esp))
      write.csv(especieserror, file="1. AFLIBER/results/ESPECIESERROR.csv" )
    }}
  else{
    write.csv(resultados, file=paste0("1. AFLIBER/results/GBIF_",esp,".csv"))}
  
  
  #malla y mapa
  especie.i <- resultados
  especie.dftot <-as.data.frame(especie.i)
  especie.df <- subset(especie.dftot, especie.dftot$decimalLongitude!="NA"|especie.dftot$decimalLatitude!="NA" )
  coordinates(especie.df)<-  c("decimalLongitude", "decimalLatitude")
  especie.locs <- SpatialPoints(especie.df)
  proj4string(especie.locs) <- "+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"
  spTransform(especie.locs, "+proj=longlat +ellps=GRS80 +no_defs")
  f <- function(x) length(gregexpr("[[:digit:]]", as.character(x))[[1]]) 
  num.occ <- nrow(especie.df[1])
  for (i in 1:num.occ){
    especie.df$precx[i]<- f(especie.df$decimalLongitude[i])
    especie.df$precy[i]<- f(especie.df$decimalLatitude [i])} 
  especie.precisos <- subset(especie.df, "precx">3|"precy">4 )
  proj4string(especie.precisos) <- "+proj=longlat +ellps=GRS80 +no_defs"
  
  #seleccionar las cuadr?culas en las que haya un punto
  
  over <- sp::over(especie.precisos, mallatr)
  presencias <- over$MGRS_10km[!is.na(over$MGRS_10km)]
  presencias.sindupl <- presencias[!duplicated(presencias)]
  #guardarlo como un csv individual
  for (i in 1:nrow(over)){
    if(exists("cuadrs"))
    {
      c<- data.frame("ESP"=esp, "MGRS_10km"=over$MGRS_10km[i])
      cuadrs <- rbind.fill (cuadrs, c)}
    else
    {
      cuadrs<- data.frame("ESP"=esp, "MGRS_10km"=over$MGRS_10km[i])
      cuad.PIB <- cuadrs[!is.na(cuadrs$MGRS_10km),]
    }}
  write.csv(cuad.PIB, file=paste0("1. AFLIBER/results/Cuadriculas_",esp,".csv"))
  write.csv(cuad.PIB[(!duplicated(cuad.PIB$MGRS_10km)),], file=paste0("1. AFLIBER/results/Cuadriculas SIN DUPLICADOS_",esp,".csv"))
  
 
  #limpiar
  rm(busq1,busq2,busq.j, c, especie.df, especie.dftot, especie.i, especie.locs, especie.precisos, 
     over, tabla.especie, resultados, especieserror, esp, keys, i,j, l, df1, df2, especiesTOTAL, 
     cuadrs, cuad.PIB)
  
}




