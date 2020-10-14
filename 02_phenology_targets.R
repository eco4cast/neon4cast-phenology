source("downloadPhenoCam.R")
##Selected Sites for Challenge
siteIDs <- c("NEON.D01.HARV.DP1.00033","NEON.D01.BART.DP1.00033","NEON.D02.SCBI.DP1.00033",
             "NEON.D05.STEI.DP1.00033","NEON.D06.UKFS.DP1.00033","NEON.D07.GRSM.DP1.00033",
             "NEON.D08.DELA.DP1.00033","NEON.D11.CLBJ.DP1.00033")

allData <- data.frame(matrix(nrow=0,ncol=5))
colnames(allData) <- c("ID","Target","Time","Value","Uncertainty")

for(i in 1:length(siteIDs)){
  siteName <- siteIDs[i]
  URL <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_1000_1day.csv",sep="")
  phenoData <- download.phenocam(URL=URL)
  subPhenoData <- cbind(rep(siteName,nrow(phenoData)),rep("gcc_mean",nrow(phenoData)),
                        subset(phenoData,select=c(date,gcc_mean,gcc_std)))
  colnames(subPhenoData) <- c("ID","Target","Time","Value","Uncertainty")
  allData <- rbind(allData,subPhenoData)
}

write.table(file="phenology-targets.csv",allData,col.names = TRUE,sep=",",row.names = FALSE)

