#Note: To do the calculations comparing 2 fragmentation methods
#took about 10 minutes of computer time.

#Set working directory to location of the data files
#setwd("C:\\Work\\Consulting\\MDR\\Other MDR Issues\\CSC\\FragmentOntologies")
setwd("C:/Work/git/NIH/scaffoldanalytics/stats")
#####################
rm(list = ls()) #Erase anything in R's working memory

DDData <- read.table("DDframes.txt", header = T, sep = "\t")
RGData <- read.table("RGD.txt", header = T, sep = "\t")

#Create a list of compounds shared between two datasets.
#Scores will not make sense if we include non-common compounds

DDCompounds <- unique(DDData$COMPOUND_ID)
RGCompounds <- unique(RGData$COMPOUND_ID)

Compounds <- intersect(DDCompounds, RGCompounds)
# Compounds <- as.numeric(as.character(intersect(DDCompounds, RGCompounds)))
N <- length(Compounds)

#Calculate PI's and common proportions for each compound
  
  #Proportions by compound

PropByCompound <- data.frame(CompoundID = rep(NA, N),
							 FragA = rep(NA, N),
							 FragB = rep(NA, N),
							              Ca = rep(NA, N),
                             Cb = rep(NA, N),
                             IntAB = rep(NA, N),
                             UnionAB = rep(NA, N),
                             CommonProp = rep(NA, N),
                             PIa = rep(NA, N),
                             PIb = rep(NA, N),
                             PIaU = rep(NA, N),
                             PIbU = rep(NA, N),
                            FragEffA = rep(NA, N),
                            FragEffB = rep(NA, N))

for (index in 1:N){
  
  PropByCompound$CompoundID[index] <- Compounds[index]
  
  MethodABelongsTo <- DDData[DDData$COMPOUND_ID == Compounds[index],
                            "StrucUniqueID"]
  MethodBBelongsTo <- RGData[RGData$COMPOUND_ID == Compounds[index],
                            "SCAFFOLD_ID"]
							
  PropByCompound$FragA[index] <- length(unique(MethodABelongsTo))
  PropByCompound$FragB[index] <- length(unique(MethodBBelongsTo))
  
  MethodACompoundCluster <- unique(DDData[DDData$StrucUniqueID %in% MethodABelongsTo,
                                   "COMPOUND_ID"])
  MethodBCompoundCluster <- unique(RGData[RGData$SCAFFOLD_ID %in% MethodBBelongsTo,
                                   "COMPOUND_ID"])
  
  PropByCompound$Ca[index] <- length(MethodACompoundCluster)
  PropByCompound$Cb[index] <- length(MethodBCompoundCluster)
  PropByCompound$IntAB[index] <- length(intersect(MethodACompoundCluster,
                                                  MethodBCompoundCluster))
  PropByCompound$UnionAB[index] <- length(union(MethodACompoundCluster,
                                                    MethodBCompoundCluster))
  PropByCompound$CommonProp[index] <- PropByCompound$IntAB[index]/PropByCompound$UnionAB[index]
  PropByCompound$PIa[index] <- PropByCompound$Ca[index]/PropByCompound$UnionAB[index]
  PropByCompound$PIb[index] <- PropByCompound$Cb[index]/PropByCompound$UnionAB[index]
  PropByCompound$PIaU[index] <- 1 - PropByCompound$Cb[index]/PropByCompound$UnionAB[index]
  PropByCompound$PIbU[index] <- 1 - PropByCompound$Ca[index]/PropByCompound$UnionAB[index]
  PropByCompound$FragEffA[index] <- PropByCompound$Ca[index]/PropByCompound$FragA[index]
  PropByCompound$FragEffB[index] <- PropByCompound$Cb[index]/PropByCompound$FragB[index]
  
}

#Create output -- averages, quantiles, and histograms

ACP <- mean(PropByCompound$CommonProp, na.rm = T)
APIaU <- mean(PropByCompound$PIaU, na.rm = T)
APIbU <- mean(PropByCompound$PIbU, na.rm = T)
AFragEffA <- mean(PropByCompound$FragEffA, na.rm = T)
AFragEffB <- mean(PropByCompound$FragEffB, na.rm = T)

CP95 <- quantile(PropByCompound$CommonProp, c(0.05,0.25,0.5,0.75, 0.95))
PIaU95 <- quantile(PropByCompound$PIaU, na.rm = T, c(0.05,0.25,0.5,0.75, 0.95))
PIbU95 <- quantile(PropByCompound$PIbU, na.rm = T, c(0.05,0.25,0.5,0.75, 0.95))
FragEffA95 <- quantile(PropByCompound$FragEffA, na.rm = T, c(0.05,0.25,0.5,0.75, 0.95))
FragEffB95 <- quantile(PropByCompound$FragEffB, na.rm = T, c(0.05,0.25,0.5,0.75, 0.95))

ACP
APIaU
APIbU
AFragEffA
AFragEffB

CP95
PIaU95
PIbU95
AFragEffA95
AFragEffB95

write.table(PropByCompound,"PropByCompound.txt",sep="\t",row.name=F,col.name=T)
#########################################################################################

################# Plot #####################
attach(PropByCompound)

hist(FragA,main="FragA")
hist(FragB,main="FragB")
hist(CommonProp,main="CommonProp")

plot(CommonProp~UnionAB)
plot(CommonProp~IntAB)
plot(Ca~Cb)

plot(UnionAB, IntAB)
plot(FragA + FragB, CommonProp)

detach(PropByCompound)




