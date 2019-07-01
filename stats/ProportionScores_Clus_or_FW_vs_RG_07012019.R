#Note: To do the calculations comparing 2 fragmentation methods
#took about 10 minutes of computer time.

#Set working directory to location of the data files

#setwd("C:\\Work\\Consulting\\MDR\\Other MDR Issues\\CSC\\FragmentOntologies")
setwd("C:/Work/git/NIH/scaffoldanalytics/stats") # Deepak's folder
#####################
rm(list = ls()) #Erase anything in R's working memory

DDData <- read.table("DDframes.txt", header = T, sep = "\t")
RGData <- read.table("RGD.txt", header = T, sep = "\t")
CLData <- read.table("CLinkClusters.txt", header = T, sep = "\t") 

#Create a list of compounds shared between two datasets.
#Scores will not make sense if we include non-common compounds

DDCompounds <- unique(DDData$COMPOUND_ID)
RGCompounds <- unique(RGData$COMPOUND_ID)
CLCompounds <- unique(CLData$COMPOUND_ID)  

# Compounds <- intersect(DDCompounds, RGCompounds)
# if comparing 3 methods pairwise, make sure we use the same set of compounds 
Compounds <- intersect(intersect(RGCompounds, DDCompounds), CLCompounds)

# Compounds <- as.numeric(as.character(intersect(DDCompounds, RGCompounds)))
N <- length(Compounds)

#Calculate PI's and common proportions for each compound
  
#Proportions by compound
statcompare <- function(ScafSetA, ScafSetB, 
                        ScafColA="StrucUniqueID", ScafColB="SCAFFOLD_ID", CompoundSet) {
  N <- length(CompoundSet)
  PropByCompound <- data.frame(CompoundID = rep(NA, N),
							                 FragA = rep(NA, N), FragB = rep(NA, N),
							                 Ca = rep(NA, N), Cb = rep(NA, N),
                               IntAB = rep(NA, N), UnionAB = rep(NA, N),
                               CommonProp = rep(NA, N),
                               PIa = rep(NA, N), PIb = rep(NA, N),
                               PIaU = rep(NA, N), PIbU = rep(NA, N),
                               FragEffA = rep(NA, N), FragEffB = rep(NA, N))

    for (index in 1:N){
  
      PropByCompound$CompoundID[index] <- CompoundSet[index]
  
      MethodABelongsTo <- ScafSetA[ScafSetA$COMPOUND_ID == Compounds[index],
                                   ScafColA] # StrucUniqueID or CLink
    
      MethodBBelongsTo <- ScafSetB[ScafSetB$COMPOUND_ID == Compounds[index],
                                   ScafColB]   # SCAFFOLD_ID
    
      PropByCompound$FragA[index] <- length(unique(MethodABelongsTo))
      PropByCompound$FragB[index] <- length(unique(MethodBBelongsTo))
    
      MethodACompoundCluster <- unique(ScafSetA[ScafSetA[,ScafColA] %in% MethodABelongsTo,
                                                "COMPOUND_ID"])
  
      MethodBCompoundCluster <- unique(ScafSetB[ScafSetB[,ScafColB] %in% MethodBBelongsTo,
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
  } # for index
  return(PropByCompound) 
} # function

# call function - A is FW, B is RGT, D is CLink
PropByCompound_AB <- statcompare(ScafSetA = DDData , ScafSetB = RGData, 
                                 ScafColA="StrucUniqueID", ScafColB="SCAFFOLD_ID", 
                                 CompoundSet = Compounds)
PropByCompound_DB <- statcompare(ScafSetA = CLData , ScafSetB = RGData, 
                                 ScafColA="CLink", ScafColB="SCAFFOLD_ID", 
                                 CompoundSet = Compounds)

#Create output -- averages, quantiles, and histograms
summarize_prop <- function(PropByCompound) {
  ACP <- mean(PropByCompound$CommonProp, na.rm = T)
  APIa <- mean(PropByCompound$PIa, na.rm = T)
  APIb <- mean(PropByCompound$PIb, na.rm = T)
  APIaU <- mean(PropByCompound$PIaU, na.rm = T)
  APIbU <- mean(PropByCompound$PIbU, na.rm = T)
  AFragA <- mean(PropByCompound$FragA, na.rm = T)
  AFragB <- mean(PropByCompound$FragB, na.rm = T)
  AFragEffA <- mean(PropByCompound$FragEffA, na.rm = T)
  AFragEffB <- mean(PropByCompound$FragEffB, na.rm = T)
  ACa <- mean(PropByCompound$Ca, na.rm = T)
  ACb <- mean(PropByCompound$Cb, na.rm = T)

  CP90 <- quantile(PropByCompound$CommonProp, c(0.1,0.5,0.9))
  PIa90 <- quantile(PropByCompound$PIa, na.rm = T , c(0.1,0.5,0.9))
  PIb90 <-  quantile(PropByCompound$PIb, na.rm = T , c(0.1,0.5,0.9)) 
  PIaU90 <- quantile(PropByCompound$PIaU, na.rm = T , c(0.1,0.5,0.9))
  PIbU90 <- quantile(PropByCompound$PIbU, na.rm = T, c(0.1,0.5,0.9))
  FragA90 <- quantile(PropByCompound$FragA, na.rm = T, c(0.1,0.5,0.9))
  FragB90 <- quantile(PropByCompound$FragB, na.rm = T, c(0.1,0.5,0.9))
  FragEffA90 <- quantile(PropByCompound$FragEffA, na.rm = T, c(0.1,0.5,0.9))
  FragEffB90 <- quantile(PropByCompound$FragEffB, na.rm = T, c(0.1,0.5,0.9))
  Ca90 <- quantile(PropByCompound$Ca, na.rm = T, c(0.1,0.5,0.9))
  Cb90 <- quantile(PropByCompound$Cb, na.rm = T, c(0.1,0.5,0.9))

  ret <- list(ACP=ACP, CP90=CP90, 
              APIa=APIa, PIa90=PIa90, APIb=APIb, PIb90=PIb90, 
              APIaU=APIaU, PIaU90=PIaU90, APIbU=APIbU, PIbU90=PIbU90, 
              AFragA=AFragA, FragA90=FragA90, AFragB=AFragB, FragB90=FragB90, 
              AFragEffA=AFragEffA, FragEffA90=FragEffA90, AFragEffB=AFragEffB, FragEffB90=FragEffB90, 
              ACa=ACa, Ca90=Ca90, ACb=ACb, Cb90=Cb90)

return(ret)
}

######### create summaries #######
SumProp_AB <- summarize_prop(PropByCompound_AB)
SumProp_DB <- summarize_prop(PropByCompound_DB)

# write output to CSV
write.table(PropByCompound_AB,"PropByCompound_FW_RGD.txt",sep="\t",row.name=F,col.name=T)
write.table(PropByCompound_DB,"PropByCompound_Cluster_RGD.txt",sep="\t",row.name=F,col.name=T)
#########################################################################################

################# Plot #####################
createPlots <- function(PropByCompound) {
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
}

createPlots(PropByCompound_AB)
createPlots(PropByCompound_DB)