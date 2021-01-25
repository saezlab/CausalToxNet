# Driver CausalToxNet Rat In Vitro

# Clean workspace
rm(list=ls());cat("\014");if(length(dev.list())>0){dev.off()}

# Set working directory
setwd("CausalToxNet/CARNIVAL_Scripts/RLS")

# TG-GATEs datasets to include
Compounds <- c("acetaminophen","colchicine","captopril","methapyrilene","nitrofurantoin")
# Compounds <- c("APAP","COL","CAP","MP","NFT")
# TimePoints <- c(2,8,24)
# TimePoints <- c(0.125,0.25,0.375,1)
TimePoints <- c(3,6,9,24)
Doses <- c('Low','Middle','High')

# Load library and settings
# library(CARNIVAL)
DRTconf <- "ABC"
TopDRT <- 50
options("scipen" = 5) # to print the full number (not as exponent - needed for writing IP constraints)

# Load resources
# CurrentPath <- getwd()
# dir.create("Resources")
# setwd(paste0(CurrentPath,"/Resources"))
# file.copy(from=system.file("Ex2_network_SBV_Omnipath.sif",package="CARNIVAL"),to=getwd(),overwrite=TRUE) # retrieve network file
# file.copy(from=system.file("HUMAN_9606_idmapping_onlyGeneName.dat",package="CARNIVAL"),to=getwd(),overwrite=TRUE) # retrieve network file
# setwd(CurrentPath)
# netFile <- "Resources/Liver_Specific_Signed_Omnipath_Network_CARNIVAL_Rat.tsv" # required

# Now use human network with rat-liver expression after humanising rat data
netFile <- "Resources/Liver_Specific_Signed_Omnipath_Network_CARNIVAL_HumanisedRat.tsv" # required

library(doParallel)
argsJob=commandArgs(trailingOnly = TRUE)
counter_cp <- as.numeric(argsJob[1])
counter_tp <- as.numeric(argsJob[2]) 
counter_do <- as.numeric(argsJob[3]) 

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! # 
counter_cp=1;counter_tp=1;counter_do=1 # For local test
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! # 

# if (file.exists(paste0("Resources/inputs_CARNIVAL_DRT_",DRTconf,"_FULL_Cutoff1p5/",Compounds[counter_cp],"_",Doses[counter_do],"_",TimePoints[counter_tp],"a_logFC_measCont_Cutoff1p5_Rat.txt"))) {
if (file.exists(paste0("Resources/inputs_CARNIVAL_DRT_",DRTconf,"_FULL_Cutoff1p5_Rat_Liver/",Compounds[counter_cp],"_",TimePoints[counter_tp],"hr_",Doses[counter_do],"__measCont_Cutoff1p5_Rat_Liver.txt"))) {
    
  print("")
  print("===========================================")
  print(paste0("Optimising: ",Compounds[counter_cp]," [",counter_cp,"/",length(Compounds),"] - ",TimePoints[counter_tp]," Hour [",counter_tp,"/",length(TimePoints),"] - ",Doses[counter_do]," [",counter_do,"/",length(Doses),"]"))
  print("===========================================")
  print("")
  
  # Specify the path to input files
  measFile <- paste0("Resources/inputs_CARNIVAL_DRT_",DRTconf,"_FULL_Cutoff1p5_Rat_Liver/",Compounds[counter_cp],"_",TimePoints[counter_tp],"hr_",Doses[counter_do],"__measCont_Cutoff1p5_Rat_Liver.txt") # required
  weightFile <- paste0("Resources/inputs_CARNIVAL_PGN_FULL_Rat_Liver/",Compounds[counter_cp],"_",TimePoints[counter_tp],"hr_",Doses[counter_do],"_scores_Rat_Liver.txt") # optional; if not, set to 'NULL'
  # inputFile <- paste0("../../GEX2NET_Paper/Rscript/Clean_CARNIVAL_Rpackage/Resources/Drug_Targets/",Compounds[counter_cp],"_drug_targets.tsv") # optional; if not, set to 'NULL'
  Result_dir <- paste0("CARNIVAL_Results_RatInVivo_",Compounds[counter_cp],"_",TimePoints[counter_tp],"hr_",Doses[counter_do])
  
  
  # Select Top-TF and convert GeneSymbol to UniprotID 
  IDmap <- read.table("Resources/HUMAN_9606_idmapping_onlyGeneName.dat",sep = "\t",stringsAsFactors = F)
  Meas <- read.table(measFile,sep = "\t",stringsAsFactors = F)
  Meas <- Meas[,order(abs(as.numeric(Meas[2,])),decreasing = T)]
  if (ncol(Meas)>TopDRT) {Meas <- Meas[,1:TopDRT]}
  
  MeasMapped <- NULL
  MeasUnmapped <- NULL
  for (counter in 1:ncol(Meas)) {
    if (length(which(Meas[1,counter]==IDmap[,3]))>0) {
      for (counter2 in 1:length(which(Meas[1,counter]==IDmap[,3]))) {
        MeasMapped <- cbind(MeasMapped,matrix(data = c(IDmap[which(Meas[1,counter]==IDmap[,3])[counter2],1],Meas[2,counter]),nrow = 2,ncol = 1))
      }
    } else {
      MeasUnmapped <- c(MeasUnmapped,Meas[1,counter])
    }
  }
  
  ifelse (test = length(MeasUnmapped)>0,yes = print(paste0("The following Gene Symbols couldn't be mapped:",MeasUnmapped)),no = print("All Gene Symbols were mapped succesfully"))
  
  NewMeasFileName <- paste0("Resources/DRT_",DRTconf,"_",Compounds[counter_cp],"_",TimePoints[counter_tp],"Day_",Doses[counter_do],"_Cutoff1p5_Rat_Liver_single.tsv")
  write.table(x = MeasMapped,file = NewMeasFileName,quote = F,sep = "\t",col.names = F,row.names = F)
  
  # Run the CARNIVAL pipeline
  # library(CARNIVAL)
  source("runCARNIVAL_CTN.R")
  source("writeLPFile_CTN.R")
  source("CARNIVAL_R_Rat/Source_CARNIVAL_Rat.R")
  # runCARNIVAL_CTN(CplexPath="/net/data.isilon/ag-saez/bq_ptrai/CR_Model/CARNIVAL/R/cplex",
  runCARNIVAL_CTN(CplexPath="~/Applications/IBM/ILOG/CPLEX_Studio1271/cplex/bin/x86-64_osx/cplex",
                  netFile=netFile,
                  measFile=NewMeasFileName,
                  inputFile=NULL,
                  weightFile=weightFile,
                  Result_dir=paste0(Result_dir,"_FI"),
                  CARNIVAL_example=NULL,
                  timelimit = 3600,
                  # timelimit = 30,
                  inverseCR = T,
                  repIndex1 = counter_cp,
                  repIndex2 = counter_tp,
                  repIndex3 = counter_do
  )
  
} else {
  print(paste0("Measurement file doesn't exist: ",Compounds[counter_cp]," [",counter_cp,"/",length(Compounds),"] - ",TimePoints[counter_tp],"Day [",counter_tp,"/",length(TimePoints),"] - ",Doses[counter_do]," [",counter_do,"/",length(Doses),"]"))
}

