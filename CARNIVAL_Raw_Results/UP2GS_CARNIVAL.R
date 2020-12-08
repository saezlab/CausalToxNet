# Gene Symbol Mapping of CARNIVAL resutls

# List all directories
AllDIR <- list.dirs()
AllDIR <- AllDIR[-1] # remove the current directory from the list

IDmap <- read.table(file = system.file("HUMAN_9606_idmapping_onlyGeneName.dat",package="CARNIVAL"),header = F,sep = "\t",stringsAsFactors = F)

for (counter_DIR in 1:length(AllDIR)) {
  print(paste0("Mapping Uniprot to Gene Symbol: ",counter_DIR,"/",length(AllDIR)))
  Unmapped <- NULL
  
  # Common SIF
  
  if (file.exists(paste0(AllDIR[counter_DIR],"/weightedModel_1.txt"))) {
    
    SIF <- read.table(file = paste0(AllDIR[counter_DIR],"/weightedModel_1.txt"),header = T,sep = "\t",stringsAsFactors = F)

    for (counter in 1:nrow(SIF)) {
      if (length(IDmap[which(IDmap[,1] == SIF[counter,1]),3])>0) {
        SIF[counter,1] <- IDmap[which(IDmap[,1] == SIF[counter,1]),3][1]
      } else {
        # print(paste0("The Uniprot ID: ",res[[1]][counter,1]," couldn't be mapped"))
        Unmapped <- c(Unmapped,SIF[counter,1])
      }
      if (length(IDmap[which(IDmap[,1] == SIF[counter,3]),3])>0) {
        SIF[counter,3] <- IDmap[which(IDmap[,1] == SIF[counter,3]),3][1]
      } else {
        # print(paste0("The Uniprot ID: ",res[[1]][counter,3]," couldn't be mapped"))
        Unmapped <- c(Unmapped,SIF[counter,3])
      }
    }
    
    InputNodes <- setdiff(sort(unique(SIF[,1])),sort(unique(SIF[,3])))
    
    write.table(x = SIF,file =  paste0(AllDIR[counter_DIR],"/weightedModel_1_GeneSymbol.txt"),quote = F,col.names = T,row.names = F,sep = "\t")
  }
  
  if (file.exists(paste0(AllDIR[counter_DIR],"/nodesAttributes_1.txt"))) {
    
    NodeAttr <- read.table(file = paste0(AllDIR[counter_DIR],"/nodesAttributes_1.txt"),header = T,sep = "\t",stringsAsFactors = F)
    
    # Common node activity
    for (counter in 1:nrow(NodeAttr)) {
      if (length(IDmap[which(IDmap[,1] == NodeAttr[counter,1]),3])>0) {
        NodeAttr[counter,1] <- IDmap[which(IDmap[,1] == NodeAttr[counter,1]),3][1]
      } else {
        # print(paste0("The Uniprot ID: ",res[[2]][counter,1]," couldn't be mapped"))
        Unmapped <- c(Unmapped,NodeAttr[counter,1])
      }
    }
    
    # Map input nodes
    for (counter in 1:length(InputNodes)) {
      if (NodeAttr[which(InputNodes[counter]==NodeAttr[,1]),6]!="D") { # if not "D"rug
        NodeAttr[which(InputNodes[counter]==NodeAttr[,1]),6] <- "I" # annotate as "I"nput
      }
    }
    
    write.table(x = NodeAttr,file =  paste0(AllDIR[counter_DIR],"/nodesAttributes_1_GeneSymbol.txt"),quote = F,col.names = T,row.names = F,sep = "\t")
    
  }
  
  print(paste0("The following node couldn't be mapped to gene symbol: ",paste(unique(Unmapped),collapse = " ")))
  
}
