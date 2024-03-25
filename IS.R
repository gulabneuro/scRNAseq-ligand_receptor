IS=function(ligands, receptors, cell_type1, cell_type2, lr.pairs, complete_dataset,num_null_reps){
  #ligands = dataframe of average expression of ligand genes in all cell types
  #receptors = dataframe of average expression of receptor genes in all cell types
  #cell_type1 = ligand cell type of interest
  #cell_type2 = receptor cell type of interest
  #lr.pairs = list of ligand-receptor pairs
  #threshold = value of IS for heatmap cutoff
  #num_null_reps = number of shuffles to generate null distribution for significance
  
  library(dplyr)
  
  #multiply values of cell type 1 in ligands with cell type 2 in receptors
  LRs=ligands[,cell_type1]*receptors[,cell_type2]
  
  #attach ligand-receptor label to rows
  LRs=data.frame(LRs)
  LRs=cbind(lr.pairs$Pair_Name,LRs)
  
  #attach cell type1-cell type2 label to each column
  colnames(LRs)=c('Pair_Name',paste(cell_type1, cell_type2,sep='-'))
  
  #generate null distribution - sample # of LR genes from whole dataset num_null_reps # of times
  n=dim(lr.pairs)[1]
  null_dist=randomize6(complete_dataset,cell_type1,n)
  for (i in 1:(num_null_reps-1)){
    null_dist=rbind(null_dist,randomize6(complete_dataset,cell_type1,n))
  }
  
  #set NA values to 0 and calculate IS cut off value
  null_dist=data.frame(null_dist)
  LRs[is.na(LRs)]<-0
  null_dist[is.na(null_dist)]<-0
  null_dist=data.matrix(null_dist)
  cut_off_p=0.01/(dim(lr.pairs)[1])^0.5
  cut_off_val=as.numeric(quantile(null_dist,(1-cut_off_p)))
  print(cut_off_val)
  
  
  #display only significant IS values
  LRs=distinct(LRs)
  LRs=filter(LRs,LRs[,2]>cut_off_val)
  LRs=LRs[order(LRs[2], decreasing = T),]
  
  
  return(LRs)
}


randomize6=function(complete_dataset,cell_type,n){
  #complete_dataset=average expression of all genes by subcluster in complete, merged dataset
  #returns one column, length=total # of genes in dataset, of randomly computed interaction scores
  #of all cell subtypes with ligand cell type
  rand_ligand=sample(complete_dataset[,cell_type], replace=T)
  rand_receptor=sample(complete_dataset, replace=T)
  IS=as.numeric(rand_ligand[1:n])*as.numeric(rand_receptor[1:n,1])
  return(IS)
}


IS_heatmap=function(ligands, receptors, cell_type, lr.pairs, threshold){
  #ligands = dataframe of average expression of ligand genes in all cell types
  #receptors = dataframe of average expression of receptor genes in all cell types
  #cell_type = ligand cell type of interest
  #lr.pairs = list of ligand-receptor pairs
  #threshold = value of IS for heatmap cutoff
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
  LRs=ligands[,cell_type]*receptors
  colnames(LRs)=paste(cell_type, colnames(receptors),sep='-')
  LRs=cbind(lr.pairs,LRs)
  LRs=distinct(LRs)
  LRs[is.na(LRs)]<-0
  colnames(LRs)[1]='Pair_Name'
  rownames(LRs)=LRs$Pair_Name
  LRs=LRs[,-1]
  LRs=filter_all(LRs, any_vars(.>threshold))
  pheatmap(LRs, display_numbers=T, fontsize_row=5, fontsize_col=6, fontsize_number=4, color=colorRampPalette(c("white", "red1", "red4"))(50))
}


#Example: cortex EC-cortex mural cell interaction scores

#ligands
cortex_ligands=read.csv('cortex_ligands.csv', row.names = 1)
#receptors
cortex_receptors=read.csv('cortex_receptors.csv', row.names=1)
#LR database
LRs=read.csv('LRs.csv')
#complete dataset
complete_dataset_merge=read.csv('complete_dataset.csv',row.names=1)

#calculate IS for cortex cECs and pericytes
EC_ligand=IS(cortex_ligands, cortex_receptors, 'cECs.1', 'pericytes.1',LRs,complete_dataset_merge, 1000)

#plot IS for cortex cECs ligands with pericyte and astrocyte receptors
IS_heatmap(cortex_ligands, cortex_receptors[,c('pericytes.1','cortex.astrocytes')], 'cECs.1',LRs$Pair_Name, 40)
