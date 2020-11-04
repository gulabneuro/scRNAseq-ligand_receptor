#Calculate IS and bootstrap null distribution for significance

#IS function for calculating interaction scores
#randomize5 function for bootstrapping significance

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
    
    #generate null distribution - # rows = # genes in dataset; # cols = num_null_reps
    n=dim(lr.pairs)[1]
    null_dist=randomize5(complete_dataset,n)
    for (i in 1:(num_null_reps-1)){
        null_dist=rbind(null_dist,randomize5(complete_dataset,n))
    }
    
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

# randomize5 - calculates random interaction scores of ligands and receptors between two random cell types

randomize5=function(complete_dataset,n){
  #complete_dataset=average expression of all genes by cluster in complete, merged dataset
  #returns one column, length=total # of genes in dataset, of randomly computed interaction scores
  rand_ligand=sample(complete_dataset, replace=T)
  rand_receptor=sample(complete_dataset, replace=T)
  return(as.numeric(rand_ligand[1:n,1])*as.numeric(rand_receptor[1:n,1]))
}

#Example: cortex EC-cortex mural cell interaction scores

#ligands
cortex_ligands=read.csv('cortex_ligands.csv')
#receptors
cortex_receptors=read.csv('cortex_receptors.csv')
#LR database
LRs=read.csv('LRs.csv')
#complete dataset
complete_dataset_merge=read.csv('merge2_averages.csv')

EC_ligand=IS(cortex_ligands, cortex_receptors, 'Endothelial.cells.1',
        'Mural.cells',LRs,complete_dataset_merge, 1000)

