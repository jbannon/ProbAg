library(TRONCO)
library(gtools)


create_dictionary_from_annotations = function(annotations){
  #print("****** CREATNG EDGE DICTIONARY*****")
  dict = c()
  ii=0
  dlines = c("gene\t num")
  for(j in annotations[,"event"]){
    dict = c(dict,ii)
    dlines = c(dlines,paste0(j,"\t",ii))
    ii=ii+1
  }
  fileConn <- file("../swapspace/node_numbering.txt")
  writeLines(dlines, fileConn)
  close(fileConn)
}

write_edge_list = function(annotations, adj_mat,fit,cancer,group){
 if (!file.exists("../swapspace/node_numbering.txt")){
  create_dictionary_from_annotations(annotations)
 }
  gene_id_map = read.delim(file = "../swapspace/node_numbering.txt", header = TRUE, sep = "\t", dec = ".")
  n_nodes = dim(annotations)[1]
  edge_list_numeric = c()
  edge_list = c()
  edge_count = 0


  for(i in 1:dim(adj_mat)[1]){
    for(j in 1:dim(adj_mat)[2]){
      if(adj_mat[i,j]==1){
        edge_count = edge_count+1
     src = annotations[rownames(adj_mat)[i],"event"]
     src_num = gene_id_map[,"num"][which(gene_id_map["gene"]==src)]
     dst = annotations[colnames(adj_mat)[j],"event"]
     dst_num = gene_id_map[,"num"][which(gene_id_map["gene"]==dst)]

     edge_list = c(edge_list,paste0(src,"\t",dst))
     edge_list_numeric = c(edge_list_numeric,paste0(src_num,"\t",dst_num))
   }
   }
  }
  edge_list = c(fit,edge_count,n_nodes,edge_list)
  edge_list_numeric = c(fit,edge_count,n_nodes,edge_list_numeric)
  fileConn <- file(paste0("../swapspace/edge_list.txt"))
  writeLines(edge_list, fileConn)
  close(fileConn)
  fileConn <- file(paste0("../swapspace/edge_list_numeric.txt"))
  writeLines(edge_list_numeric, fileConn)
  close(fileConn)
}



fit_sbcn = function(dataset){
  # meta_info = scan("../swapspace/sbcn_info.txt", character(), quote = "")
  select_max = T
  model = import.genotypes(dataset)

  model.fit = tronco.capri(model,boot.seed = 12345, nboot = 5,silent=T,pvalue=0.5)
  gene_list = model.fit$annotations

  if(select_max){
    n=nrow(model.fit$model$capri_aic$adj.matrix$adj.matrix.fit)
    adj_mat = matrix(0,n,n)
    aic =model.fit$model$capri_aic$adj.matrix$adj.matrix.fit
    bic = model.fit$model$capri_bic$adj.matrix$adj.matrix.fit
    for(i in seq(n)){
      for(j in seq(n)){
        adj_mat[i,j] =max(aic[i,j],bic[i,j])
      }
    }

  }else{
  if (model.fit$model$capri_aic$logLik >= model.fit$model$capri_bic$logLik){
    adj_mat = model.fit$model$capri_aic$adj.matrix$adj.matrix.fit
    fit="aic"
  } else {
    adj_mat = model.fit$model$capri_bic$adj.matrix$adj.matrix.fit
    fit = "bic"
  }
}


  write_edge_list(gene_list, adj_mat,fit,meta_info[2],meta_info[3])

}




args = commandArgs(trailingOnly=TRUE)

main = function(args){
  set.seed(1245)
  df = read.csv(args[1],row.names=1)
  fit_sbcn(df)
}
main(args)
