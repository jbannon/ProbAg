library(TRONCO)




args = commandArgs(trailingOnly=TRUE)
# arg1 is file
# arg2 is bool to plot
# arg3 is plot info

main = function(args){
  #print("*****CONSOLIDATIING****")
  data = read.csv(args[1],row.names=1)
  model = import.genotypes(data)
  consolidated_info = consolidate.data(model)
  indist = consolidated_info$indistinguishable
  ones = consolidated_info$ones
  zeros = consolidated_info$zeroes
  if(length(indist)>0){
    issues = indist[[1]]
    issues  = issues[,2]
    issues = unname(issues)
    keep  = issues[1]
    toss = issues[2:length(issues)]
    fileConn <- file("./swapspace/keep.txt")
    writeLines(keep, fileConn)
    close(fileConn)

    fileConn <- file("./swapspace/remove.txt")
    writeLines(toss, fileConn)
    close(fileConn)


  }

}

main(args)
