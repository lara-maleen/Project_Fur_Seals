library(readr)

tabletofile <- function(table.print,filename){
  tofile(paste(colnames(table.print),collapse="\t"),"\n",filename=filename)
  apply(table.print,1,FUN = function(x) tofile(paste(x,collapse="\t"),"\n",filename = filename))
  # print(...,file=filename,append=TRUE,sep="")
}

tofile <- function(...,filename){
  cat(...,file=filename,append=TRUE,sep="")
}

filetofile <- function(fromfile,tofile){
  hdr <- paste(rep("#",nchar(fromfile)+4),collapse="")
  tofile(hdr,"\n# ",fromfile," #\n",hdr,"\n\n",filename=tofile)
  tofile(read_file(fromfile),"\n\n",filename=tofile)
}

construct_log <- function(sessionname,files,sims){
  filename <- paste(sessionname,"overview.log",sep="")
  cat("",file=filename) # create empty log-file (without appending)
  
  tofile("Simulation started ",date()," by ",Sys.info()['user'],"\n========================================================\n\n",filename = filename)  
  tofile("Parameter values\n-----------------\n",filename=filename)
  tabletofile(sims,filename=filename)
  tofile("\n\n\n",filename=filename)
  tofile("Relevant files\n-----------------\n",filename=filename)
  for(i in files){
    filetofile(fromfile=i,tofile=filename)
  }
}

# construct_log("bla",files=c("matrix_model.R"))
# getwd()
