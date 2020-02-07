library(readr)

tofile <- function(...,filename){
  cat(...,file=filename,append=TRUE,sep="")
}

filetofile <- function(fromfile,tofile){
  hdr <- paste(rep("#",nchar(fromfile)+4),collapse="")
  tofile(hdr,"\n# ",fromfile," #\n",hdr,"\n\n",filename=tofile)
  tofile(read_file(fromfile),filename=tofile)
}

construct_log <- function(sessionname){
  filename <- paste(sessionname,".log",sep="")
  cat("",file=filename)
  
  tofile("Simulation started",date(),"\n===============================================\n\n",filename = filename)  
  tofile("Parameter values\n-----------------\n",filename=filename)
  tofile("Relevant files\n-----------------\n",filename=filename)

  filetofile(fromfile="matrix_model.R",tofile=filename)
}

construct_log("bla")
getwd()
