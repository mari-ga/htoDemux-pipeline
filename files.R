
library(argparser, quietly=TRUE)


# Create a parser
p <- arg_parser("Parameters for HTODemux")
#Output paths
p <- add_argument(p, "--htoDemuxOut",help="Path to file where the results of htoDemux will be saved", default = NULL)
p <- add_argument(p, "--nameOutputFile",help="Name for the file containing the output of HTODemux hashtag, no need to add extension, txt per default", default = "result.txt")
p <- add_argument(p, "--graphs",help="Path to folder where the graphs produced from the HTODemux function can be saved", default = NULL)

argv <- parse_args(p)


create_files <- function(name, path) {
  path_complete <- paste(path, name,".txt",sep="")
  print(path_complete)
  if (file.exists(path_complete)) {
    cat("The file already exists... sorry baby x")
    return(-1)
  } else {
    cat("Created new file with results")
    file.create(path_complete)
    return(path_complete)
  }
}

a <-create_files(argv$nameOutputFile, argv$htoDemuxOut)

print(a)
