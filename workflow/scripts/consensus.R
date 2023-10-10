suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("tidyr"))
suppressMessages(library("magrittr"))
suppressMessages(library("jsonlite"))

option_list = list(
  make_option(c("-p", "--option"), type="character", default=NULL, help="", metavar="")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}

params = opt$option # This is the parameter that specifies whether you want largest coverage, or lowest error possible

jsonFiles = list.files("data/cutadapt/",recursive=TRUE,pattern="trimParameters.json")
# jsonFiles = list.files("~/Desktop/v13_V45/cutadapt/",recursive=TRUE,pattern="trimParameters.json",full.names=TRUE)
coverageFiles = list.files("data",recursive=TRUE,pattern=params,all.files=TRUE,full.names=TRUE)
# coverageFiles = list.files("~/Desktop/v13_V45/cutadapt/",recursive=TRUE,pattern=params,all.files=TRUE,full.names=TRUE)

myJSONTable = jsonlite::fromJSON(jsonFiles[1])


if (params == "highest_coverage"){
  
}
