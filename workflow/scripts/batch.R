#!usr/bin/env Rscript

# Function to read the first header from a FASTQ file
read_first_fastq_header = function(file_path) {
  con = file(file_path, "r")
  first_header = readLines(con, n = 1)
  close(con)
  return(first_header)
}

# Function to extract and parse header parts
parse_header = function(header) {
  parts = strsplit(header, ":")[[1]]
  return(paste(parts[1:4], collapse = ":"))
}

fastq_dir = "/data/" # Directory containing your FASTQ files
fastq_files = list.files(fastq_dir, pattern = "\\.fastq$", full.names = TRUE) # Get a list of FASTQ files
fastq_files = fastq_files[grep("R1", fastq_files)] # Keep only R1 files
sample_info = read.csv("/data/files_info.csv")

# Initialize a vector to store parsed headers
all_headers = data.frame(filename=NA, batchID=NA)

# Process each FASTQ file
for (file in fastq_files) {
  header = read_first_fastq_header(file)
  parsed_header = parse_header(header)
  all_headers = rbind(data.frame(filename=file, batchID=parsed_header), all_headers)
}

# Count unique headers
theCount = table(all_headers$batchID)
theCount = data.frame(theCount)
theCount = data.frame(sequenceID=theCount$Var1, batch=paste0("Batch", 1:nrow(theCount)))

# Merge all_headers with theCount
all_headers = merge(all_headers, theCount, by.x="batchID", by.y="sequenceID")

# Keep only the file names and not the path of the files
all_headers$filename = sub(".*/", "", all_headers$filename)

# Rename columns
colnames(all_headers) = c("batchID", "fastq1", "Batch_ID")
all_headers$batchID = NULL

# Remove Batch_ID from sample_info if it exists
sample_info$Batch_ID = NULL

# Join all_headers with sample_info to add Batch_ID
sample_info = merge(sample_info, all_headers, by.x="fastq1", by.y="fastq1", all.x=TRUE)

# Output to a file
write.csv(sample_info, file = "data/files_info.csv", row.names = FALSE, quote = FALSE)
