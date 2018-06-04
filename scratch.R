library(Rsubread)
library(Rsamtools)

setwd("/Users/anthonyfederico/Work/johnson/PathoScope")

merge.samfiles <- function(outfile, samfiles) {
  outfile <- file(outfile, 'w') 
  for (fn in samfiles) {
    infile <- file(fn, 'r') 
    while ( TRUE ) {
      line = readLines(infile, n=1)
      if ( length(line) == 0 ) {
        break
      }
      if (substr(line, 1, 1) == "@") {
        writeLines(line, outfile)
      } else {
        break
      }
    }
  }
  for (fn in samfiles) {
    infile <- file(fn, 'r') 
    while ( TRUE ) {
      line = readLines(infile, n=1)
      if ( length(line) == 0 ) {
        break
      }
      if (substr(line, 1, 1) != "@") {
        writeLines(line, outfile)
      }
    }
  }
  close(outfile)
}

pathoalign <- function(index,
                       output_file,
                       readfile1="",
                       readfile2="",
                       ...) {
    # Single-end alignment
    if (readfile2 == "") {
      Rsubread::align(index=index, 
                      output_format="SAM",
                      readfile1=readfile1,
                      output_file=output_file,
                      ...)
    # Paired-end alignment
    } else {
      Rsubread::align(index=index, 
                      output_format="SAM",
                      readfile1=readfile1,
                      readfile1=readfile2,
                      output_file=output_file,
                      ...)
    }
}

pathomap <- function(outDir = ".",
                     indexDir = ".",
                     numThreads = 8,
                     outAlignFile = "",
                     inReadFilePair1 = "",
                     inReadFilePair2 = "" , 
                     targetRefFiles = c(),
                     filterRefFiles = c(),
                     targetIndexPrefixes = c(),
                     filterIndexPrefixes = c(),                  
                     memory=8000) {

  # Create index directory if it doesn't exist  
  dir.create(indexDir, showWarnings=FALSE)
  
  # Build the indices
  #for (ref in c(targetRefFiles, filterRefFiles)) {
  #  Rsubread::buildindex(basename=file.path(indexDir, basename(ref)), 
  #                       reference=ref, 
  #                       indexSplit=TRUE, 
  #                       memory=memory)
  #}
  
  # Align reads to target indicies
  dir.create(file.path(outDir, "temp"), showWarnings=FALSE)
  
  sams.target = c()
  for (ref in targetRefFiles) {
    pathoalign(index=file.path(indexDir, basename(ref)), 
               output_file=file.path(outDir, paste("temp/alignment.",  basename(ref), ".SAM", sep="")),
               readfile1=inReadFilePair1, 
               readfile2=inReadFilePair2,
               phredOffset=64)
    
    # Save the path of each bam for merging
    sams.target = c(sams.target, file.path(outDir, paste("temp/alignment.",  basename(ref), ".SAM", sep="")))
  }
  
  merge.samfiles(outfile=file.path(outDir, "temp/targets.merged.SAM"), sams.target)
}

#### CALL PATHOMAP
pathomap(outDir = "/Users/anthonyfederico/Work/johnson/PathoScope",
         indexDir = "/Users/anthonyfederico/Work/johnson/PathoScope/pathoindices",
         inReadFilePair1 = reads,
         targetRefFiles = target.refs,
         filterRefFiles = filter.refs,
         memory = 1
)


merged.sam = "/Users/anthonyfederico/Work/johnson/PathoScope/temp/targets.merged.SAM"

get.samAlignScore = function(samline) {
  l = unlist(strsplit(samline, split="\t"))
  usemapq = TRUE
  readl = as.double(nchar(l[10]))
  for (i in 12:length(l)) {
    if (usemapq & substr(l[i], 1, 5) == "AS:i:") {
      score = as.numeric(substr(l[i], 6, nchar(l[i])))
      useMapq = FALSE
    } else if (substr(l[i], 1, 5) == "YS:i:") {
      score = score+as.numeric(substr(l[i], 6, nchar(l[i])))
      readl = 2*readl
      break
    }
  }
  if (usemapq) {
    score = NULL
  } else {
    score = score+readl
  }
  return(score)
}


infile <- file(merged.sam, 'r') 
while ( TRUE ) {
  line = readLines(infile, n=1)
  if ( length(line) == 0 ) {
    break
  }
  if (substr(line, 1, 1) != "@") {
    score = get.samAlignScore(line)
    print(score)
  }
}
close(infile)
#STOP







samfiles = c("/Users/anthonyfederico/Work/johnson/PathoScope/temp/alignment.target_ref_1.fa.SAM",
             "/Users/anthonyfederico/Work/johnson/PathoScope/temp/alignment.target_ref_2.fa.SAM",
             "/Users/anthonyfederico/Work/johnson/PathoScope/temp/alignment.target_ref_3.fa.SAM")

merge.samfiles(outfile="/Users/anthonyfederico/Work/johnson/PathoScope/temp/alignment.merged.SAM", samfiles)

refs <- "/Users/anthonyfederico/Work/johnson/PathoScope/data/refs"
refs.target <- file.path(refs, "target_refs")
refs.filter <- file.path(refs, "filter_refs")

target.refs <- c(file.path(refs.target, "target_ref_1.fa"),
                 file.path(refs.target, "target_ref_2.fa"),
                 file.path(refs.target, "target_ref_3.fa"))

filter.refs <- c(file.path(refs.filter, "filter_ref_1.fa"),
                 file.path(refs.filter, "filter_ref_2.fa"),
                 file.path(refs.filter, "filter_ref_3.fa"))

reads <- "/Users/anthonyfederico/Work/johnson/PathoScope/data/srdat/reads.txt.gz"



# Use sed -i -e 's/|/;/g' ref.fa 
# to make a temp copy that removes the | delimeter
