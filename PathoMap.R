library(Rsubread)

setwd("/Users/anthonyfederico/Work/johnson/PathoScope")
source("/Users/anthonyfederico/Work/johnson/PathoScope/helpers.R")

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
  for (ref in c(targetRefFiles, filterRefFiles)) {
    Rsubread::buildindex(basename=file.path(indexDir, basename(ref)), 
                         reference=ref, 
                         indexSplit=TRUE, 
                         memory=memory)
  }
  
  # Align reads to filter indices
  dir.create(file.path(outDir, "temp"), showWarnings=FALSE)
  sams.filter = c()
  for (ref in filterRefFiles) {
    pathoalign(index=file.path(indexDir, basename(ref)), 
               output_file=file.path(outDir, paste("temp/alignment.",  basename(ref), ".SAM", sep="")),
               output_format="SAM",
               readfile1=inReadFilePair1, 
               readfile2=inReadFilePair2,
               phredOffset=64)
    sams.filter = c(sams.filter, file.path(outDir, paste("temp/alignment.",  basename(ref), ".SAM", sep="")))
  }
  for (ind in filterIndexPrefixes) {
    pathoalign(index=ind, 
               output_file=file.path(outDir, paste("temp/alignment.",  basename(ind), ".SAM", sep="")),
               output_format="SAM",
               readfile1=inReadFilePair1, 
               readfile2=inReadFilePair2,
               phredOffset=64)
    sams.filter = c(sams.filter, file.path(outDir, paste("temp/alignment.",  basename(ref), ".SAM", sep="")))
  }

  # Merge filter sam files
  merge.samfiles(samfiles=sams.filter,
                 outfile=file.path(outDir, "temp/filters.merged.SAM"))
  
  # Keep reads not mapping to filter indices
  filter.unmapped(samfile=file.path(outDir, "temp/filters.merged.SAM"),
                  outfile=file.path(outDir, "temp/filtered.SAM"))

  
  # Re-align umapped reads in merged sam to target indices
  sams.target = c()
  for (ref in targetRefFiles) {
    pathoalign(index=file.path(indexDir, basename(ref)), 
               output_file=file.path(outDir, paste("temp/alignment.",  basename(ref), ".SAM", sep="")),
               output_format="SAM",
               input_format="SAM",
               readfile1=file.path(outDir, "temp/filtered.SAM"), 
               phredOffset=64)
    sams.target = c(sams.target, file.path(outDir, paste("temp/alignment.",  basename(ref), ".SAM", sep="")))
  }
  for (ind in targetIndexPrefixes) {
    pathoalign(index=ind, 
               output_file=file.path(outDir, paste("temp/alignment.",  basename(ind), ".SAM", sep="")),
               output_format="SAM",
               input_format="SAM",
               readfile1=file.path(outDir, "temp/filtered.SAM"), 
               phredOffset=64)
    sams.target = c(sams.target, file.path(outDir, paste("temp/alignment.",  basename(ref), ".SAM", sep="")))
  }
 
  # Merge target sam files
  merge.samfiles(samfiles=sams.target,
                 outfile=file.path(outDir, "targets.merged.SAM"))   
}

refs <- "/Users/anthonyfederico/Work/johnson/PathoScope/data/refs"
inds <- "/Users/anthonyfederico/Work/johnson/PathoScope/pathoindices"
refs.target <- file.path(refs, "target_refs")
refs.filter <- file.path(refs, "filter_refs")

target.refs <- c(file.path(refs.target, "target_ref_1.fa"),
                 file.path(refs.target, "target_ref_2.fa"))

target.inds <- c(file.path(inds, "target_ref_3.fa"))

filter.refs <- c(file.path(refs.filter, "filter_ref_1.fa"))

filter.inds <- c(file.path(inds, "filter_ref_2.fa"),
                 file.path(inds, "filter_ref_3.fa"))

reads <- "/Users/anthonyfederico/Work/johnson/PathoScope/data/srdat/reads.txt.gz"

#### CALL PATHOMAP
pathomap(outDir = "/Users/anthonyfederico/Work/johnson/PathoScope",
         indexDir = "/Users/anthonyfederico/Work/johnson/PathoScope/pathoindices",
         inReadFilePair1 = reads,
         targetRefFiles = target.refs,
         filterRefFiles = filter.refs,
         targetIndexPrefixes = target.inds,
         filterIndexPrefixes = filter.inds,
         memory = 1
)