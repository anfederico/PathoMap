library(Rsubread)

pathoalign <- function(index,
                       output_file,
                       output_format,
                       readfile1="",
                       readfile2="",
                       ...) {
  # Single-end alignment
  if (readfile2 == "") {
    Rsubread::align(index=index, 
                    output_format=output_format,
                    readfile1=readfile1,
                    output_file=output_file,
                    ...)
    # Paired-end alignment
  } else {
    Rsubread::align(index=index, 
                    output_format=output_format,
                    readfile1=readfile1,
                    readfile1=readfile2,
                    output_file=output_file,
                    ...)
  }
}

merge.samfiles <- function(samfiles, outfile) {
  outfile <- file(outfile, 'w') 
  on.exit(close(outfile))
  for (fn in samfiles) {
    infile <- file(fn, 'r') 
    while ( TRUE ) {
      line <- readLines(infile, n=1)
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
      line <- readLines(infile, n=1)
      if ( length(line) == 0 ) {
        break
      }
      if (substr(line, 1, 1) != "@") {
        writeLines(line, outfile)
      }
    }
  }
}

filter.unmapped <- function(samfile, outfile) {
  infile <- file(samfile, 'r')
  outfile <- file(outfile, 'w') 
  on.exit(close(infile))
  on.exit(close(outfile))
  while ( TRUE ) {
    line <- readLines(infile, n=1)
    if ( length(line) == 0 ) {
      break
    }
    if (substr(line, 1, 1) == "@") {
      writeLines(line, outfile)
    } else {
      fields <- strsplit(line, "\t")[[1]]
      qname <- fields[1]
      rname <- fields[3]
      if (rname == "*") {
        writeLines(line, outfile)
      }
    }
  }
}

