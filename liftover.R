# Performs liftover of coordinates between hg19 and hg38

liftover_function = function(coordinates, liftover){
  
  liftover_coordinates = coordinates

  try({
    chain = import.chain(paste0("/path/to/liftover/",liftover,".over.chain"))
  
    coordinates_split = unlist(strsplit(coordinates,split=":|-"))
    chrom = paste0("chr",coordinates_split[1])
    start = as.numeric(coordinates_split[2])
    end = as.numeric(coordinates_split[3])
    
    original_ranges = GRanges(seqnames = chrom, ranges = IRanges(start=start, end=end))
    
    liftover_ranges = liftOver(original_ranges, chain)
    
    liftover_chrom = gsub(as.character(seqnames(liftover_ranges)[[1]]), patter="chr", replacement="")
    liftover_start = start(liftover_ranges)[[1]]
    liftover_end = end(liftover_ranges)[[1]]
    
    liftover_coordinates = paste0(liftover_chrom,":",liftover_start,"-",liftover_end)
  })
  
  liftover_coordinates
}
