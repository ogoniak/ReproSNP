annotate_dbsnp = function(row, output, dbsnp){
  
  vcf_id = as.character(row["ID"])
  chr_name = strsplit(vcf_id, split=":")[[1]][1]
  start = as.numeric(strsplit(vcf_id, split=":")[[1]][2])
  ref = as.character(row["REF"])
  alt = as.character(row["ALT"])
  end = start + nchar(ref) - 1
  
  annotation = getBM(attributes = c("refsnp_id","clinical_significance","allele","minor_allele_freq","minor_allele"),
                     filters = c("chr_name","start","end"),
                     values = list(chr_name,start,end),
                     mart = dbsnp)
  
  n_entries = length(annotation$refsnp_id)
  
  matching_id = NA
  
  for (i in 1:n_entries){
    
    minor_allele = gsub(as.character(annotation$minor_allele[i]), pattern="TRUE", replacement="T")
    
    if (is.na(minor_allele) & n_entries == 1){
      matching_id = annotation$refsnp_id[i]
    } else if (!is.na(minor_allele) & (minor_allele == alt | minor_allele == ref)){
      matching_id = annotation$refsnp_id[i]
    }
  }

  matching_annotation = annotation[annotation$refsnp_id == matching_id,]
  refsnp_ids = matching_annotation$refsnp_id
  
  urls = c()
  for (refsnp_id in refsnp_ids){
    url = paste0("<a href=https://www.ncbi.nlm.nih.gov/snp/?term=", refsnp_id, " target=_blank>", refsnp_id,"</a>")
    urls = c(urls, url)
  }

  MAF = round(as.numeric(matching_annotation$minor_allele_freq), 4)
  MAF = format(MAF, scientific = FALSE)
  
  return(c(rsID = paste(refsnp_ids, collapse="; "),
           dbSNP_rsID = paste(urls, collapse="; "),
           ClinicalSignificance = paste(matching_annotation$clinical_significance, collapse="; "),
           dbSNP_MAF = MAF,
           MA = gsub(paste(matching_annotation$minor_allele, collapse="; "), pattern="TRUE", replacement="T")
  ))
}



