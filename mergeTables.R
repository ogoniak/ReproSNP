# Function to merge acount and gcount output files

mergeTables = function(tmp_dir,n_cohorts){

  universal_header=c("CHROM", "ID", "REF", "ALT")
  AF_header=c("ALT_CTS", "OBS_CT")
  GT_header=c("HOM_REF_CT", "HET_REF_ALT_CTS", "TWO_ALT_GENO_CTS", "HAP_REF_CT", "HAP_ALT_CTS", "MISSING_CT")

  types=c("AF","GT")
  
  i = 0
  
  for (type in types){
    
    for (cohort in 1:n_cohorts){
      
      i = i+1
      
      if (type == "AF"){
        type_header = paste0("C",cohort,"_",AF_header)
      } else if (type == "GT"){
        type_header = paste0("C",cohort,"_",GT_header)
      }
      
      header=c(universal_header,type_header)

      table = read.table(paste0(tmp_dir,"temporary_",type,"_cohort",cohort,".",tolower(substr(type,start=1,stop=1)),"count"), col.names=header, sep="\t")
      
      if (i == 1){
        merged_table = table
      } else {
        
        if (all(merged_table$ID == table$ID)){
          merged_table = cbind(merged_table, table[,type_header])
        } else {
          message("ERROR: IDs are not identical.")
        }
        
      }
    }
  }
  
  merged_table
  
}