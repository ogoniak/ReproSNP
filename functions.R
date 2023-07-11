# OTHER FUNCTIONS

link_database = function(ID, url) {
  
  ID = unlist(str_split(ID, pattern=":"))
  CHROM = ID[1]; if(CHROM=="23"){CHROM="X"} else if (CHROM=="24"){CHROM="Y"}
  BP = ID[2]
  
  link_name = paste0(CHROM,":",BP)
  
  url = gsub(url, pattern="___CHROM___", replacement=CHROM, fixed=T)
  url = gsub(url, pattern="___BP___", replacement=BP, fixed=T)
  
  paste0('<a href="', URLdecode(url),'" target="_blank">', link_name, '</a>')
}


filter_cohort = function(phenotable, dir, cohort_name, cohort_Sex, cohort_Generation, cohort_Number_of_children, cohort_Reported_infertility, cohort_include_NA){
  
  cohort = phenotable[
    phenotable$Sex %in% cohort_Sex &
      phenotable$Generation %in% cohort_Generation &
      phenotable$Reported_infertility %in% cohort_Reported_infertility &
      ((phenotable$Number_of_children >= cohort_Number_of_children[1] &
          phenotable$Number_of_children <= cohort_Number_of_children[2]) |
         phenotable$Number_of_children == cohort_include_NA)
    ,]
  
  IDs = data.frame("FID" = cohort[,1],"IID" = cohort[,1])
  write.table(IDs, file=paste0(dir,"/temporary_IDs_cohort",cohort_name,".txt"), col.names=F, row.names=F, quote=F)
}

add_AF = function(table){
  table$C1_ALT_AF = round(table$C1_ALT_CTS / table$C1_OBS_CT, 4)
  table$C2_ALT_AF = round(table$C2_ALT_CTS / table$C2_OBS_CT, 4)
  table
}

split_ID = function(ID) {
  ID = unlist(str_split(ID, pattern=":"))
  BP = as.numeric(ID[2])
}
