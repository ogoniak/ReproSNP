# Fisher's exact test for allele frequencies

fisher = function(dataframe){
  
  dataframe$P_FISHER_VALUE = ""
  dataframe$P_FISHER = ""
  dataframe$P_ADJUSTED = ""
  
  for (row in 1:dim(dataframe)[1]){
    alt_C2 = dataframe[row,names(dataframe)=="C1_ALT_CTS"]
    total_C2 = dataframe[row,names(dataframe)=="C1_OBS_CT"]
    ref_C2 = total_C2 - alt_C2
    
    alt_C1= dataframe[row,names(dataframe)=="C1_ALT_CTS"]
    total_C1 = dataframe[row,names(dataframe)=="C1_OBS_CT"]
    ref_C1 = total_C1 - alt_C1
    
    dat = data.frame(
      "C2" = c(ref_C2,alt_C2),
      "C1" = c(ref_C1,alt_C1),
      row.names = c("REF", "ALT"),
      stringsAsFactors = FALSE
    )
    test = fisher.test(dat)
    p_value = test$p.value
    
    dataframe[row,names(dataframe)=="P_FISHER_VALUE"] = format(p_value, scientific=FALSE)
    if (p_value < 0.001){
      dataframe[row,names(dataframe)=="P_FISHER"] = "<0.001"
    } else {
      dataframe[row,names(dataframe)=="P_FISHER"] = round(p_value,3)
    }
  }
  
  dataframe$P_ADJUSTED = p.adjust(dataframe$P_FISHER_VALUE, method="bonferroni", n=length(dataframe$P_FISHER_VALUE))
  
  dataframe = dataframe[order(dataframe$P_FISHER_VALUE),]
  dataframe
}
