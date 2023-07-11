# Filters bgen files to selected genomic range and selected cohorts
# Calculates AF and GT

commandline = function(...){system(paste(...))}
commandline_test = function(...){message(paste(...))}

bgen2vcf = function(chr,start,end,cohorts){

  ## TOOLS
  plink2 = "/home/ReproSNP/TOOLS/plink2"
  bgenix = "/home/ReproSNP/TOOLS/bgenix"
  
  ## DIRECTORIES
  bgen_dir = "/home/UKB-WES/ukb23159/ukb23159-bgen/"
  sample_dir = "/home/UKB-WES/ukb23159/ReproSNP-IDs-sample/"
  tmp = "/home/ReproSNP/TMP/"
  
  ## FILES
  bgen = paste0(bgen_dir,"ukb23159_c",chr,"_b0_v1.bgen")
  sample = paste0(sample_dir,"ukb23159_c",chr,"_ReproSNP-IDs.sample")
  
  if(chr=="X"){
    chr="23"
  }
  range=paste0(chr,":",start,"-",end)
  
  tmp_bgen = paste0(tmp,"temporary.bgen")
  
  # Filter bgen to selected chromosome range
  commandline(bgenix,
              "-g", bgen,
              "-incl-range", range,
              ">", tmp_bgen)

  AF_output_files = c()
  
  sex_info = paste0(tmp,"temporary_sex_info.txt")
  commandline("echo '#IID SEX' >", sex_info, "&&",
              "cat", sample,
              "| awk 'NR > 2 { print }' | awk '{ print$1\"_\"$2, $4 }' >>", sex_info)
  
  for (cohort in 1:cohorts){
    
    # Set IDs file
    IDs = paste0(tmp,"temporary_IDs_cohort", cohort, ".txt")
    
    # Filter to selected pIDs
    filtered_vcf = paste0(tmp,"temporary_ID-filtered_cohort",cohort)
    commandline(plink2,
                "--bgen", tmp_bgen, "ref-first",
                "--sample", sample,
                "--keep", IDs,
                "--recode vcf",
                "--out", filtered_vcf)

    # Calculate allele counts
    allele_frequencies = paste0(tmp,"temporary_AF_cohort",cohort)
    commandline(plink2,
                "--vcf", paste0(filtered_vcf, ".vcf"),
                "--freq counts",
                "--update-sex", sex_info,
                "--out", allele_frequencies)
    
    # Calculate genotype counts
    genotypes = paste0(tmp,"temporary_GT_cohort",cohort)
    commandline(plink2,
                "--vcf", paste0(filtered_vcf, ".vcf"),
                "--geno-counts",
                "--update-sex", sex_info,
                "--out", genotypes)
    
    AF_output_files = c(AF_output_files, paste0(allele_frequencies, ".acount")) ##afreq
  }
  
  AF_output_files
  
}





