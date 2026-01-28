
library(dplyr)
library(data.table)
library(qqman)
library(argparser)

parser <- arg_parser("Process phenotype data")
parser <- add_argument(parser, "cohort", help = "cohort name (e.g., Health06)")
parser <- add_argument(parser, "sex", help = "Sex (Female or Male)")
args <- parse_args(parser)

# cohort <- "Health06/"
# sex <- "Male/"
cohort <- args$cohort
sex <- args$sex

input_dir <- paste0("/projects/magic-AUDIT/data/regenie/output_perTrait/", cohort, "/", sex, "/filter_maf/")
output_dir <- paste0("/projects/magic-AUDIT/data/regenie/output_perTrait/", cohort, "/", sex, "/QC_plot/")
dir.create(output_dir, recursive = TRUE)

file_list <- list.files(path = input_dir,
                        full.name =F)

qc_plot <- function(file){
  sum_stats <- paste0(input_dir, file)
  output <- paste0(output_dir, file)
  print("Processing file:")
  print(sum_stats)
  gwasResults <- fread(input =sum_stats, sep = " ", header = T, stringsAsFactors = F)
  print(paste('Number of SNP filtering by INFO >=0.3: ',nrow(gwasResults)))
  print(head(gwasResults))
  ## modify the data
  gwasResults <- gwasResults %>%
    dplyr::mutate(P = pchisq(CHISQ, 1, lower.tail = FALSE)) # %>%
    # dplyr::rename(SNP = SNP, BP = GENPOS, CHR = CHROM, ID = ID) %>%
    # dplyr::select(SNP, CHR, BP, P) %>%
    # dplyr::mutate(across(c(CHR, BP, P), as.numeric))
  print("plotting manhattan")
  threshold <- 5e-8
  lead_snps <- gwasResults$ID[gwasResults$P < threshold]
  png(paste0(output,".manhattan.png"), width = 1600, height = 800, type="cairo")
  #manhattan(gwasResults)
  manhattan(gwasResults, chr = "CHROM", bp = "GENPOS", snp = "ID", p = "P", 
            main = file, highlight = lead_snps)
  #manhattan(gwasResults, highlight = lead_snps, annotatePval = threshold)
  dev.off()
  ## QQ plot calculations
  p <- gwasResults$P[!is.na(gwasResults$P)]
  n <- length(p)
  x2obs <- qchisq(p,1,lower.tail=FALSE)
  x2exp <- qchisq((1:n - 0.5)/n,1,lower.tail=FALSE)
  lambda <- median(x2obs)/median(x2exp) # calculates your lambda value to check for inflation
  print("plotting QQ")
  #png(snakemake@output[["qq"]], type ="cairo")
  png(paste0(output,".qq.png"), width = 800, height = 800, type="cairo")
  qq(gwasResults$P, 
     main = paste("Lambda:", round(lambda, 3),basename(file)))
  dev.off()
  cat("Lambda for", file, "is:", lambda, "\n")
}

for(i in 1:length(file_list)){
  qc_plot(file_list[i])
}
