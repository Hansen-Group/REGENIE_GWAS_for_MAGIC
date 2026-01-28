library(argparser)

parser <- arg_parser("Process phenotype data")
parser <- add_argument(parser, "cohort", help="cohort name (e.g., Inter99)")
args <- parse_args(parser)

input_file <- paste0("/projects/magic-AUDIT/data/magic_phenotypes_2024_combine/", args$cohort, "/", args$cohort, "_MAGIC_Male")
output_dir <- paste0(getwd(), "/", args$cohort)
dir.create(output_dir, showWarnings = FALSE)

table <- read.csv(input_file)
table[is.na(table)] <- "NA"
table = cbind(table$particid,table)
colnames(table)[1] <- "FID"
colnames(table)[2] <- "IID"
cov_col = c("FID","IID","age","sex","bmi","age2", "fasting_glc", "C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")
not_pheno_col = c("Cohort","age","sex","bmi","age2","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")
cov_table   = table[cov_col]
pheno_table = table[!colnames(table) %in% not_pheno_col]
id_table = table[c("FID","IID")]
write.table(cov_table, file = paste0(output_dir, "/",args$cohort, "_MAGIC_Male_covariates.txt"), row.names = FALSE, sep=' ', quote=FALSE)
write.table(pheno_table, file = paste0(output_dir, "/",args$cohort, "_MAGIC_Male_phenotype.txt"), row.names = FALSE, sep=' ', quote=FALSE)
write.table(id_table, file = paste0(output_dir, "/",args$cohort, "_MAGIC_Male_id.txt"), row.names = FALSE, col.names = FALSE, sep=' ', quote=FALSE)
# female
input_file <- paste0("/projects/magic-AUDIT/data/magic_phenotypes_2024_combine/", args$cohort, "/", args$cohort, "_MAGIC_Female")
table <- read.csv(input_file)
table[is.na(table)] <- "NA"
table = cbind(table$particid,table)
colnames(table)[1] <- "FID"
colnames(table)[2] <- "IID"
cov_col = c("FID","IID","age","sex","bmi","age2", "fasting_glc","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")
not_pheno_col = c("Cohort","age","sex","bmi","age2","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")
cov_table   = table[cov_col]
pheno_table = table[!colnames(table) %in% not_pheno_col]
id_table = table[c("FID","IID")]
write.table(cov_table, file = paste0(output_dir, "/",args$cohort, "_MAGIC_Female_covariates.txt"), row.names = FALSE, sep=' ', quote=FALSE)
write.table(pheno_table, file = paste0(output_dir, "/",args$cohort, "_MAGIC_Female_phenotype.txt"), row.names = FALSE, sep=' ', quote=FALSE)
write.table(id_table, file = paste0(output_dir, "/",args$cohort, "_MAGIC_Female_id.txt"), row.names = FALSE, col.names = FALSE, sep=' ', quote=FALSE)
