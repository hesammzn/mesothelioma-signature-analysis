library(sigminer)
library(data.table)


#Loading the dataset and reference table
cn_matrix <- fread("C:\\Users\\hesam\\Desktop\\irp5\\277\\Sig_CNV\\Clonal\\277.CNV48.matrix_Clonal.tsv", sep = "\t")
rownames(cn_matrix) <- cn_matrix[[1]]  
cn_matrix[[1]] <- NULL
cosmic_sig <- fread("C:\\Users\\hesam\\Desktop\\CNV\\Sigminer\\panConusig\\COSMIC_v3.4_CN_GRCh37.txt")
cn_matrix <- as.data.frame(cn_matrix)


#Converting to data.frame
cosmic_sig <- as.data.frame(cosmic_sig)

#Setting rownames from the 'Type' column
rownames(cosmic_sig) <- cosmic_sig$Type

#Dropping the 'Type' column
cosmic_sig$Type <- NULL

#Converting to matrix
cosmic_sig <- as.matrix(cosmic_sig)
cn_matrix <- as.matrix(cn_matrix)

#Fitting using Sigminer
fit_result <- sig_fit(cn_matrix, sig = cosmic_sig)


output_dir <- ("C:/Users/hesam/Desktop/irp5/277/Sig_CNV/Clonal/")
write.csv(fit_result, file = file.path(output_dir, "Fit_result_Clonal_CNV.csv"))