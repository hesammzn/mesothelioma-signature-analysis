# Sourcing the required script from panConusig
source("C:\\Users\\hesam\\Desktop\\CNV\\Sigminer\\panConusig\\panConusig_scripts\\R\\setup_CNsigs.R")


# Loading CNV table and reference table
data <- read.csv("C:\\Users\\hesam\\Desktop\\CNV\\Sigminer\\panConusig\\Clonal.csv")
sigs <-read.csv("C:\\Users\\hesam\\Desktop\\CNV\\Sigminer\\panConusig\\COSMIC_v3.4_CN_GRCh37.txt")

# Run the matrix generation
cn_matrix <- getMatrix(data)