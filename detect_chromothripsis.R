# Code to extract chromothripsis chromosomes
# Derk Overduin, Philippe Bors, Luuk Nolden

# Load/install required packages (and dependencies)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools, cowplot, methods, BiocGenerics, graph, S4Vectors, GenomicRanges, IRanges, MASS, ggplot2, grid, gridExtra)

# Install and load shatterseek
install_github("parklab/ShatterSeek")
library(ShatterSeek)

# Read arguments provided by calling python code
args <- commandArgs(trailingOnly = TRUE)

cn_file <- args[1]
sv_file <- args[2]

cn_raw <- read.table(cn_file,       # TXT data file indicated as string or full path to the file
           header = TRUE,           # Whether to display the header (TRUE) or not (FALSE)
           sep = "\t",              # Separator of the columns of the file
           dec = ".")               # Character used to separate decimals of the numbers in the file

# Load the Copy Number (CN) profiles 
CN_data <- CNVsegs(chrom=as.character(cn_raw$chromosome),
                   start=cn_raw$start,
                   end=cn_raw$end,
                   total_cn=cn_raw$total_cn)

# Load the Structural Variation data
sv_raw <- read.table(sv_file,         # TXT data file indicated as string or full path to the file
           header = TRUE,             # Whether to display the header (TRUE) or not (FALSE)
           sep = "\t",                # Separator of the columns of the file
           dec = ".")                 # Character used to separate decimals of the numbers in the file

SV_data <- SVs(chrom1=as.character(sv_raw$chrom1),
               pos1=as.numeric(sv_raw$start1),
               chrom2=as.character(sv_raw$chrom2),
               pos2=as.numeric(sv_raw$end2),
               SVtype=as.character(sv_raw$svclass),
               strand1=as.character(sv_raw$strand1),
               strand2=as.character(sv_raw$strand2))


# Find chromothripsis in this data
chromothripsis <- shatterseek(SV.sample=SV_data, seg.sample=CN_data)

# Retrieve the chromothripsis data
summary <- chromothripsis@chromSummary

# Write the result to a csv file
write.csv(summary,'./csv/summary.csv')

# Gracefully quit the application to let python know it is time for the next genome to process
quit()