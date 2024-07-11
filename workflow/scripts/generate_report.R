suppressWarnings(suppressMessages(library(rmarkdown)))

args <- commandArgs(trailingOnly = TRUE)
root_dir <- args[1]
annot <- args[2]
conv <- args[3]
motif_method <- args[4]
rmd_file <- args[5]

output_file <- "final_report.html"

output_path <- file.path(root_dir, output_file)

rmarkdown::render(args[5], 
                  params = list(root_dir = root_dir,
                                annot = annot,
                                conv = conv,
                                motif_method = motif_method),
                  output_file = output_path)