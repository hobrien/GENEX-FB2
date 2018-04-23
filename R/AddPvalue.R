suppressMessages(library(dplyr))
suppressMessages(library(readr))

suppressMessages(library(argparser))

# parse inputs
p <- arg_parser("Selects the lowest p-value for each SNP")
p <- add_argument(p, "input", help="")
p <- add_argument(p, "output", help="")
args <- parse_args(p)

ph <- read_tsv(args$input)
ph <- mutate(ph, Coefficient_p = 2*pnorm(-abs(`Coefficient_z-score`))) %>%
  arrange(Coefficient_p)
write_tsv(ph, args$output)
