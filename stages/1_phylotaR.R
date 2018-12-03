source('parameters.R')
library(phylotaR)

# Setup ----
txids <- as.character(read.delim('txids.txt', header = FALSE)[ ,1])
ncbi_dr <- '/usr/bin'
wd <- file.path(results_dir, 'phylotar_results')
if (file.exists(wd)) {
  unlink(wd, recursive = TRUE)
}
dir.create(wd)
setup(wd = wd, txid = txids, ncbi_dr = ncbi_dr, v = TRUE)

# Run ----
run(wd)

# ps <- parameters_load(wd)
# txids <- txids_get(ps = ps)
# trm <- paste0(paste0('txid', ps[['txid']],'[Subtree]'), collapse = ' OR ')
# retmax = 1E4
# srch_rs <- rentrez::entrez_search(db = 'taxonomy', term = trm, retmax = retmax)
# srch_rs <- search_and_cache(func = rentrez::entrez_search,
#                             args = args, fnm = 'search', ps = ps)
