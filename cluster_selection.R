# Libs ----
library(phylotaR)

# Vars ----
wd <- 'phylotar_results'

# Input ----
all_cls <- read_phylota(wd)

sp_cls <- drop_by_rank(all_cls, rnk = 'species', n = 2,
                       choose_by = c("nncltds", 'pambgs'),
                       greatest = c(TRUE, FALSE))
smmry <- summary(sp_cls)
smmry <- smmry[smmry[['N_taxa']] > 20, ]
smmry <- smmry[smmry[['MAD']] > 0.4, ]
smmry <- smmry[order(smmry$N_taxa, decreasing = TRUE), ]

# SELECT
slctd_smmry <- smmry[1:10, ]
slctd_smmry$ID <- as.numeric(slctd_smmry$ID)
slctd <- drop_clstrs(sp_cls, as.character(slctd_smmry$ID))
write.csv(slctd_smmry, 'best_clusters.csv')

# OUTPUT
# write out top 10 clusters with most taxa
sqfls <- list.files(wd, pattern = 'sequences[0-9]+.fasta')
for (i in seq_along(slctd@cids)) {
  cid <- slctd@cids[i]
  sids <- slctd@clstrs[[cid]]@sids
  txids <- get_txids(slctd, cid = cid, rnk = 'species')
  scnms <- get_tx_slot(slctd, txids, 'scnm')
  n <- sapply(seq_along(scnms), function(x) 
    sum(scnms[x] == scnms[x:length(scnms)]))
  sq_nm <- paste0(scnms, '_', n)
  infile <- file.path(wd, paste0('sequences', i, '.fasta'))
  write_sqs(phylota = slctd, outfile = infile, sid = sids, sq_nm = sq_nm)
}

# ALIGN
alfls <- list.files(wd, pattern = 'alignment[0-9]+.fasta')
for (i in seq_along(slctd@cids)) {
  inpt <- file.path(wd, paste0('sequences', i, '.fasta'))
  otpt <- file.path(wd, paste0('alignment', i,'.fasta'))
  system(paste0('mafft --auto ', inpt, ' > ', otpt))
}
