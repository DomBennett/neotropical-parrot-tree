source('parameters.R')
source(file.path('tools', 'supermatrix.R'))

# Combine alignments into a supermatrix, run RAxML

# VARS ----
wd <- file.path(results_dir, 'phylotar_results')

# INPUT
alfls <- list.files(wd, pattern = 'alignment[0-9]+.fasta')
als <- vector('list', length = length(alfls))
for (i in seq_along(alfls)) {
  als[[i]] <- readSqs(file.path(wd, alfls[i]))
}

# DROP OVERHANGING EDGES
(sapply(als, function(x) nchar(x[[1]])))
als <- drpOvrhngs(als, ctff = 0.75)
(sapply(als, function(x) nchar(x[[1]])))
n_taxa <- sapply(als, length)

# GEN PARTITION TEXT
lngths <- sapply(als, function(x) nchar(x[[1]]))
partition(lngths, fl = file.path(wd, 'partition.txt'))

# GEN SUPERMARTIX
fllrs <- sapply(lngths, function(x) paste0(rep('-', x), collapse = ''))
all_nms <- unique(unlist(sapply(als, names)))
all_nms <- sort(all_nms)
pull <- !grepl('\\ssp\\.', all_nms)
all_nms <- all_nms[pull]
sprmtrx <- vector('list', length = length(all_nms))
names(sprmtrx) <- all_nms
for (nm in all_nms) {
  al <- ''
  for (i in seq_along(als)) {
    tmp <- als[[i]][[nm]]
    tmp <- ifelse(is.null(tmp), fllrs[[i]], tmp)
    al <- paste0(al, tmp)
  }
  sprmtrx[[nm]] <- al
}

# DROP TAXA WITH TOO MANY GAPS
ngaps <- sapply(gregexpr('-', sprmtrx), length)
pull <- ngaps < nchar(sprmtrx[[1]])
sprmtrx <- sprmtrx[pull]

# CHECK AND WRITE OUT
all(sapply(sprmtrx, nchar) == nchar(sprmtrx[[1]]))
names(sprmtrx) <- gsub('\\s', '_', names(sprmtrx))
writeSqs(sprmtrx, fl = file.path(wd, 'supermatrix.fasta'))

# RAxML ----
# Warning: partition.txt may need minor modification depending on gene type
inpt <- file.path(wd, 'supermatrix.fasta')
prttnfl <- file.path(wd, 'partition.txt')
library(outsider)
repo <- 'DomBennett/om..raxml..8.2.12.pthreads.sse3'
raxml <- module_import('raxml', repo)
raxml('-f', 'a', '-m', 'GTRGAMMA', '-T', '2', '-#', '100', '-p',
      sample(0:10000000, 1), '-x', sample(0:10000000, 1),
      '-n', 'parrots', '-s', inpt, '-q', prttnfl)

# consensus
raxml('-m', 'GTRCAT', '-J', 'MR', '-z', 'RAxML_bootstrap.parrots', '-n',
      'parrots_con')

# CLEAN-UP
file.rename('RAxML_bestTree.parrots', file.path(wd, 'best_tree.tre'))
file.rename('RAxML_bootstrap.parrots', file.path(wd, 'bootstraps.tre'))
file.rename('RAxML_MajorityRuleConsensusTree.parrots_con',
            file.path(wd, 'consensus.tre'))
file.remove(list.files(pattern = 'RAxML.*\\.parrots'))
