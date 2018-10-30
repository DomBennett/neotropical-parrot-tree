drpOvrhngs <- function(als, ctff=.75) {
  for(i in seq_along(als)) {
    gps_mtrx <- matrix(FALSE, ncol=nchar(als[[i]][[1]]),
                       nrow=length(als[[i]]))
    for(j in seq_along(als[[i]])) {
      pull <- gregexpr('-', als[[i]][[j]])[[1]]
      if(pull[1] == -1) {
        next
      }
      gps_mtrx[j, pull] <- TRUE
    }
    prp_mssng <- 1 - (colSums(gps_mtrx)/nrow(gps_mtrx))
    ovrlppng <- which(prp_mssng > ctff)
    strt <- ovrlppng[1]
    end <- ovrlppng[length(ovrlppng)]
    for(j in seq_along(als[[i]])) {
      als[[i]][[j]] <- substr(als[[i]][[j]], start=strt, stop=end)
    }
  }
  als
}

spltUp <- function(al, ctff=.75, mn_lngth=500) {
  gps_mtrx <- matrix(FALSE, ncol=nchar(al[[1]]),
                     nrow=length(al))
  for(i in seq_along(al)) {
    pull <- gregexpr('-', al[[i]])[[1]]
    gps_mtrx[i, pull] <- TRUE
  }
  prp_mssng <- 1 - (colSums(gps_mtrx)/nrow(gps_mtrx))
  strts <- which((prp_mssng[-1] > ctff) & (prp_mssng[-1*length(prp_mssng)] < ctff)) + 1
  ends <- which((prp_mssng[-1] < ctff) & (prp_mssng[-1*length(prp_mssng)] > ctff))
  lngths <- ends - strts
  pull <- lngths > mn_lngth
  strts <- strts[pull]
  ends <- ends[pull]
  new_als <- vector('list', length=length(strts))
  for(i in seq_along(strts)) {
    strt <- strts[i]
    end <- ends[i]
    new_al <- vector('list', length=length(al))
    for(j in seq_along(al)) {
      new_al[[j]] <- substr(al[[j]], start=strt, stop=end)
    }
    names(new_al) <- names(al)
    new_als[[i]] <- new_al
  }
  new_als
}


readSqs <- function(fl) {
  all_data <- readLines(fl)
  sqs <- list()
  for(i in seq_along(all_data)) {
    bit <- all_data[[i]]
    if(grepl(pattern='^>', x=bit)) {
      nm <- sub(pattern='^>', '', x=bit)
      sqs[[nm]] <- ''
    } else {
      sqs[[nm]] <- paste0(sqs[[nm]], bit)
    }
  }
  sqs
}

writeSqs <- function(sqs, fl) {
  fasta <- ''
  for(i in seq_along(sqs)) {
    sq <- sqs[[i]]
    dfln <- names(sqs)[[i]]
    fasta <- paste0(fasta, '>', dfln, '\n',
                    sq, '\n\n')
  }
  cat(fasta, file=fl)
}

partition <- function(lngths, fl) {
  gene <- strt <- 1
  prttn_txt <- ''
  for(lngth in lngths) {
    end <- lngth + strt - 1
    prttn_txt <- paste0(prttn_txt, 'DNA, gene',
                        gene, ' = ', strt, '-',
                        end, '\n')
    strt <- end + 1
    gene <- gene + 1
  }
  cat(prttn_txt, file=fl)
}

# Combine alignments into a supermatrix, run RAxML

# VARS ----
wd <- 'phylotar_results'

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

