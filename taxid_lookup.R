library(taxize)
library(ape)

tree <- read.tree('arbol_consen_loros_Myersoloparrots_Oc2018.txt')
# look-up ncbi ID
txids <- get_uid(sub(pattern = '_', replacement = ' ', x = tree$tip.label))
cat('Could not find ....')
missing_nms <- tree$tip.label[is.na(txids)]
for (nm in missing_nms) {
  cat(nm, '\n')
}
write(as.character(txids)[!is.na(txids)], file = 'txids.txt')
