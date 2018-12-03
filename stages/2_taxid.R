source('parameters.R')

jetz_tree_flpth <- file.path(data_dir, 'filogenialorossingruposext1.tre')

tree <- ape::read.tree(file = jetz_tree_flpth)[[1]]
# look-up ncbi ID
txids <- taxize::get_uid(sub(pattern = '_', replacement = ' ', x = tree$tip.label))
cat('Could not find ....')
missing_nms <- tree$tip.label[is.na(txids)]
for (nm in missing_nms) {
  cat(nm, '\n')
}
write(as.character(txids)[!is.na(txids)], file = file.path(results_dir,
                                                           'txids.txt'))
