# Libs
library(ape)

# Vars
wd <- 'phylotar_results'

# Read
trstr <- readLines(file.path(wd, 'consensus.tre'))
trstr <- gsub(':[0-9.]+', '', trstr)
trstr <- gsub('(\\[|\\])', '', trstr)
tree <- read.tree(text = trstr)

# Reduce
tp_lbls <- tree[['tip.label']]
tp_lbls <- sub('_[0-9]+', '', tp_lbls)
to_drp <- tree[['tip.label']][duplicated(tp_lbls)]
tree <- drop.tip(tree, to_drp)
tree[['tip.label']] <- sub('_[0-9]+', '', tree[['tip.label']])

# Root
tree <- unroot(tree)
outgroup <- "Aprosmictus_erythropterus"#, "Aprosmictus_jonquillaceus")
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)

# BRANCH SUPPORTS
spprt <- tree$node.label
spprt <- suppressWarnings(as.numeric(spprt))
spprt[is.na(spprt)] <- 0
nd_lbls <- rep('', length(spprt))
nd_lbls[spprt > 50] <- '*'
nd_lbls[spprt > 75] <- '**'
nd_lbls[spprt > 95] <- '***'

# Plot
png('tree.png', width = 2000, height = 2000)
par(mar = c(.1, .1, .1, .1))
plot(tree, edge.width = 4, cex = 1.5)
nodelabels(text = nd_lbls, frame = 'none', cex = 2.5, adj = -.25)
dev.off()

# Save
message(length(tree$tip.label), ' of 133')
write.tree(phy = tree, file = 'tree_for_sharing.tre')
