# Set the path to the source code
# setwd('TO_YOUR_WK_PATH')

source('BBmix.R')
source('ZicoSeq.R')

# Example
require(GUniFrac)
data(throat.otu.tab)
data(throat.meta)

comm <- t(throat.otu.tab)
meta.dat <- throat.meta

set.seed(123)

zico.obj <- ZicoSeq(meta.dat = meta.dat, comm = comm, grp.name = 'SmokingStatus', adj.name = 'Sex',
		# Filtering criterion
		prev.filter = 0.1, abund.filter = 10,  min.prop = 0, 
		# Winsorization to replace outliers
		is.winsor = TRUE, winsor.qt = 0.97,
		# Posterior sampling
		is.prior = TRUE, prior.dist = c('BetaMix'), post.method = c('sample'), post.sample.no = 25, 
		# Link functions
		link.func = list(function (x) x^0.25, function (x) x^0.5, function (x) x^0.75), stats.combine.func = max,
        # Permutation
		perm.no = 99,  strata = NULL, 
		# Multiple stage normalization
		stage.no = 6, topK = NULL, stage.fdr = 0.75, stage.max.pct = 0.50,  
		# Tree-based FDR control and family-wise error rate control
		is.fwer = FALSE, is.tree.fdr = FALSE, tree = NULL, 
		verbose = TRUE, return.comm = FALSE, return.perm.F = FALSE)

which(zico.obj$p.adj.fdr <= 0.1)
