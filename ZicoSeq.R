# 2019_02_20
# Beta-mixture to address sequence depth confounding
# Permutation-based FDR control
# Multiple-stage normalization
# Variance modulation by LIMMA - hasn't been finished
# Expectation of F statistics based on posterior sampling
# Omnibus test for combining diffrent effects
# Version 0.9 [PermuteDAA3.2]

require(matrixStats)
require(permute)
require(nlme)
require(vegan)
require(ape)
require(limma)
require(statmod)

GMPR <- function (comm, intersect.no = 4, ct.min = 2, verbose = FALSE) {
	# Computes the GMPR size factor
	#
	# Args:
	#   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
	#   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
	#   ct.min: the minimum number of counts required to calculate ratios ct.min = 5 has better results
	
	#
	# Returns:
	#   a list that contains:
	#      gmpr： the GMPR size factors for all samples; Samples with distinct sets of features will be output as NA.
	#      nss:   number of samples with significant sharing (> intersect.no) including itself
	
	# mask counts < ct.min
	comm[comm < ct.min] <- 0
	
	if (is.null(colnames(comm))) {
		colnames(comm) <- paste0('S', 1:ncol(comm))
	}
	
	if (verbose == TRUE)
		cat('Begin GMPR size factor calculation ...\n')
	
	comm.no <- numeric(ncol(comm))
	gmpr <- sapply(1:ncol(comm),  function(i) {	
				
				if (i %% 50 == 0) {
					if (verbose == TRUE)
						cat(i, '\n')
				}
				x <- comm[, i]
				# Compute the pairwise ratio
				pr <- x / comm
				# Handling of the NA, NaN, Inf
				pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
				# Counting the number of non-NA, NaN, Inf
				incl.no <- colSums(!is.na(pr))		
				# Calculate the median of PR
				pr.median <- colMedians(pr, na.rm=TRUE)
				# Record the number of samples used for calculating the GMPR
				comm.no[i] <<- sum(incl.no >= intersect.no)
				# Geometric mean of PR median
				if (comm.no[i] > 1) {
					return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
				} else {
					return(NA)
				}
			}
	)
	
	if (sum(is.na(gmpr))) {
		warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
						'\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
						'For these samples, their size factors are set to be NA! \n', 
						'You may consider removing these samples since they are potentially outliers or negative controls!\n',
						'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
	}
	if (verbose == TRUE) {
		cat('Completed!\n')
		cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
	}
	
	attr(gmpr, 'NSS') <- comm.no
	names(gmpr) <- colnames(comm)
	return(gmpr * median(colSums(comm)))
}

# Please indicate the code is copied from vegan package for portability
getPermuteMatrix <- function (perm, N, strata = NULL) {
	if (length(perm) == 1) {
		perm <- how(nperm = perm)
	}
	if (!missing(strata) && !is.null(strata)) {
		if (inherits(perm, "how") && is.null(getBlocks(perm))) 
			setBlocks(perm) <- strata
	}
	if (inherits(perm, "how")) 
		perm <- shuffleSet(N, control = perm)
	if (is.null(attr(perm, "control"))) 
		attr(perm, "control") <- structure(list(within = list(type = "supplied matrix"), 
						nperm = nrow(perm)), class = "how")
	perm
}

# Permutation-based FDR control
perm_fdr_adj <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	F0 <- F0[ord]
	perm.no <- ncol(Fp)
	Fp <- as.vector(Fp)
	Fp <- Fp[!is.na(Fp)]
	Fp <- sort(c(Fp, F0), decreasing = F)
	n <- length(Fp)
	m <- length(F0)
	FPN <- (n + 1) - match(F0, Fp) - 1:m
	p.adj.fdr <- FPN / perm.no / (1:m)
#		p.adj.fdr <- sapply(F0, function(x) sum(Fp >= 
#									x, na.rm=TRUE) / perm.no)/(1:length(F0))
	p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
}

# Permutation-based FWER control
# Conservative
perm_fwer_adj <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	F0 <- F0[ord]
	col.max <- colMaxs(Fp, na.rm=TRUE)
	p.adj.fwer <- sapply(F0, function(x) mean(col.max >= x))[order(ord)]
}


# Westfall young too slow
perm_fwer_adj2 <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	m <- length(F0)
	F0 <- F0[ord]
	Fp <- Fp[ord, , drop=FALSE]
	col.max <- Fp[m, ]
	p.adj.fwer <- sapply(m:1, function(i) {
				x <- F0[i]
				y <- Fp[i, ]
				col.max <<- ifelse(y > col.max, y, col.max)
				mean(col.max >= x)
			})
	p.adj.fwer <- rev(p.adj.fwer)
	p.adj.fwer <- pmin(1, rev(cummin(rev(p.adj.fwer))))[order(ord)]
}

# Pad back the NA values
na.pad <- function (vec, ind) {
	vec0 <- numeric(length(ind))
	vec0[!ind] <- vec
	vec0[ind] <- NA
	vec0
}

# Pad back the NA values
one.pad <- function (vec, ind) {
	vec0 <- numeric(length(ind))
	vec0[!ind] <- vec
	vec0[ind] <- 1
	vec0
}


# Transform to Z score 
Ztransform <- function(p.value, e.sign, eff.sign=TRUE, tol=1E-15) {
	p.value[p.value <= tol] <- tol
	p.value[p.value >= 1 - tol] <- 1 - tol
	if (eff.sign == TRUE) {
		e.sign[e.sign == 0] <- sample(c(-1, 1), sum(e.sign == 0), replace=T)
		z1 <- qnorm(p.value / 2)
		z2 <- qnorm(1 - p.value / 2)
		z <- ifelse(e.sign > 0, z2, z1)
	} else {
		z <- qnorm(1 - p.value)
	}
	return(z)
}

# Internal functions
# For gls with generic correlation structure
corHerit <- function(value, paras, form = ~1, fixed = TRUE) {
	# Place holder - check the validity of the parameter
	object <- value
	attr(object, "formula") <- form
	attr(object, "fixed") <- fixed
	attr(object, "paras") <- paras
	class(object) <- c("corHerit", "corStruct")
	object
}

Initialize.corHerit <- function (object, data, ...) {	
	# Place holder - check the validity of the parameter
	form <- formula(object)
	if (!is.null(getGroupsFormula(form))) {
		attr(object, "groups") <- getGroups(object, form, data = data)
		attr(object, "Dim") <- Dim(object, attr(object, "groups"))
	} else {
		attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
	}
	attr(object, "covariate") <- getCovariate(object, data = data)
	object
}


corMatrix.corHerit <- function (object, covariate = getCovariate(object), ...) {
	
	paras <- attr(object, "paras")
	p <- paras[['p']]
	I <- diag(p)
	hyper <- as.vector(object)
	
	lambda <- 1 / (1 + exp(hyper[1]))
	V <- exp(-((paras[['D']])) * (exp(hyper[2])))
	cor.m <- (1 - lambda) * V +  lambda * I
	cor.m
}

# Extract the coefficient
coef.corHerit <- function (object,  ...) {
	
	paras <- attr(object, "paras")
	coefs <- as.vector(object)	
	coef1 <- 1 / (1 + exp(coefs[1]))
	coef1 <- coef1 / (1 - coef1)
	
	coef2 <-  exp(coefs[2])
	coefs <- c(coef1, coef2)	
	names(coefs) <- paste("Hyper", 1:length(coefs), sep="")
	coefs
}

# Estimating the hyper parameters
EstHyper <- function (y, D, init.val=c(0, 0)) {
	
	obj.gls <- gls(model = y ~ 1, 
			correlation = corHerit(value=init.val, paras = list(p=length(y), D = D)))
	cc <- c(coef(obj.gls$modelStruct$corStruct), obj.gls$coefficients, obj.gls$logLik)
	cc	
}	

AdjStats <- function (y, V, k, mu, fudge=0.005) {
	p <- nrow(V)
	# Add small fudge
	V.inv <- solve(V + fudge * diag(p))
	I <- diag(p)
	y.adj  <- solve(I + k * V.inv) %*% (k * mu * rowSums(V.inv) + y)
	y.adj
}


# Compute the NFP on the null.ind
PermFDR <- function (F0, Fp, null.ind) {
	
	Fp <- as.matrix(Fp)
	pct <- sum(null.ind) / length(F0)
	Fp <- Fp[null.ind, ]
	
	ord <- order(F0, decreasing = T)
	F0 <- F0[ord]
	perm.no <- ncol(Fp)
	Fp <- as.vector(Fp)
	
	pct.non.na <- mean(!is.na(Fp))
	
	Fp <- Fp[!is.na(Fp)]
	Fp <- sort(c(Fp, F0), decreasing = F)
	
	n <- length(Fp)
	m <- length(F0)
	
	FPN <- (n + 1) - match(F0, Fp) - 1:m
	# Handle identical F0 situations
	FPN <- cummax(FPN)
	
	p.adj.fdr <- FPN / pct.non.na / pct / perm.no / (1:m)
	
	p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
}

# Compute tree-based FDR control	
TreeFDR <- function (X = NULL, Y = NULL, test.func = NULL, perm.func = NULL, test.obs = NULL, test.perm = NULL, tree = NULL, D, 
		eff.sign = TRUE,  B = 20,  q.cutoff = 0.5, alpha = 1,
		adaptive = c('Fisher', 'Overlap'), alt.FDR = c('BH', 'Permutation'), verbose = TRUE, ...) {
	
	if (is.null(X) | is.null(Y) | is.null(test.func) | is.null(perm.func)) {
		if (is.null(test.obs) | is.null(test.perm)) {
			stop('Please specify either "X, Y, test.func, perm.fuc", or "test.obs", "test.perm"!\n')
		} else {
			X <- matrix(test.obs$p.value, ncol = 1)
			rownames(X) <- names(test.obs$p.value)
		}
	}
	adaptive <- match.arg(adaptive)
	alt.FDR <- match.arg(alt.FDR)
	
	# Make sure the rows of X and tree tips are labeled to avoid error
	if ((is.null(rownames(X)) | is.null(tree$tip.label))) {
		warning('Both the data matrix and the tree should have labels (rownames, tip.label) to avoid potential errors!\n')
	} else {
		if (sum(!(rownames(X) %in% tree$tip.label))){
			stop('Some features in the data matrix are not in the tree! Please check!\n')
		} else {
			if (sum(!(tree$tip.label %in% rownames(X)))) {
				warning('The tree have more features than the data matrix! \n')
			} 
		}
	}
	
	# Patristic distance
	if (!is.null(tree)) {
		D <- (cophenetic(tree)) ^ alpha
	} else {
		D <- D ^ alpha
	}
	
	
	if (!is.null(rownames(X)) & !is.null(tree$tip.label)) {
		D <- D[rownames(X), rownames(X)]
	}
	
	if (is.null(test.obs)) {
		if (verbose)
			cat('Test on original data sets  ...\n')
		test.obs <- test.func(X, Y, ...)
	}
	
	
	if (!is.list(test.obs) | !all(c("e.sign", "p.value") %in% 
					names(test.obs))) {
		stop("test.func should return a list with names e.sign and p.valueif z.transform=TRUE! Please check!\n")
	}
	
	null.ind <- test.obs$p.value >= quantile(test.obs$p.value, 1 - q.cutoff)
	
	z.obs <- Ztransform(test.obs$p.value, test.obs$e.sign, eff.sign)	
	
	if (verbose)
		cat('Test on permuted data sets  ...\n')
	
	if (is.null(test.perm)) {
		z.perm <- z.perm2 <- matrix(NA, nrow(X), B)	
		for (i in 1:B) {
			perm.obj <- perm.func(X, Y, ...)
			if (!is.list(perm.obj) | !all(c("X", "Y") %in% names(perm.obj))) {
				stop("perm.func should return a list with names X and Y! Please check!\n")
			}
			X.perm <- perm.obj$X
			Y.perm <- perm.obj$Y
			test.perm <- test.func(X.perm, Y.perm, ...)
			
			z.perm[, i] <- z.perm2[, i] <- Ztransform(test.perm$p.value, test.perm$e.sign, eff.sign)
			z.perm2[!null.ind, i] <- z.obs[!null.ind]
		}
	} else {
		B <- ncol(test.perm$p.value)
		z.perm <- z.perm2 <- matrix(NA, length(test.obs$p.value), B)	
		for (i in 1:B) {
			z.perm[, i] <- z.perm2[, i] <- Ztransform(test.perm$p.value[, i], test.perm$e.sign[, i], eff.sign)
			z.perm2[!null.ind, i] <- z.obs[!null.ind]
		}
		
	}
	
	
	if (alt.FDR == 'Permutation') {
		if (verbose)
			cat('Perform ordinary permutation-based FDR control ...\n')
		if (eff.sign == TRUE) {
			p.adj0 <- PermFDR(abs(z.obs), abs(z.perm), rep(TRUE, length(z.obs)))
		} else {
			p.adj0 <- PermFDR(z.obs, z.perm, rep(TRUE, length(z.obs)))
		}
	}
	
	if (alt.FDR == 'BH') {
		if (verbose)
			cat('Perform ordinary BH-based FDR control ...\n')
		p.adj0 <- p.adjust(test.obs$p.value, 'fdr')
		
	}
	
	if (verbose)
		cat("Estimating hyperparameter ... \n")
	error <- try(obj <- EstHyper(y = z.obs, D = D))
	if (inherits(error, "try-error")) {
		if (verbose)
			cat('Hyperparameter estimation failed! Ordinary permutation-based FDR control will be used!\n')
		# p.adj <- p.adjust(test.obs$p.value,'fdr')
		p.adj <- p.adj0
		k <- NULL
		rho <- NULL
		z.adj <- NULL
	} else {
		if (verbose)
			cat("Structure-based adjustment ...\n")
		k <- obj[1]
		rho <- obj[2]
		mu <- obj[3]
		V <- exp(-1 * rho * D)
		z.adj <- stat.o <- AdjStats(y = z.obs, V = V, k = k, mu = mu)[, 1]
		stat.p <- AdjStats(y = z.perm2, V = V, k = k, mu = mu)
		
		if (eff.sign == TRUE) {
			stat.o <- abs(stat.o)
			stat.p <- abs(stat.p)
		}
		
		p.adj <- PermFDR(stat.o, stat.p, null.ind)			
		
		# Power loss check
		if (adaptive == 'Overlap') {
			# Magic numbers
			fdr.cutoff <- 0.2
			pct.cutoff <- 0.5
			
			ind0 <- p.adj0 <= fdr.cutoff
			ind <- p.adj <= fdr.cutoff
			
			if (sum(p.adj[ind0] <= fdr.cutoff) <  pct.cutoff * sum(ind0)) {
				
				# Over-adjustment checking
				if (verbose)
					cat('Potential over-adjustment! Alternative FDR control will be used!\n')
				p.adj <- p.adj0
				k <- NULL
				rho <- NULL
				z.adj <- NULL
			}
		} 
		
		if (adaptive == 'Fisher') {
			# These cutoffs are used emprically
			fdr.cutoff <- 0.2
			
			ind0 <- p.adj0 <= fdr.cutoff
			ind <- p.adj <= fdr.cutoff
			
			n <- nrow(X)
			test.p <- fisher.test(matrix(c(sum(ind), sum(ind0), n - sum(ind), n - sum(ind0)), 2, 2), alternative = 'less')$p.value
			
			if (test.p <= 0.05) {
				# Over-adjustment checking
				if (verbose)
					cat('Potential over-adjustment! Alternative FDR control will be used!\n')
				p.adj <- p.adj0
				k <- NULL
				rho <- NULL
				z.adj <- NULL
			} 
		}
	}
	
	if (verbose)
		cat("Done!\n")
	
	return(list(p.adj = p.adj,  p.unadj = test.obs$p.value, z.adj = z.adj, z.unadj = z.obs, k = k, rho = rho))
}



# Permutation-based differential abundance analysis
ZicoSeq <- function (
		meta.dat, comm, grp.name, adj.name = NULL, 
		prev.filter = 0.1, abund.filter = 10,  min.prop = 0, 
		is.winsor = TRUE, winsor.qt = 0.97,
		is.prior = TRUE, prior.dist = c('BetaMix', 'ZIBB'), 
		post.method = c('sample', 'mean'), post.sample.no = 25, 
		link.func = list(function (x) x^0.25, function (x) x^0.5, function (x) x^0.75),
		link.d.func = list(function (x) 0.25 * x^(-0.75), function (x) 0.5 * x^(-0.5), function (x) 0.75 * x^(-0.25)),
		variance.EB = FALSE, df.prior = 10, 
		perm.no = 99,  strata = NULL, stats.combine.func = max,
		stage.no = 6, topK = NULL, stage.fdr = 0.75, stage.max.pct = 0.50,  
		is.fwer = FALSE, is.tree.fdr = FALSE, tree = NULL, 
		verbose = TRUE, return.comm = FALSE, return.perm.F = FALSE, ...) {
	
	# Args:
	#   meta.dat: a data frame containing the sample information
	#   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
	#   grp.name: a character, variable of interest; it could be numeric or categorical; should be in "meta.dat"
	#   adj.name: a character vector, variable(s) to be adjusted; they could be numeric or categorical; should be in "meta.dat"
	#   prev.filter: features with prevalence (i.e., nonzero proportion) less than "prev.cutoff" or be filtered
	#   abund.filter: features with a total counts less than "abund.cutoff" or be filtered 
	#	is.winsor: a logical value indicating whether winsorization should be performed to replace outliers. The default is TRUE.		   
	#   winsor.qt: the winsorization quantile, above which the counts will be replaced
	#   is.prior: a logical value indicating whether to perform posterior inference based on some prior distribution on the proportion data
	#   prior.dist: prior distribution, either two-component beta-binomial mixture ("BetaMix") or zeroinflated beta-binomial ("ZIBB")
	#   post.method: method for posterior inference, either based on posterior sampling ("sample") or approximate posterior mean ("mean")
	#   post.sample.no:  the number of posterior samples if posterior sampling is used
	#   link.func:  a list of functions that connects the ratios to the covariates
	#   link.d.func: a list of the derivative function of "link.func"; only need to specifiy when "post.method" is "mean"
	#   variance.EB: a logical value indicating whehter to perform empirical Bayes based variance shrinkage
	#   df.prior: the degree of freedom of the prior inverse gamma distribution for variance shrinkage
	#   perm.no: the number of permutations; If the raw p values are of the major interest, set "perm.no" to at least 999
	#   strata:  a factor indicating the permutation strata; permutation will be confined to each stratum 
	#   stats.combine.func: function to combine the F-statistic for the omnibus test
	#   stage.no: the number of stages if multiple-stage ratio stategy is used
	#   topK: the number of dominant features that will be excluded in the initial stage ratio calculation
	#   stage.fdr: the fdr cutoff below which the features will be excluded for calculating the ratio
	#   stage.max.pct: the maximum percentage of features that will be excluded
	#   is.fwer: a logical value indicating whether the family-wise error rate control (West-Young) should be performed
	#   is.tree.fdr: a logical value indicating whether tree-based false discovery rate shuold be carried out
	#   tree:  a class of "phylo", the tree relats all the OTUs, and should have the same names in "comm"
	#   verbose: a logical value indicating whether the trace information should be printed out
	#   return.comm: a logical value indicating whether the wisorized, filtered "comm" matrix should be returned
	#   ...:  arguments passing to tree-based fdr control
	
	#
	# Returns:
	#   a list that contains:
	#      call: the call
	#      comm: the wisorized, filtered "comm" matrix
	#      filter.ind: a vector of logical values indicating which features are tested
	#      R2: a matrix of percent explained variance (number of features by number of functions)
	#      F0: a matrix of F-statistics (number of features by number of functions)
	#      RSS: a matrix of residual sum squares (number of features by number of functions)
	#      df.model, df.residual:  degree of freedoms for the model and residual space
	#      p.raw : the raw p-values based on permutations (not accurate if "perm.no" is small)
	#      p.adj.fdr: permutation-based FDR-adjusted p-values
	#      p.adj.tree.fdr: permutation-based tree FDR-adjusted p-values 
	#      p.adj.fwer: permutation-based FWER-adjusted (West-Young) p-values
	#      tree.fdr.obj: the object returned by the "TreeFDR"
	
	
	this.call <- match.call()
	prior.dist <- match.arg(prior.dist)
	post.method <- match.arg(post.method)
	
	###############################################################################
	# Winsorization to reduce the influence of outlier counts
	
	if (is.winsor == TRUE) {
		depth <- colSums(comm)
		comm <- apply(comm, 1, function (x) {
					x <- x / depth
					qt <- quantile(x, winsor.qt)
					x[x >= qt] <- qt
					round(x * depth)
				})
		comm <- t(comm)
		
	}
	
	###############################################################################
	# Filter to remove very rare taxa
	filter.ind <- rowMeans(comm != 0) >= prev.filter & rowSums(comm) >= abund.filter
	names(filter.ind) <- rownames(comm)
	
	if (verbose)  cat(sum(!filter.ind), ' features are filtered!\n')
	
	comm <- comm[filter.ind, ]
	
	sample.no <- ncol(comm)
	otu.no <- nrow(comm)
	row.names <- rownames(comm)
	depth <- colSums(comm)
	
	if (verbose) cat('The data has ', sample.no, ' samples and ', otu.no, ' features will be tested!\n' )
	
	rabund <- colMeans(t(comm) / colSums(comm))
	if (is.null(topK)) {
		topK <- round(length(rabund) * 0.4)  # Very aggressive
	}
	size.factor <- colSums(comm[order(rabund)[1 : (length(rabund) - topK)], ])
	
	###############################################################################
	# Generate samples from posterior distribution (stacking it)
	if (is.prior == TRUE) {
		
		if (verbose) cat('Fitting beta mixture ...\n')
		
		if (post.method == 'sample') {
			comm.p <- apply(comm, 1, function (x) {
						if (prior.dist == 'BetaMix') {
							err1 <- try(res <- bbmix.fit.MM(x, depth))
						} 
						
						if (prior.dist == 'ZIBB') {
							err1 <- try(res <- zibb.fit.MM(x, depth))
						} 
						
						# Handle error
						if (class(err1) != 'try-error') {
							if (prior.dist == 'BetaMix') {
								prop1.1 <- rbeta(sample.no * post.sample.no, shape1 = x + res$shape1.1, shape2 = res$shape1.2 + depth - x)
								prop1.2 <- rbeta(sample.no * post.sample.no, shape1 = x + res$shape2.1, shape2 = res$shape2.2 + depth - x)
								prop <- ifelse(runif(sample.no * post.sample.no) <= res$q1, prop1.1, prop1.2)
							}
							if (prior.dist == 'ZIBB') {
								prop1.2 <- rbeta(sample.no * post.sample.no, shape1 = x + res$shape2.1, shape2 = res$shape2.2 + depth - x)
								prop <- ifelse(runif(sample.no * post.sample.no) <= res$q1, 0, prop1.2)
							}
							
						} else {
							prop <- x / depth
							v <- var(prop)
							m <- mean(prop)
							
							a1 <- ((1 - m) / v - 1 / m) * m ^ 2
							a2 <- a1 * (1 / m - 1)
							
							if (is.na(a1) | a1 < 0) {
								# uniform prior
								prop <- rbeta(sample.no * post.sample.no, shape1 = x + 1, shape2 = otu.no + depth - x)
							} else {
								prop <- rbeta(sample.no * post.sample.no, shape1 = x + a1, shape2 = a2 + depth - x)
							}	
						}
						return(prop)
					})
			
			comm.p <- t(comm.p)
			comm.p.list <- list()
			st <- 1
			end <- sample.no
			for (i in 1 : post.sample.no) {
				comm.p.list[[i]] <- comm.p[, st : end]
				st <- st + sample.no
				end <- end + sample.no
			}
		} else {
			EV <- apply(comm, 1, function (x) {
						if (prior.dist == 'BetaMix') {
							err1 <- try(res <- bbmix.fit.MM(x, depth))
						} 
						
						if (prior.dist == 'ZIBB') {
							err1 <- try(res <- zibb.fit.MM(x, depth))
						} 
						
						# Handle error
						if (class(err1) != 'try-error') {
							if (prior.dist == 'BetaMix') {
								Ep <- res$q1 *  (x + res$shape1.1) / (depth + res$shape1.1 + res$shape1.2) +
										(1 - res$q1) * (x + res$shape2.1) / (depth + res$shape2.1 + res$shape2.2)
								Vp <- res$q1 *  ((x + res$shape1.1)^2 / (depth + res$shape1.1 + res$shape1.2)^2 + 
											(x + res$shape1.1) * (res$shape1.2 + depth - x) / (depth + res$shape1.1 + res$shape1.2)^2 /
											(depth + res$shape1.1 + res$shape1.2 + 1)) +
										(1 - res$q1) * ((x + res$shape2.1)^2 / (depth + res$shape2.1 + res$shape2.2)^2 + 
											(x + res$shape2.1) * ( res$shape2.2 + depth - x) / (depth + res$shape2.1 + res$shape2.2)^2 /
											(depth + res$shape2.1 + res$shape2.2 + 1)) - Ep^2
							}
							if (prior.dist == 'ZIBB') {
								Ep <- (1 - res$q1) * (x + res$shape2.1) / (depth + res$shape2.1 + res$shape2.2)
								Vp <- (1 - res$q1) * ((x + res$shape2.1)^2 / (depth + res$shape2.1 + res$shape2.2)^2 + 
											(x + res$shape2.1) * ( res$shape2.2 + depth - x) / (depth + res$shape2.1 + res$shape2.2)^2 /
											(depth + res$shape2.1 + res$shape2.2 + 1)) - Ep^2
								
							}
							
						} else {
							
							prop <- x / depth
							v <- var(prop)
							m <- mean(prop)
							
							a1 <- ((1 - m) / v - 1 / m) * m ^ 2
							a2 <- a1 * (1 / m - 1)
							if (is.na(a1) | a1 < 0) {
								# uniform prior
								
								Ep <- (x + 1) / (depth + otu.no)
								Vp <- (x + 1) * (depth - x + otu.no - 1) / (otu.no + depth)^2 /
										(otu.no + depth + 1)
							} else {
								Ep <- (x + a1) / (depth + a1 + a2)
								Vp <- (x + a1) * (depth - x + a2) / (depth + a1 + a2)^2 /
										(depth + a1 + a2 + 1)
							}	
						}
						return(c(Ep, Vp))
					})
			EV <- t(EV)
			
			comm.p.list <- list()
			comm.v.list <- list()
			comm.p.list[[1]] <- EV[, 1 : sample.no]
			comm.v.list[[1]] <- EV[, (sample.no + 1) : (2 * sample.no)]
			post.sample.no <- 1
		}
	} else {
		comm.p.list <- list()
		comm.p.list[[1]] <- t(t(comm) / depth)
		post.sample.no <- 1
	}
	
	# Replace zeros or extremely small values for log calculation
	for (i in 1 : post.sample.no) {
		temp <- comm.p.list[[i]] 
		temp[temp <= min.prop] <- min.prop
		comm.p.list[[i]] <- temp
	}
	
	###############################################################################
	# Covariate space (including intercept)
	if (!is.null(strata)) {
		strata <- factor(strata)
	}
	
	if (is.null(adj.name)) {
		M0 <- model.matrix(~ 1, meta.dat) 
	} else {
		data0 <- meta.dat[, c(adj.name), drop = FALSE]
		M0 <- model.matrix( ~ ., data0) 
	}
	
	data1 <- meta.dat[, c(grp.name), drop = FALSE]
	M1 <-  model.matrix( ~ ., data1)[, -1, drop = FALSE]  # No intercept
	
	M01 <- cbind(M0, M1)
	
	# QR decompostion
	qrX0 <- qr(M0, tol = 1e-07)
	Q0 <- qr.Q(qrX0)
	Q0 <- Q0[, 1:qrX0$rank, drop = FALSE]
	H0 <- (Q0 %*% t(Q0))
	
	qrX1 <- qr(M1, tol = 1e-07)
	Q1 <- qr.Q(qrX1)
	Q1 <- Q1[, 1:qrX1$rank, drop = FALSE] 
	
	qrX01 <- qr(M01, tol = 1e-07)
	Q01 <- qr.Q(qrX01)
	Q01 <- Q01[, 1:qrX01$rank, drop = FALSE] 
	
	R0 <- as.matrix(resid(lm(Q1 ~ Q0 - 1)))
	
	pX0 <- ncol(Q0)
	pX1 <- ncol(Q1)
	pX01 <- ncol(Q01)
	
	df.model <- pX01 - pX0
	df.residual <- sample.no - pX01
	
	func.no <- length(link.func)
	
	###############################################################################
	# Perform multiple stage normalization
    norm.ind <- NULL
	for (i in 1:stage.no) {
		
		if (verbose == TRUE)
			cat('Stage ', i, ' ')
		
		# Reference proportion
		divisor <- size.factor / depth
		
		# Create the giant Y matrix
		Y <- matrix(NA, sample.no, func.no * otu.no * post.sample.no)
		
		if (post.method == 'mean') {
			D <- matrix(NA, sample.no, func.no * otu.no * post.sample.no)
		}
		
		
		# Change order - func.no * otu.no
		for (k in 1 : post.sample.no) {
			for (j in 1 : func.no) {
				func <- link.func[[j]]
				comm.p <- comm.p.list[[k]]
				
				if (post.method == 'mean') {
					dfunc <- link.d.func[[j]]
					comm.v <- comm.v.list[[k]]
				}
				
				Y[, (k - 1) * func.no * otu.no + func.no * (0 : (otu.no - 1)) + j] <-
						func(t(comm.p) / divisor)  # No scaling first
				
				if (post.method == 'mean') {
					# row - sample, column - features
					D[, (k - 1) * func.no * otu.no + func.no * (0 : (otu.no - 1)) + j] <-
							(dfunc(t(comm.p) / divisor) / divisor)^2 * t(comm.v)
				}		
			}
		}
		
		
		Y <- t(Y)
		if (post.method == 'sample') {
			TSS <- rowSums(Y^2)
			MSS01 <- rowSums((Y %*% Q01)^2)
			MSS0 <- rowSums((Y %*% Q0)^2)
		} else {
			w0 <- colSums((Q0 %*% t(Q0))^2)
			w01 <- colSums((Q01 %*% t(Q01))^2)
			
			TSS <- rowSums(Y^2) + colSums(D)
			MSS01 <- rowSums((Y %*% Q01)^2) + colSums(w01 * D)
			MSS0 <- rowSums((Y %*% Q0)^2) + colSums(w0 * D)
		}
		
		
		MSS <- (MSS01 - MSS0) 
		RSS <- (TSS - MSS01) 
		
		perm.ind <- getPermuteMatrix(perm.no, sample.no, strata = strata)
		perm.no <- nrow(perm.ind)
		
		MRSSp <- sapply(1 : perm.no, function(ii) {
					if (verbose) {
						if (ii %% 10 == 0) cat('.')
					}
					
					Rp <- R0[perm.ind[ii, ], , drop = FALSE]
					
					# Project to the reisdual space
					Rp <- Rp - H0 %*% Rp
					
					qrRp <- qr(Rp, tol = 1e-07)
					Q1p <- qr.Q(qrRp)
					Q1p <- Q1p[, 1:qrRp$rank, drop = FALSE]
					
					if (post.method == 'sample') {
						MSS01p <- MSS0 + rowSums((Y %*% Q1p)^2) 
					} else {
						w1p <- colSums((Q1p %*% t(Q1p))^2)
						MSS01p <- MSS0 + rowSums((Y %*% Q1p)^2) + colSums(w1p * D)
					}
					
					MSSp <- (MSS01p - MSS0) 
					RSSp <- (TSS - MSS01p) 
					
					c(MSSp, RSSp)
					
				})
		
		
		unit <- func.no * otu.no * post.sample.no
		MSSp <- MRSSp[1 : unit, ]
		RSSp <- MRSSp[(unit + 1) : (2 * unit), ]
		
		# EB is based on the aggregated RSS 
		RSS.m <- array(RSS, c(func.no, otu.no,  post.sample.no)) 
		RSS.m <- t(apply(RSS.m, c(1, 2), mean))  # otu.no by func.no
		
		F0.m <- array((MSS / df.model) / (RSS / df.residual), c(func.no, otu.no,  post.sample.no))
		F0.m <-  t(apply(F0.m, c(1, 2), mean))  # otu.no by func.no
		
		R2.m <- array(MSS / TSS, c(func.no, otu.no, post.sample.no))
		R2.m <- t(apply(R2.m, c(1, 2), mean))  # otu.no by func.no
		
		###############################################################################	
		# Variance moderation - to be improved
		if (variance.EB) {
#			hist(RSS.m / df.residual)
#			out <- apply(RSS.m, 2, function (x) {
#						out <- limma::squeezeVar(x / df.residual, df.residual)
#						c(out$var.prior, out$df.prior)
#					})
#			var.prior <- out[1, ]
#			df.prior <- out[2, ]
#			
#			df.prior[df.prior > 1E6] <- 1E6
#			df.prior[df.prior < 1E-6] <- 1E-6
			
			var.prior <- colMeans(RSS.m / df.residual, na.rm = TRUE)
			
			F0 <- (MSS / df.model)  /  ((RSS  + df.prior * var.prior) / (df.prior + df.residual))
			Fp <- (MSSp / df.model) /  ((RSSp  + df.prior * var.prior) / (df.prior + df.residual))
			
			# Expectation of F0 and Fp
			F0 <- array(F0, c(func.no, otu.no, post.sample.no)) 
			Fp <- array(Fp, c(func.no, otu.no, post.sample.no, perm.no)) 
			
			F0 <- apply(F0, c(1, 2), mean)    # func.no * otu.no
			Fp <- apply(Fp, c(1, 2, 4), mean) # func.no * otu.no * perm.no
			
		} else {
			
			F0 <- (MSS / df.model)  /  (RSS  / df.residual)
			Fp <- (MSSp / df.model)  /  (RSSp  / df.residual)
			
			# Expectation of F0 and Fp
			F0 <- array(F0, c(func.no, otu.no, post.sample.no)) 
			Fp <- array(Fp, c(func.no, otu.no, post.sample.no, perm.no)) 
			
			F0 <- apply(F0, c(1, 2), mean)    # func.no * otu.no
			Fp <- apply(Fp, c(1, 2, 4), mean) # func.no * otu.no * perm.no
		}
		
		###############################################################################
		# Omnibus test by taking maximum
		F0 <- apply(F0, 2, stats.combine.func)
		Fp <- apply(Fp, c(2, 3), stats.combine.func)  # otu.no by perm.no
		
		if (verbose) cat('\n')
		
		if (mean(is.na(F0)) >= 0.1) {
			warning('More than 10% observed F stats have NA! Please check! \n')
		}
		
		if (mean(is.na(Fp)) >= 0.1) {
			warning('More than 10% permuted F stats have NA! Please check! \n')
		}
		
		na.ind <- is.na(F0)
		F0 <- F0[!na.ind]
		Fp <- Fp[!na.ind, ]
		rabund0 <- rabund[!na.ind] 
		
		which.nan.ind <- which(!na.ind)
		p.adj.fdr <- perm_fdr_adj(F0, Fp)
		
		
		if (i == stage.no) {
			break
		} else {
			# recalculating the size factor
			if (mean(p.adj.fdr <= stage.fdr) > stage.max.pct) {
				ind <- order(p.adj.fdr, -rabund0)[1 : round(length(p.adj.fdr) * stage.max.pct)]
			} else {
				ind <- which(p.adj.fdr < stage.fdr)
			}
			size.factor <- colSums(comm[setdiff(1:nrow(comm), which.nan.ind[ind]), ])
			norm.ind <- cbind(norm.ind, !(1:nrow(comm) %in% which.nan.ind[ind]))
		}

	}
	
	###############################################################################
	# Permutation-based false discovery rate control
	# Tree-based power improvement
	if (is.tree.fdr) {
		if (verbose) {
			cat("Perform Tree-based FDR control ...\n")
		}
		# Convert pseudo-F to pseudo-P
		if (is.null(tree)) {
			stop('Tree FDR requires a tree of class "phylo"!\n')
		}
		test.obs <- list()
		test.perm <- list()
		
		test.obs$p.value <- 1 - pf(F0, df.model, df.residual)
		test.perm$p.value <- 1 -  pf(Fp, df.model, df.residual)
		
		# Ignore the sign
		test.obs$e.sign <- sign(test.obs$p.value)
		test.perm$e.sign <- sign(test.perm$p.value)
		
		D <- cophenetic(tree)
		D <- D[rownames(comm), rownames(comm)]
		
		tree.fdr.obj <- TreeFDR(X = NULL, Y = NULL, test.obs = test.obs, test.perm = test.perm,  D = D, 
				eff.sign = FALSE, alt.FDR = 'Permutation', verbose = verbose, ...)
		p.adj.tree.fdr <- tree.fdr.obj$p.adj
		
		if (verbose) {
			cat('FWER control has not used the phylogenetic information!\n')
		}
		
	} else {
		# scaled F stat
		tree.fdr.obj <- NULL
		p.adj.tree.fdr <- NULL
	}
	
	if (is.fwer) {
		p.adj.fwer <- perm_fwer_adj2(F0, Fp)
	}
	
	p.raw <- rowMeans(cbind(Fp, F0) >= F0)
	p.raw <- na.pad(p.raw, na.ind)
	
	p.adj.fdr <- na.pad(p.adj.fdr, na.ind)
	
	names(p.raw) <- names(p.adj.fdr) <- rownames(R2.m) <- rownames(RSS.m) <- rownames(F0.m) <- row.names
	colnames(R2.m) <- colnames(F0.m) <- colnames(RSS.m) <- paste0('Func', 1 : func.no)
	
	if (is.tree.fdr)  {
		p.adj.tree.fdr <- na.pad(p.adj.tree.fdr, na.ind)
		names(p.adj.tree.fdr)  <- row.names
	}
	
	if (is.fwer) {
		p.adj.fwer <- na.pad(p.adj.fwer, na.ind)
		names(p.adj.fwer)  <- row.names
	} else {
		p.adj.fwer <- NULL
	}
	
	if (!return.comm) {
		comm <- NULL
	}
	
	if (!return.perm.F) {
		Fp <- NULL
	}
	
	if (verbose) cat('Completed!\n')
	
	return(list(call = this.call, comm = comm, filter.ind = filter.ind,
					R2 = R2.m, F0 = F0.m, Fp = Fp, RSS = RSS.m, df.model = df.model, df.residual = df.residual,
					p.raw = p.raw, p.adj.fdr = p.adj.fdr, p.adj.tree.fdr = p.adj.tree.fdr, p.adj.fwer = p.adj.fwer, 
					tree.fdr.obj = tree.fdr.obj, norm.ind = norm.ind))
	
}






