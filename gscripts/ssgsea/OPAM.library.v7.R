library(CvM2SL2Test)
library(MASS)
library(verification)

OPAM.project.dataset.RNAi <- function( 
		input.ds,            
		output.ds,
		gene.set.databases,
		gene.set.selection  = "ALL",  # "ALL" or list with names of gene sets
		sample.norm.type    = "rank",  # "rank", "log" or "log.rank"
		weight              = 0.25,
		statistic           = "area.under.RES",
		output.score.type   = "ES",  # "ES" or "NES"
		nperm               = 200,  # number of random permutations for NES case
		combine.mode        = "combine.off",  # "combine.off" do not combine *_UP and *_DN versions in 
		# a single score. "combine.replace" combine *_UP and 
		# *_DN versions in a single score that replaces the individual
		# *_UP and *_DN versions. "combine.add" combine *_UP and 
		# *_DN versions in a single score and add it but keeping 
		# the individual *_UP and *_DN versions.
                min.overlap         = 1,
		correl.type  = "rank")             # "rank", "z.score", "symm.rank"
{ #----------------------------------------------------------------------------------------
	
	# Load libraries
	library(gtools)
	library(verification)
	library(RColorBrewer)
	
	# Read input dataset
	
        dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
        m <- data.matrix(dataset$ds)
#	gene.names <- dataset$row.names   # in Ataris or hairpin gct files the gene symbols are in the descs column
      	gene.names <- dataset$descs
	gene.descs <- dataset$descs
	sample.names <- dataset$names
	Ns <- length(m[1,])
	Ng <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract input file name
	s <- length(temp[[1]])
	input.file.name <- temp[[1]][s]
	temp <- strsplit(input.file.name, split=".gct")
	input.file.prefix <-  temp[[1]][1]
	
	# Sample normalization
	
	if (sample.norm.type == "rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- 10000*m/Ng
	} else if (sample.norm.type == "log.rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- log(10000*m/Ng + exp(1))
	} else if (sample.norm.type == "log") {
		m[m < 1] <- 1
		m <- log(m + exp(1))
	}
	
	# Read gene set databases
	
	max.G <- 0
	max.N <- 0
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		max.G <- max(max.G, max(GSDB$size.G))
		max.N <- max.N +  GSDB$N.gs
	}
	N.gs <- 0
	gs <- matrix("null", nrow=max.N, ncol=max.G)
	gs.names <- vector(length=max.N, mode="character")
	gs.descs <- vector(length=max.N, mode="character")
	size.G <- vector(length=max.N, mode="numeric")
	start <- 1
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		N.gs <- GSDB$N.gs 
		gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
		gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
		size.G[start:(start + N.gs - 1)] <- GSDB$size.G
		gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
		start <- start + N.gs
	}
	N.gs <- max.N
	
	# Select desired gene sets
	
	if (gene.set.selection[1] != "ALL") {
		locs <- match(gene.set.selection, gs.names)
		N.gs <- sum(!is.na(locs))
		if(N.gs > 1) { 
                  gs <- gs[locs,]
		} else { 
                   gs <- t(as.matrix(gs[locs,]))   # Force vector to matrix if only one gene set specified
                }
		gs.names <- gs.names[locs]
 		gs.descs <- gs.descs[locs]
		size.G <- size.G[locs]
	}
	
	# Loop over gene sets
	
	score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
	for (gs.i in 1:N.gs) {
		#browser()
		gene.set <- gs[gs.i, 1:size.G[gs.i]]
		gene.overlap <- intersect(gene.set, gene.names)
		print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
                if (length(gene.overlap) < min.overlap) { 
			score.matrix[gs.i, ] <- rep(NA, Ns)
			next
		} else {
			gene.set.locs <- match(gene.overlap, gene.set)
			gene.names.locs <- match(gene.overlap, gene.names)
			msig <- m[gene.names.locs,]
			msig.names <- gene.names[gene.names.locs]
			if (output.score.type == "ES") {
				OPAM <- OPAM.Projection.RNAi(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = 1, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$ES.vector
			} else if (output.score.type == "NES") {
				OPAM <- OPAM.Projection.RNAi(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = nperm, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$NES.vector
			}
		}
	}

        
        locs <- !is.na(score.matrix[,1])
        print(paste("N.gs before overlap prunning:", N.gs))
        N.gs <- sum(locs)
        print(paste("N.gs after overlap prunning:", N.gs))
        score.matrix <- score.matrix[locs,]
        gs.names <- gs.names[locs]
        gs.descs <- gs.descs[locs]

	initial.up.entries <- 0
	final.up.entries <- 0
	initial.dn.entries <- 0
	final.dn.entries <- 0
	combined.entries <- 0
	other.entries <- 0
	
	if (combine.mode == "combine.off") {
		score.matrix.2 <- score.matrix
		gs.names.2 <- gs.names
		gs.descs.2 <- gs.descs
	} else if ((combine.mode == "combine.replace") || (combine.mode == "combine.add")) {
		score.matrix.2 <- NULL
		gs.names.2 <- NULL
		gs.descs.2 <- NULL
		k <- 1
		for (i in 1:N.gs) {
			temp <- strsplit(gs.names[i], split="_") 
			body <- paste(temp[[1]][seq(1, length(temp[[1]]) -1)], collapse="_")
			suffix <- tail(temp[[1]], 1)
			print(paste("i:", i, "gene set:", gs.names[i], "body:", body, "suffix:", suffix))
			if (suffix == "UP") {  # This is an "UP" gene set
				initial.up.entries <- initial.up.entries + 1
				target <- paste(body, "DN", sep="_")
				loc <- match(target, gs.names)            
				if (!is.na(loc)) { # found corresponding "DN" gene set: create combined entry
					score <- score.matrix[i,] - score.matrix[loc,]
					score.matrix.2 <- rbind(score.matrix.2, score)
					gs.names.2 <- c(gs.names.2, body)
					gs.descs.2 <- c(gs.descs.2, paste(gs.descs[i], "combined UP & DN"))
					combined.entries <- combined.entries + 1
					if (combine.mode == "combine.add") {  # also add the "UP entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.up.entries <- final.up.entries + 1
					}
				} else { # did not find corresponding "DN" gene set: create "UP" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.up.entries <- final.up.entries + 1
				}
			} else if (suffix == "DN") { # This is a "DN" gene set
				initial.dn.entries <- initial.dn.entries + 1
				target <- paste(body, "UP", sep="_")
				loc <- match(target, gs.names)            
				if (is.na(loc)) { # did not find corresponding "UP" gene set: create "DN" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.dn.entries <- final.dn.entries + 1
				} else { # it found corresponding "UP" gene set
					if (combine.mode == "combine.add") { # create "DN" entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.dn.entries <- final.dn.entries + 1
					}
				}
			} else { # This is neither "UP nor "DN" gene set: create individual entry
				score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
				gs.names.2 <- c(gs.names.2, gs.names[i])
				gs.descs.2 <- c(gs.descs.2, gs.descs[i])
				other.entries <- other.entries + 1
			}
		} # end for loop over gene sets
		print(paste("initial.up.entries:", initial.up.entries))
		print(paste("final.up.entries:", final.up.entries))
		print(paste("initial.dn.entries:", initial.dn.entries))
		print(paste("final.dn.entries:", final.dn.entries))
		print(paste("other.entries:", other.entries))
		print(paste("combined.entries:", combined.entries))

		print(paste("total entries:", length(score.matrix.2[,1])))
	}            
	
	V.GCT <- data.frame(score.matrix.2)
	names(V.GCT) <- sample.names
	row.names(V.GCT) <- gs.names.2
	write.gct(gct.data.frame = V.GCT, descs = gs.descs.2, filename = output.ds)  
	
} # end of OPAM.project.RNAi.dataset

OPAM.project.dataset.4 <- function( 
		input.ds,
		output.ds,
		gene.set.databases,
		gene.set.selection  = "ALL",  # "ALL" or list with names of gene sets
		sample.norm.type    = "rank",  # "rank", "log" or "log.rank"
		weight              = 0.25,
		statistic           = "area.under.RES",
		output.score.type   = "ES",  # "ES" or "NES"
		nperm               = 200,  # number of random permutations for NES case
		combine.mode        = "combine.off",  # "combine.off" do not combine *_UP and *_DN versions in 
		# a single score. "combine.replace" combine *_UP and 
		# *_DN versions in a single score that replaces the individual
		# *_UP and *_DN versions. "combine.add" combine *_UP and 
		# *_DN versions in a single score and add it but keeping 
		# the individual *_UP and *_DN versions.
                min.overlap         = 1,
		correl.type  = "rank")             # "rank", "z.score", "symm.rank"
{ #----------------------------------------------------------------------------------------
	
	# Load libraries
	library(gtools)
	library(verification)
	library(RColorBrewer)
	
	# Read input dataset
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	gene.names <- dataset$row.names
	gene.descs <- dataset$descs
	sample.names <- dataset$names
	Ns <- length(m[1,])
	Ng <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract input file name
	s <- length(temp[[1]])
	input.file.name <- temp[[1]][s]
	temp <- strsplit(input.file.name, split=".gct")
	input.file.prefix <-  temp[[1]][1]
	
	# Sample normalization
	
	if (sample.norm.type == "rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- 10000*m/Ng
	} else if (sample.norm.type == "log.rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- log(10000*m/Ng + exp(1))
	} else if (sample.norm.type == "log") {
		m[m < 1] <- 1
		m <- log(m + exp(1))
	}
	
	# Read gene set databases
	
	max.G <- 0
	max.N <- 0
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		max.G <- max(max.G, max(GSDB$size.G))
		max.N <- max.N +  GSDB$N.gs
	}
	N.gs <- 0
	gs <- matrix("null", nrow=max.N, ncol=max.G)
	gs.names <- vector(length=max.N, mode="character")
	gs.descs <- vector(length=max.N, mode="character")
	size.G <- vector(length=max.N, mode="numeric")
	start <- 1
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		N.gs <- GSDB$N.gs 
		gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
		gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
		size.G[start:(start + N.gs - 1)] <- GSDB$size.G
		gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
		start <- start + N.gs
	}
	N.gs <- max.N
	
	# Select desired gene sets
	
	if (gene.set.selection[1] != "ALL") {
		locs <- match(gene.set.selection, gs.names)
		N.gs <- sum(!is.na(locs))
		if(N.gs > 1) { 
                  gs <- gs[locs,]
		} else { 
                   gs <- t(as.matrix(gs[locs,]))   # Force vector to matrix if only one gene set specified
                }
		gs.names <- gs.names[locs]
 		gs.descs <- gs.descs[locs]
		size.G <- size.G[locs]
	}
	
	# Loop over gene sets
	
	score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
	for (gs.i in 1:N.gs) {
		#browser()
		gene.set <- gs[gs.i, 1:size.G[gs.i]]
		gene.overlap <- intersect(gene.set, gene.names)
		print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
                if (length(gene.overlap) < min.overlap) { 
			score.matrix[gs.i, ] <- rep(NA, Ns)
			next
		} else {
			gene.set.locs <- match(gene.overlap, gene.set)
			gene.names.locs <- match(gene.overlap, gene.names)
			msig <- m[gene.names.locs,]
			msig.names <- gene.names[gene.names.locs]
			if (output.score.type == "ES") {
				OPAM <- OPAM.Projection.3(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = 1, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$ES.vector
			} else if (output.score.type == "NES") {
				OPAM <- OPAM.Projection.3(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = nperm, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$NES.vector
			}
		}
	}

        
        locs <- !is.na(score.matrix[,1])
        print(paste("N.gs before overlap prunning:", N.gs))
        N.gs <- sum(locs)
        print(paste("N.gs after overlap prunning:", N.gs))
        score.matrix <- score.matrix[locs,]
        gs.names <- gs.names[locs]
        gs.descs <- gs.descs[locs]

	initial.up.entries <- 0
	final.up.entries <- 0
	initial.dn.entries <- 0
	final.dn.entries <- 0
	combined.entries <- 0
	other.entries <- 0
	
	if (combine.mode == "combine.off") {
		score.matrix.2 <- score.matrix
		gs.names.2 <- gs.names
		gs.descs.2 <- gs.descs
	} else if ((combine.mode == "combine.replace") || (combine.mode == "combine.add")) {
		score.matrix.2 <- NULL
		gs.names.2 <- NULL
		gs.descs.2 <- NULL
		k <- 1
		for (i in 1:N.gs) {
			temp <- strsplit(gs.names[i], split="_") 
			body <- paste(temp[[1]][seq(1, length(temp[[1]]) -1)], collapse="_")
			suffix <- tail(temp[[1]], 1)
			print(paste("i:", i, "gene set:", gs.names[i], "body:", body, "suffix:", suffix))
			if (suffix == "UP") {  # This is an "UP" gene set
				initial.up.entries <- initial.up.entries + 1
				target <- paste(body, "DN", sep="_")
				loc <- match(target, gs.names)            
				if (!is.na(loc)) { # found corresponding "DN" gene set: create combined entry
					score <- score.matrix[i,] - score.matrix[loc,]
					score.matrix.2 <- rbind(score.matrix.2, score)
					gs.names.2 <- c(gs.names.2, body)
					gs.descs.2 <- c(gs.descs.2, paste(gs.descs[i], "combined UP & DN"))
					combined.entries <- combined.entries + 1
					if (combine.mode == "combine.add") {  # also add the "UP entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.up.entries <- final.up.entries + 1
					}
				} else { # did not find corresponding "DN" gene set: create "UP" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.up.entries <- final.up.entries + 1
				}
			} else if (suffix == "DN") { # This is a "DN" gene set
				initial.dn.entries <- initial.dn.entries + 1
				target <- paste(body, "UP", sep="_")
				loc <- match(target, gs.names)            
				if (is.na(loc)) { # did not find corresponding "UP" gene set: create "DN" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.dn.entries <- final.dn.entries + 1
				} else { # it found corresponding "UP" gene set
					if (combine.mode == "combine.add") { # create "DN" entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.dn.entries <- final.dn.entries + 1
					}
				}
			} else { # This is neither "UP nor "DN" gene set: create individual entry
				score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
				gs.names.2 <- c(gs.names.2, gs.names[i])
				gs.descs.2 <- c(gs.descs.2, gs.descs[i])
				other.entries <- other.entries + 1
			}
		} # end for loop over gene sets
		print(paste("initial.up.entries:", initial.up.entries))
		print(paste("final.up.entries:", final.up.entries))
		print(paste("initial.dn.entries:", initial.dn.entries))
		print(paste("final.dn.entries:", final.dn.entries))
		print(paste("other.entries:", other.entries))
		print(paste("combined.entries:", combined.entries))

		print(paste("total entries:", length(score.matrix.2[,1])))
	}            
	
	V.GCT <- data.frame(score.matrix.2)
	names(V.GCT) <- sample.names
	row.names(V.GCT) <- gs.names.2
	write.gct(gct.data.frame = V.GCT, descs = gs.descs.2, filename = output.ds)  
	
} # end of OPAM.project.dataset.4

OPAM.project.dataset.5 <- function( 
		input.ds,
		output.ds,
		gene.set.databases,
		gene.set.selection  = "ALL",  # "ALL" or list with names of gene sets
		sample.norm.type    = "rank",  # "rank", "log" or "log.rank"
		weight              = 0.25,
		statistic           = "area.under.RES",
		output.score.type   = "ES",  # "ES" or "NES"
		nperm               = 200,  # number of random permutations for NES case
		combine.mode        = "combine.off",  # "combine.off" do not combine *_UP and *_DN versions in 
		# a single score. "combine.replace" combine *_UP and 
		# *_DN versions in a single score that replaces the individual
		# *_UP and *_DN versions. "combine.add" combine *_UP and 
		# *_DN versions in a single score and add it but keeping 
		# the individual *_UP and *_DN versions.
                min.overlap         = 1,
		correl.type  = "rank")             # "rank", "z.score", "symm.rank"
{ #----------------------------------------------------------------------------------------
	
	# Load libraries
	library(gtools)
	library(verification)
	library(RColorBrewer)
	
	# Read input dataset
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	gene.names <- dataset$row.names
	gene.descs <- dataset$descs
	sample.names <- dataset$names
	Ns <- length(m[1,])
	Ng <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract input file name
	s <- length(temp[[1]])
	input.file.name <- temp[[1]][s]
	temp <- strsplit(input.file.name, split=".gct")
	input.file.prefix <-  temp[[1]][1]
	
	# Sample normalization
	
	if (sample.norm.type == "rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- 10000*m/Ng
	} else if (sample.norm.type == "log.rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- log(10000*m/Ng + exp(1))
	} else if (sample.norm.type == "log") {
		m[m < 1] <- 1
		m <- log(m + exp(1))
	}
	
	# Read gene set databases
	
	max.G <- 0
	max.N <- 0
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		max.G <- max(max.G, max(GSDB$size.G))
		max.N <- max.N +  GSDB$N.gs
	}
	N.gs <- 0
	gs <- matrix("null", nrow=max.N, ncol=max.G)
	gs.names <- vector(length=max.N, mode="character")
	gs.descs <- vector(length=max.N, mode="character")
	size.G <- vector(length=max.N, mode="numeric")
	start <- 1
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		N.gs <- GSDB$N.gs 
		gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
		gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
		size.G[start:(start + N.gs - 1)] <- GSDB$size.G
		gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
		start <- start + N.gs
	}
	N.gs <- max.N
	
	# Select desired gene sets
	
	if (gene.set.selection[1] != "ALL") {
		locs <- match(gene.set.selection, gs.names)
                print(rbind(gene.set.selection, locs))
		N.gs <- sum(!is.na(locs))
		if(N.gs > 1) { 
                  gs <- gs[locs,]
		} else { 
                   gs <- t(as.matrix(gs[locs,]))   # Force vector to matrix if only one gene set specified
                }
		gs.names <- gs.names[locs]
 		gs.descs <- gs.descs[locs]
		size.G <- size.G[locs]
	}

        # Check for redundant gene sets

        tab <- as.data.frame(table(gs.names))
        ind <- order(tab[, "Freq"], decreasing=T)
        tab <- tab[ind,]
        print(tab[1:10,])
        print(paste("Total gene sets:", length(gs.names)))
        print(paste("Unique gene sets:", length(unique(gs.names))))

        # Loop over gene sets
	
	score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
	for (gs.i in 1:N.gs) {
		#browser()
		gene.set <- gs[gs.i, 1:size.G[gs.i]]
		gene.overlap <- intersect(gene.set, gene.names)
		print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
                if (length(gene.overlap) < min.overlap) { 
			score.matrix[gs.i, ] <- rep(NA, Ns)
			next
		} else {
			gene.set.locs <- match(gene.overlap, gene.set)
			gene.names.locs <- match(gene.overlap, gene.names)
			msig <- m[gene.names.locs,]
			msig.names <- gene.names[gene.names.locs]
			if (output.score.type == "ES") {
				OPAM <- OPAM.Projection.3(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = 1, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$ES.vector
			} else if (output.score.type == "NES") {
				OPAM <- OPAM.Projection.3(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = nperm, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$NES.vector
			}
		}
	}

        
        locs <- !is.na(score.matrix[,1])
        print(paste("N.gs before overlap prunning:", N.gs))
        N.gs <- sum(locs)
        print(paste("N.gs after overlap prunning:", N.gs))
        score.matrix <- score.matrix[locs,]
        gs.names <- gs.names[locs]
        gs.descs <- gs.descs[locs]

	initial.up.entries <- 0
	final.up.entries <- 0
	initial.dn.entries <- 0
	final.dn.entries <- 0
	combined.entries <- 0
	other.entries <- 0
	
	if (combine.mode == "combine.off") {
		score.matrix.2 <- score.matrix
		gs.names.2 <- gs.names
		gs.descs.2 <- gs.descs
	} else if ((combine.mode == "combine.replace") || (combine.mode == "combine.add")) {
		score.matrix.2 <- NULL
		gs.names.2 <- NULL
		gs.descs.2 <- NULL
		k <- 1
		for (i in 1:N.gs) {
			temp <- strsplit(gs.names[i], split="_") 
			body <- paste(temp[[1]][seq(1, length(temp[[1]]) -1)], collapse="_")
			suffix <- tail(temp[[1]], 1)
			print(paste("i:", i, "gene set:", gs.names[i], "body:", body, "suffix:", suffix))
			if (suffix == "UP") {  # This is an "UP" gene set
				initial.up.entries <- initial.up.entries + 1
				target <- paste(body, "DN", sep="_")
				loc <- match(target, gs.names)            
				if (!is.na(loc)) { # found corresponding "DN" gene set: create combined entry
					score <- score.matrix[i,] - score.matrix[loc,]
					score.matrix.2 <- rbind(score.matrix.2, score)
					gs.names.2 <- c(gs.names.2, body)
					gs.descs.2 <- c(gs.descs.2, paste(gs.descs[i], "combined UP & DN"))
					combined.entries <- combined.entries + 1
					if (combine.mode == "combine.add") {  # also add the "UP entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.up.entries <- final.up.entries + 1
					}
				} else { # did not find corresponding "DN" gene set: create "UP" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.up.entries <- final.up.entries + 1
				}
			} else if (suffix == "DN") { # This is a "DN" gene set
				initial.dn.entries <- initial.dn.entries + 1
				target <- paste(body, "UP", sep="_")
				loc <- match(target, gs.names)            
				if (is.na(loc)) { # did not find corresponding "UP" gene set: create "DN" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.dn.entries <- final.dn.entries + 1
				} else { # it found corresponding "UP" gene set
					if (combine.mode == "combine.add") { # create "DN" entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.dn.entries <- final.dn.entries + 1
					}
				}
			} else { # This is neither "UP nor "DN" gene set: create individual entry
				score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
				gs.names.2 <- c(gs.names.2, gs.names[i])
				gs.descs.2 <- c(gs.descs.2, gs.descs[i])
				other.entries <- other.entries + 1
			}
		} # end for loop over gene sets
		print(paste("initial.up.entries:", initial.up.entries))
		print(paste("final.up.entries:", final.up.entries))
		print(paste("initial.dn.entries:", initial.dn.entries))
		print(paste("final.dn.entries:", final.dn.entries))
		print(paste("other.entries:", other.entries))
		print(paste("combined.entries:", combined.entries))

		print(paste("total entries:", length(score.matrix.2[,1])))
	}            

        # Check for redundant gene sets

        tab <- as.data.frame(table(gs.names.2))
        ind <- order(tab[, "Freq"], decreasing=T)
        tab <- tab[ind,]
        print(tab[1:20,])
        print(paste("Total gene sets:", length(gs.names.2)))
        print(paste("Unique gene sets:", length(unique(gs.names.2))))
        
	V.GCT <- data.frame(score.matrix.2)
	names(V.GCT) <- sample.names
	row.names(V.GCT) <- gs.names.2
	write.gct(gct.data.frame = V.GCT, descs = gs.descs.2, filename = output.ds)  
	
} # end of OPAM.project.dataset.5

OPAM.Projection.2 <- function(
		data.array,
		gene.names,
		n.cols,
		n.rows,
		weight = 0,
		statistic = "Kolmogorov-Smirnov",  # "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises",
		# "Anderson-Darling", "Zhang_A", "Zhang_C", "Zhang_K",
		# "area.under.RES", or "Wilcoxon"
		gene.set,
		nperm = 200,
		correl.type  = "rank")             # "rank", "z.score", "symm.rank"
{
	
	ES.vector <- vector(length=n.cols)
	NES.vector <- vector(length=n.cols)
	p.val.vector <- vector(length=n.cols)
	correl.vector <- vector(length=n.rows, mode="numeric")
	
# Compute ES score for signatures in each sample
	
#   print("Computing GSEA.....")
	phi <- array(0, c(n.cols, nperm))
	for (sample.index in 1:n.cols) {
		gene.list <- order(data.array[, sample.index], decreasing=T)            
		
		#      print(paste("Computing observed enrichment for UP signature in sample:", sample.index, sep=" ")) 
		gene.set2 <- match(gene.set, gene.names)
		
		if (weight == 0) {
			correl.vector <- rep(1, n.rows)
		} else if (weight > 0) {
			if (correl.type == "rank") {
				correl.vector <- data.array[gene.list, sample.index]
			} else if (correl.type == "symm.rank") {
				correl.vector <- data.array[gene.list, sample.index]
				correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)], 
						correl.vector,
						correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)]) 
			} else if (correl.type == "z.score") {
				x <- data.array[gene.list, sample.index]
				correl.vector <- (x - mean(x))/sd(x)
			}
		}
		GSEA.results <- GSEA.EnrichmentScore5(gene.list=gene.list, gene.set=gene.set2,
				statistic = statistic, alpha = weight, correl.vector = correl.vector)
		ES.vector[sample.index] <- GSEA.results$ES
		
		if (nperm == 0) {
			NES.vector[sample.index] <- ES.vector[sample.index]
			p.val.vector[sample.index] <- 1
		} else {
			for (r in 1:nperm) {
				reshuffled.gene.labels <- sample(1:n.rows)
				if (weight == 0) {
					correl.vector <- rep(1, n.rows)
				} else if (weight > 0) {
					correl.vector <- data.array[reshuffled.gene.labels, sample.index]
				} 
				GSEA.results <- GSEA.EnrichmentScore5(gene.list=reshuffled.gene.labels, gene.set=gene.set2,
						statistic = statistic, alpha = weight, correl.vector = correl.vector)
				phi[sample.index, r] <- GSEA.results$ES
			}
			if (ES.vector[sample.index] >= 0) {
				pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
				if (length(pos.phi) == 0) pos.phi <- 0.5
				pos.m <- mean(pos.phi)
				NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
				s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
				p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
			} else {
				neg.phi <-  phi[sample.index, phi[sample.index, ] < 0]
				if (length(neg.phi) == 0) neg.phi <- 0.5 
				neg.m <- mean(neg.phi)
				NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
				s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
				p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
			}
		}
	}
	return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector))
	
} # end of OPAM.Projection.2

OPAM.Projection.3 <- function(
		data.array,
		gene.names,
		n.cols,
		n.rows,
		weight = 0,
		statistic = "Kolmogorov-Smirnov",  # "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises",
		# "Anderson-Darling", "Zhang_A", "Zhang_C", "Zhang_K",
		# "area.under.RES", or "Wilcoxon"
		gene.set,
		nperm = 200,
		correl.type  = "rank")             # "rank", "z.score", "symm.rank"
# Runs a 2-3x faster (2-2.5x for ES statistic and 2.5-3x faster for area.under.ES statsitic)
# version of GSEA.EnrichmentScore.5 internally that avoids overhead from the function call.
{
	
	ES.vector <- vector(length=n.cols)
	NES.vector <- vector(length=n.cols)
	p.val.vector <- vector(length=n.cols)
	correl.vector <- vector(length=n.rows, mode="numeric")
	
# Compute ES score for signatures in each sample
	
#   print("Computing GSEA.....")
	phi <- array(0, c(n.cols, nperm))
	for (sample.index in 1:n.cols) {
		gene.list <- order(data.array[, sample.index], decreasing=T)            
		
		#      print(paste("Computing observed enrichment for UP signature in sample:", sample.index, sep=" ")) 

		gene.set2 <- match(gene.set, gene.names)
		
		if (weight == 0) {
			correl.vector <- rep(1, n.rows)
		} else if (weight > 0) {
			if (correl.type == "rank") {
				correl.vector <- data.array[gene.list, sample.index]
			} else if (correl.type == "symm.rank") {
				correl.vector <- data.array[gene.list, sample.index]
				correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)], 
						correl.vector,
						correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)]) 
			} else if (correl.type == "z.score") {
				x <- data.array[gene.list, sample.index]
				correl.vector <- (x - mean(x))/sd(x)
			}
		}
		### Olga's Additions ###
#		ptm.new = proc.time()
		tag.indicator <- sign(match(gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
		no.tag.indicator <- 1 - tag.indicator 
		N <- length(gene.list) 
		Nh <- length(gene.set2) 
		Nm <-  N - Nh 
		orig.correl.vector <- correl.vector
		if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
		ind = which(tag.indicator==1)
		correl.vector <- abs(correl.vector[ind])^weight
		
		
		sum.correl = sum(correl.vector)
		up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
		gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
		down = gaps/Nm
		
		RES = cumsum(c(up,up[Nh])-down)
		valleys = RES[1:Nh]-up
		
		max.ES = max(RES)
		min.ES = min(valleys)
		
		if( statistic == "Kolmogorov-Smirnov" ){
			if( max.ES > -min.ES ){
				ES <- signif(max.ES, digits=5)
				arg.ES <- which.max(RES)
			} else{
				ES <- signif(min.ES, digits=5)
				arg.ES <- which.min(RES)
			}
		}
		
		if( statistic == "area.under.RES"){
			if( max.ES > -min.ES ){
				arg.ES <- which.max(RES)
			} else{
				arg.ES <- which.min(RES)
			}
			gaps = gaps+1
			RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
			ES = sum(RES)
		}
		GSEA.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
#		new.time <<- new.time + (proc.time() - ptm.new)
		### End Olga's Additions ###
		#GSEA.results <- GSEA.EnrichmentScore5(gene.list=gene.list, gene.set=gene.set2,
		#		statistic = statistic, alpha = weight, correl.vector = correl.vector)
		ES.vector[sample.index] <- GSEA.results$ES
		
		if (nperm == 0) {
			NES.vector[sample.index] <- ES.vector[sample.index]
			p.val.vector[sample.index] <- 1
		} else {
			for (r in 1:nperm) {
				reshuffled.gene.labels <- sample(1:n.rows)
				if (weight == 0) {
					correl.vector <- rep(1, n.rows)
				} else if (weight > 0) {
					correl.vector <- data.array[reshuffled.gene.labels, sample.index]
				} 
#				GSEA.results <- GSEA.EnrichmentScore5(gene.list=reshuffled.gene.labels, gene.set=gene.set2,
#						statistic = statistic, alpha = weight, correl.vector = correl.vector)
				### Olga's Additions ###
				tag.indicator <- sign(match(reshuffled.gene.labels, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
				no.tag.indicator <- 1 - tag.indicator 
				N <- length(reshuffled.gene.labels) 
				Nh <- length(gene.set2) 
				Nm <-  N - Nh 
#   orig.correl.vector <- correl.vector
				if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
				ind <- which(tag.indicator==1)
				correl.vector <- abs(correl.vector[ind])^weight   
				
				sum.correl <- sum(correl.vector)
				up = correl.vector/sum.correl
				gaps = (c(ind-1, N) - c(0, ind))
				down = gaps/Nm
				
				RES = cumsum(c(up,up[Nh])-down)
				valleys = RES[1:Nh]-up
				
				max.ES = max(RES)
				min.ES = min(valleys)
				
				if( statistic == "Kolmogorov-Smirnov" ){
					if( max.ES > -min.ES ){
						ES <- signif(max.ES, digits=5)
						arg.ES <- which.max(RES)
					} else{
						ES <- signif(min.ES, digits=5)
						arg.ES <- which.min(RES)
					}
				}
				
				if( statistic == "area.under.RES"){
					if( max.ES > -min.ES ){
						arg.ES <- which.max(RES)
					} else{
						arg.ES <- which.min(RES)
					}
					gaps = gaps+1
					RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
					ES = sum(RES)
				}
				
				GSEA.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
				### End Olga's Additions ###
				phi[sample.index, r] <- GSEA.results$ES
			}
			if (ES.vector[sample.index] >= 0) {
				pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
				if (length(pos.phi) == 0) pos.phi <- 0.5
				pos.m <- mean(pos.phi)
				NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
				s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
				p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
			} else {
				neg.phi <-  phi[sample.index, phi[sample.index, ] < 0]
				if (length(neg.phi) == 0) neg.phi <- 0.5 
				neg.m <- mean(neg.phi)
				NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
				s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
				p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
			}
		}
	}
	return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector))
	
} # end of OPAM.Projection.3


OPAM.Projection.RNAi <- function(
		data.array,
		gene.names,
		n.cols,
		n.rows,
		weight = 0,
		statistic = "Kolmogorov-Smirnov",  # "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises",
		# "Anderson-Darling", "Zhang_A", "Zhang_C", "Zhang_K",
		# "area.under.RES", or "Wilcoxon"
		gene.set,
		nperm = 200,
		correl.type  = "rank")             # "rank", "z.score", "symm.rank"
# Runs a 2-3x faster (2-2.5x for ES statistic and 2.5-3x faster for area.under.ES statsitic)
# version of GSEA.EnrichmentScore.5 internally that avoids overhead from the function call.
{
	
	ES.vector <- vector(length=n.cols)
	NES.vector <- vector(length=n.cols)
	p.val.vector <- vector(length=n.cols)
	correl.vector <- vector(length=n.rows, mode="numeric")
	
# Compute ES score for signatures in each sample
	
#   print("Computing GSEA.....")
	phi <- array(0, c(n.cols, nperm))
	for (sample.index in 1:n.cols) {
		gene.list <- order(data.array[, sample.index], decreasing=T)            
		
		#      print(paste("Computing observed enrichment for UP signature in sample:", sample.index, sep=" ")) 

#		gene.set2 <- match(gene.set, gene.names)
                gene.set2 <- seq(1:length(gene.names))[!is.na(match(gene.names, gene.set))]
		
		if (weight == 0) {
			correl.vector <- rep(1, n.rows)
		} else if (weight > 0) {
			if (correl.type == "rank") {
				correl.vector <- data.array[gene.list, sample.index]
			} else if (correl.type == "symm.rank") {
				correl.vector <- data.array[gene.list, sample.index]
				correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)], 
						correl.vector,
						correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)]) 
			} else if (correl.type == "z.score") {
				x <- data.array[gene.list, sample.index]
				correl.vector <- (x - mean(x))/sd(x)
			}
		}
		### Olga's Additions ###
#		ptm.new = proc.time()
		tag.indicator <- sign(match(gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
		no.tag.indicator <- 1 - tag.indicator 
		N <- length(gene.list) 
		Nh <- length(gene.set2) 
		Nm <-  N - Nh 
		orig.correl.vector <- correl.vector
		if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
		ind = which(tag.indicator==1)
		correl.vector <- abs(correl.vector[ind])^weight
		
		
		sum.correl = sum(correl.vector)
		up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
		gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
		down = gaps/Nm
		
		RES = cumsum(c(up,up[Nh])-down)
		valleys = RES[1:Nh]-up
		
		max.ES = max(RES)
		min.ES = min(valleys)
		
		if( statistic == "Kolmogorov-Smirnov" ){
			if( max.ES > -min.ES ){
				ES <- signif(max.ES, digits=5)
				arg.ES <- which.max(RES)
			} else{
				ES <- signif(min.ES, digits=5)
				arg.ES <- which.min(RES)
			}
		}
		
		if( statistic == "area.under.RES"){
			if( max.ES > -min.ES ){
				arg.ES <- which.max(RES)
			} else{
				arg.ES <- which.min(RES)
			}
			gaps = gaps+1
			RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
			ES = sum(RES)
		}
		GSEA.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
#		new.time <<- new.time + (proc.time() - ptm.new)
		### End Olga's Additions ###
		#GSEA.results <- GSEA.EnrichmentScore5(gene.list=gene.list, gene.set=gene.set2,
		#		statistic = statistic, alpha = weight, correl.vector = correl.vector)
		ES.vector[sample.index] <- GSEA.results$ES
		
		if (nperm == 0) {
			NES.vector[sample.index] <- ES.vector[sample.index]
			p.val.vector[sample.index] <- 1
		} else {
			for (r in 1:nperm) {
				reshuffled.gene.labels <- sample(1:n.rows)
				if (weight == 0) {
					correl.vector <- rep(1, n.rows)
				} else if (weight > 0) {
					correl.vector <- data.array[reshuffled.gene.labels, sample.index]
				} 
#				GSEA.results <- GSEA.EnrichmentScore5(gene.list=reshuffled.gene.labels, gene.set=gene.set2,
#						statistic = statistic, alpha = weight, correl.vector = correl.vector)
				### Olga's Additions ###
				tag.indicator <- sign(match(reshuffled.gene.labels, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
				no.tag.indicator <- 1 - tag.indicator 
				N <- length(reshuffled.gene.labels) 
				Nh <- length(gene.set2) 
				Nm <-  N - Nh 
#   orig.correl.vector <- correl.vector
				if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
				ind <- which(tag.indicator==1)
				correl.vector <- abs(correl.vector[ind])^weight   
				
				sum.correl <- sum(correl.vector)
				up = correl.vector/sum.correl
				gaps = (c(ind-1, N) - c(0, ind))
				down = gaps/Nm
				
				RES = cumsum(c(up,up[Nh])-down)
				valleys = RES[1:Nh]-up
				
				max.ES = max(RES)
				min.ES = min(valleys)
				
				if( statistic == "Kolmogorov-Smirnov" ){
					if( max.ES > -min.ES ){
						ES <- signif(max.ES, digits=5)
						arg.ES <- which.max(RES)
					} else{
						ES <- signif(min.ES, digits=5)
						arg.ES <- which.min(RES)
					}
				}
				
				if( statistic == "area.under.RES"){
					if( max.ES > -min.ES ){
						arg.ES <- which.max(RES)
					} else{
						arg.ES <- which.min(RES)
					}
					gaps = gaps+1
					RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
					ES = sum(RES)
				}
				
				GSEA.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
				### End Olga's Additions ###
				phi[sample.index, r] <- GSEA.results$ES
			}
			if (ES.vector[sample.index] >= 0) {
				pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
				if (length(pos.phi) == 0) pos.phi <- 0.5
				pos.m <- mean(pos.phi)
				NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
				s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
				p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
			} else {
				neg.phi <-  phi[sample.index, phi[sample.index, ] < 0]
				if (length(neg.phi) == 0) neg.phi <- 0.5 
				neg.m <- mean(neg.phi)
				NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
				s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
				p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
			}
		}
	}
	return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector))
	
} # end of OPAM.Projection.3

OPAM.project.dataset.2 <- function( 
		input.ds,
		output.ds,
		gene.set.databases,
		gene.set.selection  = "ALL",  # "ALL" or list with names of gene sets
		sample.norm.type    = "rank",  # "rank", "log" or "log.rank"
		weight              = 0.25,
		statistic           = "area.under.RES",
		output.score.type   = "ES",  # "ES" or "NES"
		nperm               = 200,  # number of random permutations for NES case
		combine.mode        = "combine.off",  # "combine.off" do not combine *_UP and *_DN versions in 
		# a single score. "combine.replace" combine *_UP and 
		# *_DN versions in a single score that replaces the individual
		# *_UP and *_DN versions. "combine.add" combine *_UP and 
		# *_DN versions in a single score and add it but keeping 
		# the individual *_UP and *_DN versions.
		correl.type  = "rank")             # "rank", "z.score", "symm.rank"
{ #----------------------------------------------------------------------------------------
	
	# Load libraries
	library(gtools)
	library(verification)
	library(RColorBrewer)
	
	# Read input dataset
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	gene.names <- dataset$row.names
	gene.descs <- dataset$descs
	sample.names <- dataset$names
	gene.descs <- dataset$descs
	Ns <- length(m[1,])
	Ng <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract input file name
	s <- length(temp[[1]])
	input.file.name <- temp[[1]][s]
	temp <- strsplit(input.file.name, split=".gct")
	input.file.prefix <-  temp[[1]][1]
	
	# Sample normalization
	
	if (sample.norm.type == "rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- 10000*m/Ng
	} else if (sample.norm.type == "log.rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- log(10000*m/Ng + exp(1))
	} else if (sample.norm.type == "log") {
		m[m < 1] <- 1
		m <- log(m + exp(1))
	}
	
	# Read gene set databases
	
	max.G <- 0
	max.N <- 0
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		max.G <- max(max.G, max(GSDB$size.G))
		max.N <- max.N +  GSDB$N.gs
	}
	N.gs <- 0
	gs <- matrix("null", nrow=max.N, ncol=max.G)
	gs.names <- vector(length=max.N, mode="character")
	gs.descs <- vector(length=max.N, mode="character")
	size.G <- vector(length=max.N, mode="numeric")
	start <- 1
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		N.gs <- GSDB$N.gs 
		gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
		gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
		size.G[start:(start + N.gs - 1)] <- GSDB$size.G
		gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
		start <- start + N.gs
	}
	N.gs <- max.N
	
	# Select desired gene sets
	
	if (gene.set.selection[1] != "ALL") {
		locs <- match(gene.set.selection, gs.names)
		N.gs <- sum(!is.na(locs))
		gs <- gs[locs,]
		gs.names <- gs.names[locs]
		gs.descs <- gs.descs[locs]
		size.G <- size.G[locs]
	}
	
	# Loop over gene sets
	
	score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
	for (gs.i in 1:N.gs) {
		gene.set <- gs[gs.i, 1:size.G[gs.i]]
		gene.overlap <- intersect(gene.set, gene.names)
		print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
		if (length(gene.overlap) == 0) { 
			score.matrix[gs.i, ] <- runif(Ns, min=1E-06, max=1.1E-06)
			next
		} else {
			gene.set.locs <- match(gene.overlap, gene.set)
			gene.names.locs <- match(gene.overlap, gene.names)
			msig <- m[gene.names.locs,]
			msig.names <- gene.names[gene.names.locs]
			if (output.score.type == "ES") {
				OPAM <- OPAM.Projection.2(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = 1, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$ES.vector
			} else if (output.score.type == "NES") {
				OPAM <- OPAM.Projection.2(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = nperm, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$NES.vector
			}
		}
	}
	
	initial.up.entries <- 0
	final.up.entries <- 0
	initial.dn.entries <- 0
	final.dn.entries <- 0
	combined.entries <- 0
	other.entries <- 0
	
	if (combine.mode == "combine.off") {
		score.matrix.2 <- score.matrix
		gs.names.2 <- gs.names
		gs.descs.2 <- gs.descs
	} else if ((combine.mode == "combine.replace") || (combine.mode == "combine.add")) {
		score.matrix.2 <- NULL
		gs.names.2 <- NULL
		gs.descs.2 <- NULL
		k <- 1
		for (i in 1:N.gs) {
			temp <- strsplit(gs.names[i], split="_") 
			body <- paste(temp[[1]][seq(1, length(temp[[1]]) -1)], collapse="_")
			suffix <- tail(temp[[1]], 1)
			print(paste("i:", i, "gene set:", gs.names[i], "body:", body, "suffix:", suffix))
			if (suffix == "UP") {  # This is an "UP" gene set
				initial.up.entries <- initial.up.entries + 1
				target <- paste(body, "DN", sep="_")
				loc <- match(target, gs.names)            
				if (!is.na(loc)) { # found corresponding "DN" gene set: create combined entry
					score <- score.matrix[i,] - score.matrix[loc,]
					score.matrix.2 <- rbind(score.matrix.2, score)
					gs.names.2 <- c(gs.names.2, body)
					gs.descs.2 <- c(gs.descs.2, paste(gs.descs[i], "combined UP & DN"))
					combined.entries <- combined.entries + 1
					if (combine.mode == "combine.add") {  # also add the "UP entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.up.entries <- final.up.entries + 1
					}
				} else { # did not find corresponding "DN" gene set: create "UP" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.up.entries <- final.up.entries + 1
				}
			} else if (suffix == "DN") { # This is a "DN" gene set
				initial.dn.entries <- initial.dn.entries + 1
				target <- paste(body, "UP", sep="_")
				loc <- match(target, gs.names)            
				if (is.na(loc)) { # did not find corresponding "UP" gene set: create "DN" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.dn.entries <- final.dn.entries + 1
				} else { # it found corresponding "UP" gene set
					if (combine.mode == "combine.add") { # create "DN" entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.dn.entries <- final.dn.entries + 1
					}
				}
			} else { # This is neither "UP nor "DN" gene set: create individual entry
				score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
				gs.names.2 <- c(gs.names.2, gs.names[i])
				gs.descs.2 <- c(gs.descs.2, gs.descs[i])
				other.entries <- other.entries + 1
			}
		} # end for loop over gene sets
		print(paste("initial.up.entries:", initial.up.entries))
		print(paste("final.up.entries:", final.up.entries))
		print(paste("initial.dn.entries:", initial.dn.entries))
		print(paste("final.dn.entries:", final.dn.entries))
		print(paste("other.entries:", other.entries))
		print(paste("combined.entries:", combined.entries))
		print(paste("total entries:", length(score.matrix.2[,1])))
	}            
	
	V.GCT <- data.frame(score.matrix.2)
	names(V.GCT) <- sample.names
	row.names(V.GCT) <- gs.names.2
	write.gct(gct.data.frame = V.GCT, descs = gs.descs.2, filename = output.ds)  
	
} # end of OPAM.project.dataset.2

OPAM.project.dataset.3 <- function( 
		input.ds,
		output.ds,
		gene.set.databases,
		gene.set.selection  = "ALL",  # "ALL" or list with names of gene sets
		sample.norm.type    = "rank",  # "rank", "log" or "log.rank"
		weight              = 0.25,
		statistic           = "area.under.RES",
		output.score.type   = "ES",  # "ES" or "NES"
		nperm               = 200,  # number of random permutations for NES case
		combine.mode        = "combine.off",  # "combine.off" do not combine *_UP and *_DN versions in 
		# a single score. "combine.replace" combine *_UP and 
		# *_DN versions in a single score that replaces the individual
		# *_UP and *_DN versions. "combine.add" combine *_UP and 
		# *_DN versions in a single score and add it but keeping 
		# the individual *_UP and *_DN versions.
		correl.type  = "rank")             # "rank", "z.score", "symm.rank"
{ #----------------------------------------------------------------------------------------
	
	# Load libraries
	library(gtools)
	library(verification)
	library(RColorBrewer)
	
	# Read input dataset
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	gene.names <- dataset$row.names
	gene.descs <- dataset$descs
	sample.names <- dataset$names
	Ns <- length(m[1,])
	Ng <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract input file name
	s <- length(temp[[1]])
	input.file.name <- temp[[1]][s]
	temp <- strsplit(input.file.name, split=".gct")
	input.file.prefix <-  temp[[1]][1]
	
	# Sample normalization
	
	if (sample.norm.type == "rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- 10000*m/Ng
	} else if (sample.norm.type == "log.rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- log(10000*m/Ng + exp(1))
	} else if (sample.norm.type == "log") {
		m[m < 1] <- 1
		m <- log(m + exp(1))
	}
	
	# Read gene set databases
	
	max.G <- 0
	max.N <- 0
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		max.G <- max(max.G, max(GSDB$size.G))
		max.N <- max.N +  GSDB$N.gs
	}
	N.gs <- 0
	gs <- matrix("null", nrow=max.N, ncol=max.G)
	gs.names <- vector(length=max.N, mode="character")
	gs.descs <- vector(length=max.N, mode="character")
	size.G <- vector(length=max.N, mode="numeric")
	start <- 1
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		N.gs <- GSDB$N.gs 
		gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
		gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
		size.G[start:(start + N.gs - 1)] <- GSDB$size.G
		gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
		start <- start + N.gs
	}
	N.gs <- max.N
	
	# Select desired gene sets
	
	if (gene.set.selection[1] != "ALL") {
		locs <- match(gene.set.selection, gs.names)
		N.gs <- sum(!is.na(locs))
		if(N.gs > 1) { 
                  gs <- gs[locs,]
		} else { 
                   gs <- t(as.matrix(gs[locs,]))   # Force vector to matrix if only one gene set specified
                }
		gs.names <- gs.names[locs]
 		gs.descs <- gs.descs[locs]
		size.G <- size.G[locs]
	}
	
	# Loop over gene sets
	
	score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
	for (gs.i in 1:N.gs) {
		#browser()
		gene.set <- gs[gs.i, 1:size.G[gs.i]]
		gene.overlap <- intersect(gene.set, gene.names)
		print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
		if (length(gene.overlap) == 0) { 
			score.matrix[gs.i, ] <- runif(Ns, min=1E-06, max=1.1E-06)
			next
		} else {
			gene.set.locs <- match(gene.overlap, gene.set)
			gene.names.locs <- match(gene.overlap, gene.names)
			msig <- m[gene.names.locs,]
			msig.names <- gene.names[gene.names.locs]
			if (output.score.type == "ES") {
				OPAM <- OPAM.Projection.3(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = 1, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$ES.vector
			} else if (output.score.type == "NES") {
				OPAM <- OPAM.Projection.3(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = nperm, correl.type = correl.type)
				score.matrix[gs.i,] <- OPAM$NES.vector
			}
		}
	}
	
	initial.up.entries <- 0
	final.up.entries <- 0
	initial.dn.entries <- 0
	final.dn.entries <- 0
	combined.entries <- 0
	other.entries <- 0
	
	if (combine.mode == "combine.off") {
		score.matrix.2 <- score.matrix
		gs.names.2 <- gs.names
		gs.descs.2 <- gs.descs
	} else if ((combine.mode == "combine.replace") || (combine.mode == "combine.add")) {
		score.matrix.2 <- NULL
		gs.names.2 <- NULL
		gs.descs.2 <- NULL
		k <- 1
		for (i in 1:N.gs) {
			temp <- strsplit(gs.names[i], split="_") 
			body <- paste(temp[[1]][seq(1, length(temp[[1]]) -1)], collapse="_")
			suffix <- tail(temp[[1]], 1)
			print(paste("i:", i, "gene set:", gs.names[i], "body:", body, "suffix:", suffix))
			if (suffix == "UP") {  # This is an "UP" gene set
				initial.up.entries <- initial.up.entries + 1
				target <- paste(body, "DN", sep="_")
				loc <- match(target, gs.names)            
				if (!is.na(loc)) { # found corresponding "DN" gene set: create combined entry
					score <- score.matrix[i,] - score.matrix[loc,]
					score.matrix.2 <- rbind(score.matrix.2, score)
					gs.names.2 <- c(gs.names.2, body)
					gs.descs.2 <- c(gs.descs.2, paste(gs.descs[i], "combined UP & DN"))
					combined.entries <- combined.entries + 1
					if (combine.mode == "combine.add") {  # also add the "UP entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.up.entries <- final.up.entries + 1
					}
				} else { # did not find corresponding "DN" gene set: create "UP" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.up.entries <- final.up.entries + 1
				}
			} else if (suffix == "DN") { # This is a "DN" gene set
				initial.dn.entries <- initial.dn.entries + 1
				target <- paste(body, "UP", sep="_")
				loc <- match(target, gs.names)            
				if (is.na(loc)) { # did not find corresponding "UP" gene set: create "DN" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.dn.entries <- final.dn.entries + 1
				} else { # it found corresponding "UP" gene set
					if (combine.mode == "combine.add") { # create "DN" entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.dn.entries <- final.dn.entries + 1
					}
				}
			} else { # This is neither "UP nor "DN" gene set: create individual entry
				score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
				gs.names.2 <- c(gs.names.2, gs.names[i])
				gs.descs.2 <- c(gs.descs.2, gs.descs[i])
				other.entries <- other.entries + 1
			}
		} # end for loop over gene sets
		print(paste("initial.up.entries:", initial.up.entries))
		print(paste("final.up.entries:", final.up.entries))
		print(paste("initial.dn.entries:", initial.dn.entries))
		print(paste("final.dn.entries:", final.dn.entries))
		print(paste("other.entries:", other.entries))
		print(paste("combined.entries:", combined.entries))
		print(paste("total entries:", length(score.matrix.2[,1])))
	}            
	
	V.GCT <- data.frame(score.matrix.2)
	names(V.GCT) <- sample.names
	row.names(V.GCT) <- gs.names.2
	write.gct(gct.data.frame = V.GCT, descs = gs.descs.2, filename = output.ds)  
	
} # end of OPAM.project.dataset.2


OPAM.match.projection.to.phenotypes <-  function(
		input.ds,
		input.cls,
		results.dir,
		normalize.score = T,
		normalization.type = "zero.one",  # "zero.one", "z.score" or "r.z.score"
		markers.num=5,
		user.colors = NA,
		markers.metric = "ROC",   # "ROC" or "T.TEST"
		markers.file = NULL,
		sort.phenotypes = T,
		sort.decreasing = T,    # T = decreasing, F = increasing
		sort.expression = T,
		sort.decreasing.genes = T,
		legend = T,
		char.res = 1,
		only.up = F,
		cmap.type = 3,
		show.desc = T,
		row.norm = T)
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
        dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)	
        m <- data.matrix(dataset$ds)
        model.names <- dataset$row.names
        model.descs <- dataset$descs
        n.models <- 
        Ns <- length(m[1,])

	for (i in 1:length(m[,1])) {
		if (sd(m[i,]) == 0) {
			val <- m[i, 1]
			m[i,] <- m[i,] + runif(n=Ns, min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
		}
	}
	dim(m)
	sample.names <- dataset$names

	temp <- strsplit(input.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
#   char.res <-  0.013 * n.models + 0.65
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list 
	
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors)) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
	
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
#      if (is.vector(cls.list)) {
#         cls.phen <- paste(phen.names, cls.phen, collapse="_")
#      } else {
#         for (i in 1:length(cls.phen)) {
#            for (j in 1:length(cls.phen[[i]])) {
#               cls.phen[[i]][j] <- paste(phen.names[i], cls.phen[[i]][j], collapse="_")
#            }
#         }
#      }
	} else {
		phen.names <- "NA"
	}
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		for (i in 1:n.phen) {
			n.classes[i] <- length(cls.phen[[i]])
		}
	}
	
	x <- rbind(sample.names, cls.list, cls.labels)
	print("before loop")
	print(x)
	print(cls.phen)
	print(phen.names)
	
	filename <- paste(results.dir, test.file.prefix, ".PHEN.MARKERS.", markers.metric, ".pdf", sep="")
	pdf(file=filename, height = 10, width = 10)
	
	# Loop over phenotypes
	
	for (k.phen in 1:n.phen) {
		
		
		if (is.vector(cls.labels)) {
			k.phen.labels <- cls.labels
			k.phen.list <- cls.list
		} else {
			k.phen.labels <- as.vector(cls.labels[k.phen,])
			k.phen.list <- as.vector(cls.list[k.phen,])
		}
		
		# Sort according to current phenotype
		
		if(sort.expression == T) {
			phen.index <- order(k.phen.labels, decreasing=sort.decreasing)
		} else {
			phen.index <- seq(1, length(k.phen.labels))
		}
		if (is.vector(cls.labels)) {
			cls.labels2 <- cls.labels[phen.index]
			cls.list2 <- cls.list[phen.index]
		} else {
			cls.labels2 <- cls.labels[, phen.index]
			cls.list2 <- cls.list[, phen.index]
		}
		k.phen.labels <- k.phen.labels[phen.index]
		k.phen.list <- k.phen.list[phen.index]
		sample.names2 <- sample.names[phen.index]
		m2 <- m[, phen.index]
		
		x <- rbind(sample.names2, cls.list2, cls.labels2)
		print(paste("inside loop phen=", k.phen))
		print(x)
		print(cls.phen)
		print(phen.names)
		
		# Markers for each class
		
		if (is.vector(cls.labels2)) {
			classes <- unique(cls.list2)
		} else {
			classes <- unique(cls.list2[k.phen, ])
		}
		if (length(classes) > 2) {
			k.only.up <- T
		} else {
			k.only.up <- only.up
		}
		
		if(length(classes) == 2) classes <- classes[1]
		markers <- NULL
		markers.descs <- NULL
		metric.list <- NULL
		p.val.list <- NULL
		k.class <- NULL
		for (k in classes) {
			if (is.vector(cls.labels2)) {
				bin.class <- ifelse(cls.list2 == k, 0, 1)
			} else {
				bin.class <- ifelse(cls.list2[k.phen, ] == k, 0, 1)
			}
			if (markers.metric == "T.TEST") {
				metric <- vector(length=n.models, mode="numeric")
				p.val <- vector(length=n.models, mode="numeric")
				for (i in 1:n.models) {
					temp <- split(m2[i, ], bin.class)
					x <- temp[[1]]
					y <- temp[[2]]
					metric[i] <- signif(t.test(x=x, y=y)$statistic, digits=3)
					p.val[i] <- signif(t.test(x=x, y=y)$p.value, digits=3)
				}
			} else if (markers.metric == "ROC") {
				bin.class <- ifelse(bin.class == 1, 0, 1)
				metric <- vector(length=n.models, mode="numeric")
				p.val <- vector(length=n.models, mode="numeric")
				for (i in 1:n.models) {
					m.score <- m2[i,]
					m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
					if (length(table(bin.class)) > 1) {
						perf.auc <- roc.area(bin.class, m.score.norm)
						metric[i] <- signif(perf.auc$A, digits=3)
						p.val[i] <- signif(perf.auc$p.value, digits=3)
					} else {
						metric[i] <- 1
						p.val[i] <- 1
					}
				}
			} else if (markers.metric == "MEAN.DIFF") {
				bin.class <- ifelse(bin.class == 1, 0, 1)
				metric <- vector(length=n.models, mode="numeric")
				p.val <- vector(length=n.models, mode="numeric")
				for (i in 1:n.models) {
					temp <- split(m2[i, ], bin.class)
					x <- temp[[1]]
					y <- temp[[2]]
					metric[i] <- signif(mean(x) - mean(y), digits=3)
					p.val[i] <- signif(t.test(x=x, y=y)$p.value, digits=3)
				}
			}
			
			if (is.na(sort.decreasing.genes)) {
				metric.order <- seq(1, length(metric))
			} else {
				metric.order <- order(metric, decreasing=sort.decreasing.genes)
			}
			if (only.up == TRUE) {
				k.markers.num <- ifelse(markers.num > n.models, n.models, markers.num)
				
#            if (length(classes) == 2) {
#               k.markers.num <- ifelse(markers.num > n.models, n.models, markers.num)
#            } else {
#               k.markers.num <- ifelse(length(classes)*markers.num > n.models, 
#                                               floor(n.models/length(classes)), markers.num)
#            }
				markers <- c(markers, model.names[metric.order][1:k.markers.num])
				markers.descs <- c(markers.descs, model.descs[metric.order][1:k.markers.num])
				metric.list <- c(metric.list, metric[metric.order][1:k.markers.num])
				p.val.list <- c(p.val.list, p.val[metric.order][1:k.markers.num])
				k.class <- c(k.class, rep(k, k.markers.num))
			} else {
				k.markers.num <- ifelse(length(classes)*markers.num > n.models, floor(n.models/length(classes)), 
						markers.num)
				markers <- c(markers, model.names[metric.order][1:k.markers.num],
						model.names[metric.order][(length(model.names) - k.markers.num +1):length(model.names)])
				markers.descs <- c(markers.descs, model.descs[metric.order][1:k.markers.num],
						model.descs[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)])
				metric.list <- c(metric.list, metric[metric.order][1:k.markers.num],
						metric[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)])
				p.val.list <- c(p.val.list, p.val[metric.order][1:k.markers.num],
						p.val[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)])
				k.class <- c(k.class, rep(k, k.markers.num), rep(paste("not", k), k.markers.num))
			}
		}
		
		V3 <- m2[markers,]
		print(V3)
		print(markers)
		
		if (show.desc == T) {
			model.descs2 <- paste(metric.list, p.val.list, k.class, markers.descs)
		} else {
			model.descs2 <- paste(metric.list, p.val.list)
		}
		height <- ifelse(length(markers) + n.phen >= 9, 10, (length(markers) + n.phen)*0.44 + 5)
#      char.res <-  0.0085 * length(markers) + 0.65
		
		
		# Sort markers inside each phenotype class
		
		if(sort.expression == T) {
			for (j in unique(k.phen.labels)) {
				V4 <- V3[ , k.phen.labels == j]
				sn <- sample.names2[k.phen.labels == j]
				if (is.vector(cls.labels)) {
					clab <- cls.labels2[k.phen.labels == j]
					clis <- cls.list2[k.phen.labels == j]
				} else {
					clab <- cls.labels2[, k.phen.labels == j]
					clis <- cls.list2[, k.phen.labels == j]
				}
				l.phen <- sum(k.phen.labels == j)
				if (l.phen > 1) {
					dist.matrix <- dist(t(V4))
					HC <- hclust(dist.matrix, method="complete")
					HC.order <- HC$order
					V4 <- V4[ , HC.order]
					sn <- sn[HC.order]
					if (is.vector(cls.labels2)) {
						clab <- clab[HC.order]
						clis <- clis[HC.order]
					} else {
						clab <- clab[, HC.order]
						clis <- clis[, HC.order]
					}
				}
				V3[ , k.phen.labels == j] <- V4
				sample.names2[k.phen.labels == j] <- sn
				if (is.vector(cls.labels2)) {
					cls.labels2[k.phen.labels == j] <- clab
					cls.list2[k.phen.labels == j] <- clis
				} else {
					cls.labels2[, k.phen.labels == j] <- clab
					cls.list2[, k.phen.labels == j] <- clis
				}
			}
		}
		x <- rbind(sample.names2, cls.list2, cls.labels2)
		print(paste("inside loop after in-class sort phen=", k.phen))
		print(x)
		print(cls.phen)
		print(phen.names)
		
		# Recompute cls.phen and cls.labels2 as order may have changed
		
		cls.phen2 <- list(NULL)
		if (is.vector(cls.labels2)) {
			classes <- unique(cls.list2)
			cls.phen2 <- classes
			cls.labels2 <- match(cls.list2, cls.phen2)
		} else {
			for (kk in 1:length(cls.list2[, 1])) {
				classes <- unique(cls.list2[kk,])
				cls.phen2[[kk]] <- classes
				cls.labels2[kk,] <- match(cls.list2[kk,], cls.phen2[[kk]])
			}
		}
		
		x <- rbind(sample.names2, cls.list2, cls.labels2)
		print(paste("inside loop after cls.phen renorm phen=", k.phen))
		print(cls.phen2)
		print(phen.names)
		
		
		library(gmodels)
		if (!is.vector(cls.labels2)) {
			if (sort.phenotypes == T) {
				phen.score <- vector(length=n.phen, mode="numeric")
				for (k.lab in 1:n.phen) {
					tab <- table(as.vector(cls.list2[k.lab,]), k.phen.list)
					print(tab)
#              phen.score[k.lab] <- 1 - chisq.test(tab)$p.value
#              phen.score[k.lab] <- 1 - fisher.test(tab)$p.value
					if ((length(tab[,1]) > 1) && (length(tab[1,]) > 1)) { 
						CT <- CrossTable(tab, chisq=T)
						phen.score[k.lab] <- CT$chisq$p.value
						print(phen.score[k.lab])
					} else {
						phen.score[k.lab] <- 0.50
						print(phen.score[k.lab])
					}
				}
				phen.order <- order(phen.score, decreasing= T)
				print(phen.order)
				cls.labels2 <- cls.labels2[phen.order,]
				cls.phen2 <- cls.phen2[phen.order]
				phen.names2 <- phen.names[phen.order]
				main.string <- paste(test.file.prefix, " - ", phen.names2[n.phen], markers.metric, " order")
			} else {
				phen.names2 <- phen.names
				main.string <- paste(test.file.prefix, " - ", phen.names2[k.phen], markers.metric, " order")
			}
		} else {
			phen.names2 <- phen.names[1]
			main.string <- paste(test.file.prefix, " - ", phen.names2, markers.metric, " order")
		}
		
#     windows(width=15, height=height)
		
		
		x <- rbind(sample.names2, cls.list2, cls.labels2)
		print(paste("inside loop after phen sort before figure phen=", k.phen))
		print(x)
		print(cls.phen2)
		print(phen.names2)
		
		phen.list <- unlist(cls.phen2)
		colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]
		
		print(rbind(phen.list, colors.list))
		
		if (show.desc == T) {
			markers <- paste(markers, seq(1, length(markers)), sep="_")
		}
		
		MSIG.HeatMapPlot.7(V = V3, row.names = markers,
				row.names2 = model.descs2, col.labels = cls.labels2, 
				col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names2,
				col.names = sample.names2, main = main.string, xlab="  ", ylab="  ", 
				row.norm = row.norm,  
				cmap.type = cmap.type, char.rescale = char.res,  legend=legend)
		
		V3 <- data.frame(V3)
		colnames(V3) <- sample.names2
		row.names(V3) <- markers
		
		if (!is.null(markers.file)) {
			write.gct(gct.data.frame = V3, descs = model.descs2, filename = markers.file)  
		}
		
	} # end loop over phenotypes
	
	dev.off()
	
}

OPAM.sort.projection.by.score.2 <- function(
		input.ds,
		input.cls,
		results.dir,
		normalize.score = T,
		normalization.type = "zero.one",
		model,
		target.phen = NA,
		target.class = NA,
		user.colors = NA,
		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T)
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	model.names <- dataset$row.names
	model.descs <- dataset$descs
	Ns <- length(m[1,])
	dim(m)
	sample.names <- dataset$names
	
	n.models <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
	char.res <-  0.013 * n.models + 0.65
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list 
	
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors[1])) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
	
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
	} else {
		phen.names <- "NA"
	}
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		for (i in 1:n.phen) {
			n.classes[i] <- length(cls.phen[[i]])
		}
	}
	
	filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
	pdf(file=paste(filename, ".pdf", sep=""), height = 8.5, width = 11)
#   windows(width=12, height=8)
	
	loc <- match(model, model.names)
	print(c("loc:", loc))
	s.order <- order(m[loc,], decreasing = decreasing.order)
	m2 <- m[, s.order]
	
	sample.names2 <- sample.names[s.order]
	
	if (is.vector(cls.labels)) {
		cls.labels2 <- cls.labels[s.order]
		cls.list2 <- cls.list[s.order]
	} else {
		cls.labels2 <- cls.labels[, s.order]
		cls.list2 <- cls.list[, s.order]
	}
	# Recompute cls.phen and cls.labels2 as order may have changed
	
	cls.phen2 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(cls.list2)
		cls.phen2 <- classes
		cls.labels2 <- match(cls.list2, cls.phen2)
	} else {
		for (kk in 1:length(cls.list2[, 1])) {
			classes <- unique(cls.list2[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2 <- c(cls.phen2, classes)
			cls.labels2[kk,] <- match(cls.list2[kk,], classes)
		}
	}
	
	
	correl <- cor(t(m2))[, loc]
	m.order <- order(correl, decreasing=T)
	correl2 <- correl[m.order]
	m2 <- m2[m.order,]
	model.names2 <- model.names[m.order]
	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list <- unlist(cls.phen2)
	colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]
	
	if (!is.na(target.phen)) {
		bin.class <- ifelse(cls.list2[target.phen,] == target.class, 1, 0)
	} else {
		bin.class <- ifelse(cls.list2[1,] == cls.list2[1,1], 1, 0)
	}
	
	sample.names2 <- paste(sample.names, bin.class, sep="_")
	
	print(bin.class)
	print(paste("n models:", n.models))
	for (i in 1:n.models) {
		print(paste(i, model.names2[i]))
		m.score <- m2[i,]
		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
		print(m.score.norm)
		if (length(unique(bin.class)) > 1) {
			perf.auc <- roc.area(bin.class, m.score.norm)
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
		} else {
			roc <- p.val <- "-"
		}
		print(paste("ROC=", roc, " p-val=", p.val))
		model.descs2[i] <- paste(roc, " (", p.val, ")")
	}
	
	MSIG.HeatMapPlot.7(V = m2, row.names = model.names2,
			row.names2 = model.descs2, col.labels = cls.labels2, 
			col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names,
			col.names = sample.names2, main = " ", xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=T)
	
	dev.off()
	
	if (!is.na(output.dataset)) {
		V.GCT <- data.frame(m2)
		colnames(V.GCT) <- sample.names2
		row.names(V.GCT) <- model.names2
		write.gct.2(gct.data.frame = V.GCT, descs = model.descs2, filename =output.dataset)  
	}
	
}

OPAM.sort.projection.by.score.3 <- function(
		input.ds,
		input.cls,
		results.dir,
		normalize.score = T,
		normalization.type = "zero.one",
		target.phen = NA,
		target.class = NA,
		user.colors = NA,
		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T)
# Calls MSIG.HeatMapPlot.8 and makes a plot sorted by the highest-scoring
# signatures and abnormalities (gene mutations or copy number alterations)
# i.e. doesn't require a "model" to score by as OPAM.sort.projection.by.score.2 does.
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	model.names <- dataset$row.names
#	model.descs <- dataset$descs
	Ns <- length(m[1,])
	dim(m)
	sample.names <- dataset$names
	
	n.models <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
	char.res <-  0.013 * n.models + 0.65
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list 
	
	
#	browser()
	
	
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors[1])) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
	
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
	} else {
		phen.names <- "NA"
	}
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
	print("cls.phen.colors:")
	print(cls.phen.colors)
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		for (i in 1:n.phen) {
			n.classes[i] <- length(cls.phen[[i]])
		}
	}
	
	filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
	pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 17 )
#   windows(width=12, height=8)
	
	roc.list = vector( length=n.models, mode="numeric" )
	p.val.list = vector( length=n.models, mode="numeric" )
	
	if (!is.na(target.phen)) {
		bin.class <- ifelse(cls.list[target.phen,] == target.class, 1, 0)
	} else {
		bin.class <- ifelse(cls.list[1,] == cls.list2[1,1], 1, 0)
	}
#	browser()
	model.descs2 = vector(length = n.models, mode="character")
	for (i in 1:n.models) {
		m.score <- m[i,]
		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
#		browser()
		if (length(unique(bin.class)) > 1) {
			perf.auc <- roc.area(bin.class, m.score.norm)
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
		} else {
			roc <- p.val <- "-"
		}
		print(paste("ROC=", roc, " p-val=", p.val)) 
		roc.list[i] = roc
		p.val.list[i] = p.val
		model.descs2[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
	m.order = order(roc.list, decreasing=TRUE)
	model.descs2 = model.descs2[m.order]
	loc = m.order[1]
	m2 <- m[m.order, ]
	model.names <- model.names[m.order]
	
#	loc <- match(model, model.names)
	print(c("loc:", loc))
	s.order <- order(m[loc,], decreasing = TRUE)
	m2 <- m2[, s.order]
	
	sample.names2 <- sample.names[s.order]
	
	
	
#	if (is.vector(cls.labels)) {
#		cls.labels2 <- cls.labels[s.order]
#		cls.list2 <- cls.list[s.order]
#	} else {
#		cls.labels2 <- cls.labels[, s.order]
#		cls.list2 <- cls.list[, s.order]
#	}
	# Recompute cls.phen and cls.labels2 as order may have changed
	
#	cls.phen2 <- list(NULL)
#	if (is.vector(cls.labels)) {
#		classes <- unique(cls.list2)
#		cls.phen2 <- classes
#		cls.labels2 <- match(cls.list2, cls.phen2)
#	} else {
#		for (kk in 1:length(cls.list2[, 1])) {
#			classes <- unique(cls.list2[kk,])
#			cls.phen2[[kk]] <- classes
#			cls.labels2[kk,] <- match(cls.list2[kk,], cls.phen2[[kk]])
#		}
#	}
#	browser()
	
	
	if (is.vector(cls.labels)) {
		cls.labels2 <- cls.labels[s.order]
		cls.list2 <- cls.list[s.order]
	} else {
		cls.labels2 <- cls.labels[, s.order]
		cls.list2 <- cls.list[, s.order]
	}
	#browser()
	
#	browser()
	m.score <- m2[1,]
	m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
	roc.list = vector(mode="numeric", length=n.phen)
	phen.descs = vector(mode="character", length=n.phen)
	for( i in 1:n.phen ){ 
		bin.gene = ifelse( cls.list2[i,]=="WT", 0, 1)
		if (length(unique(bin.gene)) > 1) {
			perf.auc <- roc.area(bin.gene, m.score.norm)
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
		} else {
			roc <- "-"
			p.val <- "-"
		}
		print(paste("ROC=", roc, " p-val=", p.val)) 
		roc.list[i] = roc
#		p.val.list[i] = p.val
		phen.descs[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
	g.order = c(1, 2, order(roc.list[3:n.phen], decreasing=TRUE)+2)  # skip PATHWAY.MUT and COPY.NUMBER
	phen.descs2 = phen.descs[g.order][1:40]
	cls.list2= cls.list2[g.order,][1:40,]
	phen.names = phen.names[g.order][1:40]
	
	# Recompute cls.list2 as some mutations or copy numbers may have been removed
	
	
	# Recompute cls.phen and cls.labels2 as order may have changed
	
	cls.phen2 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(cls.list2)
		cls.phen2 <- classes
		cls.labels2 <- match(cls.list2, cls.phen2)
	} else {
		for (kk in 1:length(cls.list2[, 1])) {
			classes <- unique(cls.list2[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2 <- c(cls.phen2, classes)
			cls.labels2[kk,] <- match(cls.list2[kk,], classes)
		}
	}
	cls.labels2 = cls.labels2[1:40,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list <- unlist(cls.phen2)
	
#	colors.list <- ifelse(unlist(cls.phen2) == target.class, 
#			ifelse(unlist(cls.phen2) == "DEL" | unlist(cls.phen2) == "AMP", 
#					ifelse(unlist(cls.phen2) == "DEL", cls.phen.colors[3], cls.phen.colors[4]), cls.phen.colors[1]), cls.phen.colors[2])
	colors.list = rep( "gray", length(phen.list))
	colors.list[phen.list=="MUT"] = cls.phen.colors[1]
	colors.list[phen.list=="DEL"] = cls.phen.colors[3]
	colors.list[phen.list=="AMP"] = cls.phen.colors[4]
	colors.list[phen.list=="ALT"] = cls.phen.colors[5]
	
#	roc.list = vector( length=n.models, mode="numeric" )
#	p.val.list = vector( length=n.models, mode="numeric" )
#	
#	if (!is.na(target.phen)) {
#		bin.class <- ifelse(cls.list2[target.phen,] == target.class, 1, 0)
#	} else {
#		bin.class <- ifelse(cls.list2[1,] == cls.list2[1,1], 1, 0)
#	}
	##	browser()
#	for (i in 1:n.models) {
#		m.score <- m2[i,]
#		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
	##	  browser()
#		if (length(unique(bin.class)) > 1) {
#			perf.auc <- roc.area(bin.class, m.score.norm)
#			roc <- signif(perf.auc$A, digits=3)
#			p.val <- signif(perf.auc$p.value, digits=3)
#		} else {
#			roc <- NA
#			p.val <- NA
#		}
#		print(paste("ROC=", roc, " p-val=", p.val)) 
#		roc.list[i] = roc
#		p.val.list[i] = p.val
#		model.descs2[i] <- paste(roc, " (", p.val, ")")
#	}
#	m.order = order(roc.list, decreasing=TRUE)
#	loc = which(m.order==1)
	##	correl <- cor(t(m2))[, loc]
#	print(c("loc:", loc))
#	s.order <- order(m[loc,], decreasing = decreasing.order)
#	m2 <- m[, s.order]
#	cls.phen2 = cls.phen2[s.order]
#	cls.labels2 = cls.labels2[s.order]
#
#	correl2 <- correl[m.order]
#	m2 = m2[m.order,]
#	model.names2 = model.names2[m.order]
#	model.descs2 = model.descs2[m.order]
	
	print("cls.phen2:")
	print(unlist(cls.phen2))
	
	print("cls.phen:")
	print(unlist(cls.phen))
	
	print("colors.list:")
	print(colors.list)
	
#	browser()
	
	MSIG.HeatMapPlot.8(V = m2, row.names = model.names,
			row.names2 = model.descs2, 
			col.labels = cls.labels2, 
			col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names,
			phen.names2 = phen.descs2,
			col.names = sample.names2, main = " ", xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	dev.off()
	
	if (!is.na(output.dataset)) {
		V.GCT <- m2
		colnames(V.GCT) <- sample.names2
		row.names(V.GCT) <- model.names2
		write.gct(gct.data.frame = V.GCT, descs = model.descs2, filename =output.dataset)  
	}
	
}

OPAM.sort.projection.by.score.4 <- function(
		input.ds,
		input.cls,
		tissue = "NA",
		results.dir,
		normalize.score = T,
		normalization.type = "zero.one",
		model = "NA",
		target.phen = NA,
		target.class = NA,
		user.colors = NA,
		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T,
		u.gene.names.known = "NA"
)
# Calls MSIG.HeatMapPlot.8 and makes a plot sorted by the highest-scoring
# signatures and abnormalities (gene mutations or copy number alterations)
# i.e. doesn't require a "model" to score by as OPAM.sort.projection.by.score.2 does.
# However, it *will* use "model" if it cannot calculate p-values on the gene signatures, which
# happens when every cell line has a genomic aberration.
#
# Runs 3 passes on the data:
# 1st pass: looks at the genes and copy number alterations specified by u.gene.names.known
# 2nd pass: looks at only the top abnormalities (using a p-value cutoff) from the 1st pass, and adjusts 
# the PATHWAY.MUT vector accordingly (only according to the genes, not by copy number data)
# 3rd pass: Takes the winning signature from the 2nd pass and then looks all the genes available
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	model.names <- dataset$row.names
#	model.descs <- dataset$descs
	Ns <- length(m[1,])
	dim(m)
	sample.names <- dataset$names
	
	n.models <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
	char.res <-  0.013 * n.models + 0.65
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list 
	
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors[1])) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
	
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
	} else {
		phen.names <- "NA"
	}
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
#	print("cls.phen.colors:")
#	print(cls.phen.colors)
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		for (i in 1:n.phen) {
			n.classes[i] <- length(cls.phen[[i]])
		}
	}
	print("--- Begin Pass 1 ---")
#	model.names.original = model.names
#	m.original = m
	phen.pass1 = c( "PATHWAY.MUT", u.gene.names.known)
	n.phen.pass1 = length(u.gene.names.known)+1
	ind.phen.pass1 = which( phen.names %in% phen.pass1 )
	phen.pass1 = phen.names[ind.phen.pass1]
	
	roc.list.pass1 = vector( length=n.models, mode="numeric" )
	p.val.list.pass1 = vector( length=n.models, mode="numeric" )
	
	cls.list.pass1 = cls.list[ind.phen.pass1,]
	cls.list.pass1.2 = ifelse(cls.list.pass1 == "WT", 0, 1)
	cls.labels.pass1 = cls.labels[ind.phen.pass1,]
#	browser()
	if (!is.na(target.phen) && length(phen.pass1) > 2 ) {
#		bin.class.pass1 = apply( cls.list2.pass1.2[3:n.phen,], MARGIN=2, FUN=sum)/(n.phen-2)
#		bin.class.pass1 = ( bin.class.pass1 - min(bin.class.pass1))/(max(bin.class.pass1) - min(bin.class.pass1))
#		bin.class.pass1 <- ifelse(cls.list[target.phen,] == target.class, 1, 0)
		bin.class.pass1 <- ifelse(apply( cls.list.pass1.2[2:n.phen.pass1,], MARGIN=2, FUN=sum) > 0, 1, 0)
		cls.labels.pass1[1,] = bin.class.pass1
		cls.list.pass1[1,] = ifelse(bin.class.pass1 == 1, "MUT", "WT")
#		if( length(unique(bin.class.pass1)) == 1) { 
#			cls.list.3 = ifelse( cls.list == "DEL" | cls.list == "AMP", 1, 0)
#			copy.number.pass1 = ifelse( apply(cls.list.3[3:n.phen,], MARGIN=2, FUN=sum) >= 1, "ALT", "WT")
#			copy.class.pass1 = ifelse( copy.number.pass1 == "ALT", 1, 0)
#			bin.class.pass1 = copy.class.pass1
#			print( "Calculating p-value with respect to copy number alterations")
#		}
	} else if (length(phen.pass1)==2 ) {
		bin.class.pass1 = ifelse(cls.list[2,]== "WT", 0,1)
		cls.labels.pass1[1,] = bin.class.pass1
		cls.list.pass1[1,] = ifelse(bin.class.pass1 == 1, "MUT", "WT")
	} else {
		bin.class.pass1 <- ifelse(cls.list[1,] == cls.list2[1,1], 1, 0)
	}
#	browser()
#pdf("ROCplots.pdf")
	model.descs2.pass1 = vector(length = n.models, mode="character")
	for (i in 1:n.models) {
		m.score <- m[i,]
		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
#		browser()
		if (length(unique(bin.class.pass1)) > 1) {
			perf.auc <- roc.area(bin.class.pass1, m.score.norm)
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.pass1[i] = perf.auc$A
			p.val.list.pass1[i] = perf.auc$p.value
			roc.plot(bin.class.pass1, m.score.norm)
		} else {
			roc <- p.val <- "-"
			roc.list.pass1[i] = NA
			p.val.list.pass1[i] = NA
		}
		print(paste(format(rownames(m)[i], width=30), "ROC=", roc, " p-val=", p.val)) 
		
		model.descs2.pass1[i] <- paste(roc, " (", p.val, ")")
	}
	dev.off()
#	browser()
	if( is.na(roc.list.pass1[1]) ){
		loc <- match(model, model.names)
		s.order.pass1 <- order(m[loc,], decreasing = decreasing.order)
#		loc = s.order.pass1[1]
		m2.pass1 <- m[, s.order.pass1]
		correl <- cor(t(m2.pass1))[, loc]
		m.order.pass1 <- order(correl, decreasing=T)
		m2.pass1 <- m2.pass1[m.order.pass1, ]
	} else{ 
		m.order.pass1 = order(roc.list.pass1, decreasing=TRUE, na.last=TRUE)
		m2.pass1 <- m[m.order.pass1, ]
		s.order.pass1 <- order(m2.pass1[1,], decreasing = TRUE )
		m2.pass1 <- m2.pass1[, s.order.pass1]
	}
#	m2.pass1 <- m2.pass1[m.order.pass1, ]
	model.descs2.pass1 = model.descs2.pass1[m.order.pass1]
	
	
	sample.names2.pass1 <- colnames(m2.pass1)
	model.names.pass1 <- rownames(m2.pass1)
	
	if (is.vector(cls.labels)) {
		cls.labels2.pass1 <- cls.labels.pass1[s.order.pass1]
		cls.list2.pass1 <- cls.list.pass1[s.order.pass1]
	} else {
		cls.labels2.pass1 <- cls.labels.pass1[, s.order.pass1]
		cls.list2.pass1 <- cls.list.pass1[, s.order.pass1]
	}
	#browser()
	
#	browser()
	m.score.pass1 <- m2.pass1[1,]
	m.score.norm.pass1 <- (m.score.pass1 - min(m.score.pass1))/(max(m.score.pass1) - min(m.score.pass1))
	roc.list.phen.pass1 = vector(mode="numeric", length=n.phen)
	p.val.list.phen.pass1 = vector(mode="numeric", length=n.phen)
	phen.descs.pass1 = vector(mode="character", length=n.phen)
	for( i in 1:n.phen.pass1 ){ 
		bin.gene = ifelse( cls.list2.pass1[i,]=="WT", 0, 1)
		if (length(unique(bin.gene)) > 1) {
			perf.auc <- roc.area(bin.gene, m.score.norm.pass1)
			
#			if( perf.auc$A < 0.5 ){
			##				browser()
#				roc = signif(1 - perf.auc$A, digits=3)
#				p.val = signif(1- perf.auc$p.val, digits=3)
#				abnormality = unique(cls.list2.pass1[i,])[which(unique(cls.list2.pass1[i,]) != "WT")]
			##				cls.list2.pass1[i,] = ifelse( cls.list2.pass1[i,] == "WT", abnormality, "WT" )
#				phen.names[i] = paste(phen.names[i], "-opposite", sep="")
#				roc.list.phen.pass1[i] = 1 - perf.auc$A
#				p.val.list.phen.pass1[i] = perf.auc$p.val   # Don't want to use these "opposite" genomic aberrations in Pass 2 
#															# because they make PATHWAY.MUT+COPY.NUMBER too dense
#			} else{
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.phen.pass1[i] = perf.auc$A
			p.val.list.phen.pass1[i] = perf.auc$p.val
#			}
		} else {
			roc <- "-"
			p.val <- "-"
			roc.list.phen.pass1[i] = NA
			p.val.list.phen.pass1[i] = NA
		}
		print(paste(format(phen.pass1[i], width=12), "ROC=", roc, " p-val=", p.val)) 
		
#		p.val.list[i] = p.val
		phen.descs.pass1[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
	g.order.pass1 = c(1, order(roc.list.phen.pass1[2:n.phen.pass1], decreasing=TRUE, na.last=TRUE)+1)  # keep PATHWAY.MUT as 1
	roc.list.phen.pass1 = roc.list.phen.pass1[g.order.pass1]
	p.val.list.phen.pass1 = p.val.list.phen.pass1[g.order.pass1] 
	phen.descs2.pass1 = phen.descs.pass1[g.order.pass1][1:n.phen.pass1]
	cls.list2.pass1 = cls.list2.pass1[g.order.pass1,][1:n.phen.pass1,]
	phen.names.pass1 = phen.pass1[g.order.pass1][1:n.phen.pass1]
	
	# Recompute cls.list2 as some mutations or copy numbers may have been removed
	
	
	# Recompute cls.phen and cls.labels2 as order may have changed
	
	cls.phen2.pass1 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(cls.list2.pass1)
		cls.phen2.pass1 <- classes
		cls.labels2.pass1 <- match(cls.list2.pass1, cls.phen2.pass1)
	} else {
		for (kk in 1:length(cls.list2.pass1[, 1])) {
			classes <- unique(cls.list2.pass1[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2.pass1 <- c(cls.phen2.pass1, classes)
			cls.labels2.pass1[kk,] <- match(cls.list2.pass1[kk,], classes)
		}
	}
	cls.labels2.pass1 = cls.labels2.pass1[1:n.phen.pass1,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list.pass1 <- unlist(cls.phen2.pass1)
	
#	colors.list <- ifelse(unlist(cls.phen2) == target.class, 
#			ifelse(unlist(cls.phen2) == "DEL" | unlist(cls.phen2) == "AMP", 
#					ifelse(unlist(cls.phen2) == "DEL", cls.phen.colors[3], cls.phen.colors[4]), cls.phen.colors[1]), cls.phen.colors[2])
	colors.list.pass1 = rep( "gray", length(phen.list.pass1))
	colors.list.pass1[phen.list.pass1=="MUT"] = cls.phen.colors[1]
	colors.list.pass1[phen.list.pass1=="DEL"] = cls.phen.colors[3]
	colors.list.pass1[phen.list.pass1=="AMP"] = cls.phen.colors[4]
	colors.list.pass1[phen.list.pass1=="ALT"] = cls.phen.colors[5]
	
#	print("cls.phen2:")
#	print(unlist(cls.phen2))
#	
#	print("cls.phen:")
#	print(unlist(cls.phen))
#	
#	print("colors.list:")
#	print(colors.list)
	
#	browser()
	
	filename <- paste(results.dir, test.file.prefix, ".3Passes.ROC", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
	pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 17 )
#   windows(width=12, height=8)
	MSIG.HeatMapPlot.9(V = m2.pass1, row.names = model.names.pass1,
			row.names2 = model.descs2.pass1, 
			col.labels = cls.labels2.pass1, 
			col.classes = cls.phen2.pass1, 
			phen.cmap = colors.list.pass1, phen.names = phen.names.pass1,
			phen.names2 = phen.descs2.pass1,
			col.names = sample.names2.pass1, main = paste(tissue, "- Phase 1: Known KRAS Pathway Abnormalities (ROC)"), 
			xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
#	dev.off()
#	break
	
	### Begin Pass 2 ###
	print( "--- Begin Pass 2 ---")
#	browser()
#	p.val.threshold = 0.1
	roc.threshold = 0.65
#	ind.top.pval = which(p.val.list.phen.pass1[2:n.phen.pass1] <= p.val.threshold )+1
	ind.top.roc = which(roc.list.phen.pass1[-1] >= roc.threshold)+1
	if( length(ind.top.roc) > 0){
		ind.roc.threshold = c(1, ind.top.roc)
	}else{
		roc.threshold = 0.6
		ind.top.roc = which(roc.list.phen.pass1[-1] >= roc.threshold)+1
		ind.roc.threshold = c(1, ind.top.roc)
	}
	
#	if( length(ind.top.pval) > 0 ){
#		ind.p.val.threshold = c(1, ind.top.pval) 
#	} else if( length(ind.top.pval) < 1 ) { 
#		p.val.threshold = 0.15
#		ind.top.pval = which(p.val.list.phen.pass1[2:n.phen.pass1] <= p.val.threshold )+1
#		ind.p.val.threshold = c(1, ind.top.pval) }
#	if( length(ind.top.pval) < 1){ 
#		p.val.threshold = 0.2
#		ind.top.pval = which(p.val.list.phen.pass1[2:n.phen.pass1] <= p.val.threshold )+1
#		ind.p.val.threshold = c(1, ind.top.pval)
#	}
#	if( length(ind.top.pval) < 1){ 
#		p.val.threshold = 0.25
#		ind.top.pval = which(p.val.list.phen.pass1[2:n.phen.pass1] <= p.val.threshold )+1
#		ind.p.val.threshold = c(1, ind.top.pval)
#	}
#	if( length(ind.top.pval) < 1){ 
#		p.val.threshold = 0.3
#		ind.top.pval = which(p.val.list.phen.pass1[2:n.phen.pass1] <= p.val.threshold )+1
#		ind.p.val.threshold = c(1, ind.top.pval)
#	}
#	if( length( ind.top.pval) < 1 ) {
#		ind.top = which(!is.na(p.val.list.phen.pass1[-1]))+1 
#		ind.p.val.threshold = c( 1, ind.top )
#	}
	
	
	n.phen.pass2 = length(ind.roc.threshold)
#	phen.descs.pass2 = phen.descs2.pass1[1:n.phen.pass2]
	cls.list2.pass2 = cls.list2.pass1[ind.roc.threshold,]
	phen.names.pass2 = phen.names.pass1[ind.roc.threshold]
#	phen.names.pass2[1] = "PATHWAY.MUT + COPY.NUMBER"
	cls.labels2.pass2 = cls.labels2.pass1[ind.roc.threshold,]
	
#	browser()
	cls.list2.pass2.2 = ifelse( cls.list2.pass2 == "WT", 0, 1)
	cls.list2.pass2.3 = ifelse( cls.list2.pass2 == "DEL" | cls.list2.pass2 == "AMP", 1, 0)
	if( n.phen.pass2 > 2 ){ 
		pathway.mut.pass2 = ifelse( apply(cls.list2.pass2.2[2:n.phen.pass2,], MARGIN=2, FUN=sum) >= 1, "MUT", "WT")
#		copy.number.pass2 = ifelse( apply(cls.list2.pass2.3[3:n.phen.pass2,], MARGIN=2, FUN=sum) >= 1, "ALT", "WT")
	} else{
		pathway.mut.pass2 = ifelse( cls.list2.pass2.2[2,] == 1, "MUT", "WT")
#		copy.number.pass2 = ifelse( cls.list2.pass2.3[3,] == 1, "ALT", "WT")
	}
#	browser()
	bin.class.pass2 = ifelse( pathway.mut.pass2 == "MUT", 1, 0 )
#	copy.class.pass2 = ifelse( copy.number.pass2 == "ALT", 1, 0)
#	if( length(unique(bin.class.pass2)) == 1) { 
#		bin.class.pass2 = copy.class.pass2
#		print( "Calculating p-value with respect to copy number alterations")
#	}
	cls.list2.pass2[1,] = pathway.mut.pass2
#	cls.list2.pass2[2,] = copy.number.pass2
	
	roc.list.pass2 = vector( length=n.models, mode="numeric" )
	p.val.list.pass2 = vector( length=n.models, mode="numeric" )
	
#	browser()
	model.descs2.pass2 = vector(length = n.models, mode="character")
	for (i in 1:n.models) {
		m.score <- m2.pass1[i,]
		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
#		browser()
		if (length(unique(bin.class.pass2)) > 1) {
			perf.auc <- roc.area(bin.class.pass2, m.score.norm)
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.pass2[i] = perf.auc$A
			p.val.list.pass2[i] = perf.auc$p.value
		} else {
			roc <- p.val <- "-"
			roc.list.pass2[i] = NA
			p.val.list.pass2[i] = NA
		}
		print(paste(format(rownames(m2.pass1)[i], width=30), "ROC=", roc, " p-val=", p.val)) 
		
		model.descs2.pass2[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
	m.order.pass2 = order(roc.list.pass2, decreasing=TRUE, na.last=TRUE)
	model.descs2.pass2 = model.descs2.pass2[m.order.pass2]
#	loc.pass2 = m.order.pass2[1]
	m2.pass2 <- m2.pass1[m.order.pass2, ]
	model.names.pass2 <- rownames(m2.pass2)
#	print(c("loc.pass2:", loc.pass2))
	s.order.pass2 <- order(m2.pass2[1,], decreasing = TRUE)
	m2.pass2 <- m2.pass2[, s.order.pass2]
	sample.names2.pass2 <- colnames(m2.pass2)
	
	if (is.vector(cls.labels)) {
		cls.labels2.pass2 <- cls.labels2.pass2[s.order.pass2]
		cls.list2.pass2 <- cls.list2.pass2[s.order.pass2]
	} else {
		cls.labels2.pass2 <- cls.labels2.pass2[, s.order.pass2]
		cls.list2.pass2 <- cls.list2.pass2[, s.order.pass2]
	}
#	browser()
	
#	browser()
	m.score.pass2 <- m2.pass2[1,]
	m.score.norm.pass2 <- (m.score.pass2 - min(m.score.pass2))/(max(m.score.pass2) - min(m.score.pass2))
	roc.list.phen.pass2 = vector(mode="numeric", length=n.phen.pass2)
	phen.descs.pass2 = vector(mode="character", length=n.phen.pass2)
	for( i in 1:n.phen.pass2 ){ 
		bin.gene = ifelse( cls.list2.pass2[i,]=="WT", 0, 1)
		if (length(unique(bin.gene)) > 1) {
			perf.auc <- roc.area(bin.gene, m.score.norm.pass2)
#			if( perf.auc$A < 0.5 ){
			##				browser()
#				roc = signif(1 - perf.auc$A, digits=3)
#				p.val = signif(1 - perf.auc$A, digits=3)
#				abnormality = unique(cls.list2.pass2[i,])[which(unique(cls.list2.pass2[i,]) != "WT")]
#				cls.list2.pass2 = ifelse( cls.list2.pass2[i,] == "WT", abnormality, "WT" )
#				phen.names.pass2[i] = paste(phen.names.pass2[i], "-opposite")
#			} else{
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.phen.pass2[i] = perf.auc$A
#			}
		} else {
			roc <- "-"
			p.val <- "-"
			roc.list.phen.pass2[i] = NA
		}
		print(paste(format(phen.names.pass2[i], width=12), "ROC=", roc, " p-val=", p.val)) 
#		p.val.list[i] = p.val
		phen.descs.pass2[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
	g.order.pass2 = c(1, order(roc.list.phen.pass2[2:n.phen.pass2], decreasing=TRUE, na.last=TRUE)+1)  # skip PATHWAY.MUT
	phen.descs2.pass2 = phen.descs.pass2[g.order.pass2]
	cls.list2.pass2 = cls.list2.pass2[g.order.pass2,]
	phen.names.pass2 = phen.names.pass2[g.order.pass2]
#	browser()
	# Recompute cls.list2 as some mutations or copy numbers may have been removed
	
	
	# Recompute cls.phen and cls.labels2 as order may have changed
	
	cls.phen2.pass2 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(as.vector(cls.list2.pass2))
		cls.phen2.pass2 <- classes
		cls.labels2.pass2 <- match(cls.list2.pass2, cls.phen2.pass2)
	} else {
		for (kk in 1:length(cls.list2.pass2[, 1])) {
			classes <- unique(cls.list2.pass2[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2.pass2 <- c(cls.phen2.pass2, classes)
			cls.labels2.pass2[kk,] <- match(cls.list2.pass2[kk,], classes)
		}
	}
	cls.labels2.pass2 = cls.labels2.pass2[1:n.phen.pass2,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list.pass2 <- unlist(cls.phen2.pass2)
	
#	colors.list <- ifelse(unlist(cls.phen2) == target.class, 
#			ifelse(unlist(cls.phen2) == "DEL" | unlist(cls.phen2) == "AMP", 
#					ifelse(unlist(cls.phen2) == "DEL", cls.phen.colors[3], cls.phen.colors[4]), cls.phen.colors[1]), cls.phen.colors[2])
	colors.list.pass2 = rep( "gray", length(phen.list.pass2))
	colors.list.pass2[phen.list.pass2=="MUT"] = cls.phen.colors[1]
	colors.list.pass2[phen.list.pass2=="DEL"] = cls.phen.colors[3]
	colors.list.pass2[phen.list.pass2=="AMP"] = cls.phen.colors[4]
	colors.list.pass2[phen.list.pass2=="ALT"] = cls.phen.colors[5]
	
	MSIG.HeatMapPlot.9(V = m2.pass2, row.names = model.names.pass2,
			row.names2 = model.descs2.pass2, 
			col.labels = cls.labels2.pass2, 
			col.classes = cls.phen2.pass2, 
			phen.cmap = colors.list.pass2, phen.names = phen.names.pass2,
			phen.names2 = phen.descs2.pass2,
			col.names = sample.names2.pass2, main = paste(tissue, "- Phase 2: only ROC >=", roc.threshold,"from 1st pass (ROC)"), 
			xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	### 3rd Pass ###	
	
	print( "--- Begin Pass 3 ---")
#	browser()
	m2.pass3 = m2.pass2
	model.names.pass3 = rownames(m2.pass3)
	sample.names2.pass3 = colnames(m2.pass3)
#	model.descs2.pass3 = model.descs2.pass2
	n.phen.pass3 = 40
#	phen.descs.pass2 = phen.descs2.pass1[1:n.phen.pass2]
	cls.list2.pass3 = cls.list[, s.order.pass1][, s.order.pass2]
	cls.labels2.pass3 = cls.labels[, s.order.pass1][, s.order.pass2]
	
#	browser()
	phen.names.pass3 = phen.names
	m.score.pass3 <- m2.pass3[1,]
	m.score.norm.pass3 <- (m.score.pass3 - min(m.score.pass3))/(max(m.score.pass3) - min(m.score.pass3))
	roc.list.phen.pass3 = vector(mode="numeric", length=n.phen)
	phen.descs.pass3 = vector(mode="character", length=n.phen)
	p.val.list.phen.pass3 = vector(mode="numeric", length=n.phen)
	for( i in 1:n.phen ){ 
		bin.gene = ifelse( cls.list2.pass3[i,]=="WT", 0, 1)
		if (length(unique(bin.gene)) > 1) {
			perf.auc <- roc.area(bin.gene, m.score.norm.pass3)
#			if( perf.auc$A < 0.5 ){
			##				browser()
#				roc = signif(1 - perf.auc$A, digits=3)
#				p.val = signif(1- perf.auc$p.value, digits=3)
#				abnormality = unique(cls.list2.pass3[i,])[which(unique(cls.list2.pass3[i,]) != "WT")]
#				cls.list2.pass3[i,] = ifelse( cls.list2.pass3[i,] == "WT", abnormality, "WT" )
#				phen.names.pass3[i] = paste(phen.names.pass3[i], "-opposite", sep="")
#				roc.list.phen.pass3[i] = 1-perf.auc$A
#				p.val.list.phen.pass3[i] = 1- perf.auc$p.value
#			} else{
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.phen.pass3[i] = perf.auc$A
			p.val.list.phen.pass3[i] = perf.auc$p.value
#			}
		} else {
			roc <- NA
			p.val <- NA
			roc.list.phen.pass3[i] = NA
			p.val.list.phen.pass3[i] = NA
		}
#		print(paste("ROC=", roc, " p-val=", p.val)) 
		
#		p.val.list[i] = p.val
		phen.descs.pass3[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
#	p.val.threshold = 0.1
	roc.threshold = 0.65
	len = length(which(roc.list.phen.pass3[-1:-2] >= roc.threshold)) + 2
#	len = length(which(p.val.list.phen.pass3[-1:-2] <= p.val.threshold))+2
#	if( len == 2 ){
#		p.val.threshold = 0.15
#		len = length(which(p.val.list.phen.pass3[-1:-2] <= p.val.threshold))+2
#	}
#	if( len == 2 ){
#		p.val.threshold = 0.2
#		len = length(which(p.val.list.phen.pass3[-1:-2] <= p.val.threshold))+2
#	}
	if( len > 40 ) len = 40
#	g.order.pass3.1 = c(1, 2, order(p.val.list.phen.pass3[3:n.phen], decreasing=FALSE, na.last=TRUE)+2 )
	g.order.pass3 =  c(1, 2, order(p.val.list.phen.pass3[-1:-2], decreasing=FALSE, na.last=TRUE)+2 )[1:len]  # skip PATHWAY.MUT and COPY.NUMBER
	phen.descs2.pass3 = phen.descs.pass3[g.order.pass3]
	cls.list2.pass3 = cls.list2.pass3[g.order.pass3,]
	cls.labels2.pass3 = cls.labels2.pass3[g.order.pass3,]
	phen.names.pass3 = phen.names.pass3[g.order.pass3]
	
	
	cls.list.mut = ifelse(cls.list2.pass3[-1:-2,] == "MUT", 1, 0)
	cls.list.alt = ifelse(cls.list2.pass3[-1:-2,] == "DEL" | cls.list2.pass3[-1:-2,] == "AMP", 1, 0)
	
#	browser()
	if( !is.vector(cls.list.alt) ){
		cls.list.mut.sum = apply(cls.list.mut, MARGIN=2, FUN=sum)
		cls.list.alt.sum = apply(cls.list.alt, MARGIN=2, FUN=sum)
		cls.list.mut.sum = ifelse(cls.list.mut.sum + cls.list.alt.sum > 0, 1, 0)
		cls.list2.pass3[1,] = ifelse( cls.list.mut.sum >= 1, "MUT", "WT")
		cls.list2.pass3[2,] = ifelse( cls.list.alt.sum >= 1, "ALT", "WT")
		bin.class.pass3 = cls.list.mut.sum
	} else{
		cls.list2.pass3[1,] = ifelse(cls.list.mut == 1, "MUT", "WT")
		cls.list2.pass3[2,] = ifelse(cls.list.alt == 1, "ALT", "WT")
		bin.class.pass3 = cls.list.mut
	}
	
	for( i in 1:2 ){ # Recalculate ROC and p-value for PATHWAY.MUT and COPY.NUMBER
		bin.gene = ifelse( cls.list2.pass3[i,]=="WT", 0, 1)
		if (length(unique(bin.gene)) > 1) {
			perf.auc <- roc.area(bin.gene, m.score.norm.pass3)
#			if( perf.auc$A < 0.5 ){
			##				browser()
#				roc = signif(1 - perf.auc$A, digits=3)
#				p.val = signif(1- perf.auc$p.value, digits=3)
#				abnormality = unique(cls.list2.pass3[i,])[which(unique(cls.list2.pass3[i,]) != "WT")]
#				cls.list2.pass3[i,] = ifelse( cls.list2.pass3[i,] == "WT", abnormality, "WT" )
#				phen.names.pass3[i] = paste(phen.names.pass3[i], "-opposite", sep="")
#				roc.list.phen.pass3[i] = 1-perf.auc$A
#				p.val.list.phen.pass3[i] = 1- perf.auc$p.value
#			} else{
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.phen.pass3[i] = perf.auc$A
			p.val.list.phen.pass3[i] = perf.auc$p.value
#			}
		} else {
			roc <- NA
			p.val <- NA
			roc.list.phen.pass3[i] = NA
			p.val.list.phen.pass3[i] = NA
		}
		print(paste(format(phen.names.pass3[i], width=12), "ROC=", roc, " p-val=", p.val)) 
		
#		p.val.list[i] = p.val
		phen.descs2.pass3[i] <- paste(roc, " (", p.val, ")")
	}
	
#	browser()
	model.descs2.pass3 = vector(length = n.models, mode="character")
	roc.list.pass3 = vector(length=n.models, mode="double")
	p.val.list.pass3 = vector(length=n.models, mode="double")
	for (i in 1:n.models) {
		m.score <- m2.pass3[i,]
		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
#		browser()
		if (length(unique(bin.class.pass3)) > 1) {
			perf.auc <- roc.area(bin.class.pass3, m.score.norm)
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.pass3[i] = perf.auc$A
			p.val.list.pass3[i] = perf.auc$p.value
		} else {
			roc <- p.val <- "-"
			roc.list.pass3[i] = NA
			p.val.list.pass3[i] = NA
		}
		print(paste(format(rownames(m2.pass3)[i], width=30), "ROC=", roc, " p-val=", p.val)) 
		
		model.descs2.pass3[i] <- paste(roc, " (", p.val, ")")
	}
	m.order.pass3 = order(roc.list.pass3, decreasing=TRUE)
	m2.pass3 = m2.pass3[m.order.pass3,]
	model.descs2.pass3 = model.descs2.pass3[m.order.pass3]
	s.order.pass3 = order(m2.pass3[1,], decreasing=TRUE)
	m2.pass3 = m2.pass3[,s.order.pass3]
	sample.names2.pass3 = colnames(m2.pass3)
	model.names.pass3 = rownames(m2.pass3)
	
	cls.phen2.pass3 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(as.vector(cls.list2.pass3))
		cls.phen2.pass3 <- classes
		cls.labels2.pass3 <- match(cls.list2.pass3, cls.phen2.pass3)
	} else {
#		browser()
		for (kk in 1:length(cls.list2.pass3[, 1])) {
#			browser()
			classes <- unique(cls.list2.pass3[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2.pass3 <- c(cls.phen2.pass3, classes)
			cls.labels2.pass3[kk,] <- match(cls.list2.pass3[kk,], classes)
		}
	}
#	cls.labels2.pass3 = cls.labels2.pass3[1:n.phen.pass3,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list.pass3 <- unlist(cls.phen2.pass3)
	
#	colors.list <- ifelse(unlist(cls.phen2) == target.class, 
#			ifelse(unlist(cls.phen2) == "DEL" | unlist(cls.phen2) == "AMP", 
#					ifelse(unlist(cls.phen2) == "DEL", cls.phen.colors[3], cls.phen.colors[4]), cls.phen.colors[1]), cls.phen.colors[2])
	colors.list.pass3 = rep( "gray", length(phen.list.pass3))
	colors.list.pass3[phen.list.pass3=="MUT"] = cls.phen.colors[1]
	colors.list.pass3[phen.list.pass3=="DEL"] = cls.phen.colors[3]
	colors.list.pass3[phen.list.pass3=="AMP"] = cls.phen.colors[4]
	colors.list.pass3[phen.list.pass3=="ALT"] = cls.phen.colors[5]
	phen.names.pass3[1] = "PATHWAY.MUT+COPY.NUMBER"
#	browser()
	MSIG.HeatMapPlot.9(V = m2.pass3, row.names = model.names.pass3,
			row.names2 = model.descs2.pass3, 
			col.labels = cls.labels2.pass3, 
			col.classes = cls.phen2.pass3, 
			phen.cmap = colors.list.pass3, phen.names = phen.names.pass3,
			phen.names2 = phen.descs2.pass3,
			col.names = sample.names2.pass3, main = paste(tissue, "- 3rd Pass: Top signature from 2nd pass with all genes ( ROC >=", roc.threshold, ") (ROC)"), 
			xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	dev.off()
	
	if (!is.na(output.dataset)) {
		V.GCT <- m2
		colnames(V.GCT) <- sample.names2
		row.names(V.GCT) <- model.names2
		write.gct(gct.data.frame = V.GCT, descs = model.descs2, filename =output.dataset)  
	}
	
}

OPAM.sort.projection.by.score.5 <- function(
		input.ds,
		input.cls,
		tissue = "NA",
		results.dir,
		normalize.score = T,
		normalization.type = "zero.one",
		model = "NA",
		target.phen = NA,
		target.class = NA,
		user.colors = NA,
		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T,
		u.gene.names.known = "NA"
)
# Calls MSIG.HeatMapPlot.9 and makes a plot sorted by the highest-scoring
# signatures and abnormalities (gene mutations or copy number alterations)
# i.e. doesn't require a "model" to score by as OPAM.sort.projection.by.score.2 does.
# However, it *will* use "model" if it cannot calculate p-values on the gene signatures, which
# happens when every cell line has a genomic aberration.
#
# Runs 3 passes on the data:
# 1st pass: looks at the genes and copy number alterations specified by u.gene.names.known
# 2nd pass: looks at only the top abnormalities (using a p-value cutoff) from the 1st pass, and adjusts 
# the PATHWAY.MUT vector accordingly (only according to the genes, not by copy number data)
# 3rd pass: Takes the winning signature from the 2nd pass and then looks all the genes available
#
# Very similar to OPAM.sort.projection.by.score.4, however this version uses rec.area instead of
# roc.area to calculate REC/ROC scores and p-values for PATHWAY.MUT, the vector of total genomic aberrations
# in all samples
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	model.names <- dataset$row.names
#	model.descs <- dataset$descs
	Ns <- length(m[1,])
	dim(m)
	sample.names <- dataset$names
	
	n.models <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
	char.res <-  0.013 * n.models + 0.65
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list 
	
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors[1])) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
	
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
	} else {
		phen.names <- "NA"
	}
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
#	print("cls.phen.colors:")
#	print(cls.phen.colors)
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		for (i in 1:n.phen) {
			n.classes[i] <- length(cls.phen[[i]])
		}
	}
#	model.names.original = model.names
#	m.original = m
	phen.pass1 = c( "PATHWAY.MUT", u.gene.names.known)
	n.phen.pass1 = length(u.gene.names.known)+1
	ind.phen.pass1 = which( phen.names %in% phen.pass1 )
	phen.pass1 = phen.names[ind.phen.pass1]
	
	roc.list.pass1 = vector( length=n.models, mode="numeric" )
	p.val.list.pass1 = vector( length=n.models, mode="numeric" )
	
	cls.list.pass1 = cls.list[ind.phen.pass1,]
	cls.list.pass1.2 = ifelse(cls.list.pass1 == "WT", 0, 1)
	cls.labels.pass1 = cls.labels[ind.phen.pass1,]
#	browser()
	if (!is.na(target.phen)) {
		bin.class.pass1 = apply( cls.list.pass1.2[-1,], MARGIN=2, FUN=sum)
		# Normalize bin.class.pass1
		if( length(unique(bin.class.pass1)) > 1){
			bin.class.pass1 = ( bin.class.pass1 - min(bin.class.pass1))/(max(bin.class.pass1) - min(bin.class.pass1))
		} else if ( length(unique(bin.class.pass1)) > 1){
			bin.class = rep(1, length(cls.list[1,]))
		}
#		bin.class.pass1 <- ifelse(cls.list[target.phen,] == target.class, 1, 0)
#		bin.class.pass1 <- ifelse(apply( cls.list.pass1.2[2:n.phen.pass1,], MARGIN=2, FUN=sum) > 0, 1, 0)
#		cls.labels.pass1[1,] = bin.class.pass1
		cls.list.pass1[1,] = ifelse(bin.class.pass1 > 0, "MUT", "WT")
#		if( length(unique(bin.class.pass1)) == 1) { 
#			cls.list.3 = ifelse( cls.list == "DEL" | cls.list == "AMP", 1, 0)
#			copy.number.pass1 = ifelse( apply(cls.list.3[3:n.phen,], MARGIN=2, FUN=sum) >= 1, "ALT", "WT")
#			copy.class.pass1 = ifelse( copy.number.pass1 == "ALT", 1, 0)
#			bin.class.pass1 = copy.class.pass1
#			print( "Calculating p-value with respect to copy number alterations")
#		}
	} else {
		bin.class.pass1 <- ifelse(cls.list[1,] == cls.list2[1,1], 1, 0)
	}
#	browser()
	model.descs2.pass1 = vector(length = n.models, mode="character")
	for (i in 1:n.models) {
		m.score <- m[i,]
		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
#		browser()
		if (length(unique(bin.class.pass1)) > 1) {
			perf.auc <- rec.area(bin.class.pass1, m.score.norm)
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.pass1[i] = perf.auc$A
			p.val.list.pass1[i] = perf.auc$p.value
		} else {
			roc <- p.val <- "-"
			roc.list.pass1[i] = NA
			p.val.list.pass1[i] = NA
		}
#		print(paste("REC=", roc, " p-val=", p.val)) 
		
		model.descs2.pass1[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
	if( is.na(roc.list.pass1[1]) ){
		loc <- match(model, model.names)
		s.order.pass1 <- order(m[loc,], decreasing = decreasing.order)
#		loc = s.order.pass1[1]
		m2.pass1 <- m[, s.order.pass1]
		correl <- cor(t(m2.pass1))[, loc]
		m.order.pass1 <- order(correl, decreasing=T)
		m2.pass1 <- m2.pass1[m.order.pass1, ]
	} else{ 
		m.order.pass1 = order(roc.list.pass1, decreasing=TRUE, na.last=TRUE)
#		m.order.pass1 = order(p.val.list.pass1, decreasing=FALSE, na.last=TRUE)	
		m2.pass1 <- m[m.order.pass1, ]
		s.order.pass1 <- order(m2.pass1[1,], decreasing = TRUE )
		m2.pass1 <- m2.pass1[, s.order.pass1]
	}
#	m2.pass1 <- m2.pass1[m.order.pass1, ]
	model.descs2.pass1 = model.descs2.pass1[m.order.pass1]
	sample.names2.pass1 <- colnames(m2.pass1)
	model.names.pass1 <- rownames(m2.pass1)
#	browser()
	if (is.vector(cls.labels)) {
		cls.labels2.pass1 <- cls.labels.pass1[s.order.pass1]
		cls.list2.pass1 <- cls.list.pass1[s.order.pass1]
	} else {
		cls.labels2.pass1 <- cls.labels.pass1[, s.order.pass1]
		cls.list2.pass1 <- cls.list.pass1[, s.order.pass1]
	}
	#browser()
	
#	browser()
	m.score.pass1 <- m2.pass1[1,]
	m.score.norm.pass1 <- (m.score.pass1 - min(m.score.pass1))/(max(m.score.pass1) - min(m.score.pass1))
	roc.list.phen.pass1 = vector(mode="numeric", length=n.phen)
	p.val.list.phen.pass1 = vector(mode="numeric", length=n.phen)
	phen.descs.pass1 = vector(mode="character", length=n.phen)
	for( i in 1:n.phen.pass1 ){ 
		bin.gene = ifelse( cls.list2.pass1[i,]=="WT", 0, 1)
		if (length(unique(bin.gene)) > 1) {
			perf.auc <- roc.area(bin.gene, m.score.norm.pass1)
			
#			if( perf.auc$A < 0.5 ){
			##				browser()
#				roc = signif(1 - perf.auc$A, digits=3)
#				p.val = signif(1- perf.auc$p.val, digits=3)
#				abnormality = unique(cls.list2.pass1[i,])[which(unique(cls.list2.pass1[i,]) != "WT")]
			##				cls.list2.pass1[i,] = ifelse( cls.list2.pass1[i,] == "WT", abnormality, "WT" )
#				phen.names[i] = paste(phen.names[i], "-opposite", sep="")
#				roc.list.phen.pass1[i] = 1 - perf.auc$A
#				p.val.list.phen.pass1[i] = perf.auc$p.val   # Don't want to use these "opposite" genomic aberrations in Pass 2 
#															# because they make PATHWAY.MUT+COPY.NUMBER too dense
#			} else{
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.phen.pass1[i] = perf.auc$A
			p.val.list.phen.pass1[i] = perf.auc$p.val
#			}
		} else {
			roc <- "-"
			p.val <- "-"
			roc.list.phen.pass1[i] = NA
			p.val.list.phen.pass1[i] = NA
		}
#		print(paste("ROC=", roc, " p-val=", p.val)) 
		
#		p.val.list[i] = p.val
		phen.descs.pass1[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
	g.order.pass1 = c(1, order(roc.list.phen.pass1[2:n.phen.pass1], decreasing=TRUE, na.last=TRUE)+1)  # keep PATHWAY.MUT and COPY.NUMBER as 1 and 2
	roc.list.phen.pass1 = roc.list.phen.pass1[g.order.pass1]
	p.val.list.phen.pass1 = p.val.list.phen.pass1[g.order.pass1] 
	phen.descs2.pass1 = phen.descs.pass1[g.order.pass1][1:n.phen.pass1]
	cls.list2.pass1 = cls.list2.pass1[g.order.pass1,][1:n.phen.pass1,]
	phen.names.pass1 = phen.pass1[g.order.pass1][1:n.phen.pass1]
	
	# Recompute cls.list2 as some mutations or copy numbers may have been removed
	
	
	# Recompute cls.phen and cls.labels2 as order may have changed
	
	cls.phen2.pass1 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(cls.list2.pass1)
		cls.phen2.pass1 <- classes
		cls.labels2.pass1 <- match(cls.list2.pass1, cls.phen2.pass1)
	} else {
		for (kk in 1:length(cls.list2.pass1[, 1])) {
			classes <- unique(cls.list2.pass1[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2.pass1 <- c(cls.phen2.pass1, classes)
			cls.labels2.pass1[kk,] <- match(cls.list2.pass1[kk,], classes)
		}
	}
	cls.labels2.pass1 = cls.labels2.pass1[1:n.phen.pass1,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list.pass1 <- unlist(cls.phen2.pass1)
	
#	colors.list <- ifelse(unlist(cls.phen2) == target.class, 
#			ifelse(unlist(cls.phen2) == "DEL" | unlist(cls.phen2) == "AMP", 
#					ifelse(unlist(cls.phen2) == "DEL", cls.phen.colors[3], cls.phen.colors[4]), cls.phen.colors[1]), cls.phen.colors[2])
	colors.list.pass1 = rep( "gray", length(phen.list.pass1))
	colors.list.pass1[phen.list.pass1=="MUT"] = cls.phen.colors[1]
	colors.list.pass1[phen.list.pass1=="DEL"] = cls.phen.colors[3]
	colors.list.pass1[phen.list.pass1=="AMP"] = cls.phen.colors[4]
	colors.list.pass1[phen.list.pass1=="ALT"] = cls.phen.colors[5]
	
#	print("cls.phen2:")
#	print(unlist(cls.phen2))
#	
#	print("cls.phen:")
#	print(unlist(cls.phen))
#	
#	print("colors.list:")
#	print(colors.list)
	
#	browser()
	
	filename <- paste(results.dir, test.file.prefix, ".3Passes.REC_ks", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
	pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 17 )
#   windows(width=12, height=8)
	MSIG.HeatMapPlot.9(V = m2.pass1, row.names = model.names.pass1,
			row.names2 = model.descs2.pass1, 
			col.labels = cls.labels2.pass1, 
			col.classes = cls.phen2.pass1, 
			phen.cmap = colors.list.pass1, phen.names = phen.names.pass1,
			phen.names2 = phen.descs2.pass1,
			col.names = sample.names2.pass1, main = paste(tissue, "- 1st Pass: Known KRAS Pathway Abnormalities (REC)"), 
			xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	### Begin Pass 2 ###
	print( "--- Begin Pass 2 ---")
#	browser()
	p.val.threshold = 0.1
	ind.top.pval = which(p.val.list.phen.pass1[2:n.phen.pass1] <= p.val.threshold )+1
	if( length(ind.top.pval) > 0 ){
		ind.p.val.threshold = c(1, ind.top.pval) 
	} else if( length(ind.top.pval) < 1 ) { 
		p.val.threshold = 0.15
		ind.top.pval = which(p.val.list.phen.pass1[2:n.phen.pass1] <= p.val.threshold )+1
		ind.p.val.threshold = c(1, ind.top.pval) }
	if( length(ind.top.pval) < 1){ 
		p.val.threshold = 0.2
		ind.top.pval = which(p.val.list.phen.pass1[2:n.phen.pass1] <= p.val.threshold )+1
		ind.p.val.threshold = c(1, ind.top.pval)
	}
	if( length(ind.top.pval) < 1){ 
		p.val.threshold = 0.25
		ind.top.pval = which(p.val.list.phen.pass1[2:n.phen.pass1] <= p.val.threshold )+1
		ind.p.val.threshold = c(1, ind.top.pval)
	}
	if( length(ind.top.pval) < 1){ 
		p.val.threshold = 0.3
		ind.top.pval = which(p.val.list.phen.pass1[2:n.phen.pass1] <= p.val.threshold )+1
		ind.p.val.threshold = c(1, ind.top.pval)
	}
	if( length( ind.top.pval) < 1 ) {
		ind.top = which(!is.na(p.val.list.phen.pass1[-1]))+1 
		ind.p.val.threshold = c( 1, ind.top )
	}
	
	
	n.phen.pass2 = length(ind.p.val.threshold)
#	phen.descs.pass2 = phen.descs2.pass1[1:n.phen.pass2]
	cls.list2.pass2 = cls.list2.pass1[ind.p.val.threshold,]
	phen.names.pass2 = phen.names.pass1[ind.p.val.threshold]
#	phen.names.pass2[1] = "PATHWAY.MUT + COPY.NUMBER"
	cls.labels2.pass2 = cls.labels2.pass1[ind.p.val.threshold,]
	
#	browser()
	cls.list2.pass2.2 = ifelse( cls.list2.pass2 == "WT", 0, 1)
	cls.list2.pass2.3 = ifelse( cls.list2.pass2 == "DEL" | cls.list2.pass2 == "AMP", 1, 0)
	if( n.phen.pass2 > 2 ){ 
		pathway.mut.pass2 = apply(cls.list2.pass2.2[2:n.phen.pass2,], MARGIN=2, FUN=sum)
		bin.class.pass2 = pathway.mut.pass2/length(pathway.mut.pass2)
		bin.class.pass2 = ( bin.class.pass2 - min(bin.class.pass2))/(max(bin.class.pass2) - min(bin.class.pass2))
		cls.list2.pass2[1,] = ifelse( bin.class.pass2 > 0, "MUT", "WT")
#		copy.number.pass2 = ifelse( apply(cls.list2.pass2.3[3:n.phen.pass2,], MARGIN=2, FUN=sum) >= 1, "ALT", "WT")
	} else{
		pathway.mut.pass2 = ifelse( cls.list2.pass2.2[2,] == 1, "MUT", "WT")
		bin.class.pass2 = ifelse( pathway.mut.pass2 == "MUT", 1, 0 )
		cls.list2.pass2[1,] = pathway.mut.pass2
#		copy.number.pass2 = ifelse( cls.list2.pass2.3[3,] == 1, "ALT", "WT")
	}
#	browser()
#	bin.class.pass2 = 
#	copy.class.pass2 = ifelse( copy.number.pass2 == "ALT", 1, 0)
#	if( length(unique(bin.class.pass2)) == 1) { 
#		bin.class.pass2 = copy.class.pass2
#		print( "Calculating p-value with respect to copy number alterations")
#	}
	
#	cls.list2.pass2[2,] = copy.number.pass2
	
	roc.list.pass2 = vector( length=n.models, mode="numeric" )
	p.val.list.pass2 = vector( length=n.models, mode="numeric" )
	
#	browser()
	model.descs2.pass2 = vector(length = n.models, mode="character")
	for (i in 1:n.models) {
		m.score <- m2.pass1[i,]
		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
#		browser()
		if (length(unique(bin.class.pass2)) > 1) {
			perf.auc <- rec.area(bin.class.pass2, m.score.norm)
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.pass2[i] = perf.auc$A
			p.val.list.pass2[i] = perf.auc$p.value
		} else {
			roc <- p.val <- "-"
			roc.list.pass2[i] = NA
			p.val.list.pass2[i] = NA
		}
		print(paste("REC=", roc, " p-val=", p.val)) 
		
		model.descs2.pass2[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
	m.order.pass2 = order(roc.list.pass2, decreasing=TRUE, na.last=TRUE)
#	m.order.pass2 = order(p.val.list.pass2, decreasing=FALSE, na.last=TRUE)
	model.descs2.pass2 = model.descs2.pass2[m.order.pass2]
#	loc.pass2 = m.order.pass2[1]
	m2.pass2 <- m2.pass1[m.order.pass2, ]
	model.names.pass2 <- rownames(m2.pass2)
#	print(c("loc.pass2:", loc.pass2))
	s.order.pass2 <- order(m2.pass2[1,], decreasing = TRUE)
	m2.pass2 <- m2.pass2[, s.order.pass2]
	sample.names2.pass2 <- colnames(m2.pass2)
	
	if (is.vector(cls.labels)) {
		cls.labels2.pass2 <- cls.labels2.pass2[s.order.pass2]
		cls.list2.pass2 <- cls.list2.pass2[s.order.pass2]
	} else {
		cls.labels2.pass2 <- cls.labels2.pass2[, s.order.pass2]
		cls.list2.pass2 <- cls.list2.pass2[, s.order.pass2]
	}
#	browser()
	
#	browser()
	m.score.pass2 <- m2.pass2[1,]
	m.score.norm.pass2 <- (m.score.pass2 - min(m.score.pass2))/(max(m.score.pass2) - min(m.score.pass2))
	roc.list.phen.pass2 = vector(mode="numeric", length=n.phen.pass2)
	phen.descs.pass2 = vector(mode="character", length=n.phen.pass2)
	for( i in 1:n.phen.pass2 ){ 
		bin.gene = ifelse( cls.list2.pass2[i,]=="WT", 0, 1)
		if (length(unique(bin.gene)) > 1) {
			perf.auc <- roc.area(bin.gene, m.score.norm.pass2)
#			if( perf.auc$A < 0.5 ){
			##				browser()
#				roc = signif(1 - perf.auc$A, digits=3)
#				p.val = signif(1 - perf.auc$A, digits=3)
#				abnormality = unique(cls.list2.pass2[i,])[which(unique(cls.list2.pass2[i,]) != "WT")]
#				cls.list2.pass2 = ifelse( cls.list2.pass2[i,] == "WT", abnormality, "WT" )
#				phen.names.pass2[i] = paste(phen.names.pass2[i], "-opposite")
#			} else{
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.phen.pass2[i] = perf.auc$A
#			}
		} else {
			roc <- "-"
			p.val <- "-"
			roc.list.phen.pass2[i] = NA
		}
		print(paste("ROC=", roc, " p-val=", p.val)) 
#		p.val.list[i] = p.val
		phen.descs.pass2[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
	g.order.pass2 = c(1, order(roc.list.phen.pass2[2:n.phen.pass2], decreasing=TRUE, na.last=TRUE)+1)  # skip PATHWAY.MUT
	phen.descs2.pass2 = phen.descs.pass2[g.order.pass2]
	cls.list2.pass2 = cls.list2.pass2[g.order.pass2,]
	phen.names.pass2 = phen.names.pass2[g.order.pass2]
#	browser()
	# Recompute cls.list2 as some mutations or copy numbers may have been removed
	
	
	# Recompute cls.phen and cls.labels2 as order may have changed
	
	cls.phen2.pass2 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(as.vector(cls.list2.pass2))
		cls.phen2.pass2 <- classes
		cls.labels2.pass2 <- match(cls.list2.pass2, cls.phen2.pass2)
	} else {
		for (kk in 1:length(cls.list2.pass2[, 1])) {
			classes <- unique(cls.list2.pass2[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2.pass2 <- c(cls.phen2.pass2, classes)
			cls.labels2.pass2[kk,] <- match(cls.list2.pass2[kk,], classes)
		}
	}
	cls.labels2.pass2 = cls.labels2.pass2[1:n.phen.pass2,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list.pass2 <- unlist(cls.phen2.pass2)
	
#	colors.list <- ifelse(unlist(cls.phen2) == target.class, 
#			ifelse(unlist(cls.phen2) == "DEL" | unlist(cls.phen2) == "AMP", 
#					ifelse(unlist(cls.phen2) == "DEL", cls.phen.colors[3], cls.phen.colors[4]), cls.phen.colors[1]), cls.phen.colors[2])
	colors.list.pass2 = rep( "gray", length(phen.list.pass2))
	colors.list.pass2[phen.list.pass2=="MUT"] = cls.phen.colors[1]
	colors.list.pass2[phen.list.pass2=="DEL"] = cls.phen.colors[3]
	colors.list.pass2[phen.list.pass2=="AMP"] = cls.phen.colors[4]
	colors.list.pass2[phen.list.pass2=="ALT"] = cls.phen.colors[5]
	
	MSIG.HeatMapPlot.9(V = m2.pass2, row.names = model.names.pass2,
			row.names2 = model.descs2.pass2, 
			col.labels = cls.labels2.pass2, 
			col.classes = cls.phen2.pass2, 
			phen.cmap = colors.list.pass2, phen.names = phen.names.pass2,
			phen.names2 = phen.descs2.pass2,
			col.names = sample.names2.pass2, main = paste(tissue, "- 2nd Pass: only p-values <=", p.val.threshold,"from 1st pass (REC)"), 
			xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	### 3rd Pass ###	
	
	print( "--- Begin Pass 3 ---")
#	browser()
	m2.pass3 = m2.pass2
	model.names.pass3 = rownames(m2.pass3)
	sample.names2.pass3 = colnames(m2.pass3)
#	model.descs2.pass3 = model.descs2.pass2
	n.phen.pass3 = 40
#	phen.descs.pass2 = phen.descs2.pass1[1:n.phen.pass2]
	cls.list2.pass3 = cls.list[, s.order.pass1][, s.order.pass2]
	cls.labels2.pass3 = cls.labels[, s.order.pass1][, s.order.pass2]
	
#	browser()
	phen.names.pass3 = phen.names
	m.score.pass3 <- m2.pass3[1,]
	m.score.norm.pass3 <- (m.score.pass3 - min(m.score.pass3))/(max(m.score.pass3) - min(m.score.pass3))
	roc.list.phen.pass3 = vector(mode="numeric", length=n.phen)
	phen.descs.pass3 = vector(mode="character", length=n.phen)
	p.val.list.phen.pass3 = vector(mode="numeric", length=n.phen)
	for( i in 1:n.phen ){ 
		bin.gene = ifelse( cls.list2.pass3[i,]=="WT", 0, 1)
		if (length(unique(bin.gene)) > 1) {
			perf.auc <- roc.area(bin.gene, m.score.norm.pass3)
#			if( perf.auc$A < 0.5 ){
			##				browser()
#				roc = signif(1 - perf.auc$A, digits=3)
#				p.val = signif(1- perf.auc$p.value, digits=3)
#				abnormality = unique(cls.list2.pass3[i,])[which(unique(cls.list2.pass3[i,]) != "WT")]
#				cls.list2.pass3[i,] = ifelse( cls.list2.pass3[i,] == "WT", abnormality, "WT" )
#				phen.names.pass3[i] = paste(phen.names.pass3[i], "-opposite", sep="")
#				roc.list.phen.pass3[i] = 1-perf.auc$A
#				p.val.list.phen.pass3[i] = 1- perf.auc$p.value
#			} else{
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.phen.pass3[i] = perf.auc$A
			p.val.list.phen.pass3[i] = perf.auc$p.value
#			}
		} else {
			roc <- NA
			p.val <- NA
			roc.list.phen.pass3[i] = NA
			p.val.list.phen.pass3[i] = NA
		}
#		print(paste("ROC=", roc, " p-val=", p.val)) 
		
#		p.val.list[i] = p.val
		phen.descs.pass3[i] <- paste(roc, " (", p.val, ")")
	}
#	browser()
	p.val.threshold = 0.1
	len = length(which(p.val.list.phen.pass3[-1:-2] <= p.val.threshold))+2
	if( len == 2 ){
		p.val.threshold = 0.15
		len = length(which(p.val.list.phen.pass3[-1:-2] <= p.val.threshold))+2
	}
	if( len == 2 ){
		p.val.threshold = 0.2
		len = length(which(p.val.list.phen.pass3[-1:-2] <= p.val.threshold))+2
	}
	if( len>40 ) len = 40
#	g.order.pass3.1 = c(1, 2, order(p.val.list.phen.pass3[3:n.phen], decreasing=FALSE, na.last=TRUE)+2 )
	g.order.pass3 =  c(1, 2, order(p.val.list.phen.pass3[-1:-2], decreasing=FALSE, na.last=TRUE)+2 )[1:len]  # skip PATHWAY.MUT and COPY.NUMBER
	phen.descs2.pass3 = phen.descs.pass3[g.order.pass3]
	cls.list2.pass3 = cls.list2.pass3[g.order.pass3,]
	cls.labels2.pass3 = cls.labels2.pass3[g.order.pass3,]
	phen.names.pass3 = phen.names.pass3[g.order.pass3]
	
	
	cls.list.mut = ifelse(cls.list2.pass3[-1:-2,] == "MUT", 1, 0)
	cls.list.alt = ifelse(cls.list2.pass3[-1:-2,] == "DEL" | cls.list2.pass3[-1:-2,] == "AMP", 1, 0)
	
#	browser()
	if( !is.vector(cls.list.alt) ){
		cls.list.mut.sum = apply(cls.list.mut, MARGIN=2, FUN=sum)
		cls.list.alt.sum = apply(cls.list.alt, MARGIN=2, FUN=sum)	
		bin.class.pass3 = cls.list.mut.sum + cls.list.alt.sum
		bin.class.pass3 = ( bin.class.pass3 - min(bin.class.pass3))/(max(bin.class.pass3) - min(bin.class.pass3))
		cls.list.mut.sum = ifelse(cls.list.mut.sum + cls.list.alt.sum > 0, 1, 0)
		cls.list2.pass3[1,] = ifelse( cls.list.mut.sum >= 1, "MUT", "WT")
		cls.list2.pass3[2,] = ifelse( cls.list.alt.sum >= 1, "ALT", "WT")
		
		
	} else{
		
		cls.list2.pass3[2,] = ifelse(cls.list.alt == 1, "ALT", "WT")
		bin.class.pass3 = cls.list.mut+cls.list.alt
		bin.class.pass3 = ( bin.class.pass3 - min(bin.class.pass3))/(max(bin.class.pass3) - min(bin.class.pass3))
		cls.list2.pass3[1,] = ifelse(bin.class.pass3 > 0 , "MUT", "WT")
	}
	
#	browser()
	for( i in 1:2 ){ # Recalculate ROC and p-value for PATHWAY.MUT and COPY.NUMBER
		bin.gene = ifelse( cls.list2.pass3[i,]=="WT", 0, 1)
		if (length(unique(bin.gene)) > 1) {
			perf.auc <- rec.area(bin.gene, m.score.norm.pass3)
#			if( perf.auc$A < 0.5 ){
			##				browser()
#				roc = signif(1 - perf.auc$A, digits=3)
#				p.val = signif(1- perf.auc$p.value, digits=3)
#				abnormality = unique(cls.list2.pass3[i,])[which(unique(cls.list2.pass3[i,]) != "WT")]
#				cls.list2.pass3[i,] = ifelse( cls.list2.pass3[i,] == "WT", abnormality, "WT" )
#				phen.names.pass3[i] = paste(phen.names.pass3[i], "-opposite", sep="")
#				roc.list.phen.pass3[i] = 1-perf.auc$A
#				p.val.list.phen.pass3[i] = 1- perf.auc$p.value
#			} else{
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.phen.pass3[i] = perf.auc$A
			p.val.list.phen.pass3[i] = perf.auc$p.value
#			}
		} else {
			roc <- NA
			p.val <- NA
			roc.list.phen.pass3[i] = NA
			p.val.list.phen.pass3[i] = NA
		}
		print(paste("ROC=", roc, " p-val=", p.val)) 
		
#		p.val.list[i] = p.val
		phen.descs2.pass3[i] <- paste(roc, " (", p.val, ")")
	}
	
#	browser()
	model.descs2.pass3 = vector(length = n.models, mode="character")
	for (i in 1:n.models) {
		m.score <- m2.pass3[i,]
		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
#		browser()
		if (length(unique(bin.class.pass3)) > 1) {
			perf.auc <- rec.area(bin.class.pass3, m.score.norm)
			roc <- signif(perf.auc$A, digits=3)
			p.val <- signif(perf.auc$p.value, digits=3)
			roc.list.pass2[i] = perf.auc$A
			p.val.list.pass2[i] = perf.auc$p.value
		} else {
			roc <- p.val <- "-"
			roc.list.pass2[i] = NA
			p.val.list.pass2[i] = NA
		}
		print(paste("REC=", roc, " p-val=", p.val)) 
		
		model.descs2.pass3[i] <- paste(roc, " (", p.val, ")")
	}
	
	cls.phen2.pass3 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(as.vector(cls.list2.pass3))
		cls.phen2.pass3 <- classes
		cls.labels2.pass3 <- match(cls.list2.pass3, cls.phen2.pass3)
	} else {
#		browser()
		for (kk in 1:length(cls.list2.pass3[, 1])) {
#			browser()
			classes <- unique(cls.list2.pass3[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2.pass3 <- c(cls.phen2.pass3, classes)
			cls.labels2.pass3[kk,] <- match(cls.list2.pass3[kk,], classes)
		}
	}
#	cls.labels2.pass3 = cls.labels2.pass3[1:n.phen.pass3,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list.pass3 <- unlist(cls.phen2.pass3)
	
#	colors.list <- ifelse(unlist(cls.phen2) == target.class, 
#			ifelse(unlist(cls.phen2) == "DEL" | unlist(cls.phen2) == "AMP", 
#					ifelse(unlist(cls.phen2) == "DEL", cls.phen.colors[3], cls.phen.colors[4]), cls.phen.colors[1]), cls.phen.colors[2])
	colors.list.pass3 = rep( "gray", length(phen.list.pass3))
	colors.list.pass3[phen.list.pass3=="MUT"] = cls.phen.colors[1]
	colors.list.pass3[phen.list.pass3=="DEL"] = cls.phen.colors[3]
	colors.list.pass3[phen.list.pass3=="AMP"] = cls.phen.colors[4]
	colors.list.pass3[phen.list.pass3=="ALT"] = cls.phen.colors[5]
	phen.names.pass3[1] = "PATHWAY.MUT+COPY.NUMBER"
#	browser()
	MSIG.HeatMapPlot.9(V = m2.pass3, row.names = model.names.pass3,
			row.names2 = model.descs2.pass3, 
			col.labels = cls.labels2.pass3, 
			col.classes = cls.phen2.pass3, 
			phen.cmap = colors.list.pass3, phen.names = phen.names.pass3,
			phen.names2 = phen.descs2.pass3,
			col.names = sample.names2.pass3, main = paste(tissue, "- 3rd Pass: Top signature from 2nd pass with all genes ( p-value <=", p.val.threshold, ") (REC)"), 
			xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	dev.off()
	
	if (!is.na(output.dataset)) {
		V.GCT <- m2
		colnames(V.GCT) <- sample.names2
		row.names(V.GCT) <- model.names2
		write.gct(gct.data.frame = V.GCT, descs = model.descs2, filename =output.dataset)  
	}
	
}

OPAM.sort.projection.by.score.6 <- function(
		input.ds,
		input.cls,
		tissue = "NA",
		results.dir,
		normalize.score = T,
		normalization.type = "zero.one",
		model = "NA",
		target.phen = NA,
		target.class = NA,
		user.colors = NA,
		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T,
		u.gene.names.known = "NA"
)
# Calls MSIG.HeatMapPlot.9 and makes a plot sorted by the highest-scoring
# signatures and abnormalities (gene mutations or copy number alterations)
# i.e. doesn't require a "model" to score by as OPAM.sort.projection.by.score.2 does.
# However, it *will* use "model" if it cannot calculate p-values on the gene signatures, which
# happens when every cell line has a genomic aberration.
#
# Runs 3 passes on the data:
# 1st pass: looks at the genes and copy number alterations specified by u.gene.names.known
# 2nd pass: looks at only the top abnormalities (using a p-value cutoff) from the 1st pass, and adjusts 
# the PATHWAY.MUT vector accordingly (only according to the genes, not by copy number data)
# 3rd pass: Takes the winning signature from the 2nd pass and then looks all the genes available
#
# Very similar to OPAM.sort.projection.by.score.4, however this version uses mutual.inf instead of
# roc.area to calculate mutual information scores and p-values for PATHWAY.MUT, the vector of total genomic aberrations
# in all samples
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	model.names <- dataset$row.names
#	model.descs <- dataset$descs
	Ns <- length(m[1,])
	dim(m)
	sample.names <- dataset$names
	
	n.models <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
	char.res <-  0.013 * n.models + 0.65
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list 
	
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors[1])) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
	
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
	} else {
		phen.names <- "NA"
	}
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
#	print("cls.phen.colors:")
#	print(cls.phen.colors)
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		for (i in 1:n.phen) {
			n.classes[i] <- length(cls.phen[[i]])
		}
	}
	print("--- Begin Pass 1 ---")
#	model.names.original = model.names
#	m.original = m
	phen.pass1 = c( "PATHWAY.MUT", u.gene.names.known)
	n.phen.pass1 = length(u.gene.names.known)+1
	ind.phen.pass1 = which( phen.names %in% phen.pass1 )
	phen.pass1 = phen.names[ind.phen.pass1]
	phen.pass1[1] = "SUMMARY"
	
	MI.list.pass1 = vector( length=n.models, mode="numeric" )
#	p.val.list.pass1 = vector( length=n.models, mode="numeric" )
	
	cls.list.pass1 = cls.list[ind.phen.pass1,]
	cls.list.pass1.2 = ifelse(cls.list.pass1 == "WT", 0, 1)
	cls.labels.pass1 = cls.labels[ind.phen.pass1,]
	browser()
	if (!is.na(target.phen)) {
		bin.class.pass1 = apply( cls.list.pass1.2[-1,], MARGIN=2, FUN=sum)
		# Normalize bin.class.pass1
		if( length(unique(bin.class.pass1)) > 1){
			bin.class.pass1 = ( bin.class.pass1 - min(bin.class.pass1))/(max(bin.class.pass1) - min(bin.class.pass1))
		} else if ( length(unique(bin.class.pass1)) == 1){
			bin.class = rep(1, length(cls.list[1,]))
		}
#		bin.class.pass1 <- ifelse(cls.list[target.phen,] == target.class, 1, 0)
#		bin.class.pass1 <- ifelse(apply( cls.list.pass1.2[2:n.phen.pass1,], MARGIN=2, FUN=sum) > 0, 1, 0)
#		cls.labels.pass1[1,] = bin.class.pass1
		cls.list.pass1[1,] = ifelse(bin.class.pass1 > 0, "MUT", "WT")
#		if( length(unique(bin.class.pass1)) == 1) { 
#			cls.list.3 = ifelse( cls.list == "DEL" | cls.list == "AMP", 1, 0)
#			copy.number.pass1 = ifelse( apply(cls.list.3[3:n.phen,], MARGIN=2, FUN=sum) >= 1, "ALT", "WT")
#			copy.class.pass1 = ifelse( copy.number.pass1 == "ALT", 1, 0)
#			bin.class.pass1 = copy.class.pass1
#			print( "Calculating p-value with respect to copy number alterations")
#		}
	} else {
		bin.class.pass1 <- ifelse(cls.list[1,] == cls.list2[1,1], 1, 0)
	}
#	browser()
#	MI.ref.models.pass1 = mutual.inf.2(bin.class.pass1, bin.class.pass1)$MI
#	print(paste("MI.ref.models.pass1 =", MI.ref.models.pass1))
#	browser()
	model.descs2.pass1 = vector(length = n.models, mode="character")
	for (i in 1:n.models) {
		m.score <- m[i,]
		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
#		browser()
		if (length(unique(bin.class.pass1)) > 1) {
#			browser()
			MI <- (mutual.inf.2(bin.class.pass1, m.score.norm)$MI)#	/MI.ref.models.pass1
#			roc <- signif(perf.auc$A, digits=3)
#			p.val <- signif(perf.auc$p.value, digits=3)
			MI.list.pass1[i] = MI
			MI.signif <- signif(MI, digits=3)
#			p.val.list.pass1[i] = perf.auc$p.value
		} else {
			MI.signif <- "-"
			MI.list.pass1[i] = NA
#			p.val.list.pass1[i] = NA
		}
#		browser()
#		print(paste("REC=", roc, " p-val=", p.val)) 
		print(paste( format(rownames(m)[i], width=30), "mutual.inf =", MI.signif))
#		browser()
		model.descs2.pass1[i] <- paste(MI.signif)
	}
#	browser()
	if( is.na(MI.list.pass1[1]) ){
		loc <- match(model, model.names)
		s.order.pass1 <- order(m[loc,], decreasing = decreasing.order)
#		loc = s.order.pass1[1]
#	s.order.pass1 = 1:Ns
		m2.pass1 <- m[, s.order.pass1]
		correl <- cor(t(m2.pass1))[, loc]
		m.order.pass1 <- order(correl, decreasing=T)
#	m.order.pass1 = 1:n.models
		m2.pass1 <- m2.pass1[m.order.pass1, ]
	} else{ 
		m.order.pass1 = order(MI.list.pass1, decreasing=TRUE, na.last=TRUE)
#m.order.pass1 = 1:n.models
		m2.pass1 <- m[m.order.pass1, ]
		s.order.pass1 <- order(m2.pass1[1,], decreasing = TRUE )
#		s.order.pass1 = 1:Ns
		m2.pass1 <- m2.pass1[, s.order.pass1]
	}
	bin.class.pass1 = bin.class.pass1[s.order.pass1]
	m2.pass1 <- m2.pass1[m.order.pass1, ]
	model.descs2.pass1 = model.descs2.pass1[m.order.pass1]
	sample.names2.pass1 <- colnames(m2.pass1)
	model.names.pass1 <- rownames(m2.pass1)
#	browser()
	if (is.vector(cls.labels)) {
		cls.labels2.pass1 <- cls.labels.pass1[s.order.pass1]
		cls.list2.pass1 <- cls.list.pass1[s.order.pass1]
	} else {
		cls.labels2.pass1 <- cls.labels.pass1[, s.order.pass1]
		cls.list2.pass1 <- cls.list.pass1[, s.order.pass1]
	}
	#browser()
	
#	pathway.name <- "KRAS_ALL_UP"
#	pathway <- m[1,]
#	pathway0 <- ifelse(pathway < median(pathway), 0, 1) # disctretized version
	
#	MI.ref.genes.pass1 <- mutual.inf.2(m[1,], m[1,])$MI
	
#	browser()
	m.score.pass1 <- m2.pass1[1,]
	m.score.norm.pass1 <- (m.score.pass1 - min(m.score.pass1))/(max(m.score.pass1) - min(m.score.pass1))
#	m.score.pass1 = ifelse( m.score.pass1 < median(m.score.pass1), -1, 1)   # discretized version
#	MI.ref.genes.pass1 <- mutual.inf.2(m.score.norm.pass1, m.score.norm.pass1)$MI
#	print(paste("MI.ref.genes.pass1 =", MI.ref.genes.pass1))
	MI.list.phen.pass1 = vector(mode="numeric", length=n.phen.pass1)
#	p.val.list.phen.pass1 = vector(mode="numeric", length=n.phen)
	phen.descs.pass1 = vector(mode="character", length=n.phen.pass1)
	
	if( length(unique(bin.class.pass1)) > 1){
		MI <-(mutual.inf.2(bin.class.pass1, m.score.norm.pass1)$MI)#/MI.ref.genes.pass1
		MI.signif <- signif(MI, digits=3)
		MI.list.phen.pass1[1] = MI
	} else{
		MI.signif <- "-"
		MI.list.phen.pass1[1] = NA
	}
	print(paste(format(phen.pass1[1], width=12), "mutual.inf =", MI.signif))
	phen.descs.pass1[1] <- paste(MI.signif)
	
#	print(m.score.pass1)
	for( i in 2:n.phen.pass1 ){ 
#		browser()
		bin.gene = ifelse( cls.list2.pass1[i,]=="WT", 0, 1)  
		# add random noise so the quantile calculation in mutual.inf doesn't return 0
		if (length(unique(bin.gene)) > 1) {
#			print(bin.gene)
#			browser()
			MI <- (mutual.inf.2(bin.gene, m.score.norm.pass1)$MI)#/MI.ref.genes.pass1
			MI.signif <- signif(MI, digits=3)
			
#			if( perf.auc$A < 0.5 ){
			##				browser()
#				roc = signif(1 - perf.auc$A, digits=3)
#				p.val = signif(1- perf.auc$p.val, digits=3)
#				abnormality = unique(cls.list2.pass1[i,])[which(unique(cls.list2.pass1[i,]) != "WT")]
			##				cls.list2.pass1[i,] = ifelse( cls.list2.pass1[i,] == "WT", abnormality, "WT" )
#				phen.names[i] = paste(phen.names[i], "-opposite", sep="")
#				roc.list.phen.pass1[i] = 1 - perf.auc$A
#				p.val.list.phen.pass1[i] = perf.auc$p.val   # Don't want to use these "opposite" genomic aberrations in Pass 2 
#															# because they make PATHWAY.MUT+COPY.NUMBER too dense
#			} else{
#			roc <- signif(perf.auc$A, digits=3)
#			p.val <- signif(perf.auc$p.value, digits=3)
			MI.list.phen.pass1[i] = MI
#			p.val.list.phen.pass1[i] = perf.auc$p.val
#			}
		} else {
			MI.signif <- "-"
#			p.val <- "-"
			MI.list.phen.pass1[i] = NA
#			p.val.list.phen.pass1[i] = NA
		}
#		browser()
#		print(paste("ROC=", roc, " p-val=", p.val)) 
		print(paste(format(phen.pass1[i], width=12), "mutual.inf =", MI.signif))
#		p.val.list[i] = p.val
		phen.descs.pass1[i] <- paste(MI.signif)
	}
#	browser()
#g.order.pass1 = 1:n.phen.pass1
	g.order.pass1 = c(1, order(MI.list.phen.pass1[2:n.phen.pass1], decreasing=TRUE, na.last=TRUE)+1)  # keep PATHWAY.MUT as 1
#	MI.list.phen.pass1 = MI.list.phen.pass1[g.order.pass1]
#	p.val.list.phen.pass1 = p.val.list.phen.pass1[g.order.pass1] 
	MI.list.phen.pass1 = MI.list.phen.pass1[g.order.pass1]
	phen.descs2.pass1 = phen.descs.pass1[g.order.pass1]#[1:n.phen.pass1]
	cls.list2.pass1 = cls.list2.pass1[g.order.pass1,]#[1:n.phen.pass1,]
	phen.names.pass1 = phen.pass1[g.order.pass1]#[1:n.phen.pass1]
	
	# Recompute cls.list2 as some mutations or copy numbers may have been removed
	
	
	# Recompute cls.phen and cls.labels2 as order may have changed
	
	cls.phen2.pass1 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(cls.list2.pass1)
		cls.phen2.pass1 <- classes
		cls.labels2.pass1 <- match(cls.list2.pass1, cls.phen2.pass1)
	} else {
		for (kk in 1:length(cls.list2.pass1[, 1])) {
			classes <- unique(cls.list2.pass1[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2.pass1 <- c(cls.phen2.pass1, classes)
			cls.labels2.pass1[kk,] <- match(cls.list2.pass1[kk,], classes)
		}
	}
	cls.labels2.pass1 = cls.labels2.pass1[1:n.phen.pass1,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list.pass1 <- unlist(cls.phen2.pass1)
	colors.list.pass1 = rep( "gray", length(phen.list.pass1))
	colors.list.pass1[phen.list.pass1=="MUT"] = cls.phen.colors[1]
	colors.list.pass1[phen.list.pass1=="DEL"] = cls.phen.colors[3]
	colors.list.pass1[phen.list.pass1=="AMP"] = cls.phen.colors[4]
	colors.list.pass1[phen.list.pass1=="ALT"] = cls.phen.colors[5]
#	browser()
#	colors.list.pass1[1,] = grey(bin.class.pass1)
#	print("cls.phen2:")
#	print(unlist(cls.phen2))
#	
#	print("cls.phen:")
#	print(unlist(cls.phen))
#	
#	print("colors.list:")
#	print(colors.list)
	
	browser()
	
	filename <- paste(results.dir, test.file.prefix, ".Phase1-2.MI", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
	pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 17 )
#   windows(width=12, height=8)
	MSIG.HeatMapPlot.10(V = m2.pass1, 
			pathway.mut = bin.class.pass1,
			row.names = model.names.pass1,
			row.names2 = model.descs2.pass1, 
			col.labels = cls.labels2.pass1, 
			col.classes = cls.phen2.pass1, 
			phen.cmap = colors.list.pass1, 
			phen.names = phen.names.pass1,
			phen.names2 = phen.descs2.pass1,
			col.names = sample.names2.pass1, 
			main = paste(tissue, "- Phase 1: Known KRAS Pathway Abnormalities (MI)"), 
			xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	### Begin Pass 2 ###
	print( "--- Begin Phase 2 ---")
#	browser()
	MI.thresholds = c(0.2, 0.1, 0.08, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01, 0)
#	MI.threshold = 0.03
	ind.top.MI = vector(mode="integer")
	MI.i = 1
	while( length(ind.top.MI) < 1)
	{
		MI.i = MI.i + 1
		ind.top.MI = which( MI.list.phen.pass1[-1] >= MI.thresholds[MI.i] ) + 1
		
	}
	ind.MI.threshold = c(1, ind.top.MI)
#	ind.top.MI = which(MI.list.phen.pass1[-1] >= MI.threshold.vector[1] )+1
#	if( length(ind.top.MI) > 0 ){
#		ind.MI.threshold = c(1, ind.top.MI) 
#	} 
#	if( length(ind.top.MI) < 1 ){
#		MI.threshold = 0.025
#		ind.top.MI = which(MI.list.phen.pass1[-1] >= MI.threshold )+1
#		ind.MI.threshold = c(1, ind.top.MI) 
#	} 
#	if( length(ind.top.MI) < 1 ){
#		MI.threshold = 0.02
#		ind.top.MI = which(MI.list.phen.pass1[-1] >= MI.threshold )+1
#		ind.MI.threshold = c(1, ind.top.MI) 
#	} 
#	if( length(ind.top.MI) < 1 ){
#		MI.threshold = 0.015
#		ind.top.MI = which(MI.list.phen.pass1[-1] >= MI.threshold )+1
#		ind.MI.threshold = c(1, ind.top.MI) 
#	} 
#	if( length(ind.top.MI) < 1 ){
#		MI.threshold = 0.01
#		ind.top.MI = which(MI.list.phen.pass1[-1] >= MI.threshold )+1
#		ind.MI.threshold = c(1, ind.top.MI) 
#	} 
#	if( length(ind.top.MI) < 1 ) { 
#		MI.threshold = 0
#		ind.top.MI = which(MI.list.phen.pass1[-1] > 0  )+1
#		ind.MI.threshold = c(1, ind.top.MI) }
#	if( length(ind.top.MI) < 1){ 
#		MI.threshold = 0.2
#		ind.top.MI = which(MI.list.phen.pass1[2:n.phen.pass1] <= MI.threshold )+1
#		ind.MI.threshold = c(1, ind.top.MI)
#	}
#	if( length(ind.top.MI) < 1){ 
#		MI.threshold = 0.25
#		ind.top.MI = which(MI.list.phen.pass1[2:n.phen.pass1] <= MI.threshold )+1
#		ind.MI.threshold = c(1, ind.top.MI)
#	}
#	if( length(ind.top.MI) < 1){ 
#		MI.threshold = 0.3
#		ind.top.MI = which(MI.list.phen.pass1[2:n.phen.pass1] <= MI.threshold )+1
#		ind.MI.threshold = c(1, ind.top.MI)
#	}
#	if( length( ind.top.MI) < 1 ) {
#		ind.top = which(!is.na(MI.list.phen.pass1[-1]))+1 
#		ind.MI.threshold = c( 1, ind.top )
#	}
	
	
	n.phen.pass2 = length(ind.MI.threshold)
#	phen.descs.pass2 = phen.descs2.pass1[1:n.phen.pass2]
	cls.list2.pass2 = cls.list2.pass1[ind.MI.threshold,]
	phen.names.pass2 = phen.names.pass1[ind.MI.threshold]
#	phen.names.pass2[1] = "PATHWAY.MUT + COPY.NUMBER"
	cls.labels2.pass2 = cls.labels2.pass1[ind.MI.threshold,]
	
#	browser()
	cls.list2.pass2.2 = ifelse( cls.list2.pass2 == "WT", 0, 1)
	cls.list2.pass2.3 = ifelse( cls.list2.pass2 == "DEL" | cls.list2.pass2 == "AMP", 1, 0)
	if( n.phen.pass2 > 2 ){ 
		pathway.mut.pass2 = apply(cls.list2.pass2.2[2:n.phen.pass2,], MARGIN=2, FUN=sum)
		bin.class.pass2 = pathway.mut.pass2/length(pathway.mut.pass2)
		bin.class.pass2 = ( bin.class.pass2 - min(bin.class.pass2))/(max(bin.class.pass2) - min(bin.class.pass2))
		bin.class.pass2.noisy = bin.class.pass2
		cls.list2.pass2[1,] = ifelse( bin.class.pass2 > 0, "MUT", "WT")
#		copy.number.pass2 = ifelse( apply(cls.list2.pass2.3[3:n.phen.pass2,], MARGIN=2, FUN=sum) >= 1, "ALT", "WT")
	} else{
		pathway.mut.pass2 = ifelse( cls.list2.pass2.2[2,] == 1, "MUT", "WT")
		bin.class.pass2 = ifelse( pathway.mut.pass2 == "MUT", 1, 0 ) #+ runif(Ns, min=-.05, max=.05)
#		bin.class.pass2.noisy = bin.class.pass2 + runif(Ns, min=-.05, max=.05)
#		bin.class.pass2.noisy = ( bin.class.pass2.noisy - min(bin.class.pass2.noisy))/(max(bin.class.pass2.noisy) - min(bin.class.pass2.noisy))
		cls.list2.pass2[1,] = pathway.mut.pass2
#		copy.number.pass2 = ifelse( cls.list2.pass2.3[3,] == 1, "ALT", "WT")
	}
#	browser()
#	bin.class.pass2 = 
#	copy.class.pass2 = ifelse( copy.number.pass2 == "ALT", 1, 0)
#	if( length(unique(bin.class.pass2)) == 1) { 
#		bin.class.pass2 = copy.class.pass2
#		print( "Calculating p-value with respect to copy number alterations")
#	}
	
#	cls.list2.pass2[2,] = copy.number.pass2
	
	MI.list.pass2 = vector( length=n.models, mode="numeric" )
#	MI.ref.models.pass2 = mutual.inf.2(bin.class.pass2, bin.class.pass2)$MI
#	print(paste("MI.ref.models.pass2 =", MI.ref.models.pass2))
#	p.val.list.pass2 = vector( length=n.models, mode="numeric" )
	
#	browser()
	model.descs2.pass2 = vector(length = n.models, mode="character")
	for (i in 1:n.models) {
		m.score <- m2.pass1[i,]
		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
#		browser()
		if (length(unique(bin.class.pass2)) > 1) {
			MI <- (mutual.inf.2(bin.class.pass2, m.score.norm)$MI)#/MI.ref.models.pass2
			MI.signif <- signif(MI, digits=3)
#			p.val <- signif(perf.auc$p.value, digits=3)
			MI.list.pass2[i] = MI
#			p.val.list.pass2[i] = perf.auc$p.value
		} else {
			MI.signif <- "-"
			MI.list.pass2[i] = NA
#			p.val.list.pass2[i] = NA
		}
		print(paste(format(rownames(m2.pass1)[i], width=30),"mutual.inf =",  MI.signif)) 
		
		model.descs2.pass2[i] <- paste(MI.signif)
	}
#	browser()
	m.order.pass2 = order(MI.list.pass2, decreasing=TRUE, na.last=TRUE)
#	m.order.pass2 = order(p.val.list.pass2, decreasing=FALSE, na.last=TRUE)
	model.descs2.pass2 = model.descs2.pass2[m.order.pass2]
#	loc.pass2 = m.order.pass2[1]
	m2.pass2 <- m2.pass1[m.order.pass2, ]
	model.names.pass2 <- rownames(m2.pass2)
#	print(c("loc.pass2:", loc.pass2))
	s.order.pass2 <- order(m2.pass2[1,], decreasing = TRUE)
	m2.pass2 <- m2.pass2[, s.order.pass2]
	bin.class.pass2 = bin.class.pass2[s.order.pass2]
#	bin.class.pass2.noisy = bin.class.pass2.noisy[s.order.pass2]
	sample.names2.pass2 <- colnames(m2.pass2)
	
	if (is.vector(cls.labels)) {
		cls.labels2.pass2 <- cls.labels2.pass2[s.order.pass2]
		cls.list2.pass2 <- cls.list2.pass2[s.order.pass2]
	} else {
		cls.labels2.pass2 <- cls.labels2.pass2[, s.order.pass2]
		cls.list2.pass2 <- cls.list2.pass2[, s.order.pass2]
	}
#	browser()
	
#	browser()
	m.score.pass2 <- m2.pass2[1,]
	m.score.norm.pass2 <- (m.score.pass2 - min(m.score.pass2))/(max(m.score.pass2) - min(m.score.pass2))
#	MI.ref.genes.pass2 = mutual.inf.2(m.score.norm.pass2, m.score.norm.pass2)$MI
#	print(paste("MI.ref.genes.pass2 =", MI.ref.genes.pass2))	
	MI.list.phen.pass2 = vector(mode="numeric", length=n.phen.pass2)
	phen.descs.pass2 = vector(mode="character", length=n.phen.pass2)
	
	if( length(unique(bin.class.pass2)) > 1){
		MI <- (mutual.inf.2(bin.class.pass2, m.score.norm.pass2)$MI)#/MI.ref.genes.pass2
		MI.signif <- signif(MI, digits=3)
		MI.list.phen.pass1[1] = MI
	} else{
		MI.signif <- "-"
		MI <- NA
		MI.list.phen.pass2[1] = MI
	}
	print(paste(format(phen.names.pass2[1], width=12), "mutual.inf =", MI.signif))
	phen.descs.pass2[1] <- paste(MI.signif)
	
	if( n.phen.pass2 == 2 ){
		phen.descs.pass2[2] <- paste(MI.signif)
		MI.list.phen.pass2[2] = MI
		g.order.pass2 = c(1,2)
	} else{
		for( i in 2:n.phen.pass2 ){ 
			bin.gene = ifelse( cls.list2.pass2[i,]=="WT", 0, 1)  
			if (length(unique(bin.gene)) > 1) {
				MI <- (mutual.inf.2(bin.gene, m.score.norm.pass2)$MI)#/MI.ref.genes.pass2
#			if( perf.auc$A < 0.5 ){
				##				browser()
#				roc = signif(1 - perf.auc$A, digits=3)
#				p.val = signif(1 - perf.auc$A, digits=3)
#				abnormality = unique(cls.list2.pass2[i,])[which(unique(cls.list2.pass2[i,]) != "WT")]
#				cls.list2.pass2 = ifelse( cls.list2.pass2[i,] == "WT", abnormality, "WT" )
#				phen.names.pass2[i] = paste(phen.names.pass2[i], "-opposite")
#			} else{
				MI.signif <- signif(MI, digits=3)
#			p.val <- signif(perf.auc$p.value, digits=3)
				MI.list.phen.pass2[i] = MI
#			}
			} else {
				MI <- "-"
#			p.val <- "-"
				MI.list.phen.pass2[i] = NA
			}
			print(paste(format(phen.names.pass2[i], width=12),"mutual.inf =", MI.signif))
#		p.val.list[i] = p.val
			phen.descs.pass2[i] <- paste(MI.signif)
		}
		
#	browser()
		g.order.pass2 = c(1, order(MI.list.phen.pass2[2:n.phen.pass2], decreasing=TRUE, na.last=TRUE)+1)  # skip PATHWAY.MUT
	}
	phen.descs2.pass2 = phen.descs.pass2[g.order.pass2]
	cls.list2.pass2 = cls.list2.pass2[g.order.pass2,]
	phen.names.pass2 = phen.names.pass2[g.order.pass2]
#	browser()
	# Recompute cls.list2 as some mutations or copy numbers may have been removed
	
	
	# Recompute cls.phen and cls.labels2 as order may have changed
	
	cls.phen2.pass2 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(as.vector(cls.list2.pass2))
		cls.phen2.pass2 <- classes
		cls.labels2.pass2 <- match(cls.list2.pass2, cls.phen2.pass2)
	} else {
		for (kk in 1:length(cls.list2.pass2[, 1])) {
			classes <- unique(cls.list2.pass2[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2.pass2 <- c(cls.phen2.pass2, classes)
			cls.labels2.pass2[kk,] <- match(cls.list2.pass2[kk,], classes)
		}
	}
	cls.labels2.pass2 = cls.labels2.pass2[1:n.phen.pass2,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list.pass2 <- unlist(cls.phen2.pass2)
	
#	colors.list <- ifelse(unlist(cls.phen2) == target.class, 
#			ifelse(unlist(cls.phen2) == "DEL" | unlist(cls.phen2) == "AMP", 
#					ifelse(unlist(cls.phen2) == "DEL", cls.phen.colors[3], cls.phen.colors[4]), cls.phen.colors[1]), cls.phen.colors[2])
	colors.list.pass2 = rep( "gray", length(phen.list.pass2))
	colors.list.pass2[phen.list.pass2=="MUT"] = cls.phen.colors[1]
	colors.list.pass2[phen.list.pass2=="DEL"] = cls.phen.colors[3]
	colors.list.pass2[phen.list.pass2=="AMP"] = cls.phen.colors[4]
	colors.list.pass2[phen.list.pass2=="ALT"] = cls.phen.colors[5]
	
	MSIG.HeatMapPlot.10(V = m2.pass2, 
			pathway.mut = bin.class.pass2,
			row.names = model.names.pass2,
			row.names2 = model.descs2.pass2, 
			col.labels = cls.labels2.pass2, 
			col.classes = cls.phen2.pass2, 
			phen.cmap = colors.list.pass2, phen.names = phen.names.pass2,
			phen.names2 = phen.descs2.pass2,
			col.names = sample.names2.pass2, main = paste(tissue, "- Phase 2: only mutual information >=", MI.thresholds[MI.i],"from Phase 1 (MI)"), 
			xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	### 3rd Pass ###	
	
#	print( "--- Begin Pass 3 ---")
	##	browser()
#	m2.pass3 = m2.pass2
#	model.names.pass3 = rownames(m2.pass3)
#	sample.names2.pass3 = colnames(m2.pass3)
	##	model.descs2.pass3 = model.descs2.pass2
#	n.phen.pass3 = 40
	##	phen.descs.pass2 = phen.descs2.pass1[1:n.phen.pass2]
#	cls.list2.pass3 = cls.list[, s.order.pass1][, s.order.pass2]
#	cls.labels2.pass3 = cls.labels[, s.order.pass1][, s.order.pass2]
#	
	##	browser()
#	phen.names.pass3 = phen.names
#	m.score.pass3 <- m2.pass3[1,]
#	m.score.norm.pass3 <- (m.score.pass3 - min(m.score.pass3))/(max(m.score.pass3) - min(m.score.pass3))
#	MI.list.phen.pass3 = vector(mode="numeric", length=n.phen)
#	phen.descs.pass3 = vector(mode="character", length=n.phen)
	##	p.val.list.phen.pass3 = vector(mode="numeric", length=n.phen)
#	for( i in 3:n.phen ){ 
#		bin.gene = ifelse( cls.list2.pass3[i,]=="WT", 0, 1)  
#		if (length(unique(bin.gene)) > 1) {
#			(MI <- mutual.inf.2(bin.gene #+ runif(Ns, min=-.01, max=.01)
#						, m.score.norm.pass3)$MI)/MI.ref
	##			if( perf.auc$A < 0.5 ){
#			##				browser()
	##				roc = signif(1 - perf.auc$A, digits=3)
	##				p.val = signif(1- perf.auc$p.value, digits=3)
	##				abnormality = unique(cls.list2.pass3[i,])[which(unique(cls.list2.pass3[i,]) != "WT")]
	##				cls.list2.pass3[i,] = ifelse( cls.list2.pass3[i,] == "WT", abnormality, "WT" )
	##				phen.names.pass3[i] = paste(phen.names.pass3[i], "-opposite", sep="")
	##				roc.list.phen.pass3[i] = 1-perf.auc$A
	##				p.val.list.phen.pass3[i] = 1- perf.auc$p.value
	##			} else{
#			MI.signif <- signif(MI, digits=3)
	##			p.val <- signif(perf.auc$p.value, digits=3)
#			MI.list.phen.pass3[i] = MI
	##			p.val.list.phen.pass3[i] = perf.auc$p.value
	##			}
#		} else {
#			MI.signif <- NA
	##			p.val <- NA
#			MI.list.phen.pass3[i] = NA
	##			p.val.list.phen.pass3[i] = NA
#		}
	##		print(paste("ROC=", roc, " p-val=", p.val)) 
#		
	##		p.val.list[i] = p.val
#		phen.descs.pass3[i] <- paste(MI.signif)
#	}
	##	browser()
	##	MI.threshold = 0.20
	##	len = length(which(MI.list.phen.pass3[-1:-2] >= MI.threshold))+2
	##	if( len>40 ) 
#	len=40
#	ind.u = match(order(unique(MI.list.phen.pass3[-1:-2]), decreasing=FALSE, na.last=TRUE), MI.list.phen.pass3[-1:-2])
	##	if( len == 2 ){
	##		MI.threshold = 0.15
	##		len = length(which(MI.list.phen.pass3[-1:-2] >= MI.threshold))+2
	##	}
	##	if( len == 2 ){
	##		MI.threshold = 0.2
	##		len = length(which(MI.list.phen.pass3[-1:-2] >= MI.threshold))+2
	##	}
	##	g.order.pass3.1 = c(1, 2, order(p.val.list.phen.pass3[3:n.phen], decreasing=FALSE, na.last=TRUE)+2 )
#	g.order.pass3 =  c(1, 2, order(MI.list.phen.pass3[-1:-2], decreasing=FALSE, na.last=TRUE)+2 )[1:len]  # skip PATHWAY.MUT and COPY.NUMBER
#	phen.descs2.pass3 = phen.descs.pass3[g.order.pass3]
#	cls.list2.pass3 = cls.list2.pass3[g.order.pass3,]
#	cls.labels2.pass3 = cls.labels2.pass3[g.order.pass3,]
#	phen.names.pass3 = phen.names.pass3[g.order.pass3]
#	
#	
#	cls.list.mut = ifelse(cls.list2.pass3[-1:-2,] == "MUT", 1, 0)
#	cls.list.alt = ifelse(cls.list2.pass3[-1:-2,] == "DEL" | cls.list2.pass3[-1:-2,] == "AMP", 1, 0)
#	
	##	browser()
#	if( !is.vector(cls.list.alt) ){
#		cls.list.mut.sum = apply(cls.list.mut, MARGIN=2, FUN=sum)
#		cls.list.alt.sum = apply(cls.list.alt, MARGIN=2, FUN=sum)	
#		bin.class.pass3 = cls.list.mut.sum + cls.list.alt.sum
#		bin.class.pass3 = ( bin.class.pass3 - min(bin.class.pass3))/(max(bin.class.pass3) - min(bin.class.pass3))
#		cls.list.mut.sum = ifelse(cls.list.mut.sum + cls.list.alt.sum > 0, 1, 0)
#		cls.list2.pass3[1,] = ifelse( cls.list.mut.sum >= 1, "MUT", "WT")
#		cls.list2.pass3[2,] = ifelse( cls.list.alt.sum >= 1, "ALT", "WT")
#		
#		
#	} else{
#		
#		cls.list2.pass3[2,] = ifelse(cls.list.alt == 1, "ALT", "WT")
#		bin.class.pass3 = cls.list.mut+cls.list.alt #+ runif(Ns, min=-.1, max=.1)
#		bin.class.pass3 = ( bin.class.pass3 - min(bin.class.pass3))/(max(bin.class.pass3) - min(bin.class.pass3))
#		cls.list2.pass3[1,] = ifelse(bin.class.pass3 > 0 , "MUT", "WT")
#	}
#	
	##	browser()
#	
#	if( length(unique(bin.class.pass3)) > 1){
#		MI <- mutual.inf.2(bin.class.pass3, m.score.norm.pass3)$MI
#		MI.signif <- signif(MI, digits=3)
#		MI.list.phen.pass3[1] = MI
#	} else{
#		MI.signif <- "-"
#		MI.list.phen.pass3[1] = NA
#	}
#	print(paste(format(phen.names.pass3[1], width=12), "mutual.inf =", MI.signif))
#	phen.descs2.pass3[1] <- paste(MI.signif)
#	for( i in 2 ){ # Recalculate MI for PATHWAY.MUT and COPY.NUMBER
#		bin.gene = ifelse( cls.list2.pass3[i,]=="WT", 0, 1)
#		if (length(unique(bin.gene)) > 1) {
#			MI <- (mutual.inf.2(bin.gene #+ runif(Ns, min=-.01, max=.01)
#						, m.score.norm.pass3)$MI)/MI.ref
	##			if( perf.auc$A < 0.5 ){
#			##				browser()
	##				roc = signif(1 - perf.auc$A, digits=3)
	##				p.val = signif(1- perf.auc$p.value, digits=3)
	##				abnormality = unique(cls.list2.pass3[i,])[which(unique(cls.list2.pass3[i,]) != "WT")]
	##				cls.list2.pass3[i,] = ifelse( cls.list2.pass3[i,] == "WT", abnormality, "WT" )
	##				phen.names.pass3[i] = paste(phen.names.pass3[i], "-opposite", sep="")
	##				roc.list.phen.pass3[i] = 1-perf.auc$A
	##				p.val.list.phen.pass3[i] = 1- perf.auc$p.value
	##			} else{
#			MI.signif <- signif(MI, digits=3)
	##			p.val <- signif(perf.auc$p.value, digits=3)
#			MI.list.phen.pass3[i] = MI
	##			p.val.list.phen.pass3[i] = perf.auc$p.value
	##			}
#		} else {
#			MI <- NA
	##			p.val <- NA
#			MI.list.phen.pass3[i] = NA
	##			p.val.list.phen.pass3[i] = NA
#		}
#		print(paste(format(phen.names.pass3[i], width=12), "mutual.inf =", MI.signif))
#		
	##		p.val.list[i] = p.val
#		phen.descs2.pass3[i] <- paste(MI.signif)
#	}
#	
	##	browser()
#	model.descs2.pass3 = vector(length = n.models, mode="character")
#	MI.list.pass3 = vector( length=n.models, mode="character")
#	for (i in 1:n.models) {
#		m.score <- m2.pass3[i,]
#		m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
	##		browser()
#		if (length(unique(bin.class.pass3)) > 1) {
#			MI <- mutual.inf.2(bin.class.pass3, m.score.norm)$MI
#			MI.signif <- signif(MI, digits=3)
	##			p.val <- signif(perf.auc$p.value, digits=3)
#			MI.list.pass3[i] = MI
	##			p.val.list.pass2[i] = perf.auc$p.value
#		} else {
#			MI.signif <- "-"
#			MI.list.pass3[i] = NA
	##			p.val.list.pass2[i] = NA
#		}
#		print(paste(format(rownames(m2.pass3)[i], width=30), "mutual.inf =", MI.signif)) 
#		
#		model.descs2.pass3[i] <- paste(MI.signif)
#	}
#	m.order.pass3 = order(MI.list.pass3, na.last=TRUE, decreasing=TRUE)
#	m2.pass3 = m2.pass3[m.order.pass3,]
#	model.descs2.pass3 = model.descs2.pass3[m.order.pass3]
#	s.order.pass3 = order(m2.pass3, decreasing=TRUE)
#	m2.pass3 = m2.pass3[,s.order.pass3]
#	bin.class.pass3 = bin.class.pass1[s.order.pass3]
#	model.names.pass3 = rownames(m2.pass3)
#	sample.names2.pass3 = colnames(m2.pass3)
#	
#	cls.phen2.pass3 <- NULL
#	if (is.vector(cls.labels)) {
#		classes <- unique(as.vector(cls.list2.pass3))
#		cls.phen2.pass3 <- classes
#		cls.labels2.pass3 <- match(cls.list2.pass3, cls.phen2.pass3)
#	} else {
	##		browser()
#		for (kk in 1:length(cls.list2.pass3[, 1])) {
	##			browser()
#			classes <- unique(cls.list2.pass3[kk,])
	##            cls.phen2[[kk]] <- classes
#			cls.phen2.pass3 <- c(cls.phen2.pass3, classes)
#			cls.labels2.pass3[kk,] <- match(cls.list2.pass3[kk,], classes)
#		}
#	}
	##	cls.labels2.pass3 = cls.labels2.pass3[1:n.phen.pass3,]
#	
#	
	##	browser()
	##	correl <- cor(t(m2))[, loc]
	##	m.order <- order(correl, decreasing=decreasing.order)
	##	correl2 <- correl[m.order]
#	
	##	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
#	phen.list.pass3 <- unlist(cls.phen2.pass3)
#	
	##	colors.list <- ifelse(unlist(cls.phen2) == target.class, 
	##			ifelse(unlist(cls.phen2) == "DEL" | unlist(cls.phen2) == "AMP", 
	##					ifelse(unlist(cls.phen2) == "DEL", cls.phen.colors[3], cls.phen.colors[4]), cls.phen.colors[1]), cls.phen.colors[2])
#	colors.list.pass3 = rep( "gray", length(phen.list.pass3))
#	colors.list.pass3[phen.list.pass3=="MUT"] = cls.phen.colors[1]
#	colors.list.pass3[phen.list.pass3=="DEL"] = cls.phen.colors[3]
#	colors.list.pass3[phen.list.pass3=="AMP"] = cls.phen.colors[4]
#	colors.list.pass3[phen.list.pass3=="ALT"] = cls.phen.colors[5]
#	phen.names.pass3[1] = "PATHWAY.MUT+COPY.NUMBER"
	##	browser()
#	MSIG.HeatMapPlot.10(V = m2.pass3, 
#			pathway.mut = bin.class.pass3,
#			row.names = model.names.pass3,
#			row.names2 = model.descs2.pass3, 
#			col.labels = cls.labels2.pass3, 
#			col.classes = cls.phen2.pass3, 
#			phen.cmap = colors.list.pass3, phen.names = phen.names.pass3,
#			phen.names2 = phen.descs2.pass3,
#			col.names = sample.names2.pass3, main = paste(tissue, "- Phase 3: Top signature from Phase 2 with all genes (Top 40 genes) (MI)"), 
#			xlab="  ", ylab="  ", row.norm = row.norm,  
#			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	dev.off()
	
	if (!is.na(output.dataset)) {
		V.GCT <- m2
		colnames(V.GCT) <- sample.names2
		row.names(V.GCT) <- model.names2
		write.gct(gct.data.frame = V.GCT, descs = model.descs2, filename =output.dataset)  
	}
	
}

OPAM.sort.projection.by.score.7 <- function(
#		input.ds,
		signatures = "NA",
		input.all.pathways.ds,
		input.cls,
		tissue = "NA",
		results.dir,
		normalize.score = T,
		normalization.type = "zero.one",
		model = "NA",
		target.phen = NA,
		target.class = NA,
		user.colors = NA,
		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T,
		u.gene.names.known = "NA",
		add.amp.del = FALSE,
		#n.random.signatures = 10,
		multiple.tissues = FALSE,
		cls.has.chrom.locs = FALSE,
		file.suffix = "",
		skip.iterative = FALSE,
		add.mut = FALSE,
		n.iter = 5,
		pdf.height = 11,
		pdf.width = 17,
		do.mRMR = F,
		skip.step2 = FALSE,
		todd.version = FALSE
)
# Calls MSIG.HeatMapPlot.9 and makes a plot sorted by the highest-scoring
# signatures and abnormalities (gene mutations or copy number alterations)
# i.e. doesn't require a "model" to score by as OPAM.sort.projection.by.score.2 does.
# However, it *will* use "model" if it cannot calculate p-values on the gene signatures, which
# happens when every cell line has a genomic aberration.
#
# Runs 3 passes on the data:
# 1st pass: looks at the genes and copy number alterations specified by u.gene.names.known
# 2nd pass: looks at only the top abnormalities (using a p-value cutoff) from the 1st pass, and adjusts 
# the PATHWAY.MUT vector accordingly (only according to the genes, not by copy number data)
# 3rd pass: Takes the winning signature from the 2nd pass and then looks all the genes available
#
# Very similar to OPAM.sort.projection.by.score.4, however this version uses mutual.inf instead of
# roc.area to calculate mutual information scores and p-values for PATHWAY.MUT, the vector of total genomic aberrations
# in all samples
#
# Differs from OPAM.sort.projection.by.score.6 by requiring the gct file of expression in 
# all pathways by the input tissue ("input.all.pathways.ds")
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
	dataset.all <- MSIG.Gct2Frame( filename = input.all.pathways.ds)
	m.all <- data.matrix(dataset.all$ds)
	model.names.all <- dataset.all$row.names
	Ns = length(m.all[1,])
	sample.names = dataset.all$names
	
	if( multiple.tissues ){
		tissue.type <- vector(length=Ns, mode="character")
#		temp = strsplit(sample.names, split="_")
		for (k in 1:Ns) {
			temp <- strsplit(sample.names[k], split="_") 
			tissue.type[k] <- paste(temp[[1]][2:length(temp[[1]])], collapse="_")
		}
		tissue.names = unique(tissue.type)
		tissue.labels = match(tissue.type, tissue.names)
	} else{
		tissue.names = tissue
		tissue.labels = rep(1, Ns)
	}
	
	if( is.na(signatures[1]) ){
		stop("Must provide a vector of signature names to evaluate, or specify 'ALL'")
	}
	
	## Remove "Commented out" signatures (with # at beginning of name)
	if( length(grep("^#", signatures)) > 0){
		signatures = signatures[-grep("^#", signatures)]
	}
	if( signatures[1] == "ALL"){
		model.names = model.names.all
		m = m.all
		model.descs = dataset.all$descs
	} else{
		model.names = signatures
		model.ind = match(signatures, model.names.all)
		m = m.all[model.ind,]
		model.descs = dataset.all$descs[model.ind]
		if( length(model.ind) == 1 ){
			m = t(as.matrix(m))
			rownames(m) = model.names
		}
		rm(list=c("m.all", "dataset.all"))
#	browser()
	}
	n.models <- length(m[,1])
	temp <- strsplit(input.all.pathways.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
	char.res <-  0.013 * n.models + 0.65
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:n.models) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list 
	
	if(cls.list[1,1] == "0" | cls.list[1,1] == "1"){
		cls.list = ifelse(cls.list=="1", "MUT", "WT")
	}
	
	#browser()
	
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors[1])) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
	
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
	} else {
		phen.names <- "NA"
	}
	
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
#	print("cls.phen.colors:")
#	print(cls.phen.colors)
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		n.classes = unlist(lapply(cls.phen, length))
#		for (i in 1:n.phen) {
#			n.classes[i] <- length(cls.phen[[i]])
#		}
	}
	pdf.options(height=pdf.height, width=pdf.width, colormodel="rgb", bg="transparent")
	phen.names = c("SUMMARY", phen.names)
	cls.list = rbind(rep("WT", length=length(cls.list[1,])), 
			#rep("WT", length=length(cls.list[1,])), 
			cls.list)
	cls.labels = rbind(rep(1, length=length(cls.labels[1,])), 
			#rep(1, length=length(cls.labels[1,])), 
			cls.labels)
	if( !todd.version ){
		print("--- Begin Pass 1 ---")
#	browser()
		
		## Remove "Commented out" gene names (with # at beginning of name)
		if( length(grep("^#", u.gene.names.known)) > 0){
			u.gene.names.known = u.gene.names.known[-grep("^#", u.gene.names.known)]
		}
		if ( add.amp.del ){
			u.gene.names.known = c( u.gene.names.known, paste(u.gene.names.known, "_AMP", sep=""), 
					paste(u.gene.names.known, "_DEL", sep="") )
		}
		if(add.mut){
			u.gene.names.known = paste(u.gene.names.known, "_MUT", sep="")
		}
		
		phen.pass1 = c( u.gene.names.known )
		
		
		## Find chromosomal locations of genes specified
		## See "if( find.chromosomal.locations)" for more transparent code
		if( cls.has.chrom.locs ){
			library(org.Hs.eg.db)
			phen.pass1.split = strsplit(phen.pass1, split="_")
			phen.pass1.noampdel = unlist( lapply(phen.pass1.split, function(x) x[1]))
			phen.pass1.egIDs = mget(phen.pass1.noampdel, org.Hs.egALIAS2EG, ifnotfound=NA)
#	phen.pass1.no.chrom.loc = which(lapply(lapply(phen.pass1.egIDs, is.na), sum) > 0)
			phen.pass1.locs.list = lapply(phen.pass1.egIDs, mget, org.Hs.egMAP)
			phen.pass1.locs = vector(mode="character", length=length(phen.pass1))
#	phen.pass1.locs[phen.pass1.no.chrom.loc] = "NA"
			phen.pass1.locs = unlist(lapply(phen.pass1.locs.list, function(x) paste(unlist(x), collapse="_")))
			phen.pass1.w.locs = paste(phen.pass1.noampdel, ".", phen.pass1.locs, sep="")
			phen.pass1.ampdel.suffix = unlist(lapply(phen.pass1.split, function(x) x[2]))
			#phen.pass1.no.suffix = which(is.na(phen.pass1.ampdel.suffix))
			#phen.pass1[phen.pass1.no.suffix] = phen.pass1.w.locs[phen.pass1.no.suffix]
			#phen.pass1[-phen.pass1.no.suffix] = paste(phen.pass1.w.locs[-phen.pass1.no.suffix], "_", phen.pass1.ampdel.suffix[-phen.pass1.no.suffix], sep="")
			phen.pass1 = paste(phen.pass1.w.locs, "_", phen.pass1.ampdel.suffix, sep="")
		}
		
		
		## Was originally immediately after "ind.phen.pass1 = ..." but since now want to find the chromosomal
		## locations of the genes, have to first find the indices of the genes specified at the onset of the
		## program, THEN find all the chromosomal locations
		
		#browser()
		
		
		phen.pass1 = c("SUMMARY", phen.pass1)
		ind.phen.pass1 = which( phen.names %in% phen.pass1 )
#	phen.names[1] = "SUMMARY"
		phen.pass1 = phen.names[ind.phen.pass1]
		
		n.phen.pass1 = length(phen.pass1)
		
		MI.list.pass1 = vector( length=n.models, mode="numeric" )
		
		
		
		cls.list.pass1 = cls.list[ind.phen.pass1,]
		cls.list.pass1.2 = ifelse(cls.list.pass1 == "WT", 0, 1)
		cls.labels.pass1 = cls.labels[ind.phen.pass1,]
		
#	browser()
		if (!is.na(target.phen)) {
			if( length( phen.pass1) > 2 ){
				bin.class.pass1 = apply( cls.list.pass1.2[-1,], MARGIN=2, FUN=sum)
			} else{ bin.class.pass1 = cls.list.pass1.2[2,] }
			# Normalize bin.class.pass1
			if( length(unique(bin.class.pass1)) > 1){
				bin.class.pass1 = normalize(bin.class.pass1) #( bin.class.pass1 - min(bin.class.pass1))/(max(bin.class.pass1) - min(bin.class.pass1))
			} else if ( length(unique(bin.class.pass1)) == 1){
				bin.class = rep(1, length(cls.list[1,]))
			}
			cls.list.pass1[1,] = ifelse(bin.class.pass1 > 0, "MUT", "WT")
		} else {
			bin.class.pass1 <- ifelse(cls.list[1,] == cls.list2[1,1], 1, 0)
		}
		#browser()
		
		### Make initial heatmap ###
		cls.phen2.pass1 <- NULL
		if (is.vector(cls.labels)) {
			classes <- unique(cls.list.pass1)
			cls.phen2.pass1 <- classes
			cls.labels.pass1 <- match(cls.list.pass1, cls.phen2.pass1)
		} else {
			for (kk in 1:length(cls.list.pass1[, 1])) {
				classes <- unique(cls.list.pass1[kk,])
#            cls.phen2[[kk]] <- classes
				cls.phen2.pass1 <- c(cls.phen2.pass1, classes)
				cls.labels.pass1[kk,] <- match(cls.list.pass1[kk,], classes)
			}
		}
		cls.labels.pass1 = cls.labels.pass1[1:length(phen.pass1),]
		
		phen.list.pass1 <- unlist(cls.phen2.pass1)
		
		colors.list = rep( "gray", length(phen.list.pass1))
		colors.list[phen.list.pass1=="MUT"] = cls.phen.colors[1]
		filename <- paste(results.dir, test.file.prefix, file.suffix, ".Step0", sep="")
		pdf(file=paste(filename, ".pdf", sep=""), height = pdf.height, width = pdf.width )
#quartz(height = 11, width = 17)
		if( multiple.tissues ){
			#browser()
			#quartz(height = 11, width = 17)
			MSIG.HeatMapPlot.10.multiple.tissues(V = m, 
					pathway.mut = bin.class.pass1,
					row.names = model.names,
					col.labels = cls.labels.pass1, 
					col.classes = cls.phen2.pass1, 
					phen.cmap = colors.list, 
					phen.names = phen.pass1,
					col.names = sample.names, 
					main = paste(tissue, "- Initial Heatmap ('Step 0')"), 
					xlab="  ", ylab="  ", row.norm = row.norm,  
					cmap.type = cmap.type, char.rescale = char.rescale,  legend=F,
					tissue.names = tissue.names,
					tissue.labels = tissue.labels)
		} else{
			MSIG.HeatMapPlot.10(V = m, 
					pathway.mut = bin.class.pass1,
					row.names = model.names,
					col.labels = cls.labels.pass1, 
					col.classes = cls.phen2.pass1, 
					phen.cmap = colors.list, 
					phen.names = phen.pass1,
					col.names = sample.names, 
					main = paste(tissue, "- Initial Heatmap ('Step 0')"), 
					xlab="  ", ylab="  ", row.norm = row.norm,  
					cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
		}
		dev.off()
		
		
		model.descs2.pass1 = vector(length = n.models, mode="character")
		
		if( length(unique(bin.class.pass1)) > 1 ){
			
			if( n.models > 1 ){
				MI.results = mutual.inf.3.v2(bin.class.pass1, m, 
						target.vector.name="SUMMARY", 
						tissue=tissue)
				MI.list.pass1  = MI.results$MI
				model.descs2.pass1 <- sapply(MI.results$MI, FUN=signif, 3)
				m.order.pass1 = order(MI.list.pass1, decreasing=TRUE, na.last=TRUE)
				m2.pass1 <- m[m.order.pass1, ]
				s.order.pass1 <- order(m2.pass1[1,], decreasing = TRUE )
				m2.pass1 <- m2.pass1[, s.order.pass1]
			} else{ 
				#browser()
				MI.ref = mutual.inf.2(bin.class.pass1, bin.class.pass1)
				MI.list.pass1 = MI.results = 
						mutual.inf.2(bin.class.pass1, m[1,])/MI.ref
				model.descs2.pass1 <- signif(MI.results, digits=3)
				m2.pass1 <- m ; m.order.pass1 = 1
				s.order.pass1 <- order(m2.pass1[1,], decreasing = TRUE )
				m2.pass1 <- t(as.matrix(m[, s.order.pass1]))
				rownames(m2.pass1) = model.names
			}
		} else{ 
			MI.list.pass1 = rep(NA, n.models)
			FDR.list.pass1 = rep(NA, n.models)
			model.descs2.pass1 = rep(" - (FDR = - )", n.models)
			if( n.models > 1 ){
				loc <- match(model, model.names)
				s.order.pass1 <- order(m[loc,], decreasing = decreasing.order)
				m2.pass1 <- m[, s.order.pass1]
				correl <- cor(t(m2.pass1))[, loc]
				m.order.pass1 <- order(correl, decreasing=T)
				m2.pass1 <- m2.pass1[m.order.pass1, ]
			} else{ 
				m2.pass1 <- t(as.matrix(m[, s.order.pass1]))
				rownames(m2.pass1) = model.names
				m.order.pass1 = 1 }
		}
		
		MI.list.pass1 = MI.list.pass1[m.order.pass1]
		bin.class.pass1 = bin.class.pass1[s.order.pass1]
		
		model.descs2.pass1.all = model.descs2.pass1
		model.descs2.pass1 = model.descs2.pass1[m.order.pass1]
		sample.names2.pass1 <- colnames(m2.pass1)
		model.names.pass1 <- rownames(m2.pass1)
		print(matrix(c(model.names.pass1, model.descs2.pass1), ncol=2), quote=F)
		if (is.vector(cls.labels)) {
			cls.labels2.pass1 <- cls.labels.pass1[s.order.pass1]
			cls.list2.pass1 <- cls.list.pass1[s.order.pass1]
		} else {
			cls.labels2.pass1 <- cls.labels.pass1[, s.order.pass1]
			cls.list2.pass1 <- cls.list.pass1[, s.order.pass1]
		}
		tissue.labels.pass1 = tissue.labels[s.order.pass1]
		sample.names2 <- colnames(m2.pass1)
		winning.model.ind.pass1 = which(model.names.pass1[1] == rownames(m2.pass1))
		
		MI.list.phen.pass1 = vector(mode="numeric", length=n.phen.pass1)
		phen.descs.pass1 = vector(mode="character", length=n.phen.pass1)
		
		if( length(unique(bin.class.pass1)) > 1){
			MI.signif <- signif(MI.list.pass1[1], digits=3)
			MI.list.phen.pass1[1] = MI.list.pass1[1]
			phen.descs.pass1[1] = model.descs2.pass1[1]
		} else{
			MI.signif <- "-"
			MI.list.phen.pass1[1] = NA
		}
		print(paste(format(phen.pass1[1], width=12), "mutual.inf =", MI.signif
				))
		print(proc.time()-t1)
		print(date())
		phen.descs.pass1[1] <- paste(MI.signif,
				sep="")
		
#	browser()
		if( n.phen.pass1 > 2 ){
			bin.gene.matrix = ifelse(cls.list2.pass1[-1,]=="WT", 0, 1)
			MI.results = mutual.inf.3.v2(
					m2.pass1[winning.model.ind.pass1,],
					bin.gene.matrix)
			MI.list.phen.pass1[-1] = MI.results$MI
			phen.descs.pass1[-1] = sapply(MI.results$MI, FUN=signif, 3)
			
			g.order.pass1 = c(1, order(MI.list.phen.pass1[-1], decreasing=TRUE, na.last=TRUE)+1)
			MI.list.phen.pass1 = MI.list.phen.pass1[g.order.pass1]
			phen.descs2.pass1 = phen.descs.pass1[g.order.pass1]
			cls.list2.pass1 = cls.list2.pass1[g.order.pass1,]
			phen.names.pass1 = phen.pass1[g.order.pass1]
		} else{
#		bin.gene.matrix = ifelse(cls.list2.pass1[-1,]=="WT", 0, 1)
#		MI.ref = mutual.inf.2(m2.pass1[winning.model.ind.pass1,],
#				m2.pass1[winning.model.ind.pass1,])
#		MI.results = mutual.inf.2(
#				m2.pass1[winning.model.ind.pass1,],
#				bin.gene.matrix)/MI.ref
#		MI.list.phen.pass1[-1] = MI.results
			phen.descs.pass1[-1] = phen.descs.pass1[1]
			
			g.order.pass1 = 1:2
#	MI.list.phen.pass1 = MI.list.phen.pass1[g.order.pass1]
			phen.descs2.pass1 = phen.descs.pass1
#	cls.list2.pass1 = cls.list2.pass1[g.order.pass1,]
			phen.names.pass1 = phen.pass1
		}
		print(matrix(c(phen.names.pass1, phen.descs2.pass1), ncol=2), quote=F)
		print(proc.time()-t1)
		print(date())
		
		# Recompute cls.list2 as some mutations or copy numbers may have been removed
		# Recompute cls.phen and cls.labels2 as order may have changed
		cls.phen2.pass1 <- NULL
		if (is.vector(cls.labels)) {
			classes <- unique(cls.list2.pass1)
			cls.phen2.pass1 <- classes
			cls.labels2.pass1 <- match(cls.list2.pass1, cls.phen2.pass1)
		} else {
			for (kk in 1:length(cls.list2.pass1[, 1])) {
				classes <- unique(cls.list2.pass1[kk,])
#            cls.phen2[[kk]] <- classes
				cls.phen2.pass1 <- c(cls.phen2.pass1, classes)
				cls.labels2.pass1[kk,] <- match(cls.list2.pass1[kk,], classes)
			}
		}
		cls.labels2.pass1 = cls.labels2.pass1[1:n.phen.pass1,]
		
		phen.list.pass1 <- unlist(cls.phen2.pass1)
		colors.list.pass1 = rep( "gray", length(phen.list.pass1))
		colors.list.pass1[phen.list.pass1=="MUT"] = cls.phen.colors[1]
		colors.list.pass1[phen.list.pass1=="DEL"] = cls.phen.colors[3]
		colors.list.pass1[phen.list.pass1=="AMP"] = cls.phen.colors[4]
		colors.list.pass1[phen.list.pass1=="ALT"] = cls.phen.colors[5]
		
		filename <- paste(results.dir, test.file.prefix, file.suffix, ".Step1.MI-HXY", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
		pdf(file=paste(filename, ".pdf", sep=""), height = pdf.height, width = pdf.width )
		#browser()
#   windows(width=12, height=8)
		if( multiple.tissues ){
			MSIG.HeatMapPlot.10.multiple.tissues(V = m2.pass1, 
					pathway.mut = bin.class.pass1,
					row.names = model.names.pass1,
					row.names2 = model.descs2.pass1, 
					col.labels = cls.labels2.pass1, 
					col.classes = cls.phen2.pass1, 
					phen.cmap = colors.list.pass1, 
					phen.names = phen.names.pass1,
					phen.names2 = c(" ", phen.descs2.pass1),
					col.names = sample.names2.pass1, 
					main = paste(tissue, "- Step 1: Known KRAS Pathway Abnormalities (MI)"), 
					xlab="  ", ylab="  ", row.norm = row.norm,  
					cmap.type = cmap.type, char.rescale = char.rescale,  legend=F,
					tissue.names = tissue.names,
					tissue.labels = tissue.labels)
		} else{
			MSIG.HeatMapPlot.10(V = m2.pass1, 
					pathway.mut = bin.class.pass1,
					row.names = model.names.pass1,
					row.names2 = model.descs2.pass1, 
					col.labels = cls.labels2.pass1, 
					col.classes = cls.phen2.pass1, 
					phen.cmap = colors.list.pass1, 
					phen.names = phen.names.pass1,
					phen.names2 = phen.descs2.pass1,
					col.names = sample.names2.pass1, 
					main = paste(tissue, "- Step 1: Known KRAS Pathway Abnormalities (MI)"), 
					xlab="  ", ylab="  ", row.norm = row.norm,  
					cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
		}
		dev.off()
#	stop("Don't do Step 2!")
		### Begin Pass 2 ###
		if( n.phen.pass1 == 2 || skip.step2 == T )
			print( "--- Begin Step 2 ---")
		if( is.na(MI.list.pass1[1]) || is.na(MI.list.phen.pass1[1]) ){
			dev.off()
			return()
		}
		MI.thresholds = c(0.2, 0.1, 0.08, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01, 0)
		#	MI.threshold = 0.03
		ind.top.MI = vector(mode="integer")
		MI.i = 0
#	FDR.i = 0
		#	browser()
		while( length(ind.top.MI) < 1)
		{
			MI.i = MI.i + 1
#		FDR.i = FDR.i + 1
			if( MI.i > length(MI.thresholds)){
				dev.off()
				print("Selected genomic aberrations do not have 
								positive mutual information with a low enough false 
								discovery rate with the selected pathways")
				return()
			}
			ind.top.MI = which( MI.list.phen.pass1[-1] >= MI.thresholds[MI.i]) +1 #& MI.list.phen.pass1[-1] > 0 
#		) + 1
			
		}
		ind.MI.threshold = c(1, ind.top.MI)
		
		
		n.phen.pass2 = length(ind.MI.threshold)
		cls.list2.pass2 = cls.list2.pass1[ind.MI.threshold,]
		phen.names.pass2 = phen.names.pass1[ind.MI.threshold]
		cls.labels2.pass2 = cls.labels2.pass1[ind.MI.threshold,]
		
		#	browser()
		cls.list2.pass2.2 = ifelse( cls.list2.pass2 == "WT", 0, 1)
		cls.list2.pass2.3 = ifelse( cls.list2.pass2 == "DEL" | cls.list2.pass2 == "AMP", 1, 0)
		if( n.phen.pass2 > 2 ){ 
			pathway.mut.pass2 = apply(cls.list2.pass2.2[2:n.phen.pass2,], MARGIN=2, FUN=sum)
			bin.class.pass2 = pathway.mut.pass2/length(pathway.mut.pass2)
			bin.class.pass2 = ( bin.class.pass2 - min(bin.class.pass2))/(max(bin.class.pass2) - min(bin.class.pass2))
			cls.list2.pass2[1,] = ifelse( bin.class.pass2 > 0, "MUT", "WT")
		} else{
			pathway.mut.pass2 = ifelse( cls.list2.pass2.2[2,] == 1, "MUT", "WT")
			bin.class.pass2 = ifelse( pathway.mut.pass2 == "MUT", 1, 0 ) #+ runif(Ns, min=-.05, max=.05)
			cls.list2.pass2[1,] = pathway.mut.pass2
		}
		
		MI.list.pass2 = vector( length=n.models, mode="double" )
		#browser()
		#### Print Step 2's initial heatmap ###
		cls.phen2.pass1.5 <- NULL
		cls.labels2.pass1.5 = cls.labels2.pass2
		if (is.vector(cls.labels)) {
			classes <- unique(cls.list2.pass1)
			cls.phen2.pass1.5 <- classes
			cls.labels2.pass1.5 <- match(cls.list2.pass1, cls.phen2.pass1.5)
		} else {
			for (kk in 1:length(cls.list2.pass2[, 1])) {
				classes <- unique(cls.list2.pass1[kk,])
#            cls.phen2[[kk]] <- classes
				cls.phen2.pass1.5 <- c(cls.phen2.pass1.5, classes)
				cls.labels2.pass1.5[kk,] <- match(cls.list2.pass2[kk,], classes)
			}
		}
#	cls.labels2.pass1.5 = cls.labels2.pass1.5[1:n.phen.pass2,]
		
		phen.list.pass1.5 <- unlist(cls.phen2.pass1.5)
		colors.list.pass1.5 = rep( "gray", n.phen.pass1)
		colors.list.pass1.5[phen.list.pass1.5=="MUT"] = cls.phen.colors[1]
		filename <- paste(results.dir, test.file.prefix, file.suffix, ".Step1.5.Heatmap", sep="")
		pdf(file=paste(filename, ".pdf", sep=""), height = pdf.height, width = pdf.width )
		if( multiple.tissues ){
			MSIG.HeatMapPlot.10.multiple.tissues(V = m2.pass1, 
					pathway.mut = bin.class.pass2,
					row.names = model.names.pass1,
					col.labels = cls.labels2.pass1.5, 
					col.classes = cls.phen2.pass1.5, 
					phen.cmap = colors.list.pass1.5, 
					phen.names = phen.names.pass1[ind.MI.threshold],
					col.names = sample.names2.pass1, 
					main = paste(tissue, "- Post-Step 1, Pre-Step 2 (Step 1.5)"), 
					xlab="  ", ylab="  ", row.norm = row.norm,  
					cmap.type = cmap.type, char.rescale = char.rescale,  legend=F,
					tissue.names = tissue.names,
					tissue.labels = tissue.labels)
		} else{
			MSIG.HeatMapPlot.10(V = m2.pass1, 
					pathway.mut = bin.class.pass2,
					row.names = model.names.pass1,
					col.labels = cls.labels2.pass1.5, 
					col.classes = cls.phen2.pass1.5, 
					phen.cmap = colors.list.pass1.5, 
					phen.names = phen.names.pass1[ind.MI.threshold],
					col.names = sample.names2.pass1, 
					main = paste(tissue, "- Post-Step 1, Pre-Step 2 (Step 1.5)"), 
					xlab="  ", ylab="  ", row.norm = row.norm,  
					cmap.type = cmap.type, char.rescale = char.rescale,  legend=F) }
		dev.off()
		
		
		#	browser()
		#pdf(file=paste(tissue, n.randomizations, "randomizations.Step2", "pdf", sep="."))
		model.descs2.pass2 = vector(length = n.models, mode="character")
		if( length(unique(bin.class.pass2)) > 1 ){
			
			if( n.models > 1 ){
				MI.results = mutual.inf.3.v2(bin.class.pass2, m2.pass1, 
						target.vector.name="SUMMARY", 
				)
				MI.list.pass2  = MI.results$MI
				model.descs2.pass2 = sapply(MI.results$MI, signif, 3)
				m.order.pass2 = order(MI.list.pass2, decreasing=TRUE, na.last=TRUE)
				m2.pass2 <- m2.pass1[m.order.pass2, ]
				s.order.pass2 <- order(m2.pass2[1,], decreasing = TRUE )
				m2.pass2 <- m2.pass2[, s.order.pass2]
			} else{
				MI.ref = mutual.inf.2(bin.class.pass2, bin.class.pass2)
				MI.list.pass2 = MI.results = 
						mutual.inf.2(bin.class.pass2, 
								m2.pass1[1,])/MI.ref
				#MI.list.pass2  = MI.results$MI
				model.descs2.pass2 = signif(MI.list.pass2, digits=3)
				m.order.pass2 = 1 #order(MI.list.pass2, decreasing=TRUE, na.last=TRUE)
				m2.pass2 <- m2.pass1#[m.order.pass2, ]
				s.order.pass2 <- order(m2.pass2[1,], decreasing = TRUE )
				m2.pass2 <- t(as.matrix(m2.pass2[, s.order.pass2]))
				rownames(m2.pass2) = model.names.pass1
			}
		} else{ 
			MI.list.pass2 = rep(NA, n.models)
			model.descs2.pass2 = rep(" - ", n.models)
			if( n.models > 1 ){
				loc <- match(model, model.names)
				s.order.pass2 <- order(m2.pass1[loc,], decreasing = decreasing.order)
				m2.pass2 <- m2.pass1[, s.order.pass2]
				correl <- cor(t(m2.pass2))[, loc]
				m.order.pass2 <- order(correl, decreasing=T)
				m2.pass2 <- m2.pass2[m.order.pass2, ]
			} else{
				#loc <- match(model, model.names)
				s.order.pass2 <- order(m2.pass1[1,], decreasing = decreasing.order)
				m2.pass2 <- t(as.matrix(m2.pass1[, s.order.pass2]))
				rownames(m2.pass2) = model.names.pass1
				m.order.pass2 = 1
#		correl <- cor(t(m2.pass2))[, loc]
#		m.order.pass2 <- order(correl, decreasing=T)
#		m2.pass2 <- m2.pass2[m.order.pass2, ]
			}
		}
		
		MI.list.pass2 = MI.list.pass2[m.order.pass2]
		bin.class.pass2 = bin.class.pass2[s.order.pass2]
		tissue.labels = tissue.labels[s.order.pass2]
		model.descs2.pass2 = model.descs2.pass2[m.order.pass2]
		sample.names2.pass2 <- colnames(m2.pass2)
		model.names.pass2 <- rownames(m2.pass2)
		print(matrix(c(model.names.pass2, model.descs2.pass2), ncol=2), quote=F)
#	browser()
		if (is.vector(cls.labels)) {
			cls.labels2.pass2 <- cls.labels2.pass2[s.order.pass2]
			cls.list2.pass2 <- cls.list2.pass2[s.order.pass2]
		} else {
			cls.labels2.pass2 <- cls.labels2.pass2[, s.order.pass2]
			cls.list2.pass2 <- cls.list2.pass2[, s.order.pass2]
		}
		
		tissue.labels.pass2 = tissue.labels.pass1[s.order.pass2]
		sample.names2 <- colnames(m2.pass2)
		
		winning.model.ind.pass2 = which(model.names.pass2[1] == rownames(m2.pass2))
		MI.list.phen.pass2 = vector(mode="numeric", length=n.phen.pass2)
		phen.descs.pass2 = vector(mode="character", length=n.phen.pass2)
		
		if( length(unique(bin.class.pass2)) > 1){
			MI.signif <- signif(MI.list.pass2[1], digits=3)
			MI.list.phen.pass2[1] = MI.list.pass2[1]
		} else{
			MI.signif <- "-"
			MI.list.phen.pass2[1] = NA
		}
		print(paste(format(phen.names.pass2[1], width=12), "mutual.inf =", MI.signif#, "  FDR =", FDR.signif
				))
		print(proc.time()-t1)
		print(date())
		phen.descs.pass2[1] <- paste(MI.signif) #, " (FDR = ", FDR.signif, ")", sep="")
		
		if( n.phen.pass2 == 2 ){
			phen.descs.pass2[2] <- phen.descs.pass2[1] 
#		FDR.list.phen.pass2[2] = FDR.list.phen.pass2[1] 
			MI.list.phen.pass2[2] = MI.list.phen.pass2[1] 
			g.order.pass2 = c(1,2)
			print(paste(format(phen.names.pass2[2], width=12), "mutual.inf =", MI.signif)) #, "  FDR =", FDR.signif))
		} else{
			bin.gene.matrix = ifelse(cls.list2.pass2[-1,]=="WT", 0, 1)
			n.aberrations = apply(bin.gene.matrix, MARGIN=1, FUN=sum)
			
			MI.results = mutual.inf.3.v2( 
					m2.pass2[winning.model.ind.pass2,],
					bin.gene.matrix,
					target.vector.name=phen.pass2,
#				n.randomizations = n.randomizations
			)
			
			MI.list.phen.pass2[-1] = MI.results$MI
			phen.descs.pass2[-1] = sapply(MI.results$MI, signif, 3)
			ind.zeros = which(n.aberrations==0) + 1
			MI.list.phen.pass2[ind.zeros] = NA
#		FDR.list.phen.pass2[ind.zeros] = NA
			phen.descs.pass2[ind.zeros] = " - "
			g.order.pass2 = c(1, order(MI.list.phen.pass2[-1], decreasing=TRUE, na.last=TRUE)+1)  # skip PATHWAY.MUT
		}
		#dev.off()
		phen.descs2.pass2 = phen.descs.pass2[g.order.pass2]
		cls.list2.pass2 = cls.list2.pass2[g.order.pass2,]
		phen.names.pass2 = phen.names.pass2[g.order.pass2]
		#	browser()
		# Recompute cls.list2 as some mutations or copy numbers may have been removed
		print(matrix(c(phen.names.pass2, phen.descs2.pass2), ncol=2), quote=F)
		
		
		# Recompute cls.phen and cls.labels2 as order may have changed
		
		cls.phen2.pass2 <- NULL
		if (is.vector(cls.labels)) {
			classes <- unique(as.vector(cls.list2.pass2))
			cls.phen2.pass2 <- classes
			cls.labels2.pass2 <- match(cls.list2.pass2, cls.phen2.pass2)
		} else {
			for (kk in 1:length(cls.list2.pass2[, 1])) {
				classes <- unique(cls.list2.pass2[kk,])
				#            cls.phen2[[kk]] <- classes
				cls.phen2.pass2 <- c(cls.phen2.pass2, classes)
				cls.labels2.pass2[kk,] <- match(cls.list2.pass2[kk,], classes)
			}
		}
		cls.labels2.pass2 = cls.labels2.pass2[1:n.phen.pass2,]
		
		phen.list.pass2 <- unlist(cls.phen2.pass2)
		colors.list.pass2 = rep( "gray", length(phen.list.pass2))
		colors.list.pass2[phen.list.pass2=="MUT"] = cls.phen.colors[1]
		colors.list.pass2[phen.list.pass2=="DEL"] = cls.phen.colors[3]
		colors.list.pass2[phen.list.pass2=="AMP"] = cls.phen.colors[4]
		colors.list.pass2[phen.list.pass2=="ALT"] = cls.phen.colors[5]
		
		filename <- paste(results.dir, test.file.prefix, file.suffix, ".Step2.MI-HXY", sep="")
		pdf(file=paste(filename, ".pdf", sep=""), height = pdf.height, width = pdf.width )
		if( multiple.tissues ){
			MSIG.HeatMapPlot.10.multiple.tissues(V = m2.pass2, 
					pathway.mut = bin.class.pass2,
					row.names = model.names.pass2,
					row.names2 = model.descs2.pass2, 
					col.labels = cls.labels2.pass2, 
					col.classes = cls.phen2.pass2, 
					phen.cmap = colors.list.pass2, phen.names = phen.names.pass2,
					phen.names2 = c(" ", phen.descs2.pass2),
					col.names = sample.names2.pass2, main = paste(tissue, 
							"- Step 2: only MI >=", MI.thresholds[MI.i],"from Step 1 (MI)"), 
					xlab="  ", ylab="  ", row.norm = row.norm,  
					cmap.type = cmap.type, char.rescale = char.rescale,  legend=F,
					tissue.names = tissue.names,
					tissue.labels = tissue.labels)
		} else{
			MSIG.HeatMapPlot.10(V = m2.pass2, 
					pathway.mut = bin.class.pass2,
					row.names = model.names.pass2,
					row.names2 = model.descs2.pass2, 
					col.labels = cls.labels2.pass2, 
					col.classes = cls.phen2.pass2, 
					phen.cmap = colors.list.pass2, phen.names = phen.names.pass2,
					phen.names2 = phen.descs2.pass2,
					col.names = sample.names2.pass2, main = paste(tissue, 
							"- Step 2: only MI >=", MI.thresholds[MI.i],"from Step 1 (MI)"), 
					xlab="  ", ylab="  ", row.norm = row.norm,  
					cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
		}
		dev.off()
	} else{ print("'todd.version' on -- skipping Steps 1 and 2 and simply 'filling in' from scratch")}
	### 3rd Pass ###	
	
	print( "--- Begin Pass 3 (Iterative Method)---")
	print("2 in explained vector = previous explained vector   1 = new additions")
	if( todd.version ){
		model.names.pass2 = rownames(m)
		m2.pass3 = m2.pass2 = m
		cls.list2.pass3 = cls.list2.pass2 = cls.list
		cls.labels2.pass3 = cls.labels2.pass2 = cls.labels
		phen.names.pass2 = phen.names
		file.suffix = paste(file.suffix, "_todd.version", sep="")
	} else{
		m2.pass3 = m2.pass2
		cls.list2.pass3 = cls.list[, s.order.pass1][, s.order.pass2]
		cls.labels2.pass3 = cls.labels[,s.order.pass1][,s.order.pass2]
	}
	model.names.pass3 = rownames(m2.pass3)
	sample.names2.pass3 = colnames(m2.pass3)
	n.phen.pass3 = 40
	
	top.genes.ind = NULL
	top.genes.names = NULL
	top.genes.vectors = NULL
	top.genes.MI = NULL
	top.diffs = NULL
	explained.vectors = NULL
	bin.gene.matrix.3 =  ifelse(cls.list2.pass3[-1,]=="WT", 0, 1)
	mid.point <- which.min(abs(m2.pass2[1,] - quantile(m2.pass2[1,], 0.5)))
	grey.and.black = c("#C0C0C0", "#000000")
	pathway.name = model.names.pass2[1]
	MI.ref = mutual.inf.2(m2.pass2[1,], m2.pass2[1,])
	
	mycol <- vector(length=512, mode = "numeric")
	for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
	for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
	mycol <- rev(mycol)
	
	explained.initial = ifelse(
			cls.list2.pass2[1,] =="WT", 0, 1)
	explained.vectors = rbind(explained.vectors, explained.initial)
	explained = explained.initial
	explained.MI.initial = mutual.inf.2(explained, m2.pass2[1,])/MI.ref
	print(paste("explained.MI.initial =", explained.MI.initial))
	print(explained)
	
	cex.axis = 1
	ncolors <- length(mycol)
	
	if(!skip.iterative){
		samples.without.mut = ifelse(
				cls.list2.pass2[2,] =="WT", 1, 0)
		#browser()
		
		wo.mut.or.blue = ifelse(c(samples.without.mut[1:mid.point], 
						rep(1, length=(Ns - mid.point)) )==1, TRUE, FALSE)
		wo.mut.and.red = ifelse(c(samples.without.mut[1:mid.point], 
						rep(0, length=(Ns - mid.point)) )==1, TRUE, FALSE)
		
		
		pdf(file=paste(results.dir, test.file.prefix, file.suffix, ".Step3_iterative.pdf", sep=""), 
				height=8.5, width=11)
		#par(mar = c(1, 15, 1, 5))
		
		
		## If we had naively searched the space without removing the explained cell lines
		MI.results = mutual.inf.3.v2(m2.pass2[1,], bin.gene.matrix.3)
		#browser()
		MI.order = order(MI.results$MI, decreasing=TRUE, na.last=TRUE)+1
		top10.names = c( #paste(c(phen.names.pass2[-1], "(from Step 2)"), collapse="  " ), 
				phen.names[MI.order[1:40]])
		top10.MI = c( #signif(explained.MI.initial, digits=4), 
				signif(MI.results$MI[MI.order[1:40]-1], digits=4))
		top10.labels = rbind(#explained+1, 
				cls.labels2.pass3[MI.order[1:40],])
		par(mar = c(1, 15, 1, 5))
		nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(1, 8), FALSE)
		max.v <- max(max(m2.pass2[1,]), -min(m2.pass2[1,]))
		V1 <- c( (ncolors/2)*normalize(m2.pass2[1,1:mid.point]) + ncolors/2, 
				(ncolors/2)*normalize(m2.pass2[1,(mid.point+1):Ns]))
		image(1:length(m2.pass2[1,]), 1:1, as.matrix(V1), 
				zlim = c(0, ncolors), col=mycol, axes=FALSE, 
				main="Naive Step 3 without exclusion of Step 2 aberrations", sub = "", xlab= "", ylab="")
		axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		V1 <- top10.labels
		V1 <- apply(V1, MARGIN=2, FUN=rev)      
#		max.v <- max(max(V1), -min(V1))
#		V1 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		image(1:length(m2.pass2[1,]), 1:dim(V1)[1], t(V1), 
				zlim = c(0, length(grey.and.black)), col=grey.and.black, axes=FALSE, 
				main="", #paste("step:", i), 
				sub = "", xlab= "", ylab="")
		axis(2, at=1:dim(V1)[1], labels=rev(top10.names), adj= 0.5, tick=FALSE, las = 1, cex.axis=1, 
				font.axis=1, line=-1)
		axis(4, at=1:dim(V1)[1], labels=rev(top10.MI), adj= 0.5, tick=FALSE, las = 1, cex.axis=1, 
				font.axis=1, line=-1)
		
		if( todd.version ){
			print(top10.labels[1:5,])
			browser()
			explained.prev = top10.labels[1,]-1
			explained.MI.prev = top10.MI[1]
		} else{
			explained.prev = explained
			explained.MI.prev = explained.MI.initial
		}
		
		for( i in 1:n.iter){
			print(paste("iteration:", i))
			MI.results = mutual.inf.3.v2( 
					m2.pass2[1,wo.mut.or.blue], 
					bin.gene.matrix.3[,wo.mut.or.blue] )
			MI.order = order(MI.results$MI, decreasing=TRUE, na.last=TRUE)+1
			top10.names = phen.names[MI.order[1:10]]
			top10.MI = MI.results$MI[MI.order[1:10]-1]
			top10.labels = cls.labels.pass3[MI.order[1:10],wo.mut.or.blue]	
			
			top.genes.ind = c(top.genes.ind, MI.order[1] )
			num.redundant = length(which(MI.results$MI == MI.results$MI[MI.order[1]-1]))-1
			top.genes.names = c( top.genes.names, paste(phen.names[MI.order[1]], "+", num.redundant, 
							ifelse(num.redundant==1, "other", "others")))
			
			mut = bin.gene.matrix.3[(MI.order[1]-1),]
			explained = ifelse(mut+explained.prev>0, 1, 0)
			explained.MI = mutual.inf.2( m2.pass2[1,wo.mut.or.blue], explained[wo.mut.or.blue])/MI.ref
			MI.diff = explained.MI - explained.MI.prev
			print(paste("Explained.MI = ", explained.MI, 
							"  MI.diff = ", ifelse(MI.diff<0, "-", "+"), 
							signif(abs(MI.diff), digits=4), sep=""))
			explained.vectors = rbind(explained.vectors, explained)
			print(2*explained.prev + mut)
			top.diffs = c(top.diffs, MI.diff)
			top.genes.vectors = rbind(top.genes.vectors, mut)
			top.genes.MI = c(top.genes.MI, paste(signif(MI.results$MI[MI.order[1]-1], digits=4), 
							sep="" ))
			
			
			par(mar = c(1, 12, 1, 12))
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(1, 8), FALSE)
			max.v <- max(max(m2.pass2[1,wo.mut.or.blue]), -min(m2.pass2[1,wo.mut.or.blue]))
			V1 <- c( (ncolors/2)*normalize(m2.pass2[1,1:mid.point]) + ncolors/2, 
					(ncolors/2)*normalize(m2.pass2[1,(mid.point+1):Ns]))[wo.mut.or.blue]
			
			image(1:length(m2.pass2[1,wo.mut.or.blue]), 1:1, as.matrix(V1), 
					zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
			axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
			V1 <- rbind( explained[wo.mut.or.blue]+1, top10.labels)
			V1 <- apply(V1, MARGIN=2, FUN=rev)      
			image(1:length(m2.pass2[1,wo.mut.or.blue]), 1:dim(V1)[1], t(V1), 
					zlim = c(0, length(grey.and.black)), col=grey.and.black, 
					axes=FALSE, main=paste("iteration:", i), sub = "", xlab= "", ylab="")
			axis(2, at=1:dim(V1)[1], labels=rev(c("explained with top result", top10.names)), 
					adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
			axis(4, at=1:dim(V1)[1], labels=rev(c( paste(signif(explained.MI, digits=4), 
											sep="" ), 
									signif(top10.MI,digits=4))), 
					adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
			
			samples.without.mut[wo.mut.or.blue] = samples.without.mut[wo.mut.or.blue] - mut[wo.mut.or.blue]
			wo.mut.or.blue = (samples.without.mut | m2.pass2[1,] <= median(m2.pass2[1,]))
			wo.mut.and.red = (samples.without.mut & m2.pass2[1,] > median(m2.pass2[1,]))
			explained.MI.prev = explained.MI
			explained.prev = ifelse(mut+explained.prev>0, 1, 0)
			print(proc.time()-t1)
			print(date())
		}
		
		explained = ifelse(apply(rbind(ifelse(cls.list2.pass2[1,]=="WT", 0,1), 
								top.genes.vectors), MARGIN=2, FUN=sum)>=1, 1, 0)
		explained.MI = mutual.inf.2(m2.pass2[1,], explained)/MI.ref
		top.genes.MI = signif(mutual.inf.3.v2(m2.pass2[1,], top.genes.vectors)$MI, digits=4)
		par(mar = c(1, 12, 1, 12))
		nf <- nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(1, 8), respect = FALSE)
		max.v <- max(max(m2.pass2[1,]), -min(m2.pass2[1,]))
		V1 <- c( (ncolors/2)*normalize(m2.pass2[1,1:mid.point]) + ncolors/2, 
				(ncolors/2)*normalize(m2.pass2[1,(mid.point+1):Ns]))
		image(1:length(m2.pass2[1,]), 1:1, as.matrix(V1), 
				zlim = c(0, ncolors), col=mycol, axes=FALSE, 
				main="Final results from iterations (removing cell lines)", sub = "", xlab= "", ylab="")
		axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		V1 <- rbind(ifelse(cls.list2.pass2[1,]=="WT", 0,1), top.genes.vectors, explained)+1
		V1 <- apply(V1, MARGIN=2, FUN=rev)      
		image(1:length(m2.pass2[1,]), 1:dim(V1)[1], t(V1), 
				zlim = c(0, length(grey.and.black)), col=grey.and.black, axes=FALSE, 
				sub = "", xlab= "", ylab="")
		axis(2, at=1:dim(V1)[1], labels=rev(c(paste(phen.names.pass2[-1], collapse=" "), 
								top.genes.names, "explained")), adj= 0.5, tick=FALSE, 
				las = 1, cex.axis=1, font.axis=1, line=-1)
		axis(4, at=1:dim(V1)[1], labels=rev( c(signif(explained.MI.initial, digits=4), 
								paste(top.genes.MI, sep=""), 
								signif(explained.MI,digits=4))), adj= 0.5, tick=FALSE, 
				las = 1, cex.axis=1, font.axis=1, line=-1)
		
		explained.vectors.MI = mutual.inf.3.v2(m2.pass2[1,], explained.vectors)$MI
		MI.diffs = explained.vectors.MI[-1] - explained.vectors.MI[1:n.iter]
		par(mar = c(1, 12, 1, 12))
		nf <- nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(1, 8), respect = FALSE)
		max.v <- max(max(m2.pass2[1,]), -min(m2.pass2[1,]))
		V1 <- c( (ncolors/2)*normalize(m2.pass2[1,1:mid.point]) + ncolors/2, 
				(ncolors/2)*normalize(m2.pass2[1,(mid.point+1):Ns]))
		image(1:length(m2.pass2[1,]), 1:1, as.matrix(V1), 
				zlim = c(0, ncolors), col=mycol, axes=FALSE, main="Final results from iterations (removing cell lines)", sub = "", xlab= "", ylab="")
		axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		V1 = apply(explained.vectors+1, MARGIN=2, FUN=rev)
		image(1:length(m2.pass2[1,]), 1:dim(V1)[1], t(V1), 
				zlim = c(0, length(grey.and.black)), col=grey.and.black, axes=FALSE, 
#			main=paste("step:", i), 
				sub = "", xlab= "", ylab="")
		left.labels = c("INITIAL cumulative", paste("cumulative, iter: ", 1:(n.iter-1), " ", sep=""),
				paste("FINAL cumulative, iter: ", n.iter, " ",sep=""))
		right.labels = c( paste(" ", signif(explained.vectors.MI[1],digits=4), sep=""), 
				paste( " ", signif(explained.vectors.MI[2:(n.iter+1)], digits=4), 
						" (",ifelse(MI.diffs < 0, "-", "+"), 
						signif(abs(MI.diffs), digits=4),")", sep="") 
		)
		axis(2, at=1:dim(V1)[1], labels=rev(left.labels), adj= 0.5, tick=FALSE, las = 1, cex.axis=1, 
				font.axis=1, line=-1)
		axis(4, at=1:dim(V1)[1], labels=rev(right.labels), adj= 0.5, tick=FALSE, 
				las = 1, cex.axis=1, font.axis=1, line=-1)
		
		cls.labels2.pass3 = rbind(cls.labels2.pass2[1,], top.genes.vectors+1)
		cls.list2.pass3 = rbind(cls.list2.pass2[1,], ifelse(top.genes.vectors==0, "WT", "MUT"))
		
		MI.results = mutual.inf.3.v2( explained,
				m2.pass2)  ## must subtract 1 from the indices because bin.gene.matrix.3 
		## doesn't include SUMMARY
		#target.vector.name=phen.pass1[ind.master],
#			n.randomizations = n.randomizations)
#	g.order.pass3.top40
		phen.descs.pass3 = ifelse( is.nan(MI.results$MI), " - ", signif(MI.results$MI, digits=3))
#	phen.descs2.pass3 = phen.descs.pass3[g.order.pass3]
		print(proc.time()-t1)
		print(date())
		
		dev.off()
	}else{ print("skipping iterative method!")}
	
	if( do.mRMR == T){
		print("--- Begin Step 3 (min redundancy Max Relevance) ---")
		
		explained.MI.initial = mutual.inf.2(explained.initial, m2.pass2[1,])/MI.ref
		print(paste("explained.MI.initial =", explained.MI.initial))
		print(explained.initial)
		
		relevance  = mutual.inf.3.v2( m2.pass2[1,], bin.gene.matrix.3, pos.and.neg=T)$MI
		redundancy = mutual.inf.3.v2( explained.initial, bin.gene.matrix.3, pos.and.neg=F)$MI
		print(proc.time()-t1)
		print(date())
		
		#browser()
		MI.D = relevance - redundancy # Mutual Information Difference
		MI.D.string = "MI.D=rel-red"
		
		MI.D.order = order(MI.D, decreasing=TRUE, na.last=TRUE)
		top.MI.D = MI.D[MI.D.order[1]]
		top.ind.MI.D = which(MI.D == top.MI.D)
		top.gene.MI.D = paste(phen.names[ top.ind.MI.D[1]+2 ], "+", length(top.ind.MI.D)-1, "others")
		## Plot and iterate with MI.D first
		print("Plot and iterate with MI.D first")
		#quartz(height=8.5, width=11)
		pdf(file=paste(results.dir, test.file.prefix, file.suffix, 
						".Step3_", MI.D.string, ".pdf", sep=""), height=8.5, width=11)
		explained.with.top.MI.D = ifelse(explained.initial + bin.gene.matrix.3[MI.D.order[1],]>=1, 1, 0)
		explained.with.top.MI.D.MI = mutual.inf.2(explained.with.top.MI.D, m2.pass2[1,])/MI.ref
		MI.diff = explained.with.top.MI.D.MI - explained.MI.initial
#	top.diffs = c(top.diffs, MI.diff)
		print(paste("Explained.MI = ", explained.with.top.MI.D.MI, 
						"  MI.diff = ", ifelse(MI.diff<0, "-", "+"), signif(abs(MI.diff), digits=4), sep=""))
		print(ifelse(cls.list2.pass2[1,]=="WT", 0,2) + bin.gene.matrix.3[MI.D.order[1],])
		top10.names = c("explained with top result ", paste(phen.names[ MI.D.order[1:10]+2 ], " "))
		top10.MI = c( paste(" MI = ", 
						signif(explained.with.top.MI.D.MI, digits=4), 
						"  diff:", ifelse(MI.diff<0, "-", "+"), signif(abs(MI.diff), digits=4), sep=""), 
				paste(" MI:", signif(relevance[MI.D.order[1:10]], digits=4), 
						"   MI.D:", signif(MI.D[MI.D.order[1:10]], digits=4), sep=""))
#	top.genes.MI = c(top.genes.MI, paste(" ", signif(top.MI.D, digits=4), 
#					#" (", ifelse(MI.diff<0, "-", "+"), signif(abs(MI.diff), digits=4), ")", 
#					sep="") )
#	top.genes.names = c(top.genes.names,
#			paste(phen.names[ top.ind.MI.D[1]+2 ], "+", length(top.ind.MI.D)-1, "others "))
#	top.genes.vectors = rbind(bin.gene.matrix.3[MI.D.order[1],])
		top10.labels = rbind( ifelse(cls.list2.pass2[1,]=="WT", 0,1), 
				explained.with.top.MI.D, bin.gene.matrix.3[MI.D.order[1:10],]) + 1
		par(mar = c(1, 12, 1, 12))
		#par(mar = c(1, 10, 1, 5))
		nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(1, 5), FALSE)
		max.v <- max(max(m2.pass2[1,]), -min(m2.pass2[1,]))
		V1 <- c( (ncolors/2)*normalize(m2.pass2[1,1:mid.point]) + ncolors/2, 
				(ncolors/2)*normalize(m2.pass2[1,(mid.point+1):Ns]))
		image(1:length(m2.pass2[1,]), 1:1, as.matrix(V1), 
				zlim = c(0, ncolors), col=mycol, axes=FALSE, main=MI.D.string, sub = "", xlab= "", ylab="")
		axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		V1 <- top10.labels
		V1 <- apply(V1, MARGIN=2, FUN=rev)      
		image(1:length(m2.pass2[1,]), 1:dim(V1)[1], t(V1), 
				zlim = c(0, length(grey.and.black)), col=grey.and.black, axes=FALSE, 
				main=paste("iteration:", 0), sub = "", xlab= "", ylab="")
		axis(2, at=1:dim(V1)[1], labels=rev(c(paste(phen.names.pass2[-1], collapse=" "), top10.names)), 
				adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		axis(4, at=1:dim(V1)[1], labels=rev(c(phen.descs.pass2[1], top10.MI)), 
				adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)	
		dev.off()
		
		#browser()
		## Plot MI.Q next
#	browser()
		print("Plot MI.Q next")
		MI.Q = (relevance)/(redundancy)   # Mutual Information Quotient
		MI.Q.string = "MI.Q=(rel)|(red)"
		
		MI.Q.order = order(MI.Q, decreasing=TRUE, na.last=TRUE)
		top.MI.Q = MI.Q[MI.Q.order[1]]
		top.ind.MI.Q = which(MI.Q == top.MI.Q)
		top.gene.MI.Q = paste(phen.names[ top.ind.MI.Q[1]+2 ], "+", length(top.ind.MI.Q)-1, "others")
		pdf(file=paste(results.dir, test.file.prefix, 
						file.suffix, ".Step3", MI.Q.string, ".pdf", sep=""), height=8.5, width=11)
		explained.with.top.MI.Q = ifelse(explained.initial + bin.gene.matrix.3[MI.Q.order[1],]>=1, 1, 0)
		explained.with.top.MI.Q.MI = mutual.inf.2(explained.with.top.MI.Q, m2.pass2[1,])/MI.ref
		MI.diff = explained.with.top.MI.Q.MI - explained.MI.initial
		#top.diffs = c(top.diffs, MI.diff)
		print(paste("Explained.MI =", explained.with.top.MI.Q.MI, 
						"  MI.diff =", ifelse(MI.diff<0, "-", "+"), signif(abs(MI.diff), digits=4), sep=""))
		print(ifelse(cls.list2.pass2[1,]=="WT", 0,2) + bin.gene.matrix.3[MI.Q.order[1],])
		top10.names = c("explained with top result", paste(phen.names[ MI.Q.order[1:10]+2 ], " "))
		top10.MI = c( paste(" MI:", 
						signif(explained.with.top.MI.Q.MI, digits=4), 
						"  diff:", ifelse(MI.diff<0, "-", "+"), signif(abs(MI.diff), digits=4), sep=""), 
				paste(" MI:", signif(relevance[MI.Q.order[1:10]], digits=4), 
						"   MI.Q:", signif(MI.Q[MI.Q.order[1:10]], digits=4), sep=""))
		top10.labels = rbind( ifelse(cls.list2.pass2[1,]=="WT", 0,1), 
				explained.with.top.MI.Q, bin.gene.matrix.3[MI.Q.order[1:10],]) +1
		#top.genes.MI = c(top.genes.MI, paste(" ", signif(top.MI.Q, digits=4), 
		#" (", ifelse(MI.diff<0, "-", "+"), signif(abs(MI.diff), digits=4), ")", 
#					sep=""))
		#top.genes.names = c(top.genes.names,
#			paste(phen.names[ top.ind.MI.Q[1]+2 ], "+", length(top.ind.MI.Q)-1, "others "))
		#top.genes.vectors = rbind(bin.gene.matrix.3[MI.Q.order[1],])
		#quartz(height=8.5, width=11)
		par(mar = c(1, 12, 1, 12))
		#par(mar = c(1, 10, 1, 5))
		nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(1, 5), FALSE)
		max.v <- max(max(m2.pass2[1,]), -min(m2.pass2[1,]))
		V1 <- c( (ncolors/2)*normalize(m2.pass2[1,1:mid.point]) + ncolors/2, (ncolors/2)*normalize(m2.pass2[1,(mid.point+1):Ns]))
		image(1:length(m2.pass2[1,]), 1:1, as.matrix(V1), 
				zlim = c(0, ncolors), col=mycol, axes=FALSE, main=MI.Q.string, sub = "", xlab= "", ylab="")
		axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		V1 <- top10.labels
		V1 <- apply(V1, MARGIN=2, FUN=rev)      
		image(1:length(m2.pass2[1,]), 1:dim(V1)[1], t(V1), 
				zlim = c(0, length(grey.and.black)), col=grey.and.black, axes=FALSE, main=paste("iteration:", 0), sub = "", xlab= "", ylab="")
		axis(2, at=1:dim(V1)[1], labels=rev(c(paste(phen.names.pass2[-1], collapse=" "), top10.names)), adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		axis(4, at=1:dim(V1)[1], labels=rev(c(phen.descs.pass2[1], top10.MI)), adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		dev.off()
	} else{ print("Skipped min redundancy Max Relevance!") }
	
	if (!is.na(output.dataset)) {
#		V.GCT <- m.all
		print("Figure out why model.descs2.pass1 does not correspond to the rows of m.all")
		browser()
#		colnames(V.GCT) <- sample.names2
#		row.names(V.GCT) <- model.names.pass1
		write.gct(gct.data.frame = m.all, descs = model.descs2.pass1, filename =output.dataset)  
	}
	
}


OPAM.sort.projection.by.score.8 <- function(
#		input.ds,
#		signatures = "NA",
		input.all.pathways.ds,
		input.cls,
		tissue = "NA",
		results.dir,
		normalize.score = T,
#		normalization.type = "zero.one",
		model = "NA",
		target.phen = NA,
		target.class = NA,
		user.colors = NA,
		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T,
		u.gene.names.known = "NA",
		n.randomizations = 10
)
# Calls MSIG.HeatMapPlot.9 and makes a plot sorted by the highest-scoring
# signatures and abnormalities (gene mutations or copy number alterations)
# i.e. doesn't require a "model" to score by as OPAM.sort.projection.by.score.2 does.
# However, it *will* use "model" if it cannot calculate p-values on the gene signatures, which
# happens when every cell line has a genomic aberration.
#
# Runs 3 passes on the data:
# 1st pass: looks at the genes and copy number alterations specified by u.gene.names.known
# 2nd pass: looks at only the top abnormalities (using a p-value cutoff) from the 1st pass, and adjusts 
# the PATHWAY.MUT vector accordingly (only according to the genes, not by copy number data)
# 3rd pass: Takes the winning signature from the 2nd pass and then looks all the genes available
#
# Very similar to OPAM.sort.projection.by.score.4, however this version uses mutual.inf instead of
# roc.area to calculate mutual information scores and p-values for PATHWAY.MUT, the vector of total genomic aberrations
# in all samples
#
# Differs from OPAM.sort.projection.by.score.6 by requiring the gct file of expression in 
# all pathways by the input tissue ("input.all.pathways.ds")
#
# Differs from OPAM.sort.projection.by.score.7 by finding the top enriched pathways that differentiate according to phenotype
# from testing all the pathways in "input.all.pathways.ds." Does not require signatures to be known a priori.
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
#	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
#	m <- data.matrix(dataset$ds)
#	model.names <- dataset$row.names
	##	model.descs <- dataset$descs
#	Ns <- length(m[1,])
#	dim(m)
#	sample.names <- dataset$names
	
	dataset.all <- MSIG.Gct2Frame( filename = input.all.pathways.ds)
	m.all <- data.matrix(dataset.all$ds)#[1:30,]
	#model.names <- dataset.all$row.names#[1:30]
	model.names <- make.unique(dataset.all$descs)
	m.all <- na.omit(t(apply(m.all, MARGIN=1, FUN=normalize)))
	Ns = length(m.all[1,])
	sample.names = dataset.all$names
	
#	if( signatures == "NA" ){
#		stop("Must provide a vector of signature names to evaluate, or specify 'ALL'")
#	}
#	if( signatures == "ALL"){
#		model.names = model.names.all
#		m = m.all
#		model.descs = dataset.all$descs
#	} else{
#		model.names = signatures
#		model.ind = match(signatures, model.names.all)
#		m = m.all[model.ind,]
#		model.descs = dataset.all$descs[model.ind]
	##	browser()
#	}
#	model.names = model.names.all
#	m = m.all
	model.descs = dataset.all$descs#[1:30]
	n.models <- length(m.all[,1])
	temp <- strsplit(input.all.pathways.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
	char.res <-  0.013 * n.models + 0.65
	
	# normalize scores
	
#	if (normalize.score == T) {
#		if (normalization.type == "zero.one") {
#			for (i in 1:n.models) {
#				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
#			}
#		} else if (normalization.type == "z.score") {
#			for (i in 1:n.models) {
#				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
#			}
#		} else if (normalization.type == "r.z.score") {
#			for (i in 1:n.models) {
#				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
#			}
#		}         
#	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list
#	browser()
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors[1])) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
#	browser()
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
	} else if( !is.null(CLS$phen.list)){
		phen.names = CLS$phen.list
	} else {
		phen.names <- "NA"
	}
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
#	print("cls.phen.colors:")
#	print(cls.phen.colors)
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		for (i in 1:n.phen) {
			n.classes[i] <- length(cls.phen[[i]])
		}
	}
	print("--- Begin Pass 1 ---")
#	model.names.original = model.names
#	m.original = m
	phen.pass1 = u.gene.names.known
	n.phen.pass1 = length(u.gene.names.known)
	ind.phen.pass1 = which( phen.names %in% phen.pass1 )
	phen.pass1 = phen.names[ind.phen.pass1]
#	phen.pass1 = c( "SUMMARY", phen.pass1)
	
#	browser()
	
	MI.list.pass1 = vector( length=n.models, mode="numeric" )
	FDR.list.pass1 = vector( length=n.models, mode="numeric" )
	
	if( !is.vector(cls.labels)){
		cls.list.pass1 = cls.list[ind.phen.pass1,]
		cls.labels.pass1 = cls.labels[ind.phen.pass1,]
	} else{ 
		cls.list.pass1 = cls.list
		cls.labels.pass1 = cls.labels
	}
	cls.list.pass1.2 = t(as.matrix(ifelse(cls.list.pass1 == "WT", 0, 1)))
	
#	browser()
	if (!is.na(target.phen)) {
		if( is.vector(cls.list.pass1.2)){ bin.class.pass1 = cls.list.pass1.2 
		} else { bin.class.pass1 = apply( cls.list.pass1.2, MARGIN=2, FUN=sum) }
		# Normalize bin.class.pass1
		if( length(unique(bin.class.pass1)) > 1){
			bin.class.pass1 = ( 
						bin.class.pass1 - min(bin.class.pass1))/(max(bin.class.pass1) - min(bin.class.pass1))
		} else if ( length(unique(bin.class.pass1)) == 1){
			bin.class = rep(1, length(cls.list[1,]))
		}
		if( is.vector( cls.list.pass1) ){
			cls.list.pass1 = ifelse(bin.class.pass1 > 0, "MUT", "WT")
		} else{	cls.list.pass1[1,] = ifelse(bin.class.pass1 > 0, "MUT", "WT") }
	} else {
#		browser()
		bin.class.pass1 <- ifelse(cls.list == cls.phen[1], 0, 1)
	}
#	browser()
#	MI.ref.models.pass1 = mutual.inf.2(bin.class.pass1, bin.class.pass1)$MI
#	print(paste("MI.ref.models.pass1 =", MI.ref.models.pass1))
#	browser()
	model.descs2.pass1 = vector(length = n.models, mode="character")
	pdf(file=paste(tissue, test.file.name, ".Phase1", "pdf", sep="."))
	#browser()
	skipped.indices = 21:(n.models)
#	browser()
	if( length(unique(bin.class.pass1)) > 1 ){
#		signature.ind = which(rownames(m.all) %in% signatures)
		MI.results = mutual.inf.3.v2(bin.class.pass1, m.all) #signature.indices = 1:n.models, )
		MI.list.pass1  = MI.results$MI
#		FDR.list.pass1 = MI.results$FDR
#		browser()
		for (i in 1:n.models) {
			MI.signif <- signif(MI.list.pass1[i], digits=3)
#			FDR.signif <- signif(FDR.list.pass1[i], digits=3)
			model.descs2.pass1[i] <- paste(MI.signif, sep="")
		}
#		browser()
		#		m.order.pass1 = order(MI.list.pass1, decreasing=TRUE, na.last=TRUE)
#		m.order.pass1 = order(MI.list.pass1, decreasing=TRUE, na.last=TRUE)
		m.order.pass1 = order(MI.list.pass1, decreasing=FALSE, na.last=TRUE)
		m.order.pass1.top10 = m.order.pass1[-skipped.indices]
#		m.order.pass1 = 1:n.models
		m2.pass1 <- m.all[m.order.pass1, ]
		s.order.pass1 <- order(m2.pass1[1,], decreasing = FALSE )
#		s.order.pass1 = 1:Ns
		m2.pass1 <- m2.pass1[-skipped.indices, s.order.pass1]
#		m2.pass1.top10 = m2.pass1[-skipped.indices,]
	} else{ 
		MI.list.pass1 = rep(NA, n.models)
		#	FDR.list.pass1 = rep(NA, n.models)
		model.descs2.pass1 = rep(" - ", n.models)
		loc <- match(model, model.names)
		s.order.pass1 <- order(m.all[loc,], decreasing = decreasing.order)
#		loc = s.order.pass1[1]
#		s.order.pass1 = 1:Ns
		m2.pass1 <- m.all[, s.order.pass1]
		correl <- cor(t(m2.pass1))[, loc]
		m.order.pass1 <- order(correl, decreasing=T)
#		m.order.pass1 = 1:n.models
		m2.pass1 <- m2.pass1[m.order.pass1, ]
	}
#	browser()
	
	dev.off()
#	skipped.indices = 11:(n.models-10)
#	browser()
	
	MI.list.pass1.top10 = MI.list.pass1[m.order.pass1.top10]
	MI.list.pass1 = MI.list.pass1[m.order.pass1]
#	FDR.list.pass1.top10 = FDR.list.pass1[m.order.pass1.top10]
#	FDR.list.pass1 = FDR.list.pass1[m.order.pass1]
	
	bin.class.pass1 = bin.class.pass1[s.order.pass1]
#	m2.pass1 <- m2.pass1[m.order.pass1, ]
	model.descs2.pass1.top10 = model.descs2.pass1[m.order.pass1.top10]	
	model.descs2.pass1 = model.descs2.pass1[m.order.pass1]
	
	sample.names2.pass1 <- colnames(m2.pass1)
	model.names.pass1.top10 <- rownames(m2.pass1)
	print(matrix(c(model.names.pass1.top10, model.descs2.pass1.top10), ncol=2), quote=F)
#	browser()
	if (is.vector(cls.labels)) {
		cls.labels2.pass1 <- cls.labels.pass1[s.order.pass1]
		cls.list2.pass1 <- cls.list.pass1[s.order.pass1]
	} else {
		cls.labels2.pass1 <- cls.labels.pass1[, s.order.pass1]
		cls.list2.pass1 <- cls.list.pass1[, s.order.pass1]
	}
	m.all = m.all[, s.order.pass1]
#	browser()
	winning.model.ind.pass1 = which(model.names.pass1.top10[1] == rownames(m.all))
	
#	pathway.name <- "KRAS_ALL_UP"
#	pathway <- m[1,]
#	pathway0 <- ifelse(pathway < median(pathway), 0, 1) # disctretized version
	
#	MI.ref.genes.pass1 <- mutual.inf.2(m[1,], m[1,])$MI
	
#	browser()
#	m.score.pass1 <- m2.pass1[1,]
#	m.score.norm.pass1 <- (m.score.pass1 - min(m.score.pass1))/(max(m.score.pass1) - min(m.score.pass1))
#	m.score.pass1 = ifelse( m.score.pass1 < median(m.score.pass1), -1, 1)   # discretized version
#	MI.ref.genes.pass1 <- mutual.inf.2(m.score.norm.pass1, m.score.norm.pass1)$MI
#	print(paste("MI.ref.genes.pass1 =", MI.ref.genes.pass1))
	MI.list.phen.pass1 = vector(mode="numeric", length=n.phen.pass1)
#	FDR.list.phen.pass1 = vector(mode="numeric", length=n.phen.pass1)
	phen.descs2.pass1 = vector(mode="character", length=n.phen.pass1)
	
	if( length(unique(bin.class.pass1)) > 1){
#		MI.results <-(mutual.inf.3(bin.class.pass1, m.all, 
#							winning.model.ind.pass1, gene.target.name = phen.pass1[1]))#/MI.ref.genes.pass1
		MI.signif <- signif(MI.list.pass1[1], digits=3)
		MI.list.phen.pass1[1] = MI.list.pass1[1]
#		FDR.signif <- signif(FDR.list.pass1[1], digits=3)
#		FDR.list.phen.pass1[1] = FDR.list.pass1[1]
	} else{
		MI.signif <- "-"
		MI.list.phen.pass1[1] = NA
#		FDR.signif <- "- "
#		FDR.list.phen.pass1[1] = NA
	}
	print(paste(format(phen.pass1[1], width=12), "mutual.inf =", MI.signif))
	phen.descs2.pass1[1] <- paste(MI.signif)
#	browser()
	if( n.phen >= 2 ){
		bin.gene.matrix = ifelse(cls.list2.pass1[-1,]=="WT", 0, 1)
		n.aberrations = apply(bin.gene.matrix, MARGIN=1, FUN=sum)
		u.n.aberrations = unique(n.aberrations[n.aberrations != 0])
		for( i in 1:length(u.n.aberrations)){
			ind.without.SUMMARY = which(n.aberrations == u.n.aberrations[i])
			ind.master = ind.without.SUMMARY + 1
#		browser()
#		bin.gene.matrix.temp = bin.gene.matrix[ind.without.SUMMARY,]
			
			MI.results = mutual.inf.3.v2(bin.gene.matrix[ind.without.SUMMARY,], 
					m.all, winning.model.ind.pass1, gene.target.name=phen.pass1[ind.master],
					n.randomizations = n.randomizations)
			
			MI.list.phen.pass1[ind.master] = MI.results$MI
			#		FDR.list.phen.pass1[ind.master] = MI.results$FDR
			for( j in 1:length(ind.master)){
				phen.descs.pass1[ind.master[j]] = 
						paste( signif(MI.results$MI[j], digits=3), 
								sep="")
			}
		}
		ind.zeros = which(n.aberrations==0) + 1
		MI.list.phen.pass1[ind.zeros] = NA
#		FDR.list.phen.pass1[ind.zeros] = NA
		phen.descs.pass1[ind.zeros] = " - "
	}
	
	phen.names.pass1 = phen.pass1#[g.order.pass1]#[1:n.phen.pass1]
#	browser(text="Figure out how to print phen.descs.pass1 and phen.names.pass1 in a nice table")
	
	# Recompute cls.list2 as some mutations or copy numbers may have been removed
	
	
	# Recompute cls.phen and cls.labels2 as order may have changed
	
	cls.phen2.pass1 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(cls.list2.pass1)
		cls.phen2.pass1 <- classes
		cls.labels2.pass1 <- match(cls.list2.pass1, cls.phen2.pass1)
	} else {
		for (kk in 1:length(cls.list2.pass1[, 1])) {
			classes <- unique(cls.list2.pass1[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2.pass1 <- c(cls.phen2.pass1, classes)
			cls.labels2.pass1[kk,] <- match(cls.list2.pass1[kk,], classes)
		}
	}
#	cls.labels2.pass1 = cls.labels2.pass1[1:n.phen.pass1,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list.pass1 <- unlist(cls.phen2.pass1)
	colors.list.pass1 = rep( "gray", length(phen.list.pass1))
	colors.list.pass1[phen.list.pass1=="MUT"] = cls.phen.colors[1]
	colors.list.pass1[phen.list.pass1=="DEL"] = cls.phen.colors[3]
	colors.list.pass1[phen.list.pass1=="AMP"] = cls.phen.colors[4]
	colors.list.pass1[phen.list.pass1=="ALT"] = cls.phen.colors[5]
#	browser()
#	colors.list.pass1[1,] = grey(bin.class.pass1)
#	print("cls.phen2:")
#	print(unlist(cls.phen2))
#	
#	print("cls.phen:")
#	print(unlist(cls.phen))
#	
#	print("colors.list:")
#	print(colors.list)
	
#	browser()
	
	filename <- paste(results.dir, test.file.prefix, ".Phase1.MI|HXY", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
	pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 17 )
#	browser()
#   windows(width=12, height=8)
	MSIG.HeatMapPlot.9(V = m2.pass1, 
#			pathway.mut = bin.class.pass1,
			row.names = model.names.pass1.top10,
			row.names2 = model.descs2.pass1.top10, 
			col.labels = cls.labels2.pass1, 
			col.classes = cls.phen2.pass1, 
			phen.cmap = colors.list.pass1, 
			phen.names = phen.names.pass1,
			phen.names2 = phen.descs2.pass1,
			col.names = sample.names2.pass1, 
			main = paste(tissue, test.file.prefix), 
			xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	dev.off()
#	browser()
	if (!is.na(output.dataset)) {
#		V.GCT <- m.all
#		colnames(V.GCT) <- sample.names2
#		row.names(V.GCT) <- model.names2
		write.gct(gct.data.frame = m.all, descs = model.descs2.pass1, filename = paste(output.dataset, ".gct", sep=""))
		write.cls.2( class.v = cls.labels2.pass1, phen = cls.phen, filename = paste(output.dataset, ".cls", sep=""))
	}
	
}




MSIG.HeatMapPlot.11 <- function(
		## For Plotting expression heatmap only!! (No phenotypes)
		V, 
		row.names = "NA",
		row.names2 = "NA", 
		col.labels = "NA",
		col.labels2 = "NA", 
		col.classes = "NA", 
		phen.cmap = NULL, 
		col.names = "NA",
		phen.names = "NA", 
		phen.names2 = "NA",
		main = " ", 
		sub = " ", 
		xlab=" ", 
		ylab=" ",
		row.norm = TRUE,
		char.rescale = 0.85,                               
		cmap.type = 1,   # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
		max.v = "NA",
		legend = T)
{
#
# Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	n.rows <- length(V[,1])
	n.cols <- length(V[1,])
	V1 <- matrix(0, nrow=n.rows, ncol=n.cols)
	
#       if ((cmap.type == 5) | (cmap.type == 3)) {
	if (cmap.type == 5) {
		row.norm <- F
	}
	
	if (row.norm == TRUE) {
		row.mean <- apply(V, MARGIN=1, FUN=mean)
		row.sd <- apply(V, MARGIN=1, FUN=sd)
		row.n <- length(V[,1])
		for (i in 1:n.rows) {
			if (row.sd[i] == 0) {
				V1[i,] <- 0
			} else {
				V1[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
			}
			V1[i,] <- ifelse(V1[i,] < -4, -4, V1[i,])
			V1[i,] <- ifelse(V1[i,] > 4, 4, V1[i,])
		}
	} else {
		V1 <- V
	}
	
	if (cmap.type == 1) { 
		mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
				"#FF9DB0", "#FF7080", 
				"#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage,
		# pre-gene cluster, original pinkogram color map
	} else if (cmap.type == 2) {
		violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
		mycol <- rev(violet.palette(20))
		
#          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
#                        "#BCBDDC","#A8A6CF",
#                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
#                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
	} else if (cmap.type == 6) {
		mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
				"#7A0177", "#49006A")
	} else if (cmap.type == 7) {
		mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
				"#A63603", "#7F2704")
	} else if (cmap.type == 8) {
		mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
				"#006D2C", "#00441B")
	} else if (cmap.type == 9) {
		mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
				"#08519C", "#08306B")
	} else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
		mycol <- vector(length=512, mode = "numeric")
		
		for (k in 1:256) {
			mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
		}
		for (k in 257:512) {
			mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
		}
		mycol <- rev(mycol)
	}
	
	ncolors <- length(mycol)
	
	if (cmap.type == 5) {
		if (max.v == "NA") {
			max.v <- max(max(V1), -min(V1))
		}
		V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		
	} else {
		V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
	}
	
	if (col.labels[1] == "NA") {      
		heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
		heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
		tot.cols <- ncolors
		browser()
		if (legend == T) {
			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(8, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap[1:length(col.classes)])
		image(1:n.cols, 1:n.rows, t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
		n.rows.phen <- 0
	} else {
		tot.cols <- ncolors
		if (is.vector(col.labels)) {
			heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			n.rows.phen <- 1
			heatm[n.rows + 1,] <- tot.cols + col.labels
			cols.row <- length(unique(col.labels))
			tot.cols <- tot.cols + cols.row
			phen.cmap <- phen.cmap[1:cols.row]
		} else {
			n.rows.phen <- length(col.labels[,1])
			cols.row <- vector(length=n.rows.phen, mode = "numeric")
			heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
				heatm[k,] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1,]
				cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1,]))
				tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
#                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))
				
			}
			phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
		}
		if (legend == T) {
#              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
			nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(8, 4, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(5, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap)
		image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
	}
	
# Add lines to separate phenotypes or subgroups
	
	if (col.labels2[1] != "NA") {
		groups <-  split(col.labels2, col.labels2)
		len.vec <- lapply(groups, length)
		plot.div <- c(0.51, cumsum(len.vec) + 0.5)
		for (i in plot.div) {
			lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
		}
		lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48), lwd = 2,
				cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50), lwd = 2,
				cex = 0.9, col = "black")
	}
	if (row.names[1] != "NA") {
#		browser()
		numC <- nchar(row.names)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names[i] <- substr(row.names[i], 1, 40)
			row.names[i] <- paste(row.names[i], " ", sep="")
		}
		if (phen.names[1] == "NA") {
			head.names <- paste("Class", seq(n.rows.phen, 1, -1))
		} else {
			head.names <- as.character(rev(phen.names))
		}
		row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
#            print(paste("n.rows:", n.rows))
#            print(paste("Phen names:", phen.names))
#            print(paste("Head names:", head.names))
#            print(paste("Row names:", row.names))
		axis(2, at=1:(n.rows + n.rows.phen), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char,
				font.axis=2, line=-1)
	}
	
	if (row.names2[1] != "NA") {
#		browser()
		numC <- nchar(row.names2)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names2[i] <- substr(row.names2[i], 1, 40)
			row.names2[i] <- paste(" ", row.names2[i], sep="")
			
		}
		for( i in 1:n.rows.phen ){
			phen.names2[i] <- substr(phen.names2[i], 1, 40)
			phen.names2[i] <- paste( " ", phen.names2[i], sep="")
		}
		
		row.names2 <- rev(row.names2)
		phen.names2 <- rev(phen.names2)
		axis(4, at=1:(n.rows + n.rows.phen), labels=c(row.names2, phen.names2), adj= 0.5, tick=FALSE, las = 1, 
				cex.axis=size.row.char, font.axis=2, line=-1)
	}
	
	if (col.names[1] != "NA") {
		size.col.char <- char.rescale*20/(n.cols + 25)
		axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
	}
	
	# Phenotype Legend 
	
#      print("--------------------------------------------------------------------------------------------")
	if (legend == T) {
		leg.txt <- NULL
		p.vec <- NULL
		c.vec <- NULL
		c2.vec <- NULL
		ind <- 1
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		for (i in 1:n.rows.phen) {  
			if (is.vector(col.labels)) {
				phen.v <- as.character(col.classes)
			} else {
				phen.v <- as.character(col.classes[[i]])
			}
			p.name <- paste(as.character(rev(head.names)[i]), ":   ", sep="")
			leg.txt <- c(p.name, phen.v)  
			p.vec <-  rep(22, cols.row[i] + 1)
			c.vec <-  c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
			c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
			ind <- ind + cols.row[i]
			offset <- 0.07
			legend(x=0, y= 1 - offset*i, 
					horiz = T, x.intersp = 0.5, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, 
					pt.bg = c.vec, col = c2.vec, cex = 1.20, pt.cex=1.75)
		}
	}
	
	# Color map legend
	
#       print(c("range V=", range(V)))
#       print(c("range V1=", range(V1)))
#       print(c("range V2=", range(V2)))
	
	par(mar = c(2, 12, 2, 12))
	num.v <- 20
	range.v <- range(V2)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
	image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col=mycol, axes=FALSE,
			main=" ", sub = " ", xlab= ylab, ylab=xlab)
	range.v <- range(V1)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
#          print(c("heatm.v2=", heatm.v2))
	axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.5*char.rescale, font.axis=1)
	
	return()
	
}

MSIG.HeatMapPlot.9 <- function(
		V, 
		row.names = "NA",
		row.names2 = "NA", 
		col.labels = "NA",
		col.labels2 = "NA", 
		col.classes = "NA", 
		phen.cmap = "NA", 
		col.names = "NA",
		phen.names = "NA", 
		phen.names2 = "NA",
		main = " ", 
		sub = " ", 
		xlab=" ", 
		ylab=" ",
		row.norm = TRUE,
		char.rescale = 0.85,                               
		cmap.type = 1,   # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
		max.v = "NA",
		legend = F)
{
#
# Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
# Doesn't plot the spectrum on the bottom
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	n.rows <- length(V[,1])
	n.cols <- length(V[1,])
	V1 <- matrix(0, nrow=n.rows, ncol=n.cols)
	
#       if ((cmap.type == 5) | (cmap.type == 3)) {
	if (cmap.type == 5) {
		row.norm <- F
	}
	
	if (row.norm == TRUE) {
		row.mean <- apply(V, MARGIN=1, FUN=mean)
		row.sd <- apply(V, MARGIN=1, FUN=sd)
		row.n <- length(V[,1])
		for (i in 1:n.rows) {
			if (row.sd[i] == 0) {
				V1[i,] <- 0
			} else {
				V1[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
			}
			V1[i,] <- ifelse(V1[i,] < -4, -4, V1[i,])
			V1[i,] <- ifelse(V1[i,] > 4, 4, V1[i,])
		}
	} else {
		V1 <- V
	}
	
	if (cmap.type == 1) { 
		mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
				"#FF9DB0", "#FF7080", 
				"#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage,
		# pre-gene cluster, original pinkogram color map
	} else if (cmap.type == 2) {
		violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
		mycol <- rev(violet.palette(20))
		
#          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
#                        "#BCBDDC","#A8A6CF",
#                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
#                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
	} else if (cmap.type == 6) {
		mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
				"#7A0177", "#49006A")
	} else if (cmap.type == 7) {
		mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
				"#A63603", "#7F2704")
	} else if (cmap.type == 8) {
		mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
				"#006D2C", "#00441B")
	} else if (cmap.type == 9) {
		mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
				"#08519C", "#08306B")
	} else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
		mycol <- vector(length=512, mode = "numeric")
		
		for (k in 1:256) {
			mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
		}
		for (k in 257:512) {
			mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
		}
		mycol <- rev(mycol)
	}
	
	ncolors <- length(mycol)
	
	if (cmap.type == 5) {
		if (max.v == "NA") {
			max.v <- max(max(V1), -min(V1))
		}
		V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		
	} else {
		V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
	}
	
	if (col.labels[1] == "NA") {      
		heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
		heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
		tot.cols <- ncolors
		if (legend == T) {
			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(8, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap[1:length(col.classes)])
		image(1:n.cols, 1:n.rows, t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
		n.rows.phen <- 0
	} else {
		tot.cols <- ncolors
		if (is.vector(col.labels)) {
			heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			n.rows.phen <- 1
			heatm[n.rows + 1,] <- tot.cols + col.labels
			cols.row <- length(unique(col.labels))
			tot.cols <- tot.cols + cols.row
			phen.cmap <- phen.cmap[1:cols.row]
		} else {
			n.rows.phen <- length(col.labels[,1])
			cols.row <- vector(length=n.rows.phen, mode = "numeric")
			heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
				heatm[k,] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1,]
				cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1,]))
				tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
#                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))
				
			}
			phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
		}
		if (legend == T) {
#              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
			nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(8, 4, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(5, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap)
		image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
	}
	
# Add lines to separate phenotypes or subgroups
	
	if (col.labels2[1] != "NA") {
		groups <-  split(col.labels2, col.labels2)
		len.vec <- lapply(groups, length)
		plot.div <- c(0.51, cumsum(len.vec) + 0.5)
		for (i in plot.div) {
			lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
		}
		lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48), lwd = 2,
				cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50), lwd = 2,
				cex = 0.9, col = "black")
	}
	if (row.names[1] != "NA") {
#		browser()
		numC <- nchar(row.names)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names[i] <- substr(row.names[i], 1, 40)
			row.names[i] <- paste(row.names[i], " ", sep="")
		}
		if (phen.names[1] == "NA") {
			head.names <- paste("Class", seq(n.rows.phen, 1, -1))
		} else {
			head.names <- as.character(rev(phen.names))
		}
		row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
#            print(paste("n.rows:", n.rows))
#            print(paste("Phen names:", phen.names))
#            print(paste("Head names:", head.names))
#            print(paste("Row names:", row.names))
		axis(2, at=1:(n.rows + n.rows.phen), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char,
				font.axis=2, line=-1)
	}
	
	if (row.names2[1] != "NA") {
#		browser()
		numC <- nchar(row.names2)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names2[i] <- substr(row.names2[i], 1, 40)
			row.names2[i] <- paste(" ", row.names2[i], sep="")
			
		}
		for( i in 1:n.rows.phen ){
			phen.names2[i] <- substr(phen.names2[i], 1, 40)
			phen.names2[i] <- paste( " ", phen.names2[i], sep="")
		}
		
		row.names2 <- rev(row.names2)
		phen.names2 <- rev(phen.names2)
		axis(4, at=1:(n.rows + n.rows.phen), labels=c(row.names2, phen.names2), adj= 0.5, tick=FALSE, las = 1, 
				cex.axis=size.row.char, font.axis=2, line=-1)
	}
	
	if (col.names[1] != "NA") {
		size.col.char <- char.rescale*20/(n.cols + 25)
		axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
	}
	
	# Phenotype Legend 
	
#      print("--------------------------------------------------------------------------------------------")
	if (legend == T) {
		leg.txt <- NULL
		p.vec <- NULL
		c.vec <- NULL
		c2.vec <- NULL
		ind <- 1
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		for (i in 1:n.rows.phen) {  
			if (is.vector(col.labels)) {
				phen.v <- as.character(col.classes)
			} else {
				phen.v <- as.character(col.classes[[i]])
			}
			p.name <- paste(as.character(rev(head.names)[i]), ":   ", sep="")
			leg.txt <- c(p.name, phen.v)  
			p.vec <-  rep(22, cols.row[i] + 1)
			c.vec <-  c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
			c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
			ind <- ind + cols.row[i]
			offset <- 0.07
			legend(x=0, y= 1 - offset*i, 
					horiz = T, x.intersp = 0.5, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, 
					pt.bg = c.vec, col = c2.vec, cex = 1.20, pt.cex=1.75)
		}
	}
	
	# Color map legend
	
#       print(c("range V=", range(V)))
#       print(c("range V1=", range(V1)))
#       print(c("range V2=", range(V2)))
	
#	par(mar = c(2, 12, 2, 12))
#	num.v <- 20
#	range.v <- range(V2)
#	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
#	heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
#	image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col=mycol, axes=FALSE,
#			main=" ", sub = " ", xlab= ylab, ylab=xlab)
#	range.v <- range(V1)
#	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
#	heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
	##          print(c("heatm.v2=", heatm.v2))
#	axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.5*char.rescale, font.axis=1)
	
	return()
	
}

MSIG.HeatMapPlot.10<- function(
		V, 
		pathway.mut,
		row.names = "NA",
		row.names2 = "NA", 
		col.labels = "NA",
		col.labels2 = "NA", 
		col.classes = "NA", 
		phen.cmap = "NA", 
		col.names = "NA",
		phen.names = "NA", 
		phen.names2 = "NA",
		main = " ", 
		sub = " ", 
		xlab=" ", 
		ylab=" ",
		row.norm = TRUE,
		char.rescale = 0.85,                               
		cmap.type = 1,   # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
		max.v = "NA",
		legend = T,
		tissue.names = "NA",
		tissue.labels = NA)
{
#
# Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
# Doesn't plot the spectrum on the bottom
#
# Plots PATHWAY.MUT as a continuous vector in a greyscale spectrum
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	n.rows <- length(V[,1])
	n.cols <- length(V[1,])
	V1 <- matrix(0, nrow=n.rows, ncol=n.cols)
	set3 = brewer.pal(12, "Set3")
	accent = brewer.pal(8, "Accent")
	
#       if ((cmap.type == 5) | (cmap.type == 3)) {
	if (cmap.type == 5) {
		row.norm <- F
	}
	
	if (row.norm == TRUE) {
		row.mean <- apply(V, MARGIN=1, FUN=mean)
		row.sd <- apply(V, MARGIN=1, FUN=sd)
		row.n <- length(V[,1])
		for (i in 1:n.rows) {
			if (row.sd[i] == 0) {
				V1[i,] <- 0
			} else {
				V1[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
			}
			V1[i,] <- ifelse(V1[i,] < -4, -4, V1[i,])
			V1[i,] <- ifelse(V1[i,] > 4, 4, V1[i,])
		}
	} else {
		V1 <- V
	}
	
	if (cmap.type == 1) { 
		mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
				"#FF9DB0", "#FF7080", 
				"#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage,
		# pre-gene cluster, original pinkogram color map
	} else if (cmap.type == 2) {
		violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
		mycol <- rev(violet.palette(20))
		
#          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
#                        "#BCBDDC","#A8A6CF",
#                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
#                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
	} else if (cmap.type == 6) {
		mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
				"#7A0177", "#49006A")
	} else if (cmap.type == 7) {
		mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
				"#A63603", "#7F2704")
	} else if (cmap.type == 8) {
		mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
				"#006D2C", "#00441B")
	} else if (cmap.type == 9) {
		mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
				"#08519C", "#08306B")
	} else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
		mycol <- vector(length=512, mode = "numeric")
		
		for (k in 1:256) {
			mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
		}
		for (k in 257:512) {
			mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
		}
		mycol <- rev(mycol)
	}
#	browser()
	ncolors <- length(mycol)
	pathway.mut = (-(pathway.mut*.749 + 0.251 - 1))
#	image(1:n.cols, 1, as.matrix(pathway.mut), col=gray(n.cols:0/n.cols))
	if (cmap.type == 5) {
		if (max.v == "NA") {
			max.v <- max(max(V1), -min(V1))
		}
		V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		
	} else {
		V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
	}
	
	if (col.labels[1] == "NA") {      
		heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
		heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
		tot.cols <- ncolors
		if (legend == T) {
			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(8, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap[1:length(col.classes)])
		image(1:n.cols, 1:n.rows, t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
		n.rows.phen <- 0
	} else {
		tot.cols <- ncolors
		if (is.vector(col.labels)) {
			heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			n.rows.phen <- 1
			heatm[n.rows + 1,] <- tot.cols + col.labels
			cols.row <- length(unique(col.labels))
			tot.cols <- tot.cols + cols.row
			phen.cmap <- phen.cmap[1:cols.row]
			pathway.mut.grey = grey(pathway.mut)
			u.pathway.mut.grey = unique(pathway.mut.grey)
			heatm[n.rows + n.rows.phen,] = match(pathway.mut.grey, u.pathway.mut.grey) + tot.cols
		} else {
			n.rows.phen <- length(col.labels[,1])
			cols.row <- vector(length=n.rows.phen, mode = "numeric")
			heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
#			heatm[n.rows+n.rows.phen,] = t(as.matrix(gray(pathway.mut)))
			for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
				heatm[k,] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1,]
				cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1,]))
				tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
#                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))
				
			}
#			browser()
			pathway.mut.grey = grey(pathway.mut)
			u.pathway.mut.grey = unique(pathway.mut.grey)
			heatm[n.rows + n.rows.phen,] = match(pathway.mut.grey, u.pathway.mut.grey) + tot.cols
			tot.cols = tot.cols + length(u.pathway.mut.grey)
			phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
		}
#		image(as.matrix(grey(pathway.mut)))
		if (legend == T) {
#              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
			nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(8, 4, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(5, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
#		browser()
		mycol <- c(mycol, phen.cmap, u.pathway.mut.grey)
#		browser()
		image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
	}
	
# Add lines to separate phenotypes or subgroups
	
	if (col.labels2[1] != "NA") {
		groups <-  split(col.labels2, col.labels2)
		len.vec <- lapply(groups, length)
		plot.div <- c(0.51, cumsum(len.vec) + 0.5)
		for (i in plot.div) {
			lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
		}
		lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48), lwd = 2,
				cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50), lwd = 2,
				cex = 0.9, col = "black")
	}
	if (row.names[1] != "NA") {
#		browser()
		numC <- nchar(row.names)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names[i] <- substr(row.names[i], 1, 40)
			row.names[i] <- paste(row.names[i], " ", sep="")
		}
		if (phen.names[1] == "NA") {
			head.names <- paste("Class", seq(n.rows.phen, 1, -1))
		} else {
			head.names <- as.character(rev(phen.names))
		}
		row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
#            print(paste("n.rows:", n.rows))
#            print(paste("Phen names:", phen.names))
#            print(paste("Head names:", head.names))
#            print(paste("Row names:", row.names))
#		browser()
		axis(2, at=1:(n.rows + n.rows.phen), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char,
				font.axis=2, line=-1)
	}
	
	if (row.names2[1] != "NA") {
#		browser()
		numC <- nchar(row.names2)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names2[i] <- substr(row.names2[i], 1, 40)
			row.names2[i] <- paste(" ", row.names2[i], sep="")
			
		}
		for( i in 1:n.rows.phen ){
			phen.names2[i] <- substr(phen.names2[i], 1, 40)
			phen.names2[i] <- paste( " ", phen.names2[i], sep="")
		}
		
		row.names2 <- rev(row.names2)
		phen.names2 <- rev(phen.names2)
		axis(4, at=1:(n.rows + n.rows.phen), labels=c(row.names2, phen.names2), adj= 0.5, tick=FALSE, las = 1, 
				cex.axis=size.row.char, font.axis=2, line=-1)
	}
	
	if (col.names[1] != "NA") {
		size.col.char <- char.rescale*20/(n.cols + 25)
		axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
	}
	
	# Phenotype Legend 
	
#      print("--------------------------------------------------------------------------------------------")
	if (legend == T) {
		leg.txt <- NULL
		p.vec <- NULL
		c.vec <- NULL
		c2.vec <- NULL
		ind <- 1
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		for (i in 1:n.rows.phen) {  
#			browser()
			if (is.vector(col.labels)) {
				phen.v <- as.character(col.classes)
			} else {
				phen.v <- as.character(col.classes[[i]])
			}
			p.name <- paste(as.character(rev(head.names)[i]), ":   ", sep="")
			leg.txt <- c(p.name, phen.v)  
			p.vec <-  rep(22, cols.row[i] + 1)
			c.vec <-  c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
			c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
			ind <- ind + cols.row[i]
			offset <- 0.07
			legend(x=0, y= 1 - offset*i, 
					horiz = T, x.intersp = 0.5, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, 
					pt.bg = c.vec, col = c2.vec, cex = 1.20, pt.cex=1.75)
		}
	}
	
	# Color map legend
	
#       print(c("range V=", range(V)))
#       print(c("range V1=", range(V1)))
#       print(c("range V2=", range(V2)))
	
	par(mar = c(2, 12, 2, 12))
	num.v <- 20
	range.v <- range(V2)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
	image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col=mycol, axes=FALSE,
			main=" ", sub = " ", xlab= ylab, ylab=xlab)
	range.v <- range(V1)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
	#          print(c("heatm.v2=", heatm.v2))
	axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.5*char.rescale, font.axis=1)
	
	return()
	
}

MSIG.HeatMapPlot.10.multiple.tissues <- function(
		V, 
		pathway.mut,
		row.names = "NA",
		row.names2 = "NA", 
		col.labels = "NA",
		col.labels2 = "NA", 
		col.classes = "NA", 
		phen.cmap = "NA", 
		col.names = "NA",
		phen.names = "NA", 
		phen.names2 = "NA",
		main = " ", 
		sub = " ", 
		xlab=" ", 
		ylab=" ",
		row.norm = TRUE,
		char.rescale = 0.85,                               
		cmap.type = 1,   # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
		max.v = "NA",
		legend = T,
		tissue.names = "NA",
		tissue.labels = NA)
{
#
# Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
# Doesn't plot the spectrum on the bottom
#
# Plots PATHWAY.MUT as a continuous vector in a greyscale spectrum
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	library(RColorBrewer)
	
	n.tissues = length(tissue.names)
	
	n.rows <- length(V[,1])
	n.cols <- length(V[1,])
	V1 <- matrix(0, nrow=n.rows, ncol=n.cols)
	
	
#       if ((cmap.type == 5) | (cmap.type == 3)) {
	if (cmap.type == 5) {
		row.norm <- F
	}
	
	if (row.norm == TRUE) {
		row.mean <- apply(V, MARGIN=1, FUN=mean)
		row.sd <- apply(V, MARGIN=1, FUN=sd)
		row.n <- length(V[,1])
		for (i in 1:n.rows) {
			if (row.sd[i] == 0) {
				V1[i,] <- 0
			} else {
				V1[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
			}
			V1[i,] <- ifelse(V1[i,] < -4, -4, V1[i,])
			V1[i,] <- ifelse(V1[i,] > 4, 4, V1[i,])
		}
	} else {
		V1 <- V
	}
	
	if (cmap.type == 1) { 
		mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
				"#FF9DB0", "#FF7080", 
				"#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage,
		# pre-gene cluster, original pinkogram color map
	} else if (cmap.type == 2) {
		violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
		mycol <- rev(violet.palette(20))
		
#          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
#                        "#BCBDDC","#A8A6CF",
#                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
#                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
	} else if (cmap.type == 6) {
		mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
				"#7A0177", "#49006A")
	} else if (cmap.type == 7) {
		mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
				"#A63603", "#7F2704")
	} else if (cmap.type == 8) {
		mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
				"#006D2C", "#00441B")
	} else if (cmap.type == 9) {
		mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
				"#08519C", "#08306B")
	} else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
		mycol <- vector(length=512, mode = "numeric")
		
		for (k in 1:256) {
			mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
		}
		for (k in 257:512) {
			mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
		}
		mycol <- rev(mycol)
	}
#	browser()
	ncolors <- length(mycol)
	pathway.mut = (-(pathway.mut*.749 + 0.251 - 1))
#	image(1:n.cols, 1, as.matrix(pathway.mut), col=gray(n.cols:0/n.cols))
	if (cmap.type == 5) {
		if (max.v == "NA") {
			max.v <- max(max(V1), -min(V1))
		}
		V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		
	} else {
		V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
	}
	
	if (col.labels[1] == "NA") {      
		heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
		heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
		tot.cols <- ncolors
		if (legend == T) {
			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(8, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap[1:length(col.classes)])
		image(1:n.cols, 1:n.rows, t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
		n.rows.phen <- 0
	} else {
		tot.cols <- ncolors
		if (is.vector(col.labels)) {
			heatm <- matrix(0, nrow = n.rows + 2, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			n.rows.phen <- 1
			heatm[n.rows + 1,] <- tot.cols + col.labels
			cols.row <- length(unique(col.labels))
			tot.cols <- tot.cols + cols.row
			phen.cmap <- phen.cmap[1:cols.row]
			pathway.mut.grey = grey(pathway.mut)
			u.pathway.mut.grey = unique(pathway.mut.grey)
			heatm[n.rows + n.rows.phen,] = match(pathway.mut.grey, u.pathway.mut.grey) + tot.cols
		} else {
			n.rows.phen <- length(col.labels[,1])
			cols.row <- vector(length=n.rows.phen, mode = "numeric")
			heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
#			heatm[n.rows+n.rows.phen,] = t(as.matrix(gray(pathway.mut)))
			for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
				heatm[k,] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1,]
				cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1,]))
				tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
#                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))
				
			}
#			browser()
			pathway.mut.grey = grey(pathway.mut)
			u.pathway.mut.grey = unique(pathway.mut.grey)
			heatm[n.rows + n.rows.phen,] = match(pathway.mut.grey, u.pathway.mut.grey) + tot.cols
			tot.cols = tot.cols + length(u.pathway.mut.grey)
			phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
		}
#		image(as.matrix(grey(pathway.mut)))
		if (legend == T) {
#              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
			nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(8, 4, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(5, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		#browser()
		
		mycol <- c(mycol, phen.cmap, u.pathway.mut.grey)
		if( length(tissue.names) > 1 ){
			#browser()
			tissue.colors = c(brewer.pal(12, "Set3"), brewer.pal(12,"Paired"))[1:length(tissue.names)]
			#row.names = c(row.names, "Tissue Types")
			mycol <- c(mycol, tissue.colors)
			n.rows.phen = n.rows.phen + 1
			heatm = rbind(heatm, (tissue.labels + tot.cols))
			tot.cols = tot.cols + length(tissue.colors)
			
		}
		#browser()
		
		
		#browser()
		image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
	}
	
# Add lines to separate phenotypes or subgroups
	
	if (col.labels2[1] != "NA") {
		groups <-  split(col.labels2, col.labels2)
		len.vec <- lapply(groups, length)
		plot.div <- c(0.51, cumsum(len.vec) + 0.5)
		for (i in plot.div) {
			lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
		}
		lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48), lwd = 2,
				cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50), lwd = 2,
				cex = 0.9, col = "black")
	}
	if (row.names[1] != "NA") {
#		browser()
		numC <- nchar(row.names)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names[i] <- substr(row.names[i], 1, 40)
			row.names[i] <- paste(row.names[i], " ", sep="")
		}
		if (phen.names[1] == "NA") {
			head.names <- paste("Class", seq(n.rows.phen, 1, -1))
		} else {
			head.names <- as.character(rev(phen.names))
		}
		row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
		if( length(tissue.names) > 1){ row.names = c(row.names, "Tissue Types")}
#            print(paste("n.rows:", n.rows))
#            print(paste("Phen names:", phen.names))
#            print(paste("Head names:", head.names))
#            print(paste("Row names:", row.names))
#		browser()
		axis(2, at=1:(n.rows + n.rows.phen), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char,
				font.axis=2, line=-1)
	}
	
	if (row.names2[1] != "NA") {
#		browser()
		numC <- nchar(row.names2)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names2[i] <- substr(row.names2[i], 1, 40)
			row.names2[i] <- paste(" ", row.names2[i], sep="")
			
		}
		for( i in 1:n.rows.phen ){
			phen.names2[i] <- substr(phen.names2[i], 1, 40)
			phen.names2[i] <- paste( " ", phen.names2[i], sep="")
		}
		
		row.names2 <- rev(row.names2)
		phen.names2 <- rev(phen.names2)
		axis(4, at=1:(n.rows + n.rows.phen), labels=c(row.names2, phen.names2), adj= 0.5, tick=FALSE, las = 1, 
				cex.axis=size.row.char, font.axis=2, line=-1)
	}
	
	if (col.names[1] != "NA") {
		size.col.char <- char.rescale*20/(n.cols + 25)
		axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
	}
	
	# Phenotype Legend 
	
#      print("--------------------------------------------------------------------------------------------")
	if (legend == T) {
		leg.txt <- NULL
		p.vec <- NULL
		c.vec <- NULL
		c2.vec <- NULL
		ind <- 1
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		for (i in 1:(n.rows.phen-1)) {  
#			browser()
			if (is.vector(col.labels)) {
				phen.v <- as.character(col.classes)
			} else {
				phen.v <- as.character(col.classes[[i]])
			}
			p.name <- paste(as.character(rev(head.names)[i]), ":   ", sep="")
			leg.txt <- c(p.name, phen.v)  
			p.vec <-  rep(22, cols.row[i] + 1)
			c.vec <-  c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
			c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
			ind <- ind + cols.row[i]
			offset <- 0.07
			legend(x=0, y= 1 - offset*i, 
					horiz = T, x.intersp = 0.5, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, 
					pt.bg = c.vec, col = c2.vec, cex = 1.20, pt.cex=1.75)
		}
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		legend(x=0, y= 10, horiz = T, x.intersp = 0.5, legend=tissue.names, bty="n", xjust=0, yjust= 1, 
				fill = tissue.colors, cex = 1.20, pt.cex=1.75, ncol=1)
	}
	#browser()
	## Tissue Legend
	if(length(tissue.names)>1){
		#browser()
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		legend(x=0, y= 1, #horiz = T, x.intersp = 0.5, y.intersp=.25, 
				legend=tissue.names, bty="n", xjust=0, yjust= 1, 
				fill = tissue.colors, #cex = 1.20, pt.cex=1.75, 
				ncol=4)
	}
	# Color map legend
	
#       print(c("range V=", range(V)))
#       print(c("range V1=", range(V1)))
#       print(c("range V2=", range(V2)))
	
	if(legend==TRUE){
		par(mar = c(2, 12, 2, 12))
		num.v <- 20
		range.v <- range(V2)
		incr <-  (range.v[1] - range.v[2])/(num.v - 1)
		heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
		image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=" ", sub = " ", xlab= ylab, ylab=xlab)
		range.v <- range(V1)
		incr <-  (range.v[1] - range.v[2])/(num.v - 1)
		heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
		#          print(c("heatm.v2=", heatm.v2))
		axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.5*char.rescale, font.axis=1)
	}
	
	return()
	
}


MSIG.Gct2Frame <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, blank.lines.skip=T, comment.char="", as.is=T, na.strings = "")
	descs <- ds[,1]
	ds <- ds[-1]
	row.names <- row.names(ds)
	names <- names(ds)
	return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}

Read.GeneSets.db <- function(
		gs.db,
		thres.min = 2,
		thres.max = 2000,
		gene.names = NULL)
{
	
	temp <- readLines(gs.db)
	max.Ng <- length(temp)
	temp.size.G <- vector(length = max.Ng, mode = "numeric") 
	for (i in 1:max.Ng) {
		temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
	}
	max.size.G <- max(temp.size.G)      
	gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
	temp.names <- vector(length = max.Ng, mode = "character")
	temp.desc <- vector(length = max.Ng, mode = "character")
	gs.count <- 1
	for (i in 1:max.Ng) {
		gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
		gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
		gene.set.name <- gs.line[1] 
		gene.set.desc <- gs.line[2] 
		gene.set.tags <- vector(length = gene.set.size, mode = "character")
		for (j in 1:gene.set.size) {
			gene.set.tags[j] <- gs.line[j + 2]
		}
		if (is.null(gene.names)) {
			existing.set <- rep(TRUE, length(gene.set.tags))
		} else {
			existing.set <- is.element(gene.set.tags, gene.names)
		}
		set.size <- length(existing.set[existing.set == T])
		if ((set.size < thres.min) || (set.size > thres.max)) next
		temp.size.G[gs.count] <- set.size
		gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
		temp.names[gs.count] <- gene.set.name
		temp.desc[gs.count] <- gene.set.desc
		gs.count <- gs.count + 1
	}
	Ng <- gs.count - 1
	gs.names <- vector(length = Ng, mode = "character")
	gs.desc <- vector(length = Ng, mode = "character")
	size.G <- vector(length = Ng, mode = "numeric") 
	
	gs.names <- temp.names[1:Ng]
	gs.desc <- temp.desc[1:Ng]
	size.G <- temp.size.G[1:Ng]
	
	return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = size.G, max.N.gs = max.Ng))
}

write.cls.2 <- function (class.v, phen, filename) 
{
	f <- file(filename, "w")
	n <- length(phen)
	l <- length(class.v)
	cat(l, n, "1", "\n", file = f, append = TRUE, sep = " ")
	cat("#", unlist(phen), "\n", file = f, append = TRUE, sep = " ")
	if (is.vector(class.v)) {
		class.v <- phen[class.v]
		cat(class.v, "\n", file = f, append = TRUE, sep = " ")
	} else {
		class.list <- matrix(0, nrow=length(class.v[,1]), ncol=length(class.v[1,]))
		for (i in 1:length(class.v[,1])) {
			class.list[i,] <- unlist(phen[[i]])[class.v[i,]]
			cat(class.list[i,], "\n", file = f, append = TRUE, sep = " ")
		}
	}
	close(f)
}

write.gct <- function(gct.data.frame, descs = "", filename) 
{
	f <- file(filename, "w")
	cat("#1.2", "\n", file = f, append = TRUE, sep = "")
	cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
	cat("Name", "\t", file = f, append = TRUE, sep = "")
	cat("Description", file = f, append = TRUE, sep = "")
	
	names <- names(gct.data.frame)
	cat("\t", names[1], file = f, append = TRUE, sep = "")
	
	if (length(names) > 1) {
		for (j in 2:length(names)) {
			cat("\t", names[j], file = f, append = TRUE, sep = "")
		}
	}
	cat("\n", file = f, append = TRUE, sep = "\t")
	
	oldWarn <- options(warn = -1)
	m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
	m[, 1] <- row.names(gct.data.frame)
	if (length(descs) > 1) {
		m[, 2] <- descs
	} else {
		m[, 2] <- row.names(gct.data.frame)
	}
	index <- 3
	for (i in 1:dim(gct.data.frame)[2]) {
		m[, index] <- gct.data.frame[, i]
		index <- index + 1
	}
	write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
	close(f)
	options(warn = 0)
	
}

MSIG.ReadPhenFile <- function(file = "NULL") {
#
# Reads a matrix of class vectors from a CLS file and defines phenotype and class labels vectors
#  (numeric and character) for the samples in a gene expression file (RES or GCT format)
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	cls.cont <- readLines(file)
	num.lines <- length(cls.cont)
	temp <- unlist(strsplit(cls.cont[[1]], " "))
	if (length(temp) == 3) {
		phen.names <- NULL
		col.phen <- NULL
	} else {
		l.phen.names <- match("phen.names:", temp)
		l.col.phen <- match("col.phen:", temp)
		phen.names <- temp[(l.phen.names + 1):(l.col.phen - 1)]
		col.phen <- temp[(l.col.phen + 1):length(temp)]
	}
	temp <- unlist(strsplit(cls.cont[[2]], " "))
	phen.list <- temp[2:length(temp)]
	
	for (k in 1:(num.lines - 2)) {
		temp <- unlist(strsplit(cls.cont[[k + 2]], " "))
		if (k == 1) {
			len <- length(temp)
			class.list <- matrix(0, nrow = num.lines - 2, ncol = len)
			class.v <- matrix(0, nrow = num.lines - 2, ncol = len)
			phen <- list(NULL)
		}
		class.list[k, ] <- temp
		classes <- unique(temp)
		class.v[k, ] <- match(temp, classes)
		phen[[k]] <- classes
	}
	if (num.lines == 3) {
		class.list <- as.vector(class.list)
		class.v <- as.vector(class.v)
		phen <- unlist(phen)
	}
	return(list(phen.list = phen.list, phen = phen, phen.names = phen.names, col.phen = col.phen,
					class.v = class.v, class.list = class.list))
}

MSIG.ReadPhenFile.2 <- function(file = "NULL") { 
#
# Reads a matrix of class vectors from a CLS file and defines phenotype and class labels vectors
#  (numeric and character) for the samples in a gene expression file (RES or GCT format)
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	cls.cont <- readLines(file)
	num.lines <- length(cls.cont)
	temp <- unlist(strsplit(cls.cont[[1]], " "))
	if (length(temp) == 3) {
		phen.names <- NULL
		col.phen <- NULL
	} else {
#		browser()
		l.phen.names <- match("phen.names:", temp)
		l.col.phen <- match("col.phen:", temp)
		phen.names <- temp[(l.phen.names + 1):(l.col.phen - 1)]
		col.phen <- temp[(l.col.phen + 1):length(temp)]
	}
	temp <- unlist(strsplit(cls.cont[[2]], " "))
	phen.list <- temp[2:length(temp)]
	
	phen <- NULL
	for (k in 1:(num.lines - 2)) {
		temp <- unlist(strsplit(cls.cont[[k + 2]], " "))
		if (k == 1) {
			len <- length(temp)
			class.list <- matrix(0, nrow = num.lines - 2, ncol = len)
			class.v <- matrix(0, nrow = num.lines - 2, ncol = len)
#           phen <- NULL
		}
		class.list[k, ] <- temp
		classes <- unique(temp)
		class.v[k, ] <- match(temp, classes)
#        phen[[k]] <- classes
		phen <- c(phen, classes)
	}
	if (num.lines == 3) {
		class.list <- as.vector(class.list)
		class.v <- as.vector(class.v)
#         phen <- unlist(phen)
	}
	return(list(phen.list = phen.list, phen = phen, phen.names = phen.names, col.phen = col.phen,
					class.v = class.v, class.list = class.list))
}

MSIG.Subset.Dataset.2 <- function(
		input.ds,
		input.cls = NULL,
		column.subset = "ALL",    # subset of column numbers or names (or phenotypes)
		column.sel.type = "samples",  # "samples" or "phenotype"
		row.subset = "ALL",       # subset of row numbers or names
		output.ds,
		output.cls = NULL) {
	
# start of methodology
	
	print(c("Running MSIG.Subset.Dataset... on GCT file:", input.ds))
	print(c("Running MSIG.Subset.Dataset... on CLS file:", input.cls))
	
# Read input datasets
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)
	m <- data.matrix(dataset$ds)
	gs.names <- dataset$row.names
	gs.descs <- dataset$descs
	sample.names <- dataset$names
	
# Read CLS file
	
	if (!is.null(input.cls)) {
		CLS <- MSIG.ReadPhenFile.2(file=input.cls)
		class.labels <- CLS$class.v
		class.phen <- CLS$phen
		class.list <- CLS$class.list 
	}
	
# Select desired column subset
	
	if (column.sel.type == "samples") {
		if (column.subset[1] == "ALL") {
			m2 <- m
			sample.names2 <- sample.names
			if (!is.null(input.cls)) {
				class.labels2 <- class.labels
			}
		} else {
			if (is.numeric(column.subset[1])) {
				m2 <- m[,column.subset]
				sample.names2 <- sample.names[column.subset]
				if (!is.null(input.cls)) {
					if (is.vector(class.labels)) {
						class.labels2 <- class.labels[column.subset]
					} else {
						class.labels2 <- class.labels[, column.subset]
					}
				}
			} else {
				locations <- !is.na(match(sample.names, column.subset))
				sample.names2 <- sample.names[locations]
				m2 <- m[, locations]
				if (!is.null(input.cls)) {
					if (is.vector(class.labels)) {
						class.labels2 <- class.labels[locations]
					} else {
						class.labels2 <- class.labels[, locations]
					}
				}
			}
		}
	} else if (column.sel.type == "phenotype") {
		locations <- !is.na(match(class.list, column.subset))
		sample.names2 <- sample.names[locations]
		m2 <- m[,locations]
		if (!is.null(input.cls)) {
			if (is.vector(class.labels)) {
				class.labels2 <- class.labels[locations]
			} else {
				class.labels2 <- class.labels[, locations]
			}
		}
	}
	
	if (row.subset[1] == "ALL") {
		m3 <- m2
		gs.names2 <- gs.names
		gs.descs2 <- gs.descs
	} else {
		locations <- !is.na(match(gs.names, row.subset))
		m3 <- m2[locations,]
		gs.names2 <- gs.names[locations]
		gs.descs2 <- gs.descs[locations]
	}
	
# Save datasets
	
	V <- data.frame(m3)
	names(V) <- sample.names2
	row.names(V) <- gs.names2
	write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)  
	
	if (!is.null(input.cls)) {
		write.cls.2(class.v = class.labels2, phen = class.phen, filename = output.cls) 
	}
}

OPAM.match.projection.to.pathway  <- function(
		input.ds,
		input.cls          = NA,
		results.dir,
		normalize.score    = F,
		normalization.type = "zero.one",
		pathway,
		max.n              = 10,
		user.colors        = NA,
		decreasing.order   = T,
		sort.columns       = F,
		char.rescale       = 1.25,
		cmap.type          = 3,
		row.norm           = T,
		output.dataset     = NA)
{
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	pathway.names <- dataset$row.names
	pathway.descs <- dataset$descs
	Ns <- length(m[1,])
	dim(m)
	sample.names <- dataset$names
	
	n.pathways <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
#   char.res <-  0.013 * n.pathways + 0.65
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			for (i in 1:n.pathways) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:n.pathways) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:n.pathways) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
	
	loc <- match(pathway, pathway.names)
	print(c("loc:", loc))
	if (sort.columns == T) {
		s.order <- order(m[loc,], decreasing = decreasing.order)
		m2 <- m[, s.order]
		sample.names2 <- sample.names[s.order]
	} else {
		m2 <- m
		sample.names2 <- sample.names
	}
	correl <- cor(t(m2))[, loc]
	m.order <- order(correl, decreasing=T)
	correl2 <- correl[m.order]
	m2 <- m2[m.order[1:max.n],]
	pathway.names2 <- pathway.names[m.order]
	pathway.descs2 <- signif(correl2, digits=3)
	
	if (input.cls == "NA") {
		cls.labels2 <- c(rep(0, 10), rep(1, length(sample.names2) - 10))
		cls.phen2 <- c(" ")
		colors.list <- c("white")
		phen.names2 <- "    "
	} else {
		CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
		cls.labels <- CLS$class.v
		cls.phen <- CLS$phen
		cls.list <- CLS$class.list 
		if (!is.null(CLS$phen.names)) {
			phen.names <- CLS$phen.names
		} else {
			phen.names <- "  "
		}
		if (is.vector(cls.labels)) {
			if (sort.columns == T) {
				cls.labels2 <- cls.labels[s.order]
				cls.list2 <- cls.list[s.order]
			} else {
				cls.labels2 <- cls.labels
				cls.list2 <- cls.list
			}
			n.phen <- 1
		} else {
			if (sort.columns == T) {
				cls.labels2 <- cls.labels[, s.order]
				cls.list2 <- cls.list[, s.order]
			} else {
				cls.labels2 <- cls.labels
				cls.list2 <- cls.list
			}
			n.phen <- length(cls.labels2[,1])
		}
		cls.phen2 <- list(NULL)
		if (is.vector(cls.labels2)) {
			classes <- unique(cls.list2)
			cls.phen2 <- classes
			cls.labels2 <- match(cls.list2, cls.phen2)
		} else {
			for (kk in 1:length(cls.list2[, 1])) {
				classes <- unique(cls.list2[kk,])
				cls.phen2[[kk]] <- classes
				cls.labels2[kk,] <- match(cls.list2[kk,], cls.phen2[[kk]])
			}
		}
		phen.names2 <- phen.names
		if (!is.na(user.colors[1])) {
			c.test <- user.colors
		} else {
			if (!is.null(CLS$col.phen)) {
				c.test <- CLS$col.phen
			} else {
				c.test <- c(brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Set1"),
						brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
						brewer.pal(n=8, name="BuGn"),
						brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
						brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
						brewer.pal(n=8, name="BuGn"),
						brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
						brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
						brewer.pal(n=8, name="BuGn"))
			}
		}
	}
	cls.phen.index <- unlist(cls.phen2)
	colors.list <- c.test[1:length(cls.phen.index)]
	
	filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ.TO.", pathway, sep="")
	pdf(file=paste(filename, ".pdf", sep=""), height = 8.5, width = 10.5)
	
	MSIG.HeatMapPlot.7(V = m2, row.names = pathway.names2[1:max.n],
			row.names2 = pathway.descs2[1:max.n], col.labels = cls.labels2, 
			col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names,
			col.names = sample.names2, main = " ", xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=T)
	dev.off()
	
	if (!is.na(output.dataset)) {
		V.GCT <- m2
		colnames(V.GCT) <- sample.names2
		row.names(V.GCT) <- pathway.names2
		write.gct(gct.data.frame = V.GCT, descs = pathway.descs2, filename =output.dataset)  
	}
	
}

MSIG.Define.Dataset.from.Table2 <- function(
		input.gct,
		table.txt,
		output.gct,
		output.txt = NULL,  # optional version of table with overlap (GCT & TAB) samples
		output.cls,
		prefix_entries = F)
{
# Read input dataset
	
	library(RColorBrewer)
	
	dataset1 <- MSIG.Gct2Frame(filename = input.gct)
	m <- data.matrix(dataset1$ds)
	gene.names <- dataset1$row.names
	gene.decs  <- dataset1$descs
	sample.names.gct <- dataset1$names
	Ns <- length(sample.names.gct)
	
# Read Table 
	
	tab <- read.delim(table.txt, header=T, row.names = 1, sep="\t", skip=0, blank.lines.skip=T, comment.char="", as.is=T)
	sample.names.tab <- row.names(tab)
	phen.names <- names(tab)
	
	overlap <- intersect(sample.names.tab, sample.names.gct)
	print("sample names GCT")
	print(sample.names.gct)
	print("sample names TAB")
	print(sample.names.tab)
	
	locs.gct <- match(overlap, sample.names.gct)
	print(match(sample.names.tab, sample.names.gct))
	print(match(sample.names.gct, sample.names.tab))
	
	locs.tab <- match(overlap, sample.names.tab)
	print(locs.tab)
	print(c("GCT matching set (", length(locs.gct), " samples):", sample.names.gct[locs.gct]))
	print(c("TAB matching set (", length(overlap), " samples):", sample.names.tab[locs.tab]))
	print(c("overlap set (", length(overlap), " samples):", overlap))
	
	m2 <- m[, locs.gct]
	sample.names.gct <- sample.names.gct[locs.gct]
	sample.names.tab <- sample.names.tab[locs.tab]
	
	if (!is.null(output.txt)) {
		tab2 <- tab[locs.tab,]
		sample.names.tab2 <- sample.names.tab[locs.tab]
		col.names <- paste(colnames(tab2), collapse = "\t")
		col.names <- paste("SAMPLE", col.names, sep= "\t")
		write(noquote(col.names), file = output.txt, append = F, ncolumns = length(col.names))
		write.table(tab2, file=output.txt, quote=F, col.names = F, row.names = T, append = T, sep="\t")
	}
	
	cls.table <- t(tab[locs.tab,])
	
	if (prefix_entries == TRUE) {
		for (i in 1:length(cls.table[,1])) {
#        cls.table[i,] <- paste(row.names(cls.table)[i], cls.table[i,], sep=".")
			cls.table[i,] <- paste(colnames(tab)[i], tab[,i], sep=".")
		}
	}
	
	if (!is.null(output.gct)) {      
		V <- data.frame(m2)
		names(V) <- sample.names.gct
		row.names(V) <- gene.names
		write.gct(gct.data.frame = V, descs = gene.decs, filename = output.gct)
	}
	
	class.phen <- unique(cls.table)
	n <- length(class.phen)
	l <- length(cls.table[1,])
	
	col.list <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"),
			brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"),
			brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"))
	num <- 0
	class.order.list <- NULL
	for (i in 1:length(cls.table[,1])) {
		num <- num + length(unique(cls.table[i,]))
		class.order.list <- c(class.order.list, unique(cls.table[i,]))
	}
	
	phen.names.string <- paste("phen.names:", paste(phen.names, collapse=" "), sep=" ")
	sig.col <- col.list[1:num]
	col.phen.string <- paste("col.phen:", paste(sig.col, collapse=" "), sep=" ")
	cat(paste(l, num, length(cls.table[, 1]), phen.names.string, col.phen.string, sep=" "), "\n", 
			file = output.cls, append = FALSE, sep = "")
	cat("# ", paste(class.order.list, collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	for (i in 1:length(cls.table[,1])) {
		cat(paste(cls.table[i,], collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	}
}

MSIG.Define.Dataset.from.Table.2 <- function(
		input.gct,
		table.txt,
		output.gct,
		output.cls,
		prefix_entries = F)
{
# Read input dataset
	
	library(RColorBrewer)
	
	dataset1 <- MSIG.Gct2Frame(filename = input.gct)
	m <- data.matrix(dataset1$ds)
	gene.names <- dataset1$row.names
	gene.decs  <- dataset1$descs
	sample.names.gct <- dataset1$names
	Ns <- length(sample.names.gct)
	
#	browser()
# Read Table 
	
	
	tab <- read.delim(table.txt, header=T, row.names = 1, sep="\t", skip=0, blank.lines.skip=T, comment.char="", as.is=T)
	sample.names.tab <- row.names(tab)
	phen.names <- names(tab)
	overlap <- intersect(sample.names.tab, sample.names.gct)
	if(length(overlap)==0){ return(NULL)}
#	print(overlap)
#	print("sample names GCT")
#	print(sample.names.gct)
#	print("sample names TAB")
#	print(sample.names.tab)
	
	if(length(overlap)==0){ return(NULL)}
	
	locs.gct <- match(overlap, sample.names.gct)
	print(match(sample.names.tab, sample.names.gct))
	print(match(sample.names.gct, sample.names.tab))
	locs.tab <- match(overlap, sample.names.tab)
#	print(locs.tab)
#	print(c("GCT matching set (", length(locs.gct), " samples):", sample.names.gct[locs.gct]))
#	print(c("TAB matching set (", length(overlap), " samples):", sample.names.tab[locs.tab]))
#	print(c("overlap set (", length(overlap), " samples):", overlap))
	
	
	
	m2 <- m[, locs.gct]
	sample.names.gct <- sample.names.gct[locs.gct]
	sample.names.tab <- sample.names.tab[locs.tab]
	cls.table <- t(tab[locs.tab,])
	
	if (prefix_entries == TRUE) {
		for (i in 1:length(cls.table[,1])) {
#        cls.table[i,] <- paste(row.names(cls.table)[i], cls.table[i,], sep=".")
			cls.table[i,] <- paste(colnames(tab)[i], tab[,i], sep=".")
		}
	}
	
	if (!is.null(output.gct)) {      
		V <- data.frame(m2)
		names(V) <- sample.names.gct
		row.names(V) <- gene.names
		write.gct(gct.data.frame = V, descs = gene.decs, filename = output.gct)
	}
	
	class.phen <- unique(cls.table)
	n <- length(class.phen)
	l <- length(cls.table[1,])
	
	col.list <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"),
			brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"),
			brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"))
#	num <- 0
#	class.order.list <- NULL
	class.order.list = apply(cls.table, 1, unique)
	num = sum(unlist(lapply(class.order.list, length)))
	class.order.list = unlist(class.order.list)
	
	
#	for (i in 1:length(cls.table[,1])) {
	##		num <- num + length(unique(cls.table[i,]))
#		class.order.list <- c(class.order.list, unique(cls.table[i,]))
#	}
	
	phen.names.string <- paste("phen.names:", paste(phen.names, collapse=" "), sep=" ")
	sig.col <- col.list[1:num]
	col.phen.string <- paste("col.phen:", paste(sig.col, collapse=" "), sep=" ")
	cat(paste(l, num, length(cls.table[, 1]), phen.names.string, col.phen.string, sep=" "), "\n", 
			file = output.cls, append = FALSE, sep = "")
	cat("# ", paste(class.order.list, collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	for (i in 1:length(cls.table[,1])) {
		cat(paste(cls.table[i,], collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	}
}


rec.area <- function(
		obs,
		pred,
		metric = "absolute.deviation",       # Either "squared.error" or "absolute.deviation"
#		null.distribution = "gaussian",  # Either "gaussian" [null.model = mean(obs)] or "laplacian" [null.model = median(obs)]
		interval = 0.01
){
#	browser()
	error.windows = seq(0, 1, by=interval)
	n.errors = length(error.windows)
	intervals = rep( interval, n.errors )
	n.obs = length(obs)
	n.pred = length(pred)
	if( n.obs != n.pred ){ stop( "The number of observations does not equal the number of predictions." ) }
#	if( null.distribution == "gaussian" ){ null.model = mean(obs) 
#	} else if( null.distribution == "laplacian" ){ null.model = median(obs) }
	
	if( metric == "squared.error" ){
		difference = (obs-pred)^2
		accuracy = unlist(lapply(error.windows, FUN=squared.error, difference, n.obs))
	} else if( metric == "absolute.deviation" ){
		difference = abs(obs-pred)
		accuracy = unlist(lapply(error.windows, FUN=absolute.deviation, difference, n.obs))
	}
#	plot(accuracy, type="l"); par(new=TRUE); plot(error.windows, type="l")
	
	triangle.heights = accuracy - c(0, accuracy[1:(n.errors-1)])
	triangles = triangle.heights*intervals/2
	rectangle.heights = c(0, accuracy[1:(n.errors-1)])
	rectangles = rectangle.heights*intervals
#	A = (cumsum(accuracy)*intervals)[n.errors]
	A = sum( rectangles + triangles)
	
	# Calculate p-value using Kolmogorov-Smirnov Test
#	Dn = max(accuracy-error.windows)
#	i = 1:100
#	x = sqrt( n.obs*n.pred/(n.obs+n.pred) )*Dn
#	p.value = 1 - (sqrt(2*pi)/x)*sum( exp(-(2*i - 1)^2 * pi^2/ (8*x^2)) )
#	browser()
#	pred.scrambled = sample(pred)
#	difference.scrambled = abs(pred.scrambled - obs)
#	accuracy.scrambled = unlist(lapply(error.windows, FUN=squared.error, difference.scrambled, n.obs))
#	triangle.heights = accuracy.scrambled - c(0, accuracy.scrambled[1:(n.errors-1)])
#	triangles = triangle.heights*intervals/2
#	rectangle.heights = c(0, accuracy.scrambled[1:(n.errors-1)])
#	rectangles = rectangle.heights*intervals
#	A.scrambled = sum( rectangles + triangles)
#	T2.scrambled = .5*(sum((accuracy.scrambled-error.windows)^2))
#	p.value.scrambled = cvmts.pval(T2.scrambled, n.errors, n.errors)
	
	# Calculate p-value using Cramer-Von-Mises Criterion
	T2 = .25*(sum((accuracy-error.windows)^2))  # accuracy-error.windows = integral difference between null model and REC
#	browser()
	p.value = cvmts.pval(T2, n.errors, n.errors)
	
	
#	T2.norm = (T2- min(T2))/(max(T2)-min(T2))
#	wilcox.test(T2.norm, error.)
#browser()
#	U =  2*n.errors^2*(T2 + (4*n.errors^2-1)/12*n.errors)
#	p.value.u = cvmts.pval(U, n.errors, n.errors)
	
#	print('calculating REC...')
#	browser()
#	rec.list.ccle[master.ind] <<- A
#	p.value.list.ccle[master.ind] <<- p.value
#	T2.list.ccle[master.ind] <<- T2
#	
#	
#	
#	rec.list.scrambled[master.ind] <<- A.scrambled
#	p.value.list.scrambled[master.ind] <<- p.value.scrambled
#	T2.list.scrambled[master.ind] <<- T2.scrambled
#	master.ind <<- master.ind+1
#	browser()
#	stat = ks.test(accuracy, error.windows, exact=TRUE)
	return( list(A = A, p.value = p.value, T2=T2) )
}

squared.error <- function( error, squared.difference, n ){
	return( length(which( squared.difference <= error ))/n )
}

absolute.deviation <- function( error, absolute.difference, n ) {
	return( length(which( absolute.difference <= error ))/n )
}

mutual.inf <- function(x, y, n.grid=100) {
	
	kde2d.xy <- kde2d(x, y, n = n.grid, h = c(width.SJ(x, method="dpi"), width.SJ(y, method="dpi")))
#	X <- kde2d.xy$x
#	Y <- kde2d.xy$y
	PXY <- kde2d.xy$z/sum(kde2d.xy$z)
	
	PX <- apply(PXY, MARGIN=1, sum)
	PX <- PX/sum(PX)
	PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
	
	PY <- apply(PXY, MARGIN=2, sum)
	PY <- PY/sum(PY)
	PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
	
	MI <- sum(PXY * log2(PXY/(PX*PY)))
#	browser()
	return(MI)
}

mutual.inf.2 <- function(x, y, n.grid=100, normalize.by ="HXY", # Whether to normalize by HXY, HX, or HY
		pos.and.neg = T

) {
	# x and y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	if( length(unique(x)) == 1 || length(unique(y)) == 1 ){
#		browser()
		return( NA )
	}
#	bandwidth.x = ifelse(IQR(x) == 0, bcv(x, n.grid), width.SJ(x, method="dpi"))
#	bandwidth.y = ifelse(IQR(y) == 0, bcv(y, n.grid), width.SJ(y, method="dpi"))
#	print("---")
#	print(x)
#	print(y)
	kde2d.xy <- kde2d(x, y, n = n.grid, h = c(suppressWarnings(bcv(x)), suppressWarnings(bcv(y))) )
#	X <- kde2d.xy$x
#	Y <- kde2d.xy$y
#	Z = kde2d.xy$z
	PXY <- kde2d.xy$z/sum(kde2d.xy$z)
	
	PX <- apply(PXY, MARGIN=1, sum)
	PX <- PX/sum(PX)
	HX = -sum(PX * log2(PX))
	PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
	
	PY <- apply(PXY, MARGIN=2, sum)
	PY <- PY/sum(PY)
	HY = -sum( PY * log2(PY))
	PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
	
#	browser()
	
#	MIXY = PXY * log2(PXY/(PX*PY))
#	
#	if( pos.and.neg ){
#	q1 = MIXY[1:(n.grid/2), 1:(n.grid/2)]
#	q2 = MIXY[1:(n.grid/2), (n.grid/2 + 1):n.grid]
#	q3 = MIXY[(n.grid/2+1):n.grid, 1:(n.grid/2)]
#	q4 = MIXY[(n.grid/2+1):n.grid, (n.grid/2+1):n.grid]
#	
#	# q's divide MIXY into quarters. If the sum of q2 and q3 is greater than the sum of q1 and q4, then
#	# x and y are negatively correlated.
#	# on heatmap:   q2  q4
#	#               q1  q3
#	
	## Ignore NaN's that are a result of underflow (experimentally derived)
#	MI <- ifelse( sum(q1+q4, na.rm=TRUE) < sum(q2+q3, na.rm=TRUE), 
#			-sum(MIXY, na.rm=TRUE), sum(MIXY, na.rm=TRUE))
#} else{ MI = sum(MIXY, na.rm=TRUE)}
#	MI <- ifelse( sum(q1+q4, na.rm=TRUE) < sum(q2+q3, na.rm=TRUE), -sum(q2+q3-q1-q4, na.rm=TRUE), sum(q1+q4-q2-q3, na.rm=TRUE))
	HXY <- - sum(PXY * log2(PXY), na.rm=TRUE)
	
#	HX = -sum( PX * log2(PX) )
#	HY = -sum( PY * log2(PY) )
#	MI.norm = (HX+HY)/HXY
#browser()
#	normalization.factor = 1 #ifelse(normalize.by=="HXY", HXY, ifelse(normalize.by=="HX", HX, HY))
	## browser()
	MI.norm =  ifelse(pos.and.neg, sign(cor(x, y)), 1) * ((HX + HY)/HXY - 1) #MI/normalization.factor
#	browser()
	
	return( MI.norm )#list(MI=MI, HXY=HXY))
}

#mutual.inf.2.single.gene.target <- function( signature, gene.target){
#	return(mutual.inf.2(gene.target, signature))
#}

mutual.inf.2.multiple.gene.targets <- function( signature, gene.targets ){
	return(apply(gene.targets, 
					MARGIN=1, FUN=mutual.inf.2, 
					signature ) )
}

mutual.inf.3 <- function(gene.target, signature.matrix, signature.indices, n.grid=100, gene.target.name = "",
		n.randomizations = 100, tissue = "NA") {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the indices of signatures that you are interested in using.
	# This code is used in comparing the chosen signatures to SUMMARY, the mutation
	# status of all the cell lines. 
#	browser()
	n.signatures = length(signature.matrix[,1])
	MI.vector = vector(length=n.signatures, mode="double")
#	MI.vector.rand = vector(length=length(signature.matrix[,1]), mode="double")
	gene.target.rand = t(replicate(n.randomizations, sample(gene.target)))
	MI.matrix.rand = matrix(ncol = n.signatures, nrow = n.randomizations)
#	browser()
#	for( i in 1:length(signature.matrix[,1]) ){
	MI.vector = apply(signature.matrix, MARGIN=1, FUN=mutual.inf.2, gene.target)
#		browser()
#		for( j in 1:n.randomizations ){
	MI.matrix.rand = apply(signature.matrix, 
			MARGIN=1, FUN=mutual.inf.2.multiple.gene.targets, 
			gene.target.rand)
	#mutual.inf.2(gene.target.rand[j,], signature.matrix[i,])$MI
#		}
#		browser()
#		MI.vector.rand[i] = mean(temp.MI.rand)
#	}
#	x.rand = sample(x)
#	print("Make plot of densities!! And save the output!")
#	browser()
	
	quartz()
#	if( gene.target.name =="SUMMARY"){
	temp <- density(MI.vector, adjust=1, n = 512, from=min(MI.vector), to=max(MI.vector))
	x <- temp$x
	y <- temp$y/sum(temp$y)
	
	temp.rand <- density(MI.matrix.rand, adjust=1, n = 512, from=min(MI.matrix.rand), to=max(MI.matrix.rand))
	x.rand <- temp.rand$x
	y.rand <- temp.rand$y/sum(temp.rand$y)
#		pdf(file=paste(tissue, gene.target.name, n.randomizations, "pdf", sep=".") )
#		quartz(file=paste(tissue, gene.target.name, n.randomizations, "pdf", sep="."))
	plot(x.rand, y.rand, type="l", lwd=2, xlab="MI", #xlim = c(max(min(x), 10^-5), max(x)), ylim = range(c(y, y.rand)), 
			col="red", 
			ylab = "P(MI)", main=paste(tissue, gene.target.name, n.randomizations, sep="  "))
	points(x, y, type="l", lwd=2, col = "black")
	legend("topright", c("actual gene target vs all gene sets", "randomized gene target vs all gene sets"), col=c("black", "red"), lwd=c(2,2))
	browser()
#		dev.off()
#	}
	
	MI = MI.vector[signature.indices]
	FDR = vector(length=length(MI))
#	browser()
	ranked.MI.vector = rank(-MI.vector)  # take negative so rank 1 corresponds to highest value
	ranked.MI.matrix.rand = rank(MI.matrix.rand)
	median.MI.rand = median(MI.matrix.rand)
	if( gene.target.name[1] == "EGFR_AMP" #|| gene.target.name =="TP53"
			){ browser() }
	for( i in 1:length(MI)){
		if( MI[i] > median.MI.rand ){
			rank.observed = ranked.MI.vector[signature.indices[i]]
			rank.randomized = sum(MI[i] < MI.matrix.rand)
		} else{ 
			
			rank.observed = n.signatures - ranked.MI.vector[signature.indices[i]] + 1
			rank.randomized = sum(MI[i] > MI.matrix.rand)
#			browser()
		}
		
		FDR[i] = (rank.randomized/n.randomizations)/rank.observed
#		if( MI[i] <= median.MI.rand ){ browser() }
	}
	
	
#	MI.rand.ind = which(x.rand >= MI)
	
#	MI.integral = sum(x[MI.ind]*y[MI.ind], na.rm=T)
#	MI.rand.integral = sum(x.rand[MI.ind]*y.rand[MI.ind], na.rm=T)
#	FDR = MI.rand.integral/MI.integral
#	browser()
	
	return(list(MI=MI, FDR=FDR))
}

mutual.inf.3.v2 <- function(target.vector, comparison.matrix, n.grid=100, target.vector.name = "",
		tissue = "NA", normalize.by = "HXY", pos.and.neg=T) {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the indices of signatures that you are interested in using.
	# This code is used in comparing the chosen signatures to SUMMARY, the mutation
	# status of all the cell lines. 
	
	n.signatures = length(comparison.matrix[,1])
	MI.vector = vector(length=n.signatures, mode="double")
	
	MI.ref = mutual.inf.2(target.vector, target.vector, normalize.by=normalize.by)
	print(paste("MI.ref =", MI.ref))
	MI.vector = apply(comparison.matrix, MARGIN=1, FUN=mutual.inf.2, target.vector, normalize.by=normalize.by, 
			pos.and.neg=pos.and.neg)
	MI = MI.vector/MI.ref
	FDR = rep(1, length=length(MI))
	
	return(list(MI=MI, FDR=FDR))
}



mutual.inf.4 <- function(gene.targets, signature.matrix, signature.index, n.grid=100, gene.target.name = "",
		n.randomizations = 100) {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the index of the "winning" signature that is used to compare to
	# all the genomic aberrations.
#	browser()
	if( is.vector(gene.targets)){
		gene.targets = t(as.matrix(gene.targets))
	}
	n.gene.targets = length(gene.targets[,1])
	n.signatures = length(signature.matrix[,1])
	n.samples = length(gene.targets[1,])
	MI.matrix = matrix(ncol = n.signatures, nrow = n.gene.targets)
	MI.matrix.rand = matrix(ncol = n.signatures, nrow = n.randomizations)
	gene.target.rand = t(replicate(n.randomizations, sample(gene.targets[1,])))
#	temp.MI.rand = vector(length=n.iter)
#	browser()
#	for( i in 1:length(signature.matrix[,1]) ){
#		browser()
	MI.matrix = apply(signature.matrix, MARGIN=1, 
			FUN=mutual.inf.2.multiple.gene.targets, 
			gene.targets)
	MI.matrix.rand = apply(signature.matrix, 
			MARGIN=1, FUN=mutual.inf.2.multiple.gene.targets, 
			gene.target.rand)
#		MI.matrix.rand[i,] = apply(gene.target.rand, 
#						MARGIN=1, FUN=mutual.inf.2, 
#						signature.matrix[i,]) 
#		browser()
#		MI.vector.rand[i] = mean(temp.MI.rand)
#	}
#	x.rand = sample(x)
	
#	browser()
	quartz()
	temp <- density(MI.matrix, adjust=1, n = 512, from=min(MI.matrix), to=max(MI.matrix))
	x <- temp$x
	y <- temp$y/sum(temp$y)
	
	temp.rand <- density(MI.matrix.rand, adjust=1, n = 512, from=min(MI.matrix), to=max(MI.matrix))
	x.rand <- temp.rand$x
	y.rand <- temp.rand$y/sum(temp.rand$y)
	
#	pdf(file=paste(tissue, paste(gene.target.name, collapse="-"), rownames(signature.matrix)[signature.index], n.randomizations, "pdf", sep="."))
#	quartz(file=paste(tissue, paste(gene.target.name, collapse="-"), rownames(signature.matrix)[signature.index], n.randomizations, "pdf", sep="."))
#if( gene.gene.target.name =="SUMMARY"){
	plot(x.rand, y.rand, type="l", lwd=2, xlab="MI", #xlim = c(max(min(x), 10^-5), max(x)), ylim = range(c(y, y.rand)), 
			col="red", 
			ylab = "P(MI)", main=paste(tissue, paste(gene.target.name, collapse=" "), rownames(signature.matrix)[signature.index], n.randomizations, sep="  "))
	points(x, y, type="l", lwd=2, col = "black")
	legend("topright", c("actual gene target(s) vs all gene sets", "randomized gene target vs all gene sets"), col=c("black", "red"), lwd=c(2,2))
#	if( gene.target.name[1] =="KRAS_AMP") {browser()}
	browser()
#	dev.off()
#}
#
#	browser()
#	MI = ifelse(is.matrix(MI.matrix), MI.matrix[,signature.index], MI.matrix[signature.index])
#	ranked.MI.matrix = ifelse( is.matrix(MI.matrix), apply(-MI.matrix, MARGIN=1, rank), rank(-MI.matrix))
#	MI.vector = ifelse(is.matrix(MI.matrix))
	
	
	if(is.matrix(MI.matrix)){
		MI = MI.matrix[,signature.index]
		ranked.MI.matrix =  apply(-MI.matrix, MARGIN=1, rank)
		FDR = vector(length=n.gene.targets, mode="numeric")
#		browser()
		for( i in 1:n.gene.targets){
			if( MI[i] > median(MI.matrix.rand) ){
				rank.observed = ranked.MI.matrix[signature.index, i]
				rank.randomized = sum(MI[i] < MI.matrix.rand)
			} else{ 
#				browser()
				rank.observed = n.signatures - ranked.MI.matrix[signature.index, i]
				rank.randomized = sum(MI[i] > MI.matrix.rand)
			}
			FDR[i] = (rank.randomized/n.randomizations)/rank.observed
#		browser()
		}
	} else{ 
		MI = MI.matrix[signature.index]
		ranked.MI.matrix = rank(-MI.matrix)
		if( MI > median(MI.matrix.rand)){
			rank.observed = ranked.MI.matrix[signature.index]
			rank.randomized = sum(MI <= MI.matrix.rand)
		} else{
			rank.observed = n.signatures - ranked.MI.matrix[signature.index]
			rank.randomized = sum(MI >= MI.matrix.rand)
		}
		FDR = (rank.randomized/n.randomizations)/rank.observed
	}
#	if( gene.target.name == "EGFR_AMP" #|| gene.target.name =="TP53"
#			){ browser() }
	
#	if( MI > 0 ){
#		MI.ind = which(x >= MI) 
#	} else{ MI.ind = which( x <= MI) }
#	MI.rand.ind = which(x.rand >= MI)
	
#	MI.integral = sum(x[MI.ind]*y[MI.ind], na.rm=T)
#	MI.rand.integral = sum(x.rand[MI.ind]*y.rand[MI.ind], na.rm=T)
#	FDR = MI.rand.integral/MI.integral
#	browser()
	
	return(list(MI=MI, FDR=FDR))
}

mutual.inf.4.v2 <- function(gene.targets, signature.matrix, signature.index, n.grid=100, gene.target.name = "",
		n.randomizations = 100) {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the index of the "winning" signature that is used to compare to
	# all the genomic aberrations.
#	browser()
	if( is.vector(gene.targets)){
		gene.targets = t(as.matrix(gene.targets))
	}
	n.gene.targets = length(gene.targets[,1])
	n.signatures = length(signature.matrix[,1])
	n.samples = length(gene.targets[1,])
	MI.matrix = matrix(ncol = n.signatures, nrow = n.gene.targets)
	if( n.randomizations > 0 ){
		MI.matrix.rand = matrix(ncol = n.signatures, nrow = n.randomizations)
		gene.target.rand = t(replicate(n.randomizations, sample(gene.targets[1,])))
		MI.matrix.rand = apply(signature.matrix, 
				MARGIN=1, FUN=mutual.inf.2.multiple.gene.targets, 
				gene.target.rand)
		MI.matrix.rand = normalize(MI.matrix.rand)
	}
#	temp.MI.rand = vector(length=n.iter)
#	browser()
#	for( i in 1:length(signature.matrix[,1]) ){
#		browser()
	MI.matrix = apply(signature.matrix, MARGIN=1, 
			FUN=mutual.inf.2.multiple.gene.targets, 
			gene.targets)
	MI.matrix = normalize(MI.matrix)
#		MI.matrix.rand[i,] = apply(gene.target.rand, 
#						MARGIN=1, FUN=mutual.inf.2, 
#						signature.matrix[i,]) 
#		browser()
#		MI.vector.rand[i] = mean(temp.MI.rand)
#	}
#	x.rand = sample(x)
	
#	browser()
#	quartz()
#	temp <- density(MI.matrix, adjust=1, n = 512, from=min(MI.matrix), to=max(MI.matrix))
#	x <- temp$x
#	y <- temp$y/sum(temp$y)
#	
#	temp.rand <- density(MI.matrix.rand, adjust=1, n = 512, from=min(MI.matrix), to=max(MI.matrix))
#	x.rand <- temp.rand$x
#	y.rand <- temp.rand$y/sum(temp.rand$y)
#	
	##	pdf(file=paste(tissue, paste(gene.target.name, collapse="-"), rownames(signature.matrix)[signature.index], n.randomizations, "pdf", sep="."))
	##	quartz(file=paste(tissue, paste(gene.target.name, collapse="-"), rownames(signature.matrix)[signature.index], n.randomizations, "pdf", sep="."))
	##if( gene.gene.target.name =="SUMMARY"){
#	plot(x.rand, y.rand, type="l", lwd=2, xlab="MI", #xlim = c(max(min(x), 10^-5), max(x)), ylim = range(c(y, y.rand)), 
#			col="red", 
#			ylab = "P(MI)", main=paste(tissue, paste(gene.target.name, collapse=" "), rownames(signature.matrix)[signature.index], n.randomizations, sep="  "))
#	points(x, y, type="l", lwd=2, col = "black")
#	legend("topright", c("actual gene target(s) vs all gene sets", "randomized gene target vs all gene sets"), col=c("black", "red"), lwd=c(2,2))
	##	if( gene.target.name[1] =="KRAS_AMP") {browser()}
#	browser()
#	dev.off()
#}
#
#	browser()
#	MI = ifelse(is.matrix(MI.matrix), MI.matrix[,signature.index], MI.matrix[signature.index])
#	ranked.MI.matrix = ifelse( is.matrix(MI.matrix), apply(-MI.matrix, MARGIN=1, rank), rank(-MI.matrix))
#	MI.vector = ifelse(is.matrix(MI.matrix))
#	browser()
#MI.ref = mutual.inf.2( signature.matrix )
	
	if(is.matrix(MI.matrix)){
		MI = MI.matrix[,signature.index]
		if( n.randomizations > 0 ){
			ranked.MI.matrix =  apply(-MI.matrix, MARGIN=1, rank)
			FDR = vector(length=n.gene.targets, mode="numeric")
#		browser()
			for( i in 1:n.gene.targets){
				if( MI[i] > median(MI.matrix.rand) ){
					rank.observed = ranked.MI.matrix[signature.index, i]
					rank.randomized = sum(MI[i] < MI.matrix.rand)
				} else{ 
#				browser()
					rank.observed = n.signatures - ranked.MI.matrix[signature.index, i]
					rank.randomized = sum(MI[i] > MI.matrix.rand)
				}
				FDR[i] = (rank.randomized/n.randomizations)/rank.observed
#		browser()
			}
		} else{ FDR = rep(1, length=n.gene.targets) }
	} else{ 
		MI = MI.matrix[signature.index]
		if( n.randomizations > 0 ){
			ranked.MI.matrix = rank(-MI.matrix)
			if( MI > median(MI.matrix.rand)){
				rank.observed = ranked.MI.matrix[signature.index]
				rank.randomized = sum(MI <= MI.matrix.rand)
			} else{
				rank.observed = n.signatures - ranked.MI.matrix[signature.index]
				rank.randomized = sum(MI >= MI.matrix.rand)
			}
			FDR = (rank.randomized/n.randomizations)/rank.observed
		} else{ FDR = rep(1, length=n.gene.targets) }
	}
#	if( gene.target.name == "EGFR_AMP" #|| gene.target.name =="TP53"
#			){ browser() }
	
#	if( MI > 0 ){
#		MI.ind = which(x >= MI) 
#	} else{ MI.ind = which( x <= MI) }
#	MI.rand.ind = which(x.rand >= MI)
	
#	MI.integral = sum(x[MI.ind]*y[MI.ind], na.rm=T)
#	MI.rand.integral = sum(x.rand[MI.ind]*y.rand[MI.ind], na.rm=T)
#	FDR = MI.rand.integral/MI.integral
#	browser()
	
	return(list(MI=MI, FDR=FDR))
}



mise <- function( x ) {
	n.x = length(x)
	r = seq(2,10)
	f = vector(length=n.x, mode=mode(x))
	for( i in 1:n.x ){
		f = f + (x - x[i])^2
	}
	expected.value = sum(f)/n.x
}

amise <- function( v ){
	# Reference: 
	# "Very fast optimal bandwith selection for univariate kernel density estimation"
	# Vikas Chandrakant Raykar and Ramani Duraiswami
	# [CS-TR-4774/UMIACS-TR-2005-73] June 28, 2006
	
	H4 <- function( x ){ x^4 - 6*x^2 + 3 }
	H6 <- function( x ){ x^6 - 15*x^4 + 45*x^2 - 15 } 
	
	N = length(v)
	
	# Step 1 on page 11 of reference
	sigma = mean(v)
	
	# Step 2 on p. 11
	Phi6 = sigma^(-7)*(-15/(16*sqrt(pi)))
	Phi8 = sigma^(-9)*(-105/(32*sqrt(pi)))
	
	# Step 3 on p. 11
	g1 = ( -6/(sqrt(2*pi) * Phi6 * N))^(1/7)
	g2 = ( 30/(sqrt(2*pi) * Phi8 * N))^(1/9)
	
	# Make a matrix Z where Z(i,j) = x_i - x_j
	Z = matrix(v, ncol = N, nrow = N, byrow=TRUE) - matrix(v, ncol = N, nrow = N, byrow=FALSE)
	
	Phi4 <- function( g ) 1/(N*(N-1)*sqrt(2*pi)*g^5) * sum( H4(Z/g) * exp( -(Z^2)/(2*g^2)) )
	Phi4.g1 = Phi4(g1)
	Phi6.g2 = 1/(N*(N-1)*sqrt(2*pi)*g1^5) * sum( H4(Z/g1) * exp( -(Z^2)/(2*g1^2)) )
	
	# Step 4
	Y <- function( h ){
		( (-6*sqrt(2)*Phi4.g1)/Phi6.g2)^(1/7)*h^(5/7)
	}
	
	fxn <- function( h ){
		h - ( 1/ (sqrt(2) * Phi4(Y(h)) * N))^(1/5)
	}
	newtonraphson(fxn, mean(v))
}

#parzen.window <- function(z, h){
#	Sigma = cov(z,z)
#	
#	
#}

write.cls.with.locs <- function( output.cls,
		cls.labels,
		phen.names){
	
	class.order.list = apply(cls.labels, 1, unique)
	num = sum(unlist(lapply(class.order.list, length)))
	class.order.list = unlist(class.order.list)
	
	class.phen <- unique(cls.labels)
	n <- length(class.phen)
	l <- length(cls.list[1,])
	
	phen.names.string <- paste("phen.names:", paste(phen.names, collapse=" "), "col.names:", sep=" ")
	cat(paste(l, num, length(cls.list[, 1]), phen.names.string, sep=" "), "\n", 
			file = output.cls, append = FALSE, sep = "")
	cat("# ", paste(class.order.list, collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	for (i in 1:length(cls.list[,1])) {
		cat(paste(cls.list[i,], collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	}
	
}

normalize <- function( v ){
	(v - min(v))/(max(v) - min(v))
}


mutual.inf.P <- function(x, y, n.grid=100) {
    # for definitions of mutual information and the universal metric (NMI) see the 
    # definition of "Mutual Information" in wikipedia and Thomas and Cover's book

#   kde2d.xy <- kde2d(x, y, n = n.grid, h = c(width.SJ(x, method="dpi"), width.SJ(y, method="dpi")))
   kde2d.xy <- kde2d(x, y, n = n.grid, h = c(bcv(x), bcv(y)))
   X <- kde2d.xy$x
   Y <- kde2d.xy$y  
#   PXY <- kde2d.xy$z/sum(kde2d.xy$z)
   PXY <- kde2d.xy$z/sum(kde2d.xy$z) + .Machine$double.eps

#   PX <- apply(PXY, MARGIN=1, sum)
   PX <- rowSums(PXY)
   PX <- PX/sum(PX)
   HX <- -sum(PX * log2(PX))
   PX <- matrix(PX, nrow=n.grid, ncol=n.grid)

#   PY <- apply(PXY, MARGIN=2, sum)
   PY <- colSums(PXY)
   PY <- PY/sum(PY)
   HY <- -sum(PY * log2(PY))
   PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)

   MI <- sum(PXY * log2(PXY/(PX*PY)))
   MI
   HXY <- - sum(PXY * log2(PXY))
   NMI <- sign(cor(x, y)) * ((HX + HY)/HXY - 1)  # use peason correlation the get the sign (directionality)

   return(list(MI=MI, HXY=HXY, HX=HX, HY=HY, NMI=NMI))
}

OPAM.Evaluate.Results <- function(
   input.ds,
   input.cls,
   phenotype = NULL,
   target.class = NULL,
   target.type = "discrete",
   output.txt,
   output.pdf) {

   pdf(file=output.pdf, height=8.5, width=11)
   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   model.names <- dataset$row.names
   model.descs <- dataset$descs
   Ns <- length(m[1,])
   for (i in 1:length(m[,1])) {
      if (sd(m[i,]) == 0) {
         val <- m[i, 1]
	 m[i,] <- m[i,] + runif(n=Ns, min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
      }
   }

   dim(m)
   sample.names <- dataset$names
   CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
   cls.labels <- CLS$class.v
   cls.phen <- CLS$phen
   cls.list <- CLS$class.list 

   library(verification)

   if (is.null(phenotype)) {
      phen.loc <- 1
   } else {
      phen.loc <- match(phenotype, CLS$phen.names)
   }
   if (is.vector(CLS$class.list)) {
       target.vec <- CLS$class.list
   } else {
       target.vec <- CLS$class.list[phen.loc,]
   }
   if (target.type == "continuous") {
      target <- target.vec
   } else if (target.type == "discrete") {
      target <- ifelse(target.vec == target.class, 1, 0)    
   }

   ind <- order(target)
   target <- target[ind]
   target.vec <- target.vec[ind]
   m <- m[, ind]
   sample.names[ind]
   class.v <- CLS$class.v
   if (is.vector(class.v)) {
       class.v <- class.v[ind]
   } else {
       class.v <- class.v[, ind]
   }
   annot <- MI <- AUC <- AUC.pval <- t.stat <- t.pval <- vector(length=dim(m)[1], mode="numeric")

   NMI.ref <- mutual.inf.P(x = target, y = target, n.grid=100)$NMI

   for (i in 1:dim(m)[1]) {
      feature <- m[i,]
      MI[i] <- signif(mutual.inf.P(target, feature, n.grid=100)$NMI/NMI.ref, 4)
      if (target.type == "continuous") {
         AUC[i] <- AUC.pval[i] <- t.stat[i] <- t.pval[i] <- "-"
      } else if (target.type == "discrete") {
         feature.norm <- (feature - min(feature))/(max(feature) - min(feature))
         perf.auc <- roc.area(target, feature.norm)
         AUC[i] <- ifelse(perf.auc$A < 0.5, -(1 - perf.auc$A), perf.auc$A)
         AUC[i] <- signif(AUC[i], digits=4)
         p.val <- perf.auc$p.value
         p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
         AUC.pval[i] <- signif(p.val, digits=4)
         temp <- split(feature, target)
	 x <- temp$'1'
         y <- temp$'0'
         t.stat[i] <- signif(t.test(x=x, y=y)$statistic, digits=4)
         p.val <- t.test(x=x, y=y)$p.value
         p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
         t.pval[i] <- signif(p.val, digits=4)
     }
     annot[i] <- paste(MI[i], "     ", AUC[i], " (", AUC.pval[i], ")    ", t.stat[i], " (", t.pval[i], ")", sep="")
   }

   mycol <- vector(length=512, mode = "numeric")
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   cex.axis = 1
   ncolors <- length(mycol)

   nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(4, 10), FALSE)
   par(mar = c(1, 15, 5, 15))
   max.v <- max(max(target), -min(target))
   V1 <- target
   image(1:length(target), 1:1, as.matrix(V1), zlim = c(0, 1), col=c("yellow", "purple"), axes=FALSE, main="", sub = "", xlab= "", ylab="")
   axis(2, at=1:1, labels=paste(phenotype, target.class), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
   axis(4, at=1:1, labels=paste("NMI     AUC (p-val)     t-test (p-val)", sep=""), adj= 0.5, tick=FALSE, 
       las = 1, cex.axis=0.80, font.axis=1, line=-1) 
   par(mar = c(5, 15, 1, 15))
   V1 <- m
   for (i in 1:dim(V1)[1]) V1[i,] <- (V1[i,] - mean(V1[i,]))/sd(V1[i,])
   max.v <- max(max(V1), -min(V1))
   V1 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
   V1 <- apply(V1, MARGIN=2, FUN=rev)
   image(1:dim(V1)[2], 1:dim(V1)[1], t(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
   axis(2, at=1:dim(V1)[1], labels=row.names(V1), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
   axis(4, at=1:dim(V1)[1], labels=rev(annot), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
   axis(1, at=1:dim(V1)[2], labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=0.60, font.axis=1, line=-1)

   dev.off()

   annot2 <- data.frame(cbind(MI, AUC, AUC.pval, t.stat, t.pval))
   row.names(annot2) <- row.names(m)
   write(paste(c("gene set ", noquote(colnames(annot2))), collapse="\t"), file = output.txt, append = F, ncolumns = length(colnames(annot2)))
   write.table(annot2, file=output.txt, append=T, quote=F, sep="\t", eol="\n", col.names=F, row.names=T)

}

OPAM.Evaluate.Results.2 <- function(
   input.ds,
   input.cls,
   phenotype = NULL,
   target.class = NULL,
   target.type = "discrete",
   sort.results = T,
   display.top.n = 20,
   output.txt,
   output.pdf) {

   pdf(file=output.pdf, height=8.5, width=11)
   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   N <- dim(m)[1]
   model.names <- dataset$row.names
   model.descs <- dataset$descs
   Ns <- length(as.matrix(m)[1,])
   for (i in 1:N) {
      if (sd(as.matrix(m)[i,]) == 0) {
         val <- as.matrix(m)[i, 1]
	 m[i,] <- as.matrix(m)[i,] + runif(n=Ns, min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
      }
   }

   dim(m)
   sample.names <- dataset$names
   CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
   cls.labels <- CLS$class.v
   cls.phen <- CLS$phen
   cls.list <- CLS$class.list 

   library(verification)

   if (is.null(phenotype)) {
      phen.loc <- 1
   } else {
      phen.loc <- match(phenotype, CLS$phen.names)
   }
   if (is.vector(CLS$class.list)) {
       target.vec <- CLS$class.list
   } else {
       target.vec <- CLS$class.list[phen.loc,]
   }
   if (target.type == "continuous") {
      target <- target.vec
   } else if (target.type == "discrete") {
      target <- ifelse(target.vec == target.class, 1, 0)    
   }

   ind <- order(target)
   target <- target[ind]
   target.vec <- target.vec[ind]
   m <- as.matrix(m)[, ind]
   sample.names <- sample.names[ind]
   class.v <- CLS$class.v
   if (is.vector(class.v)) {
       class.v <- class.v[ind]
   } else {
       class.v <- class.v[, ind]
   }
   annot <- MI <- AUC <- AUC.pval <- t.stat <- t.pval <- vector(length=N, mode="numeric")

   NMI.ref <- mutual.inf.P(x = target, y = target, n.grid=100)$NMI

   for (i in 1:N) {
      if (N == 1) {
         feature <- m
      } else {
         feature <- m[i,]
      }
      MI[i] <- signif(mutual.inf.P(target, feature, n.grid=100)$NMI/NMI.ref, 4)
      if (target.type == "continuous") {
         AUC[i] <- AUC.pval[i] <- t.stat[i] <- t.pval[i] <- "-"
      } else if (target.type == "discrete") {
         feature.norm <- (feature - min(feature))/(max(feature) - min(feature))
         perf.auc <- roc.area(target, feature.norm)
         AUC[i] <- ifelse(perf.auc$A < 0.5, -(1 - perf.auc$A), perf.auc$A)
         AUC[i] <- signif(AUC[i], digits=4)
         p.val <- perf.auc$p.value
         p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
         AUC.pval[i] <- signif(p.val, digits=4)
         temp <- split(feature, target)
	 x <- temp$'1'
         y <- temp$'0'
         t.stat[i] <- signif(t.test(x=x, y=y)$statistic, digits=4)
         p.val <- t.test(x=x, y=y)$p.value
         p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
         t.pval[i] <- signif(p.val, digits=4)
     }
     annot[i] <- paste(MI[i], "     ", AUC[i], " (", AUC.pval[i], ")    ", t.stat[i], " (", t.pval[i], ")", sep="")
   }

  if ((N > 1) & (sort.results == T)) {
     MI.order <- order(MI, decreasing=T)
     MI <- MI[MI.order]
     AUC <- AUC[MI.order]
     AUC.pval <- AUC.pval[MI.order]
     t.stat <- t.stat[MI.order]
     t.pval <= t.pval[MI.order]
     m <- as.matrix(m)[MI.order,]
     annot <- annot[MI.order]
   }

   mycol <- vector(length=512, mode = "numeric")
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   cex.axis = 1
   ncolors <- length(mycol)

   nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(4, 10), FALSE)
   par(mar = c(1, 15, 5, 15))
   max.v <- max(max(target), -min(target))
   V1 <- target
   image(1:length(target), 1:1, as.matrix(V1), zlim = c(0, 1), col=c("yellow", "purple"), axes=FALSE, main="", sub = "", xlab= "", ylab="")
   axis(2, at=1:1, labels=paste(phenotype, target.class), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
   axis(4, at=1:1, labels=paste("NMI     AUC (p-val)     t-test (p-val)", sep=""), adj= 0.5, tick=FALSE, 
       las = 1, cex.axis=0.80, font.axis=1, line=-1) 
   par(mar = c(5, 15, 1, 15))

   if (display.top.n > N) display.top.n <- N

   if (N == 1) {
      V1 <- m
      V1 <- (V1 - mean(V1))/sd(V1)
   } else {
      V1 <- m[1:display.top.n, ]
      for (i in 1:display.top.n) V1[i,] <- (V1[i,] - mean(V1[i,]))/sd(V1[i,])
   }

   max.v <- max(max(V1), -min(V1))
   V1 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))

   if (N > 1) {
      V1 <- apply(V1, MARGIN=2, FUN=rev)
      image(1:Ns, 1:display.top.n, t(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      axis(2, at=1:display.top.n, labels=row.names(V1), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(4, at=1:display.top.n, labels=rev(annot[1:display.top.n]), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(1, at=1:Ns, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=0.60, font.axis=1, line=-1)
   } else {
      image(1:Ns, 1:1, as.matrix(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      axis(2, at=1:1, labels=model.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(4, at=1:1, labels=annot[1], adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(1, at=1:Ns, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=0.60, font.axis=1, line=-1)
   }
   dev.off()

   annot2 <- data.frame(cbind(MI, AUC, AUC.pval, t.stat, t.pval))
   row.names(annot2) <- row.names(m)
   write(paste(c("gene set ", noquote(colnames(annot2))), collapse="\t"), file = output.txt, append = F, 
          ncolumns = length(colnames(annot2)))
   write.table(annot2, file=output.txt, append=T, quote=F, sep="\t", eol="\n", col.names=F, row.names=T)

}



OPAM.Evaluate.Results.2.1 <- function(
   input.ds,
   input.cls,
   phenotype = NULL,
   target.class = NULL,
   target.type = "discrete",
   sort.results = T,
   display.top.n = 20,
   output.txt,
   output.tiff) {

   tiff(file=output.tiff, width = 1200, height = 1000, units = "px", pointsize = 17)
   
   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   N <- dim(m)[1]
   model.names <- dataset$row.names
   model.descs <- dataset$descs
   Ns <- length(as.matrix(m)[1,])
   for (i in 1:N) {
      if (sd(as.matrix(m)[i,]) == 0) {
         val <- as.matrix(m)[i, 1]
	 m[i,] <- as.matrix(m)[i,] + runif(n=Ns, min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
      }
   }

   dim(m)
   sample.names <- dataset$names
   CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
   cls.labels <- CLS$class.v
   cls.phen <- CLS$phen
   cls.list <- CLS$class.list 

   library(verification)

   if (is.null(phenotype)) {
      phen.loc <- 1
   } else {
      phen.loc <- match(phenotype, CLS$phen.names)
   }
   if (is.vector(CLS$class.list)) {
       target.vec <- CLS$class.list
   } else {
       target.vec <- CLS$class.list[phen.loc,]
   }
   if (target.type == "continuous") {
      target <- target.vec
   } else if (target.type == "discrete") {
      target <- ifelse(target.vec == target.class, 1, 0)    
   }

   ind <- order(target)
   target <- target[ind]
   target.vec <- target.vec[ind]
   m <- as.matrix(m)[, ind]
   sample.names <- sample.names[ind]
   class.v <- CLS$class.v
   if (is.vector(class.v)) {
       class.v <- class.v[ind]
   } else {
       class.v <- class.v[, ind]
   }
   annot <- MI <- AUC <- AUC.pval <- t.stat <- t.pval <- vector(length=N, mode="numeric")

   NMI.ref <- mutual.inf.P(x = target, y = target, n.grid=100)$NMI

   for (i in 1:N) {
      if (N == 1) {
         feature <- m
      } else {
         feature <- m[i,]
      }
      MI[i] <- signif(mutual.inf.P(target, feature, n.grid=100)$NMI/NMI.ref, 2)
      if (target.type == "continuous") {
         AUC[i] <- AUC.pval[i] <- t.stat[i] <- t.pval[i] <- "-"
      } else if (target.type == "discrete") {
         feature.norm <- (feature - min(feature))/(max(feature) - min(feature))
         perf.auc <- roc.area(target, feature.norm)
         AUC[i] <- ifelse(perf.auc$A < 0.5, -(1 - perf.auc$A), perf.auc$A)
         AUC[i] <- signif(AUC[i], digits=2)
         p.val <- perf.auc$p.value
         p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
         AUC.pval[i] <- signif(p.val, digits=2)
         temp <- split(feature, target)
	 x <- temp$'1'
         y <- temp$'0'
         t.stat[i] <- signif(t.test(x=x, y=y)$statistic, digits=2)
         p.val <- t.test(x=x, y=y)$p.value
         p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
         t.pval[i] <- signif(p.val, digits=2)
     }
     annot[i] <- paste(MI[i], "     ", AUC[i], " (", AUC.pval[i], ")    ", t.stat[i], " (", t.pval[i], ")", sep="")
   }

  if ((N > 1) & (sort.results == T)) {
     MI.order <- order(MI, decreasing=T)
     MI <- MI[MI.order]
     AUC <- AUC[MI.order]
     AUC.pval <- AUC.pval[MI.order]
     t.stat <- t.stat[MI.order]
     t.pval <= t.pval[MI.order]
     m <- as.matrix(m)[MI.order,]
     annot <- annot[MI.order]
   }

   mycol <- vector(length=512, mode = "numeric")
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   cex.axis = 1
   ncolors <- length(mycol)

   nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(2, 10), FALSE)
   par(mar = c(1, 19, 4, 11))
   max.v <- max(max(target), -min(target))
   V1 <- target
   image(1:length(target), 1:1, as.matrix(V1), zlim = c(0, 1), col=c("lightgray", "black"), axes=FALSE, main="", sub = "", xlab= "", ylab="")
   axis(2, at=1:1, labels=paste(phenotype, target.class), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
   axis(4, at=1:1, labels=paste("NMI     AUC (p-val)     t-test (p-val)", sep=""), adj= 0.5, tick=FALSE, 
       las = 1, cex.axis=0.80, font.axis=1, line=-1) 
   par(mar = c(4, 19, 1, 11))

   if (display.top.n > N) display.top.n <- N

   if (N == 1) {
      V1 <- m
      V1 <- (V1 - mean(V1))/sd(V1)
   } else {
      V1 <- m[1:display.top.n, ]
      for (i in 1:display.top.n) V1[i,] <- (V1[i,] - mean(V1[i,]))/sd(V1[i,])
   }

   max.v <- max(max(V1), -min(V1))
   V1 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))

   if (N > 1) {
      V1 <- apply(V1, MARGIN=2, FUN=rev)
      image(1:Ns, 1:display.top.n, t(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      axis(2, at=1:display.top.n, labels=row.names(V1), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(4, at=1:display.top.n, labels=rev(annot[1:display.top.n]), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(1, at=1:Ns, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=0.60, font.axis=1, line=-1)
   } else {
      image(1:Ns, 1:1, as.matrix(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      axis(2, at=1:1, labels=model.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(4, at=1:1, labels=annot[1], adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(1, at=1:Ns, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=0.60, font.axis=1, line=-1)
   }
   dev.off()

   annot2 <- data.frame(cbind(MI, AUC, AUC.pval, t.stat, t.pval))
   row.names(annot2) <- row.names(m)
   write(paste(c("gene set ", noquote(colnames(annot2))), collapse="\t"), file = output.txt, append = F, 
          ncolumns = length(colnames(annot2)))
   write.table(annot2, file=output.txt, append=T, quote=F, sep="\t", eol="\n", col.names=F, row.names=T)

}



OPAM.Evaluate.Results.4 <- function(
   input.ds,
   input.cls,
   phenotype = NULL,
   target.class = NULL,
   target.type = "discrete",
   sort.phenotype = T,
   sort.results = T,
   display.top.n = 20,
   output.txt,
   output.pdf,
   cex.axis=0.7) {

   pdf(file=output.pdf, height=8.5, width=11)
   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   N <- dim(m)[1]
   model.names <- dataset$row.names
   model.descs <- dataset$descs
   Ns <- length(as.matrix(m)[1,])
   for (i in 1:N) {
      if (sd(as.matrix(m)[i,]) == 0) {
         val <- as.matrix(m)[i, 1]
	 m[i,] <- as.matrix(m)[i,] + runif(n=Ns, min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
      }
   }

   dim(m)
   sample.names <- dataset$names
   CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
   cls.labels <- CLS$class.v
   cls.phen <- CLS$phen
   cls.list <- CLS$class.list 

   library(verification)

   if (is.null(phenotype)) {
      phen.loc <- 1
   } else {
      phen.loc <- match(phenotype, CLS$phen.names)
   }
   if (is.vector(CLS$class.list)) {
       if (target.type == "discrete") {
          target <- ifelse(CLS$class.list == target.class, 1, 0)    
       } else {
          target <- as.numeric(CLS$class.v) + runif(length(CLS$class.v), min=0, max=10*.Machine$double.eps)
       }
   } else {
       if (target.type == "discrete") {
          target <- ifelse(CLS$class.list[phen.loc,] == target.class, 1, 0)    
       } else {
          target <- as.numeric(CLS$class.v[phen.loc,]) + runif(length(CLS$class.v[phen.loc,]), min=0, max=10*.Machine$double.eps)

       }
   }

   print(paste("target:", target))

   class.v <- CLS$class.v
   if (sort.phenotype == T) {
      ind <- order(target)
      target <- target[ind]
      m <- as.matrix(m)[, ind]
      sample.names <- sample.names[ind]
      if (is.vector(class.v)) {
          class.v <- class.v[ind]
      } else {
          class.v <- class.v[, ind]
      }
   }

   annot <- MI <- AUC <- AUC.pval <- t.stat <- t.pval <- vector(length=N, mode="numeric")

   NMI.ref <- mutual.inf.P(x = target, y = target, n.grid=100)$NMI

   for (i in 1:N) {
      if (N == 1) {
         feature <- m
      } else {
         feature <- m[i,]
      }
      
      MI[i] <- signif(mutual.inf.P(target, feature, n.grid=100)$NMI/NMI.ref, 4)
      if (target.type == "continuous") {
         AUC[i] <- AUC.pval[i] <- t.stat[i] <- t.pval[i] <- "-"
      } else if (target.type == "discrete") {
         feature.norm <- (feature - min(feature))/(max(feature) - min(feature))
         perf.auc <- roc.area(target, feature.norm)
         AUC[i] <- ifelse(perf.auc$A < 0.5, -(1 - perf.auc$A), perf.auc$A)
         AUC[i] <- signif(AUC[i], digits=4)
         p.val <- perf.auc$p.value
         p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
         AUC.pval[i] <- signif(p.val, digits=4)
         temp <- split(feature, target)
	 x <- temp$'1'
         y <- temp$'0'
         t.stat[i] <- signif(t.test(x=x, y=y)$statistic, digits=4)
         p.val <- t.test(x=x, y=y)$p.value
         p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
         t.pval[i] <- signif(p.val, digits=4)
     }
     annot[i] <- paste(MI[i], "     ", AUC[i], " (", AUC.pval[i], ")    ", t.stat[i], " (", t.pval[i], ")", sep="")
   }

  if ((N > 1) & (sort.results == T)) {
     MI.order <- order(MI, decreasing=T)
     MI <- MI[MI.order]
     AUC <- AUC[MI.order]
     AUC.pval <- AUC.pval[MI.order]
     t.stat <- t.stat[MI.order]
     t.pval <= t.pval[MI.order]
     m <- as.matrix(m)[MI.order,]
     annot <- annot[MI.order]
   }

   mycol <- vector(length=512, mode = "numeric")
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   ncolors <- length(mycol)

   nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(4, 10), FALSE)
   par(mar = c(1, 18, 3, 11))
   max.v <- max(max(target), -min(target))
   V1 <- target
   image(1:length(target), 1:1, as.matrix(V1), zlim = c(0, 1), col=c("yellow", "purple"), axes=FALSE, main="", sub = "", xlab= "", ylab="")
   axis(2, at=1:1, labels=paste(phenotype, target.class), adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.axis, font.axis=1, line=-1)
   axis(4, at=1:1, labels=paste("NMI     AUC (p-val)     t-test (p-val)", sep=""), adj= 0.5, tick=FALSE, 
       las = 1, cex.axis=cex.axis, font.axis=1, line=-1) 
   par(mar = c(5, 18, 1, 11))

   if (display.top.n > N) display.top.n <- N

   if (N == 1) {
      V1 <- m
      V1 <- (V1 - mean(V1))/sd(V1)
   } else {
      V1 <- m[1:display.top.n, ]
      for (i in 1:display.top.n) V1[i,] <- (V1[i,] - mean(V1[i,]))/sd(V1[i,])
   }

   max.v <- max(max(V1), -min(V1))
   V1 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))

   if (N > 1) {
      V1 <- apply(V1, MARGIN=2, FUN=rev)
      image(1:Ns, 1:display.top.n, t(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      axis(2, at=1:display.top.n, labels=row.names(V1), adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.axis, font.axis=1, line=-1)
      axis(4, at=1:display.top.n, labels=rev(annot[1:display.top.n]), adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.axis, font.axis=1, line=-1)
      axis(1, at=1:Ns, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=cex.axis, font.axis=1, line=-1)
   } else {
      image(1:Ns, 1:1, as.matrix(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      axis(2, at=1:1, labels=model.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.axis, font.axis=1, line=-1)
      axis(4, at=1:1, labels=annot[1], adj= 0.5, tick=FALSE, las = 1, cex.axis=cex.axis, font.axis=1, line=-1)
      axis(1, at=1:Ns, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=cex.axis, font.axis=1, line=-1)
   }
   dev.off()

   annot2 <- data.frame(cbind(MI, AUC, AUC.pval, t.stat, t.pval))
   row.names(annot2) <- row.names(m)
   write(paste(c("gene set ", noquote(colnames(annot2))), collapse="\t"), file = output.txt, append = F, 
          ncolumns = length(colnames(annot2)))
   write.table(annot2, file=output.txt, append=T, quote=F, sep="\t", eol="\n", col.names=F, row.names=T)

}



OPAM.Evaluate.Results.Yan <- function(
# this version has to option of using the Kernel method, its adjusted variant and the knn approach from the parmigene package
                                      
   input.ds,
   input.cls,
   phenotype = NULL,
   target.class = NULL,
   target.type = "discrete",
   MI.method = "kernel_I", # kernel_I, kernel_Adj or knn                                    
   sort.results = T,
   display.top.n = 20,
   output.txt,
   output.pdf) {

   library(parmigene)

  
   pdf(file=output.pdf, height=8.5, width=11)
   dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
   m <- data.matrix(dataset$ds)
   N <- dim(m)[1]
   model.names <- dataset$row.names
   model.descs <- dataset$descs
   Ns <- length(as.matrix(m)[1,])
   for (i in 1:N) {
      if (sd(as.matrix(m)[i,]) == 0) {
         val <- as.matrix(m)[i, 1]
	 m[i,] <- as.matrix(m)[i,] + runif(n=Ns, min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
      }
   }

   dim(m)
   sample.names <- dataset$names
   CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
   cls.labels <- CLS$class.v
   cls.phen <- CLS$phen
   cls.list <- CLS$class.list 

   library(verification)

   if (is.null(phenotype)) {
      phen.loc <- 1
   } else {
      phen.loc <- match(phenotype, CLS$phen.names)
   }
   if (is.vector(CLS$class.list)) {
       target.vec <- CLS$class.list
   } else {
       target.vec <- CLS$class.list[phen.loc,]
   }
   if (target.type == "continuous") {
      target <- target.vec
   } else if (target.type == "discrete") {
      target <- ifelse(target.vec == target.class, 1, 0)    
   }

   ind <- order(target)
   target <- target[ind]
   target.vec <- target.vec[ind]
   m <- as.matrix(m)[, ind]
   sample.names <- sample.names[ind]
   class.v <- CLS$class.v
   if (is.vector(class.v)) {
       class.v <- class.v[ind]
   } else {
       class.v <- class.v[, ind]
   }
   annot <- MI <- AUC <- AUC.pval <- t.stat <- t.pval <- vector(length=N, mode="numeric")

   if (MI.method == "kernel_I") {
      NMI.ref <- NMI(x=target, y=target, n.grid=50, make.plot=F)$NMI 
    } else if (MI.method == "kernel_Adj") {
       rho <- cor(target, target)
       adj <- log(1/(abs(rho) + 0.25)) + 0.75
       delta.param <- 1
       delta <- c(delta.param * adj * bcv(target), delta.param * adj * bcv(target))
       NMI.ref <- NMI(x=target, y=target, n.grid=50, delta = delta, make.plot=F)$NMI  
    } else { # knn method
       NMI.ref <- cor(target, target) * knnmi(target, target, k=3, noise=1e-12)/NMI(x=target, y=target, n.grid=50, make.plot=F)$HXY
    }
      
   for (i in 1:N) {
      if (N == 1) {
         feature <- m
      } else {
         feature <- m[i,]
      }

   if (MI.method == "kernel_I") {
       NMI.val <- NMI(x=target, y=feature, n.grid=50, make.plot=F)$NMI  # MI according to our kernel method
    } else if (MI.method == "kernel_Adj") {
       rho <- cor(target, feature)
       adj <- log(1/(abs(rho) + 0.25)) + 0.75
       delta.param <- 1
       delta <- c(delta.param * adj * bcv(target), delta.param * adj * bcv(feature))
       NMI.val <- NMI(x=target, y=feature, n.grid=50, delta = delta, make.plot=F)$NMI  # MI according to our kernel method
    } else { # knn method
       NMI.val<- cor(target, feature) * knnmi(target, feature, k=3, noise=1e-12)/NMI(x=target, y=feature, n.grid=50, make.plot=F)$HXY
    }
      
      MI[i] <- signif(NMI.val/NMI.ref, 4)
      if (target.type == "continuous") {
         AUC[i] <- AUC.pval[i] <- t.stat[i] <- t.pval[i] <- "-"
      } else if (target.type == "discrete") {
         feature.norm <- (feature - min(feature))/(max(feature) - min(feature))
         perf.auc <- roc.area(target, feature.norm)
         AUC[i] <- ifelse(perf.auc$A < 0.5, -(1 - perf.auc$A), perf.auc$A)
         AUC[i] <- signif(AUC[i], digits=4)
         p.val <- perf.auc$p.value
         p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
         AUC.pval[i] <- signif(p.val, digits=4)
         temp <- split(feature, target)
	 x <- temp$'1'
         y <- temp$'0'
         t.stat[i] <- signif(t.test(x=x, y=y)$statistic, digits=4)
         p.val <- t.test(x=x, y=y)$p.value
         p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
         t.pval[i] <- signif(p.val, digits=4)
     }
     annot[i] <- paste(MI[i], "     ", AUC[i], " (", AUC.pval[i], ")    ", t.stat[i], " (", t.pval[i], ")", sep="")
   }

  if ((N > 1) & (sort.results == T)) {
     MI.order <- order(MI, decreasing=T)
     MI <- MI[MI.order]
     AUC <- AUC[MI.order]
     AUC.pval <- AUC.pval[MI.order]
     t.stat <- t.stat[MI.order]
     t.pval <= t.pval[MI.order]
     m <- as.matrix(m)[MI.order,]
     annot <- annot[MI.order]
   }

   mycol <- vector(length=512, mode = "numeric")
   for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
   for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
   mycol <- rev(mycol)
   cex.axis = 1
   ncolors <- length(mycol)

   nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(4, 10), FALSE)
   par(mar = c(1, 15, 5, 15))
   max.v <- max(max(target), -min(target))
   V1 <- target
   image(1:length(target), 1:1, as.matrix(V1), zlim = c(0, 1), col=c("yellow", "purple"), axes=FALSE, main="", sub = "", xlab= "", ylab="")
   axis(2, at=1:1, labels=paste(phenotype, target.class), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
   axis(4, at=1:1, labels=paste("NMI     AUC (p-val)     t-test (p-val)", sep=""), adj= 0.5, tick=FALSE, 
       las = 1, cex.axis=0.80, font.axis=1, line=-1) 
   par(mar = c(5, 15, 1, 15))

   if (display.top.n > N) display.top.n <- N

   if (N == 1) {
      V1 <- m
      V1 <- (V1 - mean(V1))/sd(V1)
   } else {
      V1 <- m[1:display.top.n, ]
      for (i in 1:display.top.n) V1[i,] <- (V1[i,] - mean(V1[i,]))/sd(V1[i,])
   }

   max.v <- max(max(V1), -min(V1))
   V1 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))

   if (N > 1) {
      V1 <- apply(V1, MARGIN=2, FUN=rev)
      image(1:Ns, 1:display.top.n, t(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      axis(2, at=1:display.top.n, labels=row.names(V1), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(4, at=1:display.top.n, labels=rev(annot[1:display.top.n]), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(1, at=1:Ns, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=0.60, font.axis=1, line=-1)
   } else {
      image(1:Ns, 1:1, as.matrix(V1), zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
      axis(2, at=1:1, labels=model.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(4, at=1:1, labels=annot[1], adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
      axis(1, at=1:Ns, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=0.60, font.axis=1, line=-1)
   }
   dev.off()

   annot2 <- data.frame(cbind(MI, AUC, AUC.pval, t.stat, t.pval))
   row.names(annot2) <- row.names(m)
   write(paste(c("gene set ", noquote(colnames(annot2))), collapse="\t"), file = output.txt, append = F, 
          ncolumns = length(colnames(annot2)))
   write.table(annot2, file=output.txt, append=T, quote=F, sep="\t", eol="\n", col.names=F, row.names=T)

}

NMI <- function(x, y, n.grid=100, make.plot=F, delta = c(bcv(x), bcv(y))) {
    # for definitions of mutual information and the universal metric (NMI) see the 
    # definition of "Mutual Information" in wikipedia and Thomas and Cover's book

   kde2d.xy <- kde2d(x, y, n = n.grid, h = delta)
   X <- kde2d.xy$x
   Y <- kde2d.xy$y  
   PXY <- kde2d.xy$z + .Machine$double.eps
   PXY <- PXY/sum(PXY)
   PX <- rowSums(PXY)
   PX <- PX/sum(PX)
   HX <- -sum(PX * log2(PX))
   PY <- colSums(PXY)
   PY <- PY/sum(PY)
   HY <- -sum(PY * log2(PY))
   HXY <- - sum(PXY * log2(PXY))
   NMI <- sign(cor(x, y)) * ((HX + HY)/HXY - 1)  # use pearson correlation the get the sign (directionality)
   PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
   PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
   MI <- sum(PXY * log2(PXY/(PX*PY)))

   if (make.plot != F) {
      mycol <- vector(length=512, mode = "numeric")
      for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
      for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
      mycol <- rev(mycol)
      quartz(width=12, height=8)
      nf <- layout(matrix(c(1,2,3,4), 2, 2, byrow=T), c(1,1), c(1,1), TRUE)
      plot(X, PX[,1], type="l", main="X")
      plot(Y, PY[1,], type="l", main="Y")
      MIXY <- PXY * log2(PXY/(PX*PY))
      image(PXY, main = paste("P(x, y)", " pearson cor=", signif(cor(x, y), 3)), col=mycol)
      sub <- ifelse(make.plot != T, make.plot, " ")
      image(MIXY, main=paste("MI:", signif(MI, 3)), col=mycol, sub=sub)
   }
   return(list(NMI=NMI, MI=MI, HXY=HXY, HX=HX, HY=HY))
}
