logit <- function(x, min=0, max=1) {
	p <- (x-min)/(max-min)
	p <- ifelse(p <   .Machine$double.neg.eps,  .Machine$double.neg.eps,p) # fix -Inf
	p <- ifelse(p > 1-.Machine$double.neg.eps,1-.Machine$double.neg.eps,p) # fix +Inf
	log(p/(1-p))
}

inv.logit <- function(x, min=0, max=1) {
	p <- exp(x)/(1+exp(x))
	p <- ifelse(is.na(p) & !is.na(x), 1, p )                # fix problems with +Inf
	p <- ifelse(p <= exp(logit(0))/(1+exp(logit(0))), 0, p) # fix 0 rounding
	p * (max-min) + min
}


pinvr <- function(X, max.sv = min(dim(X)), tol = sqrt(.Machine$double.eps)) {
	#
	# based on suggestions of R. M. Heiberger, T. M. Hesterberg and WNV
	#
	if(length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
		stop("'X' must be a numeric or complex matrix")
	if(!is.matrix(X)) X <- as.matrix(X)
	Xsvd <- svd(X)
	if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
	Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)

	# My extension: START
	max.sv      <- min(ifelse(max.sv < 0, 1, max.sv),min(dim(X)))
	L           <- logical(length(Positive))
	L[1:max.sv] <- TRUE
	Positive    <- Positive & L
	# My extension: END

	if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
	else if(!any(Positive)) array(0, dim(X)[2L:1L])
	else Xsvd$v[, Positive, drop=FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop=FALSE]))
}

missMethy <- function(dat, min = NULL, max = NULL, max.sv = NULL, col.list = NULL) {
	out    <- dat
	NAcols <- colSums(is.na(dat)) > 0; NAcols <- which(NAcols)

	# Convert col.list, if any, from names to numbers
	if(is.character(col.list)) {
		if(is.null(colnames(dat)))
			col.list <- NULL
		else
			col.list <- which(colnames(dat) %in% col.list)
	}

	# If all the colums have a missing value we cannot do anything
	if(length(NAcols) < ncol(dat)) {
		# Columns with all NAs or a single not NA value (to be excluded: not enough information for imputation)
		NAall  <- colSums(is.na(dat))<(nrow(dat)-1); NAall <- which(NAall)
		# List of columns to impute
		NAlist <- intersect(NAcols, NAall)
		# Filter the columns to impute
		if(!is.null(col.list))
			NAlist <- intersect(NAlist,col.list)

		while(length(NAlist) != 0) {
			col_id <- NAlist[1]

			# List of rows for which col_id is NA
			row_id <- which(is.na(dat[,col_id])==TRUE)

			# Colum indexes of NA columns for all the row_id(s)
			if(length(row_id) == 1)
				tmp1 <- which(is.na(dat[row_id,])==TRUE)
			else
				tmp1 <- which(colSums(is.na(dat[row_id,])) == length(row_id))

			# Column indexes: no colum element is NA for the rows not in row_id
			tmp2   <- which(colSums(is.na(dat[-row_id,])) == 0)

			# List of colums in NAlist that are NA only for all the row_id(s)
			NAcols_rowid <- intersect(intersect(tmp1,tmp2),NAlist)

			# Extract submatrices for regression
			A <- dat[-row_id,-NAcols]
			B <- dat[-row_id,NAcols_rowid]
			C <- dat[row_id,-NAcols]

			#print(paste("Actual positions:",length(NAlist)," #SampIDs:",length(row_id),"#CpGIDs:",length(NAcols_rowid),"Dim(A):",nrow(A),ncol(A),"-",format(object.size(A), units = "auto")))

			# Updates or computes min, max range from B distribution
			min    <- min(min(B),min)
			max    <- max(max(B),max)
			# Updates or computes max.sv from A. Negative or zero value not allowed
			max.sv <- max(ifelse(is.null(max.sv),min(dim(A)),max.sv),1)

			# X <- pinvr(A,rank)%*%logit(B,min,max) (X = A^-1*logit(B))
			# P <- inv.logit(C%*%X,min,max)         (P = logit^-1(C*X))
			out[row_id,NAcols_rowid] <- inv.logit(C%*%(pinvr(A,max.sv)%*%logit(B,min,max)),min,max)

			# Update NA column list
			NAlist <- setdiff(NAlist,NAcols_rowid)
		}
	}
	return(out)
}
