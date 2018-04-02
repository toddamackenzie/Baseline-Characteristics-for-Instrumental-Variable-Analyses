
##############################################################################
# Baseline Statistics by Compliance Strata
# Comparisons of Baseline Statistics Between Compliers, Takers and Non-takers
# R: binary instrument
# X: binary treatment
# Z: covariates
##############################################################################

Proportions.CTN <- function(R, X) {
	p.T <- sum(X==1 & R == 0) / sum(R == 0)
	p.N <- sum(X==0 & R == 1) / sum(R == 1)
	p.C <- 1 - p.T - p.N
	c("C"=p.C, "T"=p.T, "N"=p.N)
}

Means.CTN <- function(R, X, Z) {
	p <- Proportions.CTN(R, X)
	mean.All <- mean(Z)
	mean.T <- mean(Z[R==0 & X==1])
	mean.N <- mean(Z[R==1 & X==0])
	mean.C <- (mean.All - p["T"]*mean.T - p["N"]*mean.N) / p["C"]
	o <- c(mean.C, mean.T, mean.N, mean.All)
	names(o) <- c("C", "T", "N", "All") 
	o
}

Variance.CTN <- function(R, X, Z) {
	p <- Proportions.CTN(R, X)
	mns <- Means.CTN(R, X, Z)
	var.All <- var(Z)
	var.T <- var(Z[R==0 & X==1])
	var.N <- var(Z[R==1 & X==0])
	var.mns <- sum(p * (mns[1:3] - mns["All"])^2) 
	var.C <- (var.All - var.mns - p["N"]*var.N - p["T"]*var.T) / p["C"]
	o <- c(var.C, var.T, var.N, var.All)
	names(o) <- c("C", "T", "N", "All") 
	o	
}

Covariance.of.Means.CTN <- function(R, X, Z) {
	p <- Proportions.CTN(R, X)
	mns <- Means.CTN(R, X, Z)
	vars <- Variance.CTN(R, X, Z)
	cov.TT <- vars["T"]/sum(!R & X)
	cov.NN <- vars["N"]/sum(R & !X)
	cov.TN <- 0
	# Derivation of the complier var and covar of the means
	n <- length(R)
	n.RX <- c("00"=sum(!R & !X),"01"=sum(!R & X),"10"=sum(R & !X),"11"=sum(R& X))
	odds.R01 <- sum(R==0)/sum(R==1)
	var.RX <-c("00"=var(Z[!R&!X]), "01"=var(Z[!R&X]), "10"=var(Z[R&!X]), "11"=var(Z[R&X]))
	cov.CC <- sum( (n.RX*var.RX) * c(1, (1/odds.R01)^2, (odds.R01)^2, 1) / (n * p["C"])^2 )
	cov.CT <- - (1/odds.R01) * vars["T"] / n
	cov.CN <- - (odds.R01) * vars["N"] / n
	o <- matrix(nrow=3, ncol=3, c(cov.CC, cov.CT, cov.CN,  cov.CT, cov.TT, cov.TN, cov.CN, cov.TN, cov.NN))
	dimnames(o)[[1]] <- dimnames(o)[[2]] <- c("C", "T", "N")
	o
}

Differences.CTN <- function(R, X, Z) {
	mns <- Means.CTN(R, X, Z)
	Cov.Mns <- Covariance.of.Means.CTN(R, X, Z)
	mn.diffs <- mns["C"] - mns[c("T","N")]
	L <- matrix(nrow=3, ncol=2, c(1,-1,0,   1,0,-1))
	Cov.diffs <- t(L) %*% Cov.Mns %*% L
	Chi.2 <- mn.diffs %*% solve(Cov.diffs) %*% mn.diffs
	p <- 1 - pchisq(Chi.2, df=2)
	se.mn.diffs <- diag(Cov.diffs)^0.5
	mn.diffs <- c(mn.diffs, mns["T"] - mns["N"] )
	var.T.vs.N <- sum(diag(Covariance.of.Means.CTN(R, X, Z))[c(2,3)])
	se.mn.diffs <- c(se.mn.diffs, sqrt(var.T.vs.N))
	Chi.1 <- (mn.diffs/se.mn.diffs)^2
	p.mn.diffs <- 1 - pchisq(Chi.1, 1)
	nms <- c("Complier vs Taker", "Complier vs Non-Taker", "Taker vs Non-Taker")
	names(mn.diffs) <- names(se.mn.diffs) <- names(p.mn.diffs) <- nms
	list(TS=Chi.2, p.overall=p, mn.diffs=mn.diffs, se.mn.diffs=se.mn.diffs, chi.mn.diffs=Chi.1, p.mn.diffs=p.mn.diffs,Cov.diffs)
}

