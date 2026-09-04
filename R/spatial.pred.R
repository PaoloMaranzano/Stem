#' @keywords internal
#' @noRd

`spatial.pred` <-
function(mu1,mu2,Sigma11,Sigma12,Sigma22,X1){
	### Nothing observed at this time point: there is nothing to condition on,
	### so the predictor is the unconditional mean and its variance the
	### unconditional one.
	if (length(X1) == 0L) {
		return(list(pred = mu2, se.pred = sqrt(diag(Sigma22))))
	}
	pred = mu2 + Sigma12 %*% solve(Sigma11) %*% (X1-mu1)
	var.pred = Sigma22 - Sigma12 %*%solve(Sigma11)%*%t(Sigma12)
	se.pred = sqrt(diag(var.pred))
	return(list(pred=pred,se.pred=se.pred))
}

