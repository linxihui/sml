#' @title Backward feature elimination based on Wald type test for Cox model
#' @param formula A model formula to start backward elimination
#' @param data A data frame
#' @param alpha Significant level. Default to 0.05
#' @param method Method to resolve tied event times pass to `coxph`. Default to 'breslow' instead of 'efron`
#' @param trace A integer for different levels of tracing. 0 = no trace, 1 = print dropped variable, 2 = print all considered models during the procedure
#' @param ... Other parameters passed `coxph`
#' @return A list with two entries: \cr
#'     $ fit: a `coxph` object of selected model
#'     $ wald.test: wald test p values for each feature and possible interaction.
#' @details Notice that model hierachical structure is maintained, i.e., iff interactions involved in the model, its main effects will also be kept in the model.
#' @export
bwe.cox <- function(formula, data, alpha = 0.05, trace = 1, method = 'breslow', ...) {
	formula <- stats::formula(terms(formula, data = data));
	fm <- formula;
	wald.ps <- 1;
	rm.v <- '';
	while(length(attr(terms(fm), 'term.labels')) > 0 && !is.null(rm.v)) {
		fit <- coxph(fm, data, method = method);
		fit$call$formula <- fm;
		switch(as.character(trace), 
			'2' = {print(fit);}, 
			'1' = {print(fm); cat(rep('-', 40), '\n', sep='')}
			); 
		wald.ps <- sapply(
			X = fit$assign, 
			FUN = function(ind) {
				chisq <- t(fit$coef[ind]) %*% solve(fit$var[ind, ind, drop=FALSE]) %*% fit$coef[ind];
				pchisq(chisq, df = length(ind), lower.tail = FALSE)
				}
			);

		# keep the hierachical structure
		wald.ps <- sort(wald.ps, decreasing = TRUE);
		isignif.v <- names(wald.ps)[wald.ps > alpha];
		fm.str <- colnames(attr(terms(fm), 'factor'));
		rm.v <- NULL;
		if (length(isignif.v) > 0) {
			for (v in isignif.v) {
				v.pos <- grep(paste0('(^|:)', v, '(:|$)'), fm.str);
				if (length(v.pos) < 2) {rm.v <- v; if(trace > 0) cat(rm.v, ' dropped\n'); break};
				}
			}
	
		if (!is.null(rm.v)) {
			fm <- update(fm, formula(paste0('.~.-', rm.v)));
			}
		}
	if (trace > 1) {cat('\n', rep('-', 40), '\n', 'Finaly Model:\n', rep('-', 40), '\n', sep=''); print(fit);}
	return(list(fit = fit, wald.test = wald.ps));
	}
