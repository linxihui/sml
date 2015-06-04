#' @title Kernel Weighted Kaplan-Meier model for survival
#' @param formula Model formula
#' @param data Data frame
#' @param x Design matrix (NO intercept)
#' @param y Reponse vector of `Surv` object 
#' @param xtest If the 'formula-data' format used, `xtest` is a data frame of the test set.
#' 	If the 'x-y' method called, `xtest` is the design matrix of the test set
#' @param ytest (optional) Survival outcome for testset
#' @param times Times to evaluate survival probabilities. Currently no used. All unique event times in training are used
#' @param type Type of estimator, either 'Kaplan-meier' or 'nelson-aalen'
#' @param k Number of nearest neighbour used. Default to `nrow(x)` which seems the best
#' @param scaled Logical value indicating if to standardize x, y
#' @param kernel String or 'kernel' object (see kernlab::gausspr)
#' @param kpar A list of Kernel parameters or 'automatic' if a radial kernel specified
#' @param ... Further argument passed to internal functions
#' @return A 'kkm' object
#' @examples
#'
#' library(survival);
#' data(pbc, package = 'randomForestSRC');
#' pbc <- na.omit(pbc);
#' i.tr <- sample(nrow(pbc), 100);
#' kkm.pred <- kkm(Surv(days, status) ~., data = pbc[i.tr, ], xtest = pbc[-i.tr, ], kernel = 'laplacedot', kpar = list(sigma = 0.1));
#' # concordance index if using 30th event time
#' survConcordance(y[-i.tr] ~ I(1- kkm.pred$test.predicted.survival[, 30]))$concordance
#'
#' @export
kkm <- function(x, ...) UseMethod('kkm');

#' @rdname kkm
#' @export 
kkm.formula <- function(
	formula, data, 
	xtest = NULL, ytest = NULL,
	...) {
	formula <- stats::formula(terms(formula, data = data));
	mf <- model.frame(formula, data);
	x <- model.matrix(formula, data);
	y <- model.response(mf);
	xint <- match("(Intercept)", colnames(x), nomatch = 0);
	if (xint > 0) x <- x[, -xint, drop = FALSE];
	rm(mf);
	if (is.null(ytest) && !is.null(xtest)) {
		resp.var <- all.names(formula[2]);
		if (length(resp.var) > 1) resp.var <- resp.var[-1];
		if (all(resp.var %in% names(xtest))) ytest <- eval(formula[[2]], envir = xtest);
		}
	if (!is.null(xtest)) {
		xtest <- model.matrix(formula[-2], xtest)[, -xint, drop = FALSE];
		}
	out <- kkm.default(x, y, xtest, ytest, ...);
	out$formula <- formula;
	out$call <- match.call();
	class(out) <- c('kkm.formula', class(out));
	return(out);
	}

#' @rdname kkm
#' @export
kkm.default <- function(
	x, y, xtest = NULL, ytest = NULL, 
	return.train.prediction = FALSE || is.null(xtest),
	scaled = TRUE, k = min(30, nrow(x)), 
	times = y[as.logical(y[, 2]), 1],
	type = c('kaplan-meier', 'nelson-aalen'),
	kernel = 'rbfdot', kpar = 'automatic'
	) {
	stopifnot(suppressMessages(require(survival)));
	stopifnot(suppressMessages(require(kernlab)));
	type <- match.arg(type);
	x.scale <- NULL;
	if (scaled) {
		x <- scale(x);
		x.scale <- attributes(x)[c('scaled:center', 'scaled:scale')];
		attributes(x)[c('scaled:center', 'scaled:scale')] <- NULL;
		}

	if (!is.null(xtest) && !is.null(x.scale)) {
		xtest <- scale(xtest, center = x.scale[['scaled:center']], scale = x.scale[['scaled:scale']]);
		}

	if (return.train.prediction) {
		if.test <- c(rep(FALSE, nrow(x)), rep(TRUE, ifelse(is.null(xtest), 0,  nrow(xtest))));
		xtest <- rbind(x, xtest);
	} else {
		if.test <- rep(TRUE, ifelse(is.null(xtest), 0, nrow(xtest)));
		}

	if(is.matrix(kernel))  {
		kern.weight <- kernel;
	} else {
		if(is.character(kernel)) {
			kernel <- match.arg(
				kernel, 
				c(
					"rbfdot","polydot","tanhdot","vanilladot",
					"laplacedot","besseldot","anovadot","splinedot"
					)
				);
			if(is.character(kpar)) {
				if((kernel %in% c("tanhdot", "vanilladot", "polydot", "besseldot", "anovadot", "splinedot")) && "automatic" == kpar ) {
					# cat(" Setting default kernel parameters ","\n");
					kpar <- list();
					}
				}
			}

		if (!is.function(kernel)) {
			if (!is.list(kpar) && is.character(kpar) && c(class(kernel) %in% c("rbfkernel", "laplacedot") || kernel %in% c("laplacedot", "rbfdot"))) {
				kp <- match.arg(kpar, "automatic");
				if("automatic" == kp) {
					kpar <- list(sigma = mean(sigest(x, scaled = FALSE)[c(1,3)]));
					}
				# cat("Using automatic sigma estimation (sigest) for RBF or laplace kernel","\n");
				}
			}

		if(!is(kernel,"kernel")) {
			if (is(kernel,"function")) {kernel <- deparse(substitute(kernel));}
			kernel <- do.call(kernel, kpar);
			}
		if(!is(kernel,"kernel")) {stop("kernel must inherit from class `kernel'");}

		kern.weight <- kernelMatrix(kernel, x, xtest)@.Data;
		}

	if (k < nrow(x)) {
		knn.mask <- apply(kern.weight, 2, function(z) {
			o <- rep(0, length(z));
			o[order(-z)[1:k]] <- 1;
			return(o);
			})
		kern.weight <- kern.weight*knn.mask;
		}

	time.order <- order(y[,1], -y[,2]);
	kern.weight <- kern.weight[time.order, ];
	y <- y[time.order]; # time sorted from smallest to largest

	is.event <- as.logical(y[, 2]);
	stopifnot(all(times > 0));
	times <- unique(sort(times));
	is.nonduplicated.event <- !duplicated(y[is.event, 1]);
	uniq.event.times <- y[is.event, 1][is.nonduplicated.event];
	is.subset <- all(times %in% uniq.event.times);
	is.exact <- identical(times, uniq.event.times);
	if (is.subset &!is.exact) {
		time.subset <- uniq.event.times %in% times;
		}

	kkm.est <- t(apply(
		X = kern.weight,
		MARGIN = 2,
		FUN = function(w) {
			event.mask <- is.event & w > 0;
			at.risk <- rev(cumsum(rev(w)))[event.mask];
			surv <- switch(type,
				'kaplan-meier' = exp(cumsum(log(1 - w[event.mask] / at.risk))),
				'nelson-aalen' = exp(-cumsum(w[event.mask.mask] / at.risk))
				);
			if (any(event.mask != is.event) || !is.subset) {
				times.ok <- y[event.mask, 1];
				c(1., surv)[sapply(times, function(z) {a <- which(z >= c(0., times.ok)); a[length(a)]})]
			} else if (!is.exact) {
				# only event time has survival prob, but many be duplicated
				surv[is.nonduplicated.event][time.subset];
			} else surv[is.nonduplicated.event];
			}
		));

	if (all(0 != times)) {
		times <- c(0, times);
		kkm.est <- cbind(1., kkm.est);
		}
	colnames(kkm.est) <- as.character(times);
	rownames(kkm.est) <- rownames(xtest);

	return(
		structure(
			list(
				test.predicted.survival = if(any(if.test)) kkm.est[if.test, ] else NULL,
				test.observed.survival = ytest,
				train.predicted.survival = kkm.est[!if.test, ],
				train.observed.survival = y,
				time.interest = times, 
				type = type,
				kernel = kernel, x.scale = x.scale,
				call = match.call()
				),
			class = 'kkm'
			)
		);
	}


#' @title Plot fitted survival curves
#' @param x `kkm` object from `kkm`
#' @param type Either 'test' or 'train'
#' @param subset Subset of samples to plot, either a vector of logical or subscript/index
#' @param xlab Name of x label
#' @param ylab Name of y label
#' @param ... Other parameters passed to `lattice::xyplot`
#' @return A `trellis` object generated by `lattice::xyplot`
#' @export
plot.kkm <- function(x, type = c('test', 'train'), subset, xlab = 'Time', ylab = 'Survival', ...) {
	type <- match.arg(type);
	surv <- paste0(type, '.predicted.survival');
	stopifnot(require(reshape2));
	stopifnot(require(lattice));
	z <- if (missing(subset)) melt(x[[surv]]) else melt(x[[surv]][subset, ]);
	xyplot(value ~ Var2, data = z, groups = Var1, type = 's', xlab = xlab, ylab = ylab, ...);
	}


if (0) {
	library(survival);
	data(pbc, package = 'randomForestSRC');
	pbc <- na.omit(pbc);
	i.tr <- sample(nrow(pbc), 100);

	x <- model.matrix(Surv(days, status) ~., data = pbc)[, -1];
	y <- Surv(pbc$days, pbc$status);

	kkm.pred1 <- kkm(x[i.tr, ], y[i.tr, ], xtest = x[-i.tr,], ytest = y[i.tr,], 
		kernel = 'laplacedot', kpar = list(sigma = 0.1));
	kkm.pred2 <- kkm(Surv(days, status) ~., data = pbc[i.tr, ], xtest = pbc[-i.tr, ], 
		k = 20L, kernel = 'laplacedot', kpar = list(sigma = 0.1));

	plot(kkm.pred1, subset = sample(nrow(xtest), 5))
	plot(kkm.pred2)

	survConcordance(y[-i.tr] ~I(1- kkm.pred1$test.predicted.survival[, 30]))$concordance

	library(randomForestSRC);
	rf <- rfsrc(Surv(days, status) ~., data = pbc[i.tr, ]);
	rf.pred <- predict(rf, pbc[-i.tr, ]);
	survConcordance(y[-i.tr] ~ rf.pred$predicted)$concordance
	survConcordance(y[-i.tr] ~ I(1-rf.pred$survival[, 30]))$concordance

	comp <- t(sapply(
		seq_along(rf.pred$time.interest), 
		function(i) {
			c(rf = survConcordance(y[-i.tr] ~ I(1-rf.pred$survival[, i]))$concordance,
				kernel = survConcordance(y[-i.tr] ~ I(1- kkm.pred$test.pred[, i]))$concordance)
			}
		));
	comp <- cbind(comp, diff=comp[, 2] - comp[, 1]);
	print(comp);
	apply(comp[, 1:2], 2, max);
	}
