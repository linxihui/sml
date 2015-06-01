#' @title Kernel Weighted Kaplan-Meier model for survival
#' @param formula Model formula
#' @param data Data frame
#' @param x Design matrix (NO intercept)
#' @param y Reponse vector of `Surv` object 
#' @param xtest If the 'formula-data' format used, `xtest` is a data frame of the test set.
#' 	If the 'x-y' method called, `xtest` is the design matrix of the test set
#' @param ytest (optional) Survival outcome for testset
#' @param times Times to evaluate survival probabilities. Currently no used. All unique event times in training are used
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
	scaled = TRUE, k = nrow(x), 
	times = unique(sort(y[as.logical(y[, 2]), 1])),
	kernel = 'rbfdot', kpar = 'automatic'
	) {
	stopifnot(suppressMessages(require(survival)));
	stopifnot(suppressMessages(require(kernlab)));
	x.scale <- NULL;
	if (scaled) {
		x <- scale(x);
		x.scale <- attributes(x)[c('scaled:center', 'scaled:scale')];
		attributes(x)[c('scaled:center', 'scaled:scale')] <- NULL;
		}
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

	if (!is.null(xtest)) {
		xtest <- scale(xtest, center = x.scale[['scaled:center']], scale = x.scale[['scaled:scale']]);
		}

	if (return.train.prediction) {
		xpred <- rbind(x, xtest);
		if.test <- c(rep(FALSE, nrow(x)), rep(TRUE, ifelse(is.null(xtest), 0,  nrow(xtest))));
	} else {
		xpred <- xtest;
		if.test <- rep(TRUE, ifelse(is.null(xtest), 0, nrow(xtest)));
		}

	kern.dist <- kernelMatrix(kernel, x, xpred)@.Data;

	if (k < nrow(x)) {
		knn.mask <- apply(kern.dist, 2, function(z) {
			o <- rep(0, length(z));
			o[order(-z)[1:k]] <- 1;
			return(o);
			})
		kern.dist <- kern.dist*knn.mask;
		}

	time.order <- order(-y[,1], y[,2]); # time sorted from largest to smallest
	is.event.sorted <- as.logical(y[time.order, 2]);
	event.sorted <- y[time.order, 1][is.event.sorted];

	output.index <- rev(!duplicated(event.sorted));
	times <- rev(y[time.order, 1][is.event.sorted])[output.index];

	kkm.est <- t(apply(
		X = kern.dist,
		MARGIN = 2,
		FUN = function(w) {
			w <- w[time.order];
			at.risk <- cumsum(w)[is.event.sorted];
			exp(cumsum(rev(log((at.risk - w[is.event.sorted]) / (at.risk)))))[output.index];
			}
		))

	colnames(kkm.est) <- as.character(times);

	return(
		structure(
			list(
				test.predicted.survival = if(any(if.test)) kkm.est[if.test, ] else NULL,
				test.observed.survival = ytest,
				train.predicted.survival = kkm.est[!if.test, ],
				train.observed.survival = y,
				time.interest = times, 
				kernel = kernel, x.scale = x.scale,
				call = match.call()
				),
			class = 'kkm'
			)
		);
	}



if (0) {
	library(survival);
	data(pbc, package = 'randomForestSRC');
	pbc <- na.omit(pbc);
	i.tr <- sample(nrow(pbc), 100);

	x <- model.matrix(Surv(days, status) ~., data = pbc)[, -1];
	y <- Surv(pbc$days, pbc$status);

	kkm.pred <- kkm(x[i.tr, ], y[i.tr, ], xtest = x[-i.tr,], ytest = y[i.tr,], kernel = 'laplacedot', kpar = list(sigma = 0.1));
	# kkm.pred <- kkm(Surv(days, status) ~., data = pbc[i.tr, ], xtest = pbc[-i.tr, ], kernel = 'laplacedot', kpar = list(sigma = 0.1));

	survConcordance(y[-i.tr] ~I(1- kkm.pred$test.predicted.survival[, 30]))$concordance

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