# implement extreme learning machine for survival, regression, classification
# TODO:
# 	Implement deep extreme learning machine
#
#' @title Extreme learning machine for survival, classification, regression (univariate / multivariate)
#' @param formula Model formula
#' @param data Training data frame
#' @param ... Further parameters passed to internal function
#' @param x Design matrix
#' @param y Survival object, factor or a vector/matrix of continuous variable
#' @param nhid Number of (random) hidden neurons
#' @param actfun Activating function
#' @param scaled If to standardize input units
#' @param family Choices of 'cox', 'binomial', 'multinomial', 'poisson', 'gaussian', 'mgaussian' 
#' @param lambda Shrinkage factor
#' @return elm object.
#' @examples
#' library(survival);
#' data(pbc, package = 'randomForestSRC');
#' pbc <- na.omit(pbc);
#' i.tr <- sample(nrow(pbc), 100);
#' 
#' elm.f <- elm(Surv(days, status) ~., data = pbc[i.tr, ], nhid = 500);
#' elm.pred <- predict(object = elm.f, newdata = pbc[-i.tr, ], type = 'link');
#' survConcordance(Surv(days, status) ~ elm.pred, data = pbc[-i.tr, ])
#' @export
elm <- function(x, ...) UseMethod('elm');

#' @rdname elm
#' @export
elm.formula <- function(formula, data = environment(formula), ...) {
	formula <- stats::formula(terms(formula, data = data));
	mf <- model.frame(formula, data);
	x <- model.matrix(formula, data);
	y <- model.response(mf);
	xint <- match("(Intercept)", colnames(x), nomatch = 0);
	if (xint > 0) x <- x[, -xint, drop = FALSE];
	rm(mf);
	out <- elm.default(x, y, ...);
	out$formula <- formula;
	out$call <- match.call();
	class(out) <- c('elm.formula', class(out));
	return(out);
	}

#' @rdname elm
#' @export
elm.default <- function(
	x, y, nhid = max(200, 2*ncol(x)), actfun = 'sig', scaled = TRUE,
	family = switch(class(y), 
		'Surv' = 'cox', 
		'factor' = ifelse(nlevels(y) > 2, 'multinomial', 'binomial'),
		'integer' = 'poisson', 
		'matrix' = 'mgaussian',
		'gaussian'), 
	lambda = ifelse(family == 'mgaussian', 0, 1.0),
	...) {
	stopifnot(suppressMessages(require(glmnet)));
	if (nhid < 1) {
		stop("ERROR: number of hidden neurons must be >= 1")
		}
	x.scale <- NULL;
	y.scale <- NULL;
	if (scaled) {
		x <- scale(x);
		x.scale <- attributes(x)[c('scaled:center', 'scaled:scale')];
		attributes(x)[c('scaled:center', 'scaled:scale')] <- NULL;
		if ('gaussian' == family) {
			y <- scale(y);
			y.scale <- attributes(y)[c('scaled:center', 'scaled:scale')];
			attributes(y)[c('scaled:center', 'scaled:scale')] <- NULL;
			if(1 == ncol(y)) y <- c(y);
			}
		}
	inpweight <- matrix(runif(nhid*ncol(x), -1, 1), ncol = nhid);
	tempH <- x %*% inpweight;
	biashid <- runif(nhid, min = -1, max = 1);
	biasMatrix <- matrix(
		rep(biashid, nrow(x)), 
		nrow = nrow(x), ncol = nhid, byrow = TRUE
		);
	tempH = tempH + biasMatrix;
	actfun2 <- actfun;
	H <- switch(actfun,
		"sig"      = 1/(1 + exp(-1 * tempH)),
		"radbas"   = exp(-1 * (tempH^2)),
		"hardlim"  = ifelse(tempH >= 0, 1,  0),
		"hardlims" = ifelse(tempH >= 0, 1, -1),
		"satlins"  = ifelse(tempH >= 1, 1, ifelse(x <= -1, -1, x)),
		"tansig"   = 2/(1 + exp(-2 * tempH)) - 1,
		"tribas"   = ifelse(tempH >= -1 & tempH <= 1, 1-abs(tempH), 0),
		"poslin"   = ifelse(tempH < 0, 0, tempH),
		"purelin"  = tempH,
		do.call(actfun2, list(tempH))
		);
	if ('cox' == family) {
		stopifnot(suppressMessages(require(survival)));
		}
	if (family %in% c('gaussian', 'mgaussian') && lambda <= 0) {
		outlayer <- ginv(H) %*% y;
	} else {
		outlayer <- glmnet(
			x = H, y = y, family = family, 
			alpha = 0, lambda = lambda, ...
			);
		}
	model = list(
		inpweight = inpweight, biashid = biashid,
		outlayer = outlayer, nhid = nhid, actfun = actfun, 
		x.scale = x.scale, y.scale = y.scale, 
		family = family, lambda = lambda 
		);
	model$call <- match.call();
	class(model) <- "elm";
	return(model);
	}

#' @param object elm object
#' @param newdata A data frame (if elm called via formula) or a design matrix (if elm called via design matrix)
#' @param type Type of output
#' @rdname elm
#' @export
predict.elm <- function (
	object, newdata, 
	type = c('response', 'probabilities', 'link', 'risk')
	) {
	type <- match.arg(type);
	if (is(object, 'elm.formula')) {
		newdata <- model.matrix(object$formula[-2], newdata);
		xint <- match("(Intercept)", colnames(newdata), nomatch = 0);
		if (xint > 0) newdata <- newdata[, -xint, drop = FALSE];
		}
	if (!is.null(object$x.scale)) {
		newdata <- scale(newdata, 
			center = object$x.scale[['scaled:center']], 
			scale = object$x.scale[['scaled:scale']]
			);
		}
	tmpHTest = newdata %*% object$inpweight;
	biasMatrix <- matrix(
		rep(object$biashid, nrow(newdata)), 
		nrow = nrow(newdata), ncol = object$nhid, byrow = TRUE
		);
	tmpHTest = tmpHTest + biasMatrix;
	actfun2 <- object$actfun;
	HTest <- switch(object$actfun,
		"sig"      = 1/(1 + exp(-1 * tmpHTest)),
		"radbas"   = exp(-1 * (tmpHTest^2)),
		"hardlim"  = ifelse(tmpHTest >= 0, 1,  0),
		"hardlims" = ifelse(tmpHTest >= 0, 1, -1),
		"satlins"  = ifelse(tmpHTest >= 1, 1, ifelse(x <= -1, -1, x)),
		"tansig"   = 2/(1 + exp(-2 * tmpHTest)) - 1,
		"tribas"   = ifelse(tmpHTest >= -1 & tmpHTest <= 1, 1-abs(tmpHTest), 0),
		"poslin"   = ifelse(tmpHTest < 0, 0, tmpHTest),
		"purelin"  = tmpHTest,
		do.call(actfun2, list(tmpHTest))
		);
	if (is(object$outlayer, 'glmnet')) {
		if (object$family %in% c('binomial', 'multinomial')) {
			if ('response' == type) {type <- 'class';}
			if ('probabilities' == type) {type <- 'response';}
			} else if ('cox' == object$family) { 
				if ('risk' == type) {type <- 'response';}
				}
		pred <- predict(object$outlayer, newx = HTest, type = type);
	} else {
		pred <- HTest %*% object$outlayer;
		}
	if (object$family %in% c('gaussian', 'mgaussian') && !is.null(object$y.scale)) {
		pred <- pred*object$y.scale[['scaled:scale']] + object$y.scale[['scaled:center']];
		}
	if(is.matrix(pred)) {rownames(pred) <- rownames(newdata);}
	return(pred);
	}
