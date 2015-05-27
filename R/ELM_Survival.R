# implement extreme learning machine for survival, regression, classification
# TODO:
# 	Implement deep extreme learning machine

elm <- function(x, ...) UseMethod('elm');

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

elm.default <- function(
	x, y, nhid = 2*ncol(x), actfun = 'sig', scaled = TRUE,
	family = switch(class(y), 
		'Surv' = 'cox', 
		'factor' = ifelse(nlevels(y) > 2, 'multinomial', 'binomial'),
		'integer' = 'poisson', 
		'gaussian'), 
	lambda = 0.001,
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
    biashid <- runif(nhid, min = -1, max = 1)
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
	glmnet.fit <- glmnet(
		x = H, y = y, family = family, 
		alpha = 0, lambda = lambda, ...
		);
    model = list(
		inpweight = inpweight, biashid = biashid,
		nhid = nhid, actfun = actfun, 
		x.scale = x.scale, y.scale = y.scale, 
		family = family, fitted = glmnet.fit
		);
    model$call <- match.call();
    class(model) <- "elm";
    return(model);
    }

# predict function for extreme learning machine

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
	if (object$family %in% c('binomial', 'multinomial')) {
		if ('response' == type) {type <- 'class';}
		if ('probabilities' == type) {type <- 'response';}
	} else if ('cox' == object$family) { 
		if ('risk' == type) {type <- 'response';}
		}
	pred <- predict(object$fitted, newx = HTest, type = type);
	if ('gaussian' == object$family && !is.null(object$y.scale)) {
		pred <- pred*object$y.scale[['scaled:scale']] + object$y.scale[['scaled:center']];
		}
	if(is.matrix(pred)) {rownames(pred) <- rownames(newdata);}
	return(pred);
	}

if (0) {

library(survival);
data(pbc, package = 'randomForestSRC');
pbc <- na.omit(pbc);
i.tr <- sample(nrow(pbc), 100);

elm.f <- elm(Surv(days, status) ~., data = pbc[i.tr, ], nhid = 500, lambda = 0.01);
elm.pred <- predict(object = elm.f, newdata = pbc[-i.tr, ], type = 'link')
cat('C-index:', 
	survConcordance(Surv(days, status) ~ elm.pred, data = pbc[-i.tr, ])$concordance, 
	'\n');

}
