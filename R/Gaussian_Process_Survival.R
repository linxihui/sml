#' @title Gaussian process for survival, regression and classification
#' @description It works for classification, regression, cox regression and poission regression
#' The first/main purpose is to implement Gaussian process into Cox's model
#' The glmnet is used to solve a general Ridge-regularized regression solver. This may be changed in the future.
#' @param x Design matrix (NO intercept)
#' @param y Reponse vector of type double, integer, factor or Surv
#' @param kernel, kpar check kernlab::gausspr
#' @param family check glmnet::glmnet
#' @return A 'gprsc' object
#' @examples 
#' ## classification
#' library(pROC);
#' library(mlbench);
#' data(Sonar);
#' i.tr <- sample(nrow(Sonar), size = 0.6*nrow(Sonar));
#' 
#' gp <- gpsrc(Class ~., data = Sonar[i.tr,], kpar = list(sigma = 0.01))
#' gp.pred <- predict(gp, Sonar[-i.tr, ], type = 'prob')
#' cat('AUC:', roc(y[-i.tr], gp.pred[, 2])$auc, '-- gpsrc.formula', '\n');
#' 
#' ## regression
#' 
#' data(BostonHousing);
#' BH <- BostonHousing;
#' i.tr <- sample(nrow(BH), 200);
#' 
#' gp2 <- gpsrc(medv ~., data = BH[i.tr, ], kpar = list(sigma = 0.1));
#' gp2.pred <- predict(gp2, BH[-i.tr, ])
#' cat('RMSE:', sqrt(mean((BH$medv[-i.tr] - gp2.pred)^2)), '-- gpsrc.formula\n');
#' 
#' ## survival
#' 
#' library(survival);
#' data(pbc, package = 'randomForestSRC');
#' pbc <- na.omit(pbc);
#' i.tr <- sample(nrow(pbc), 100);
#' 
#' gp <- gpsrc(Surv(days, status) ~., data = pbc[i.tr, ], 
#' 	kernel = 'laplacedot', kpar = list(sigma = 0.1));
#' gp.pred <- predict(gp, pbc[-i.tr, ])
#' cat('C-index:', survConcordance(y[-i.tr] ~ gp.pred)$concordance, '-- gpsrc.default\n');
#' @export
gpsrc <- function(x, ...) UseMethod('gpsrc');


#' @rdname gpsrc
#' @export gpsrc.formula
gpsrc.formula <- function(formula, data, ...) {
	formula <- stats::formula(terms(formula, data = data));
	mf <- model.frame(formula, data);
	x <- model.matrix(formula, data);
	y <- model.response(mf);
	xint <- match("(Intercept)", colnames(x), nomatch = 0);
	if (xint > 0) x <- x[, -xint, drop = FALSE];
	rm(mf);
	out <- gpsrc.default(x, y, ...);
	out$formula <- formula;
	out$call <- match.call();
	class(out) <- c('gpsrc.formula', class(out));
	return(out);
	}

gpsrc.default <- function(
	x, y, scaled = TRUE, 
	kernel = 'rbf', kpar = 'automatic',
	var = 1.0,
	family = switch(class(y), 
		'Surv' = 'cox', 
		'factor' = 'binomial', 
		'integer' = 'poisson', 
		'gaussian'), 
	...) {

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
					cat(" Setting default kernel parameters ","\n");
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
				cat("Using automatic sigma estimation (sigest) for RBF or laplace kernel","\n");
				}
			}

		if(!is(kernel,"kernel")) {
			if (is(kernel,"function")) {kernel <- deparse(substitute(kernel));}
			kernel <- do.call(kernel, kpar);
			}
		if(!is(kernel,"kernel")) {stop("kernel must inherit from class `kernel'");}

		K <- kernelMatrix(kernel, x);
		# cholesky decomposition K = R'R,  f = theta = K alpha = R' beta =>  R alpha = beta
		R <- chol(K);
		
		require(glmnet);
		lambda <- switch(family, 
			'gaussian' = var/nrow(x),
			1 / nrow(x)
			);
		beta <- as.vector(glmnet(t(R), y, alpha = 0, lambda = lambda, family = family, ...)$beta);
		alpha <- backsolve(R, beta);

		return(structure(
			list(
				alpha = alpha, 
				kernel = kernel,
				xmatrix = x,
				x.scale = x.scale,
				y.scale = y.scale,
				family = family,
				lev = if(is.factor(y)) levels(y) else NULL
				), 
			class = 'gpsrc'
			));
		}

#' @param object 'gpsrc' object from \code{gpsrc}
#' @param newdata Design matrix of test set
#' @param type Type of output
#' @rdname gpsrc
predict.gpsrc <- function(
	object, newdata, 
	type = c('response', 'probabilities', 'link', 'risk')
	) {
	type <- match.arg(type);
	if (is(object, 'gpsrc.formula')) {
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
	pred <- kernelMult(
		object$kernel, newdata, 
		object$xmatrix, as.matrix(object$alpha)
		); # need to check carefully
	if ('gaussian' == object$family && !is.null(object$y.scale)) {
		pred <- pred*object$y.scale[['scaled:scale']] + object$y.scale[['scaled:center']];
		}
	if ('binomial' == object$family) {
		if ('response' == type) {
			pred <- object$lev[as.numeric(pred > 0) + 1];
		} else if ('probabilities' == type) {
			pred <- 1 / (1 + exp(-pred));
			pred <- cbind(1 - pred, pred);
			colnames(pred) <- object$lev
			}
		}
	if ('cox' == object$family && 'risk' == type ) {
		pred <- exp(pred)
		}
	return(pred);
	}


if (0) { # don't run
### TEST ###########################################################################################

## classification

library(pROC);
library(kernlab);
library(mlbench);
data(Sonar);
i.tr <- sample(nrow(Sonar), size = 0.6*nrow(Sonar));
x <- as.matrix(Sonar[, -ncol(Sonar)]);
y <- factor(Sonar[[ncol(Sonar)]]);

enet <- cv.glmnet(x[i.tr, ], y[i.tr], alpha = 0.5, family = 'binomial')
enet.pred <- predict(enet, x[-i.tr, ], type = 'response', s = 'lambda.min')
cat('AUC:', roc(y[-i.tr], enet.pred)$auc, '-- enet(CV)', '\n');

kl.gp <- gausspr(x[i.tr, ], y[i.tr], kpar = list(sigma = 0.01));
kl.gp.pred <- predict(kl.gp, x[-i.tr, ], type = 'prob');
cat('AUC:', roc(y[-i.tr], kl.gp.pred[, 2])$auc, '-- kernlab', '\n');


gp <- gpsrc(x[i.tr, ], y[i.tr], kpar = list(sigma = 0.01));
gp.pred <- predict(gp, x[-i.tr, ], type = 'prob');
cat('AUC:', roc(y[-i.tr], gp.pred[, 2])$auc, '-- gpsrc.default', '\n');

gp2 <- gpsrc(Class ~., data = Sonar[i.tr,], kpar = list(sigma = 0.01))
gp2.pred <- predict(gp2, Sonar[-i.tr, ], type = 'prob')
cat('AUC:', roc(y[-i.tr], gp2.pred[, 2])$auc, '-- gpsrc.formula', '\n');


## regression

data(BostonHousing);
BH <- BostonHousing;
i.tr <- sample(nrow(BH), 200);

x <- model.matrix(medv ~., data = BH)[, -1];
y <- BH$medv;

ridge <- cv.glmnet(x[i.tr, ], y[i.tr], alpha = 0);
ridge.pred <- predict(ridge, x[-i.tr, ], s = 'lambda.min');
cat('RMSE:', sqrt(mean((BH$medv[-i.tr] - ridge.pred)^2)), '-- Ridge\n');

kl.gp <- gausspr(medv ~., data = BH[i.tr, ], kpar = list(sigma = 0.1));
kl.gp.pred <- predict(kl.gp, BH[-i.tr, ])
cat('RMSE:', sqrt(mean((BH$medv[-i.tr] - kl.gp.pred)^2)), '-- kernlab\n');

gp <- gpsrc(x[i.tr, ], y[i.tr], kpar = list(sigma = 0.1));
gp.pred <- predict(gp, x[-i.tr, ]);
cat('RMSE:', sqrt(mean((BH$medv[-i.tr] - gp.pred)^2)), '-- gpsrc.default\n');

gp2 <- gpsrc(medv ~., data = BH[i.tr, ], kpar = list(sigma = 0.1));
gp2.pred <- predict(gp2, BH[-i.tr, ])
cat('RMSE:', sqrt(mean((BH$medv[-i.tr] - gp2.pred)^2)), '-- gpsrc.formula\n');


## survival

library(survival);
library(randomForestSRC);
data(pbc, package = 'randomForestSRC');
pbc <- na.omit(pbc);

i.tr <- sample(nrow(pbc), 100);

x <- model.matrix(Surv(days, status) ~., data = pbc)[, -1];
y <- Surv(pbc$days, pbc$status);

cox <- coxph(Surv(days, status) ~ ., data = pbc[i.tr, ], method = 'breslow')
cox.pred <- predict(cox, pbc[-i.tr, ]);
cat('C-index:', survConcordance(y[-i.tr] ~ cox.pred)$concordance, '-- CoxPH\n');

ridge <- cv.glmnet(x[i.tr, ], y[i.tr], family = 'cox', alpha = 0);
ridge.pred <- predict(ridge, x[-i.tr, ]);
cat('C-index:', survConcordance(y[-i.tr] ~ ridge.pred)$concordance, '-- Ridge(CV)\n');

lasso <- cv.glmnet(x[i.tr, ], y[i.tr], family = 'cox', alpha = 1);
lasso.pred <- predict(lasso, x[-i.tr, ]);
cat('C-index:', survConcordance(y[-i.tr] ~ lasso.pred)$concordance, '-- Lasso(CV)\n');

gp <- gpsrc(x[i.tr, ], y[i.tr], kernel = 'laplacedot', kpar = list(sigma = 0.1));
gp.pred <- predict(gp, x[-i.tr, ])
cat('C-index:', survConcordance(y[-i.tr] ~ gp.pred)$concordance, '-- gpsrc.default\n');

rf <- rfsrc(Surv(days, status) ~., data = pbc[i.tr, ]);
rf.pred <- predict(rf, pbc[-i.tr, ]);
cat('C-index:', survConcordance(y[-i.tr] ~ rf.pred$predicted)$concordance, '-- RF\n');

}
