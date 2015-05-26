# extent the idea of random feature in extreme learning machine to survival 

library(elmNN);
library(survival);
library(doMC);
registerDoMC(cores = 10);

source('../../GroupShare/load_prostate.R');
seed <- 20150521;
dat <- load.prostate(train.subset = 'train', seed = seed);
test <- load.prostate(train.subset = 'holdout', seed = seed);

## backward selected (in SAS)
bwe.var <- c("REGION_C", "TRT2_ID", "ALP", "AST", "HB", "ALB", "LYMPH_NODES", "LIVER", 
        "ADRENAL", "PROSTATECTOMY", "ANALGESICS", "ACE_INHIBITORS", "ESTROGENS", "MHCARD");

# survival version of extreme learning machine
elmtrain.Surv <- function (x, y, lambda = 0.001, nhid, actfun, ...) {
        require(glmnet);
    if (nhid < 1) 
        stop("ERROR: number of hidden neurons must be >= 1")
    inpweight <- elmNN:::randomMatrix(nhid, ncol(x), -1, 1)
        tempH <- x %*% inpweight;
    biashid <- runif(nhid, min = -1, max = 1)
    biasMatrix <- matrix(rep(biashid, nrow(x)), nrow = nrow(x), 
        ncol = nhid, byrow = T)
    tempH = tempH + biasMatrix
        actfun2 <- actfun
        H <- switch(actfun,
                "sig"      = 1/(1 + exp(-1 * tempH)),
                "sin"      = sin(tempH),
                "radbas"   = exp(-1 * (tempH^2)),
                "hardlim"  = hardlim(tempH),
                "hardlims" = hardlims(tempH),
                "satlins"  = satlins(tempH),
                "tansig"   = 2/(1 + exp(-2 * tempH)) - 1,
                "tribas"   = tribas(tempH),
                "poslin"   = poslin(tempH),
                "purelin"  = tempH,
                do.call(actfun2, list(tempH))
                )
        # H <- as.data.frame(H);
        # data <- data.frame(surv = y, H);
        cox <- glmnet(x = H, y = y, family = 'cox', alpha = 0, lambda = lambda, ...);

    model = list(inpweight = inpweight, biashid = biashid,  nhid = nhid,
        actfun = actfun, fitted = cox)
    model$call <- match.call()
    class(model) <- "elmNN.coxph"
    model
        }


# predict function for extreme learning machine

predict.elmNN.coxph <- function (object, newdata = NULL, ...) {
        tmpHTest = newdata %*% object$inpweight
        biasMatrix <- matrix(rep(object$biashid, nrow(newdata)), nrow = nrow(newdata), 
                ncol = object$nhid, byrow = T)
        tmpHTest = tmpHTest + biasMatrix
        actfun2 <- object$actfun
        HTest <- switch(object$actfun,
                "sig"      = 1/(1 + exp(-1 * tmpHTest)),
                "sin"      = sin(tmpHTest),
                "radbas"   = exp(-1 * (tmpHTest^2)),
                "hardlim"  = hardlim(tmpHTest),
                "hardlims" = hardlims(tmpHTest),
                "satlins"  = satlins(tmpHTest),
                "tansig"   = 2/(1 + exp(-2 * tmpHTest)) - 1,
                "tribas"   = tribas(tmpHTest),
                "poslin"   = poslin(tmpHTest),
                "purelin"  = tmpHTest,
                do.call(actfun2, list(tmpHTest))
                );
        predict(object$fitted, newx = HTest, lambda = object$fitted$lambda)
    }


# all covariate
sel.var <- setdiff(names(dat), c('survDays', 'death'));
# sel.var <- bwe;
fm <- as.formula(paste('~', paste(sel.var, collapse = '+')))
x <- model.matrix(fm, data = dat)[, -1];
x <- scale(x);
y <- with(dat, Surv(survDays, death));
newdata <- model.matrix(fm, data = test)[, -1]
newdata <- scale(newdata, center = attr(x, 'scaled:center'), scale = attr(x, 'scaled:scale'));

coxELM <- elmtrain.Surv(x, y, nhid = 500, actfun = 'sig', lambda = 0.5);
cox.pred <- predict.elmNN.coxph(coxELM, newdata);
print(survConcordance(Surv(test$survDays, test$death) ~ cox.pred)$concordance)


# bagged elm cox

bagging <- function(x, y, B = 500, sampsize = nrow(x), replace = TRUE, fit, predict = predict, aggregate = function(x) Reduce('+', x) / length(x), ...) {
        f <- foreach(i = 1:B) %dopar% { 
                ind <- sample(nrow(x), size = sampsize, replace = replace);
                fit(x[ind, ], if (is.matrix(y)) y[ind, ] else y[ind], ...)
                }
        structure(list(fitted = f, predict = predict, aggregate = aggregate) , class = 'bagging')
        }
predict.bagging <- function(object, newdata, ...) {
        pred <- foreach(i = 1:length(object$fitted)) %dopar% {
                object$predict(object$fitted[[i]], newdata = newdata)
                }
        object$aggregate(pred)
        }


bag.coxELM <- bagging(
        x, y, B = 500, 
        fit = elmtrain.Surv, predict = predict.elmNN.coxph, 
        nhid = 30, actfun = 'sig', lambda = 0.01, standardize = FALSE);
bag.coxELM.pred <- predict(bag.coxELM, newdata);
print(survConcordance(Surv(test$survDays, test$death) ~ bag.coxELM.pred)$concordance)
