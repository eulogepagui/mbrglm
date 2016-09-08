## * the functions 'mbrglm' and 'mbrglm.fit' were written using as basis structure
## the code of  'brglm' and 'brglm.fit', respectively.
## * 'brglm' and 'brglm.fit' are implemented in the "R package brglm" version 0.5-9 by Ioannis Kosmidis.
## * the functions 'print.mbrglm', 'summary.mbrglm' and  'print.summary.mbrglm' are 
## modifications of  'print.glm', 'summary.glm' and  'print.summary.glm', respectively.
## Euloge C. Kenne Pagui <kenne@stat.unipd.it> [05/09/2016]

mbrglm <-
function (formula, family = binomial, data, weights, subset, 
    na.action, start = NULL, etastart, mustart, offset,
    model = TRUE, method = "mbrglm.fit",  x = FALSE,
    y = TRUE, contrasts = NULL,control.glm=glm.control()
    , control.mbrglm=mbrglm.control(),  ...)
{
   call <- match.call()
   if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
   if (is.function(family)) 
      family <- family()
   if (is.null(family$family)) 
   {
      print(family)
      stop("'family' not recognized")
   }
   mbr <- method == "mbrglm.fit"
   if (mbr & family$family != "binomial") 
      stop("families other than 'binomial' are not currently implemented")
   if (missing(data)) 
      data <- environment(formula)
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data", "subset", "weights", "na.action", 
     "etastart", "mustart", "offset"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf$drop.unused.levels <- TRUE
   mf[[1]] <- as.name("model.frame")
   mf <- eval(mf, parent.frame())
    
   switch(method, model.frame = return(mf), glm.fit = fit.proc <- glm.fit, 
     mbrglm.fit = fit.proc <- mbrglm.fit, stop("invalid 'method': ", 
         method))
   if (mbr) {
      formals(fit.proc)$control.mbrglm <- control.mbrglm
   }
    
   mt <- attr(mf, "terms")
   Y <- model.response(mf, "any")
   if (length(dim(Y)) == 1) 
   {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm)) 
         names(Y) <- nm
   }
    
   X <- if (!is.empty.model(mt)) 
      model.matrix(mt, mf, contrasts)
   else matrix(, NROW(Y), 0L) 
#      Xor <- if (!is.empty.model(mt)) 
#         model.matrix(mt, mf, contrasts)
#     else matrix(, NROW(Y), 0)
#     Xmax <- apply(abs(Xor), 2, max)
#     Xmax[Xmax == 0] <- 1
#     X <- sweep(Xor, 2, Xmax, "/")
   weights <- as.vector(model.weights(mf))
   if (!is.null(weights) && !is.numeric(weights)) 
      stop("'weights' must be a numeric vector")
    
   offset <- as.vector(model.offset(mf))
   if (!is.null(weights) && any(weights < 0)) 
      stop("negative weights not allowed")
   if (!is.null(offset)) {
      if (length(offset) == 1) 
         offset <- rep(offset, NROW(Y))
      else if (length(offset) != NROW(Y)) 
         stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
              length(offset), NROW(Y)), domain = NA)
   }
   mustart <- model.extract(mf, "mustart")
   etastart <- model.extract(mf, "etastart")
   par.names <- colnames(X)

   fit <- fit.proc(x = X, y = Y, weights = weights, start = start,
     etastart = etastart, mustart = mustart, offset = offset, 
     family = family, control = control.glm,control.mbrglm = control.mbrglm,
     intercept = attr(mt,"intercept") > 0) 
    
   fit <- c(fit, list(call = call, formula = formula, terms = mt, 
       data = data, offset = offset, method = method,  
       contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
   class(fit) <- c("mbrglm","glm","lm")
   fit
}


mbrglm.fit <-
   function (x, y, weights = rep(1, nobs), start = NULL,etastart = NULL,
             mustart = NULL, offset = rep(0, nobs), family=binomial(),
           control = glm.control(),  control.mbrglm = mbrglm.control(),intercept = TRUE)
{
   
modification = function(X,mu.eta,mu,A,B,InfoInv,weights)
{
   X<-as.matrix(X)
   n <- nrow(X)
   p <- ncol(X)
   nu_r_s_t <- nu_r_st <- array(0,c(p,p,p))
   for (r in 1:p)
   {
      nu_r_s_t[r,,] <- t(X)%*%((A^3*mu*(1-mu)*(1-2*mu)/weights^2*X[,r])*X)
      nu_r_st[r,,] <- t(X)%*%((B*mu.eta*X[,r])*X)
   }
   mod <- rep(0,p)
   out <- .C('modification',
            as.integer(p),
            as.double(InfoInv),
            as.double( nu_r_s_t),
            as.double(nu_r_st),
            mod=as.double(mod)
            )
   out$mod
}
   link <- family$link
   ## here we add the second derivative of the link function to 
   ## each element of the binomial's family  
   family$mu.eta.eta<-enrich(make.link(link),with="inverse link derivatives")$d2mu.deta
   x <- as.matrix(x)
   nobs <- NROW(y)
   nvars <- ncol(x)
   EMPTY <- nvars == 0
   if (is.null(weights)) 
      weights <- rep.int(1, nobs)
   if (is.null(offset)) 
      offset <- rep.int(0, nobs)
   dev.resids <- family$dev.resids
   variance <- family$variance
   linkinv <- family$linkinv
   if (!is.function(variance) || !is.function(linkinv)) 
      stop("'family' argument seems not to be a valid family object")
   if (EMPTY) 
   {
      return(glm.fit(x = x, y = y, weights = weights, start = start, 
         etastart = etastart, mustart = mustart, offset = offset, 
         family = family, control = control, intercept = intercept))
   }
   valideta <- family$valideta
   if (is.null(valideta)) 
      valideta <- function(eta) TRUE
   validmu <- family$validmu
   if (is.null(validmu)) 
      validmu <- function(mu) TRUE
   if (is.null(mustart)) 
   {
      eval(family$initialize)
      etastart <- family$linkfun(mustart)
      Astart <- weights
      mu.etastart <- family$mu.eta(etastart)
      if (link =="probit" | link =="cloglog" )
      {
         Astart <- (weights*mu.etastart)/(mustart*(1-mustart))
      } 
      W <-diag(Astart*mu.etastart)
      par <- solve(t(x)%*%W%*%x)%*%t(x)%*%W%*%y  
   }
   else 
   {
      mukeep <- mustart
      eval(family$initialize)
      mustart <- mukeep
      etastart <- family$linkfun(mustart)
      Astart <- weights
      mu.etastart <- family$mu.eta(etastart)
      if (link =="probit" | link =="cloglog" )
      {
         Astart <- (weights*mu.etastart)/(mustart*(1-mustart))
      } 
      W <-diag(Astart*mu.etastart)
      par <- solve(t(x)%*%W%*%x)%*%t(x)%*%W%*%y  
   }
   if (!is.null(start))
   {
      if (length(start) != nvars) 
           stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                nvars, paste(deparse(xnames), collapse = ", ")), 
                domain = NA)
      else par <-start
   } 
   # print(y)
   warnValue <- options(warn = -1)     
   prior.weights <- weights       
   nonzero.w <- which(weights != 0)
   x <- as.matrix(x[nonzero.w, ]); y <- y[nonzero.w];weights <- weights[nonzero.w]
   etastart <- etastart[nonzero.w]; mustart <- mustart[nonzero.w];
   offset <- offset[nonzero.w];
   temp.fit <- glm.fit(x = x, y = y, weights = weights, start = start, 
     etastart = etastart, mustart = mustart, offset = offset, 
     family = family, control = control, intercept = intercept)
   redundant <- is.na(temp.fit$coefficients)
   if (any(redundant)) 
   {
      x <- x[, -which(redundant), drop = FALSE]
      nvars <- nvars - sum(redundant)
   }
   temp.fit <- glm.fit(x = x, y = y, weights = weights, start = start, 
     etastart = etastart, mustart = mustart, offset = offset, 
     family = family, control = control, intercept = intercept)
      
   step <- .Machine$integer.max
   nIter <- 0
   test <- TRUE
   while ( test & (nIter < control.mbrglm$mbr.maxit))
   {
      nIter <- nIter + 1  
      eta = drop(x%*%par) + offset   
      mu <- family$linkinv(eta)
      mu.eta <- family$mu.eta(eta) 
      mu.eta.eta <- family$mu.eta.eta(eta)
      A <- (weights*mu.eta)/(mu*(1-mu))
      B <- (weights*( mu.eta.eta/(mu*(1-mu)) 
           + (mu.eta^2*(2*mu-1))/(mu*(1-mu))^2))
      info <- t(x)%*%((W.i<-A*mu.eta)*x)
      score <- drop(t(x)%*%(A*(y-mu)))
      InfoInv <- try(chol2inv(chol(info)),TRUE)
      if(failedInv <- inherits(InfoInv, "try-error")) 
      {
         warning("failed to invert the information matrix: iteration stopped prematurely")
         break
      }
      mod <- modification(x,mu.eta,mu,A,B,InfoInv,weights)
      modscore <- score + info%*%mod
      u <- (y-mu)/mu.eta
      y.adj <- x%*%(par+mod)+u
      par <- InfoInv%*%t(x)%*%(W.i*y.adj)
      if (control.mbrglm$mbr.trace) 
      {
         cat("\n")
         cat("Iteration:", nIter, "\n")
         cat("Coefficients:", drop(par),"\n") 
         cat("Modified scores:", modscore, "\n")
      }
      test <- sqrt(crossprod(drop(modscore))) > control.mbrglm$mbr.epsilon
   }
   
   rownames(par)<- colnames(x)
   options(warnValue)
   fitted <- family$linkinv(eta)
   W <- diag(weights*fitted*(1-fitted))
   h <- diag(W^{1/2}%*%x%*%chol2inv(chol(t(x)%*%W%*%x))%*%t(x)%*%W^{1/2})
   PearsonResid <- (y*weights-weights*fitted)/sqrt(weights*fitted*(1-fitted))
   StdPearsonResid <- PearsonResid /sqrt(1-h)
   ## fitted probabilities required for null deviance calculation
   fitted.null <- (1+6*sum(y*weights))/(2+6*sum(weights))
   if(link=="probit")
   {
      par.null <- qnorm(fitted.null)
   } 
   if(link=="cloglog")
   {
      par.null<- log(-log(1-fitted.null))
   }
   temp.fit$residuals <- StdPearsonResid
   temp.fit$coefficients <- drop(par)
   temp.fit$fitted.values <- fitted
   temp.fit$linear.predictors <- eta
   temp.fit$FisherInfo <- info
   temp.fit$FisherInfoInvs <- InfoInv
   temp.fit$nIter <- nIter
   temp.fit$ModifiedScores <- c(modscore)
   temp.fit$prior.weights <- prior.weights
   temp.fit$weights <- weights
   temp.fit$deviance <- sum(dev.resids(temp.fit$y, temp.fit$fitted.values,
        temp.fit$weights))
   temp.fit$null.deviance <- sum(dev.resids(temp.fit$y, rep(fitted.null,length(temp.fit$y)),
        temp.fit$weights))
   temp.fit$converged <- nIter < control.mbrglm$mbr.maxit
   if (!temp.fit$converged) 
      warning("Iteration limit reached")    
   temp.fit
}     

print.mbrglm <-
function (x, digits = max(3, getOption("digits") - 3), ...)
     {
          if (x$method == "glm.fit" | !(nPars <- length(coef(x)))) {
               class(x) <- class(x)[-match("brglm", class(x))]
               return(print(x, digits, ...))
          }
          cat("\nCall: ", deparse(x$call), "\n\n")
          if (nPars) {
               cat("Standardized Pearson residual:\n") 
               print.default(format(summary(x$residuals), digits = digits),
                    print.gap = 2, quote = FALSE)
               cat("\n")
               cat("Coefficients")
               if (is.character(co <- x$contrasts))
                    cat("  [contrasts: ", apply(cbind(names(co), co),
                         1, paste, collapse = "="), "]")
               cat(":\n")
               print.default(format(x$coefficients, digits = digits),
                    print.gap = 2, quote = FALSE)
          }
          else cat("No coefficients\n\n")
          cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
               x$df.residual, "Residual\n")
          if (nchar(mess <- naprint(x$na.action)))
               cat("  (", mess, ")\n", sep = "")
          if (!is.null(x$penalized.deviance))
               cat("Deviance:\t   ", format(round(x$deviance, digits)),
                    "\nPenalized Deviance:", format(round(x$penalized.deviance,
                         digits)), "\tAIC:", format(round(x$aic, digits)),
                    "\n")
          else cat("Deviance:\t   ", format(round(x$deviance, digits)),"\n")
               #"\tAIC:", format(round(x$aic, digits)), "\n")
          invisible(x)
     }


print.summary.mbrglm <-
function (x, digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
          signif.stars = getOption("show.signif.stars"), ...)
     {
        
          cat("\nCall:\n")
          cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
               "\n\n", sep = "")
          if (length(x$aliased) == 0) {
               cat("\nNo Coefficients\n")
          }
          else {
               df <- if ("df" %in% names(x))
                    x[["df"]]
               else NULL
               if (!is.null(df) && (nsingular <- df[3] - df[1]))
                    cat("\n Coefficients: (", nsingular, " not defined because of singularities)\n",
                         sep = "")
               else cat("\nCoefficients:\n")
               coefs <- x$coefficients
               if (!is.null(aliased <- x$aliased) && any(aliased)) {
                    cn <- names(aliased)
                    coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
                         colnames(coefs)))
                    coefs[!aliased, ] <- x$coefficients
               }
               printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                    na.print = "NA", ...)
          }
          cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
               format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null",
                    "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance",
                         "deviance")]), digits = max(5, digits + 1)), " on",
                    format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
                    1, paste, collapse = " "), sep = "")
#           cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
#                format(x$dispersion), ")\n", sep = "")
#           if (!is.null(x$penalized.deviance))
#                cat("Penalized deviance:", format(round(x$penalized.deviance,
#                     digits = max(5, digits + 1))), "\n")
          if (nchar(mess <- naprint(x$na.action)))
               cat("  (", mess, ")\n", sep = "")
#           cat("AIC: ", format(x$aic, digits = max(4, digits + 1)),
#                "\n")
          correl <- x$correlation
          if (!is.null(correl)) {
               p <- NCOL(correl)
               if (p > 1) {
                    cat("\nCorrelation of Coefficients:\n")
                    if (is.logical(symbolic.cor) && symbolic.cor) {
                         print(symnum(correl, abbr.colnames = NULL))
                    }
                    else {
                         correl <- format(round(correl, 2), nsmall = 2,
                              digits = digits)
                         correl[!lower.tri(correl)] <- ""
                         print(correl[-1, -p, drop = FALSE], quote = FALSE)
                    }
               }
          }
          cat("\n")
          invisible(x)
     }

summary.mbrglm <-
function (object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE,
          ...)
     {
          if (object$method == "glm.fit")
               return(summary.glm(object, dispersion = NULL, correlation = FALSE,
                    symbolic.cor = FALSE, ...))
          df.r <- object$df.residual
          if (is.null(dispersion))
               dispersion <- 1
          aliased <- is.na(coef(object))
          p <- object$rank
          if (p > 0) {
               p1 <- 1:p
               Qr <- object$qr
#                coef.p <- object$coefficients[Qr$pivot[p1]]
               coef.p <- object$coefficients
#                covmat.unscaled <- chol2inv(chol(object$FisherInfo))
               covmat.unscaled <- object$FisherInfoInv
               covmat <- dispersion * covmat.unscaled
               var.cf <- diag(covmat)
               s.err <- sqrt(var.cf)
               tvalue <- coef.p/s.err
               dn <- c("Estimate", "Std. Error")
               pvalue <- 2 * pnorm(-abs(tvalue))
               coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
               dimnames(coef.table) <- list(names(coef.p), c(dn, "z value",
                    "Pr(>|z|)"))
               residuals = object$residuals
               df.f <- NCOL(Qr$qr)
          }
          else {
               coef.table <- matrix(, 0, 4)
               dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
                    "t value", "Pr(>|t|)"))
               covmat.unscaled <- covmat <- matrix(, 0, 0)
               df.f <- length(aliased)
          }
#           keep <- match(c("call", "terms", "family", "deviance", "aic",
#                "contrasts", "df.residual", "null.deviance", "df.null",
#                "iter", "na.action", "penalized.deviance"), names(object),
#                0)
          keep <- match(c("call", "terms", "family", 
               "contrasts",  "nIter", "na.action", "df.residual", "null.deviance", "df.null","deviance"), names(object),
               0)
          ans <- c(object[keep], list(coefficients = coef.table, aliased = aliased,
               dispersion = dispersion, df = c(object$rank, df.r, df.f),
               cov.unscaled = covmat.unscaled, cov.scaled = covmat,residuals=residuals))
          if (correlation && p > 0) {
               dd <- sqrt(diag(covmat.unscaled))
               ans$correlation <- covmat.unscaled/outer(dd, dd)
               ans$symbolic.cor <- symbolic.cor
          }
          class(ans) <- "summary.mbrglm"
          return(ans)
     }

