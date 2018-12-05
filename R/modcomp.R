#' A Function for Running Sequential Mixed-Effects Model Comparisons
#'
#' This function allows you run mixed-effects model comparisons with bootstrapping switched on or off
#' @param bootstr should bootstrapping be used. Defaults to FALSE.
#' @keywords bootstrapping, comparisons, mixed
#' @export
#' @examples
#' modelcomp()

modelcomp <- function(origmodel, resdf=data.frame(), compdf=data.frame(), bootstr=FALSE){
  pv <- function(val) { paste( "=", sub("^(-?)0.", "\\1.", sprintf("%.3f", val))) } ## Print P-value to text in APA format
  if (bootstr == FALSE){
    compmodel <- refitML(origmodel)
    resdf <- data.frame(subset(tidy(origmodel), group=="fixed")) ## create tidy df of the model summary
    names(resdf)[names(resdf) == 'estimate'] <- 'beta' ## rename beta column
    names(resdf)[names(resdf) == 'statistic'] <- 't' ## rename t column
    t <- resdf$beta/resdf$std.error
    tpval <- round( (2*(1 - pnorm(abs(t)))), 3) ## calculate p based on t and z distribution (use for intercept only)
    resdf$beta <- ifelse(abs(resdf$beta)<0.0001, scientific(resdf$beta), round(resdf$beta, 4))
    resdf$std.error <- ifelse(abs(resdf$std.error)<0.0001, scientific(resdf$std.error), round(resdf$std.error, 4))
    resdf$t <- ifelse(abs(t)<0.01, scientific(t), round(t, 2))
    resdf$t.pz <- ifelse(tpval<0.001,"< .001", pv(tpval))
    terms <- rev(attr(terms(compmodel),"term.labels")) ## extract predictor names
    for (i in 1:length(terms)){ ## sequentially decompose the model from highest to lowest order
      newformula <- paste(". ~ . - ",terms[i],"")
      compmodel2 <- update(compmodel, as.formula(newformula))
      an <- anova(compmodel, compmodel2)
      x2 <- an$Chisq[2]
      x2 <- ifelse(x2<0.01, scientific(x2), round(x2, 2))
      pval <- an$`Pr(>Chisq)`[2]
      sig <- ifelse(pval<0.1, ifelse(pval<0.05, ifelse(pval<0.001," *** "," * ") ," . "),"")
      x2.p <- ifelse(pval<0.001,"< .001", pv(pval))
      temp <- data.frame(term=terms[i], x2,x2.p,sig)
      compdf <- rbind(compdf,temp)
      compmodel <- compmodel2
    }
    compdf <- compdf[dim(compdf)[1]:1,] ## orders the comparisons to match lmer output
    rownames(compdf) <- 1:nrow(compdf)
    compdf$sig <- as.character(compdf$sig)
    firstrow <- c(NA, NA, ifelse(resdf$t.pz[1]<0.05," *** ",""))
    resdf <- cbind(resdf, rbind(firstrow, compdf[-c(1)])) ## bind model summary with comparison results
    # effR <- r.squaredGLMM(origmodel)
    # resdf$r2m <- effR[1]
    # resdf$r2c <- effR[2]
    resdf$version <- "anova"

  } else {

    compmodel <- refitML(origmodel)
    mod.bootmer <- bootMer(origmodel, fixefBoot, nsim=rsamp, seed=sn) ## bootstrap the original model
    resdf <- data.frame(tidy(mod.bootmer))
    names(resdf)[names(resdf) == 'statistic'] <- 'beta'
    t <- resdf$beta/resdf$std.error
    resdf$beta <- ifelse(abs(resdf$beta)<0.0001, scientific(resdf$beta), round(resdf$beta, 4))
    tpval <- round( (2*(1 - pnorm(abs(t)))), 3) ## calculate p based on t and z distribution (use for intercept only)
    resdf$std.error <- ifelse(abs(resdf$std.error)<0.0001, scientific(resdf$std.error), round(resdf$std.error, 4))
    resdf$t <- ifelse(abs(t)<0.01, scientific(t), round(t, 2))
    resdf$t.pz <- ifelse(tpval<0.001,"< .001", pv(tpval))
    for(i in 1:nrow(resdf)){ ## compute confidence intervals for the model betas
      bootInt <- boot.ci(mod.bootmer, type = "norm", index = i)
      lower <- as.numeric(bootInt$normal[2])
      beta.lower <- ifelse(abs(lower)<0.0001, scientific(lower), round(lower, 4))
      upper <- as.numeric(bootInt$normal[3])
      beta.upper <- ifelse(abs(upper)<0.0001, scientific(upper), round(upper, 4))
      resdf$beta[i] <- paste(resdf$beta[i], " [", beta.lower, ", ", beta.upper, "]", sep = "")
    }
    terms <- rev(attr(terms(compmodel),"term.labels")) ## extract predictor names
    for (i in 1:length(terms)){
      newformula <- paste(". ~ . - ",terms[i],"")
      compmodel2 <- update(compmodel, as.formula(newformula))
      an <- PBmodcomp(compmodel, compmodel2, nsim=rsamp, seed=sn) ## bootstrap likelihood ratio test
      x2 <- an$test[2,1]
      x2 <- ifelse(x2<0.01, scientific(x2), round(x2, 2))
      pval <- an$test[2,3]
      sig <- ifelse(pval<0.1, ifelse(pval<0.05, ifelse(pval<0.001," *** "," * ") ," . "),"")
      x2.p <- ifelse(pval<0.001,"< .001", pv(pval))
      temp <- data.frame(term=terms[i], x2,x2.p,sig)
      compdf <- rbind(compdf,temp)
      compmodel <- compmodel2
    }
    compdf <- compdf[dim(compdf)[1]:1,] ## orders the comparisons to match lmer output
    rownames(compdf) <- 1:nrow(compdf)
    compdf$sig <- as.character(compdf$sig)
    firstrow <- c("-", "-", ifelse(resdf$t.pz[1]<0.05," *** ",""))
    resdf <- cbind(resdf, rbind(firstrow, compdf[-c(1)])) ## bind model summary with comparison results
    # effR <- r.squaredGLMM(origmodel)
    # resdf$r2m <- effR[1]
    # resdf$r2c <- effR[2]
    resdf$version <- "bootstrap"
  }
  return(resdf)
}
