#' A Function for Running Sequential Mixed-Effects Model Comparisons
#'
#' This function allows you run mixed-effects model comparisons with bootstrapping switched on or off
#' @keywords bootstrapping, comparisons, mixed
#' @export
#' @examples 
#' lmercomp()

# Function for running bootstrapped model comparisons in Mixed-Effect Models
lmercomp <- function(origmodel, para_stats=data.frame(), comp_stats=data.frame(),
                     nsamp, sn, cores){
  
  # Sub-functions for extracting fixed effects and formatting p-values
  fixefBoot <- function(fit) {return(fixef(fit))}
  pv <- function(val) { paste( "=", sub("^(-?)0.", "\\1.", sprintf("%.3f", val))) }
  
  # Model fit statistics
  lmerfit <- function(mod){
    BIC <- BIC(mod)
    AICc <- AICc(mod)
    r2 <- r.squaredGLMM(mod)
    R2m <- r2[1, 1]
    R2c <- r2[1, 2]
    fit <- as_tibble(data.frame(BIC, AICc, R2m, R2c))
    return(fit)
  }
  
  compmodel <- refitML(origmodel)
  fit_stats <- lmerfit(compmodel)
  
  # Model parameter estimates
  cluster <- makeCluster(cores, type = 'SOCK')
  registerDoParallel(cluster)
  mod.bootmer <- bootMer(origmodel, FUN=fixef, nsim=nsamp, seed=sn,
                         cl=cluster, parallel='multicore', ncpus=cores)
  stopCluster(cluster)
  para_stats <- tidy(mod.bootmer, effects="fixed")[c("term", "statistic", "std.error")]
  names(para_stats)[names(para_stats) == 'statistic'] <- 'beta'
  para_stats$t <- para_stats$beta/para_stats$std.error
  para_stats$t.pz <- 2 * (1 - pnorm(abs(para_stats$t)))
  
  para_stats$beta.lower <- NA
  para_stats$beta.upper <- NA
  for(i in 1:nrow(para_stats)){
    bootInt <- boot.ci(mod.bootmer, type="norm", index = i)
    lower <- as.numeric(bootInt$normal[2])
    upper <- as.numeric(bootInt$normal[3])
    para_stats$beta[i] <- para_stats$beta[i]
    para_stats$beta.lower[i] <- lower
    para_stats$beta.upper[i] <- upper
  }
  
  # Model comparisons via sequential decomposition
  terms <- rev(attr(terms(compmodel),"term.labels"))
  for(i in 1:length(terms)){
    newformula <- paste(". ~ . - ",terms[i],"")
    compmodel2 <- update(compmodel, as.formula(newformula))
    cluster <- makeCluster(cores, type = 'SOCK')
    registerDoParallel(cluster)
    an <- PBmodcomp(compmodel, compmodel2, nsim=nsamp, seed=sn, cl=cluster)
    stopCluster(cluster)
    x2 <- an$test[2,1]
    x2.p <- an$test[2,3]
    temp <- data.frame(term=terms[i], x2, x2.p)
    comp_stats <- rbind(comp_stats, temp)
    compmodel <- compmodel2
  }
  comp_stats <- as_tibble(comp_stats[dim(comp_stats)[1]:1,])
  firstrow <- NA
  para_stats <- as_tibble(cbind(para_stats, rbind(firstrow, comp_stats[-c(1)])))
  para_stats <- para_stats[, c("term", "beta", "beta.lower", "beta.upper",
                               "t", "t.pz", "x2", "x2.p")]
  
  # Output all stats
  all_stats <- list(model_fit=fit_stats, model_parameters=para_stats)
  return(all_stats)
}