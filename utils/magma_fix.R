.fit_otu_table <- function (otu_table, distrib, sequencing_depth="GMPR", Xv=NULL)
{
  #libsize
  if (!length(sequencing_depth) %in% c(1,dim(otu_table)[1]) ) {
    stop("sequencing_depth has an incorrect format")
  } else if (length(sequencing_depth) == dim(otu_table)[1] ) {
    libsize <- sequencing_depth
  } else if (sequencing_depth == "GMPR") {
    libsize <- as.numeric(GMPR(otu_table,intersect.no = 10))
    if (any(is.na(libsize))){libsize <- as.numeric(GMPR(otu_table,intersect.no = 5))}
    if (any(is.na(libsize))){
      libsize <- apply(otu_table,1,function(x){exp(mean(log(x+1)))})
      cat("Due to the presence of one or more atypical samples mentioned above, the size of the library can not be estimated with the GMPR procedure.\n")
      cat("The estimate is made with the geometric mean of the OTU reads of each samples, with the addition of a pseudocount.\n")
    }
  } else if (sequencing_depth == "TS") {
    libsize <- rowSums(otu_table)
  } else if (sequencing_depth == "GM.pseudocount") {
    libsize <- apply(otu_table,1,function(x){exp(mean(log(x+1)))})
  } else if (sequencing_depth == "unif") {
    libsize <- rep(1,dim(otu_table)[1])
  } else { stop("sequencing_depth has an incorrect format") }
  
  
  # no covariate
  if (is.null(Xv)){
    
    .fit_distrib <- switch(distrib, P=.fit_P, ZIP=.fit_ZIP, NB=.fit_NB, ZINB=.fit_ZINB)
    if (is.null(.fit_distrib)) stop ("distrib must be", "\"P\",", "\"NB\",", "\"ZIP\",", "\"ZINB\"")
    
    suppressWarnings(apply(X=otu_table, MARGIN=2, FUN = .fit_distrib, seqdep=libsize))
    
    
    # with covariate
  } else {
    
    if(!is.null(colnames(Xv))){
      if (colnames(Xv)[1] != "(Intercept)") {
        if (is.character(Xv)) {Xv <- factor(Xv)} else if (!is.numeric) {Xv <- droplevels(Xv)}
        Xv <- try(stats::model.matrix(~.,data=data.frame(Xv)), silent=T)}
    } else {
      if (is.character(Xv)) {Xv <- factor(Xv)} else if (!is.numeric(Xv)) {Xv <- droplevels(Xv)}
      Xv <- try(stats::model.matrix(~.,data=data.frame(Xv)), silent=T)
    }
    #error messages
    if (methods::is(Xv,"try-error")) stop("X has an incorrect format")
    if (dim(Xv)[1] != dim(otu_table)[1]) stop("X and otu_table do not have the same number of rows")
    
    
    .fit_distrib <- switch(distrib, P=.fit_P_covariate, ZIP=.fit_ZIP_covariate, NB=.fit_NB_covariate, ZINB=.fit_ZINB_covariate)
    if (is.null(.fit_distrib)) stop ("distrib must be", "\"P\",", "\"NB\",", "\"ZIP\",", "\"ZINB\"")
    
    suppressWarnings(apply(X=otu_table, MARGIN=2, FUN = .fit_distrib, seqdep=libsize, Xv=Xv))
  }
}


.ll_NB <- function(par, donnees) {
  Yobs   <- donnees[,1]
  seqdep <- donnees[,2]
  
  # parameters
  pxV   <- exp (par)
  mu    <- pxV[1] * seqdep
  theta <- pxV[2]
  
  # Maximum-likelihood
  MLE <- log ( stats::dnbinom(Yobs, size=theta, mu=mu) )
  return( - sum(MLE) )
}


.fit_NB <- function(data, seqdep)
{
  nouveau_par_init = c(mean(log(data[is.finite(data)])), 0)
  par_init0  <- c(0,0)
  mu_init    <- mean(data/seqdep)
  theta_init <- max(mu_init^2 / (stats::var(data/seqdep)-mu_init), 10^(-10))
  par_init   <- c(log(mu_init),log(theta_init))
  print(par_init)
  print("magma fix")
  res_optim <- try(stats::optim(fn=.ll_NB, par=par_init, donnees=cbind(data,seqdep), hessian = FALSE ), silent=T)
  
  if (methods::is(res_optim,"try-error")) {
    mu_init    <- mean(data/seqdep)
    theta_init <- max(mu_init^2 / (stats::var(data/seqdep)-mu_init), 10^(-10))
    par_init   <- c(log(mu_init),log(theta_init))
    res_optim <- stats::optim(fn=.ll_NB, par=par_init, donnees=cbind(data,seqdep), hessian = FALSE)
  }
  
  mu        <- exp(res_optim$par[1]) * seqdep
  theta     <- exp(res_optim$par[2])
  res <- list(mu=mu, theta=theta)
  return (res)
}


.fit_ZINB <- function(data,seqdep)
{
  if (sum(data==0) > 0 ) {
    ctrl        <- pscl::zeroinfl.control(method = "L-BFGS-B")
    ctrl$reltol <- NULL
    ctrl$factr  <- 1e-3/.Machine$double.eps
    fit  <- try(pscl::zeroinfl(data~1+offset(log(seqdep))|1, dist = "negbin", link = "logit" , control=ctrl), silent=T)
    
    if(methods::is(fit,"try-error")) {
      fit2 <- .fit_NB(data,seqdep)
      mu    <- fit2$mu
      theta <- fit2$theta
      p0    <- 0
    }else{
      mu    <- exp(fit$coefficients$count) * seqdep
      theta <- fit$theta
      p0    <- .logit(fit$coefficients$zero , inverse =T)
    }
  } else {
    fit2  <- .fit_NB(data,seqdep)
    mu    <- fit2$mu
    theta <- fit2$theta
    p0    <- 0
  }
  
  return (list(mu=mu, theta=theta, p0=p0))
}


# from VGAM package
.logit <- function (theta, bvalue = NULL, inverse = FALSE, deriv = 0, short = TRUE,
                    tag = FALSE)
{
  if (is.character(theta)) {
    string <- if (short)
      paste("logit(", theta, ")", sep = "")
    else paste("log(", .as.char.expression(theta), "/(1-",
               .as.char.expression(theta), "))", sep = "")
    if (tag)
      string <- paste("Logit:", string)
    return(string)
  }
  if (!inverse && length(bvalue)) {
    theta[theta <= 0] <- bvalue
    theta[theta >= 1] <- 1 - bvalue
  }
  if (inverse) {
    switch(as.character(deriv), `0` = stats::plogis(theta), `1` = 1/Recall(theta = theta,
                                                                           bvalue = bvalue, inverse = FALSE, deriv = deriv),
           `2` = theta * (1 - theta) * (1 - 2 * theta), `3` = (1 -
                                                                 6 * theta * (1 - theta)) * theta * (1 - theta),
           stop("argument 'deriv' unmatched"))
  }
  else {
    switch(as.character(deriv), `0` = stats::qlogis(theta), `1` = 1/(theta *
                                                                       (1 - theta)), `2` = (2 * theta - 1)/(theta * (1 -
                                                                                                                       theta))^2, `3` = 2 * (1 - 3 * theta * (1 - theta))/(theta *
                                                                                                                                                                             (1 - theta))^3, stop("argument 'deriv' unmatched"))
  }
}

magma.norm <- function (data, distrib="ZINB", seq_depth="GMPR", X=NULL)
{
  if (distrib=="empirical") {
    if (is.null(dim(data))) {
      F_emp <- stats::ecdf(as.numeric(data))
      data <- stats::qnorm((F_emp(data-1)+F_emp(data))/2) #median approximation
    } else {
      for (j in 1:dim(data)[2]){
        F_emp <- stats::ecdf(data[,j])
        data[,j] <- stats::qnorm( (F_emp(data[,j]-1)+F_emp(data[,j]))/2 )
      }
    }
    res <- data
  } else {
    F_parameters <- .fit_otu_table(otu_table=data, distrib=distrib, sequencing_depth=seq_depth, Xv=X)
    data <- .copula_norm_from_param(data=data, distrib=distrib, param=F_parameters)
    res = list(z=data, param = F_parameters)
  }
  return(res)
}


magma <- function (data, distrib="ZINB", X=NULL, data.sup=NA, method="glasso",  criterion.select="stars", lambda=NULL, seq_depth="GMPR", magma.select=TRUE, verbose=TRUE, ...)
{
  if (verbose) {
    cat("Conducting copula median transformations....")
    utils::flush.console()
  }
  
  res_norm <- magma.norm(data=data, distrib=distrib, seq_depth=seq_depth, X=X)
  z <- res_norm$z
  Fparam = res_norm$param
  
  if (verbose) {
    cat("done.\n")
    utils::flush.console()
  }
  if (!any(is.na(data.sup))){
    z <- cbind(z, magma.norm(data=data.sup, distrib="empirical"))
  }
  
  args <- list(...)
  args.huge <- list(lambda = lambda,
                    nlambda = args$nlambda,
                    lambda.min.ratio = args$lambda.min.ratio,
                    scr = args$scr,
                    scr.num = args$scr.num,
                    cov.output = args$cov.output,
                    sym = args$sym)
  args.huge.select <- list(ebic.gamma = args$ebic.gamma,
                           stars.thresh = args$stars.thresh,
                           stars.subsample.ratio = args$stars.subsample.ratio,
                           rep.num = args$rep.num)
  args.huge[sapply(args.huge, is.null)] <- NULL
  args.huge.select[sapply(args.huge.select, is.null)] <- NULL
  
  magma.res <- do.call(huge::huge, c(args.huge, list(x = z, method = method, verbose = verbose)))
  magma.res$rawdata <- data
  magma.res$param <- Fparam
  
  if (!any(is.na(data.sup))){magma.res$rawdata <- cbind(data,data.sup)}
  
  if (length(lambda)!=1){
    if (isTRUE(magma.select)){
      magma.res <- do.call(huge::huge.select, c(args.huge.select,list(est = magma.res, criterion = criterion.select, verbose = verbose)))
    } else { if (verbose) {cat("no model selection is made, plot.magma can not be used\n")}}
  } else {
    magma.res$opt.icov    <- magma.res$icov[[1]]
    magma.res$refit       <- (as.matrix(magma.res$icov[[1]])!=0)*1
    diag(magma.res$refit) <- 0
  }
  
  class(magma.res) <- "magma"
  magma.res
}

.copula_norm_from_param <- function (data, distrib, param)
{
  copula_norm_j <- switch(distrib, P=.copula_norm_j_P, ZIP=.copula_norm_j_ZIP, NB=.copula_norm_j_NB, ZINB=.copula_norm_j_ZINB)
  if (is.null(copula_norm_j)) stop ("distrib must be", "\"P\",", "\"NB\",", "\"ZIP\",", "\"ZINB\"", "\"binomial\",")
  
  for (j in 1:dim(data)[2]){
    data[,j] <- copula_norm_j (data[,j], param[[j]])
  }
  .fixInf(data)
}



.copula_norm_j_P <- function (data_j, param_j)
{
  unif <- stats::ppois(data_j -1, param_j$mu) + stats::ppois(data_j, param_j$mu)
  stats::qnorm(unif/2)
}

.copula_norm_j_ZIP <- function (data_j, param_j)
{
  unif <- .pzipois(data_j-1, lambda = param_j$mu, pstr0 = param_j$p0) + .pzipois(data_j, lambda = param_j$mu, pstr0 = param_j$p0)
  stats::qnorm(unif/2)
}

.copula_norm_j_NB <- function (data_j, param_j)
{
  unif <- stats::pnbinom(data_j-1, mu=param_j$mu, size=param_j$theta) + stats::pnbinom(data_j, mu=param_j$mu, size=param_j$theta)
  stats::qnorm(unif/2)
}

.copula_norm_j_ZINB <- function (data_j, param_j)
{
  unif <- .pzinegbin(data_j-1, munb = param_j$mu, size = param_j$theta, pstr0 = param_j$p0) + .pzinegbin(data_j, munb = param_j$mu, size = param_j$theta, pstr0 = param_j$p0)
  stats::qnorm(unif/2)
}



# from SpiecEasi
.fixInf <- function(data)
{
  # hacky way of replacing infinite values with the col max + 1
  if (any(is.infinite(as.matrix(data)))) {
    data <-  apply(data, 2, function(x) {
      if (any(is.infinite(x))) {
        x[ind<-which(is.infinite(x))] <- NA
        x[ind] <- max(x, na.rm=TRUE)+1
      }
      x
    })
  }
  as.matrix(data)
}



#from VGAM package
.pzipois <- function (q, lambda, pstr0 = 0)
{
  LLL <- max(length(pstr0), length(lambda), length(q))
  if (length(pstr0) != LLL)
    pstr0 <- rep_len(pstr0, LLL)
  if (length(lambda) != LLL)
    lambda <- rep_len(lambda, LLL)
  if (length(q) != LLL)
    q <- rep_len(q, LLL)
  ans <- stats::ppois(q, lambda)
  ans <- ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)
  deflat.limit <- -1/expm1(lambda)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN
  ans
}

#from VGAM package
.pzinegbin <- function (q, size, prob = NULL, munb = NULL, pstr0 = 0)
{
  if (length(munb)) {
    if (length(prob))
      stop("arguments 'prob' and 'munb' both specified")
    prob <- size/(size + munb)
  }
  LLL <- max(length(pstr0), length(size), length(prob), length(q))
  if (length(q) != LLL)
    q <- rep_len(q, LLL)
  if (length(size) != LLL)
    size <- rep_len(size, LLL)
  if (length(prob) != LLL)
    prob <- rep_len(prob, LLL)
  if (length(pstr0) != LLL)
    pstr0 <- rep_len(pstr0, LLL)
  ans <- stats::pnbinom(q = q, size = size, prob = prob)
  ans <- ifelse(q < 0, 0, pstr0 + (1 - pstr0) * ans)
  prob0 <- prob^size
  deflat.limit <- -prob0/(1 - prob0)
  ans[pstr0 < deflat.limit] <- NaN
  ans[pstr0 > 1] <- NaN
  ans
}
