#' Gaussian Response fitting function with warm starts
#'
lspath <- function(x, y, e, df = 5,
                   # main.effect.names, interaction.names,
                   lambda.beta, lambda.gamma,
                   weights,
                   lambda.factor,
                   nlambda.gamma,
                   nlambda.beta,
                   nlambda,
                   threshold, max.iter,
                   initialization.type,
                   center, normalize, verbose,
                   cores) {

  rm(list = ls())
  devtools::load_all()
  options(scipen = 999, digits = 4)
  source("R/sim-data.R")
  x = X; y = Y; e = E
  lambda.beta = NULL ; lambda.gamma = NULL
  threshold = 1e-7 ; max.iter = 500 ; initialization.type = "ridge";
  nlambda.gamma = 10; nlambda.beta = 10;
  nlambda = 100 ; lambda.factor = 0.001
  cores = 1;
  df = 5;
  center=TRUE; normalize=TRUE; verbose = TRUE
# ==============================================================

  y <- drop(y)
  e <- drop(e)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])

  # Expand X's
  Phi_j <- do.call(cbind, lapply(seq_len(nvars), function(j) splines::bs(x[,j], df = df)))
  main_effect_names <- paste(paste0("X", rep(seq_len(nvars), each = df)), rep(seq_len(df), times = nvars), sep = "_")
  dimnames(Phi_j)[[2]] <- main_effect_names

  # X_E x Phi_j
  XE_Phi_j <- e * Phi_j
  interaction_names <- paste(main_effect_names, "X_E", sep = ":")
  dimnames(XE_Phi_j)[[2]] <- interaction_names

  design <- cbind(Phi_j, "X_E" = e, XE_Phi_j)
  design[1:5,1:10]
  colnames(design)
  obj <- standardize(x = design, y = y, center = center, normalize = normalize)
  x <- obj$x
  y <- obj$y
  bx <- obj$bx
  by <- obj$by
  sx <- obj$sx

  nulldev <- as.numeric(crossprod(y))

  if (is.null(lambda.gamma) & is.null(lambda.beta)) {

    # the sequence needs to have beta fixed first and then iterate over
    # lambda_gamma. Ive tried it the other oway arpund and the solutions are
    # too sparse
    lamb <- rev(lambda_sequence(x, y, nlambda = nlambda.beta,
                                lambda.factor = lambda.factor))
    lambda.beta <- rep(lamb, each = nlambda.beta)
    lambda.gamma <- rep(lamb, nlambda.beta)

    lambda_gamma_list <- lapply(seq_len(length(lambda.gamma)),
                                function(i) lambda.gamma[i])

    lambda_beta_list <- lapply(seq_len(length(unlist(lambda.beta))),
                               function(i) unlist(lambda.beta)[i])

  } else {

    # convert to a list. each element corresponds to a value of lambda_gamma
    # these are already of the proper length i.e., if the user specifies
    # lambda.beta and lambda.gamma then they this will not take all possible
    # combinations of lambda.beta and lambda.gamma. It will be the first element
    # of each as a pair, and so on. This is done on purpose for use with
    # the cv.shim function which uses the same lambda sequences for each fold...
    lambda_gamma_list <- lapply(seq_len(length(lambda.gamma)),
                                function(i) lambda.gamma[i])

    lambda_beta_list <- lapply(seq_len(length(unlist(lambda.beta))),
                               function(i) unlist(lambda.beta)[i])
  }

  tuning_params_mat <- matrix(c(lambda_gamma_list, lambda_beta_list),
                              nrow = 2, ncol = nlambda, byrow = T)
  dimnames(tuning_params_mat)[[1]] <- c("lambda.gamma","lambda.beta")
  dimnames(tuning_params_mat)[[2]] <- paste0("s",seq_len(nlambda))

  # determine which of the lambda's should use the previous warm start
  # for example, if there are 100 total lambdas = 10 x 10 grid, then
  # the 11th combination should not use previous values. it should re-start
  # using the initializations
  switchWarmStart <- data.frame(X1 = cbind(seq(from = 1 + nlambda.gamma, to = nlambda, by = nlambda.gamma)))
  # this matrix is in the correct order such that the result of one tuning
  # parameter should be the warm start for the next parameter, accounting
  # for the fact that the search is being done over a grid. This should be
  # called 2-D warm starts
  # lambdaNames <- paste0("s",as.vector(out))
  lambdaNames <- dimnames(tuning_params_mat)[[2]]

  adaptive.weights <- ridge_weights(x = x, y = y,
                                    main.effect.names = c(main_effect_names,"X_E"),
                                    interaction.names = interaction_names)

  # used in the warm start strategy because we ONLY want to use warm starts
  # for each lambda_gamma for a fixed lambda_beta. For the next lambda_beta
  # we restart using the initial values
  adaptive.weights.start <- adaptive.weights
  # initialization
  betas_and_alphas <- uni_fun(variables = colnames(x), x = x, y = y,
                              include.intercept = F,
                              type = initialization.type)

  # this converts the alphas to gammas
  uni_start <- convert(betas_and_alphas, main.effect.names = c(main_effect_names,"X_E"),
                       interaction.names = interaction_names)

  # used in the warm start strategy because we ONLY want to use warm starts
  # for each lambda_gamma for a fixed lambda_beta. For the next lambda_beta
  # we restart using the initial values
  uni_start_iteration1 <- uni_start

  # for all the x_tilde in zero_x_tilde, return the following matrix with 0 for each coefficient
  # this is like a place holder.
  coef_zero_gamma_matrix <- matrix(data = 0,
                                   nrow = length(interaction_names),
                                   ncol = 1,
                                   dimnames = list(interaction_names))

  # index data.frame to figure out which j < j'
  index <- data.frame(c(main_effect_names,"X_E"), seq_along(c(main_effect_names,"X_E")),
                      stringsAsFactors = F)
  colnames(index) <- c("main.effect.names","index")

  # matrix to store results of betas and alphas on standardized scale
  coefficientMat <- matrix(nrow = length(c(main_effect_names,"X_E",interaction_names)),
                           ncol = nlambda,
                           dimnames = list(c(main_effect_names,"X_E", interaction_names),
                                           lambdaNames))

  betaMat <- matrix(nrow = length(c(main_effect_names,"X_E")), ncol = nlambda,
                    dimnames = list(c(main_effect_names,"X_E"),
                                    lambdaNames))

  gammaMat <- matrix(nrow = length(interaction_names), ncol = nlambda,
                     dimnames = list(c(interaction_names),
                                     lambdaNames))

  outPrint <- matrix(NA, nrow = nlambda, ncol = 6,
                     dimnames = list(lambdaNames,
                                     c("dfBeta","dfAlpha","deviance",
                                       "percentDev",
                                       "lambdaBeta", "lambdaGamma")))

  # trying to implement this one at a time, so that
  pb <- progress::progress_bar$new(
    format = "  fitting over all pairs of tuning parameters [:bar] :percent eta: :eta",
    total = 100, clear = FALSE, width= 90)
  pb$tick(0)

  for (LAMBDA in lambdaNames) {
    # (LAMBDA <- lambdaNames[1])
    lambdaIndex <- which(LAMBDA==lambdaNames)
    lambda_beta <- tuning_params_mat["lambda.beta",LAMBDA][[1]]
    lambda_gamma <- tuning_params_mat["lambda.gamma",LAMBDA][[1]]

    if (verbose) {
      message(paste("Index:",LAMBDA, ", lambda_beta:", lambda_beta, ", lambda_gamma:", lambda_gamma))
    }

    beta_hat_previous <- uni_start[c(main_effect_names,"X_E"), , drop = F]

    gamma_hat_previous <- uni_start[interaction_names, , drop = F]

    # store likelihood values at each iteration in a matrix Q
    # rows: iteration number
    Q <- matrix(nrow = max.iter + 1, ncol = 1)

    # store the value of the likelihood at the 0th iteration
    Q[1,1] <- Q_theta(x = x, y = y,
                      beta = beta_hat_previous,
                      gamma = gamma_hat_previous,
                      lambda.beta = lambda_beta,
                      lambda.gamma = lambda_gamma,
                      weights = adaptive.weights,
                      main.effect.names = c(main_effect_names,"X_E"),
                      interaction.names = interaction_names)

    m <- 1 # iteration counter
    delta <- 1 # threshold initialization
    # to see which lambdas have converged: 0=not converged, 1=converged
    converged <- 0

    while (threshold < delta && m < max.iter){

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update gamma (interaction parameter)
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      # this is a nsubjects x lambda matrix for each tuning parameter stored in a list
      # each element of the list corresponds to a tuning parameter
      # need to keep y_tilde_list and x_tilde_list of length nlambda

      y_tilde <- y - x[,c(main_effect_names,"X_E"), drop = F] %*% beta_hat_previous

      x_tilde <- xtilde(interaction.names = interaction_names,
                        data.main.effects = x[,c(main_effect_names,"X_E"), drop = F],
                        beta.main.effects = beta_hat_previous)

      # indices of the x_tilde matrices that have all 0 columns
      zero_x_tilde <- is.null(colnames(check_col_0(x_tilde)))

      # this will store the results but will be shorter than nlambda
      gamma_hat_next <- if (zero_x_tilde) coef_zero_gamma_matrix else
        as.matrix(coef(glmnet::glmnet(
          x = x_tilde,
          y = y_tilde,
          penalty.factor = adaptive.weights[interaction.names,,drop=F],
          lambda = lambda_gamma,
          standardize = F, intercept = F))[-1,,drop = F])

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update beta (main effect parameter) step 4 of algortihm in Choi et al
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      beta_hat_next <- beta_hat_previous

      for (j in main.effect.names) {

        #j="x1"
        # determine the main effects not in j
        j_prime_not_in_j <- setdiff(main.effect.names,j)

        # for (notconverged in not_converged) {
        y_tilde_2 <- y -
          x[,j_prime_not_in_j, drop = F] %*%
          beta_hat_next[j_prime_not_in_j, , drop = F] -
          rowSums(
            xtilde_mod(beta.main.effects = beta_hat_next[j_prime_not_in_j, , drop = F],
                       gamma.interaction.effects = gamma_hat_next,
                       interaction.names = interaction.names[-grep(paste0("\\b",j,"\\b"), interaction.names)],
                       data.main.effects = x[,j_prime_not_in_j, drop = F]))

        # j' less than j
        j.prime.less <- index[which(index[,"index"] < index[which(index$main.effect.names == j),2]),
                              "main.effect.names"]

        # need to make sure paste(j.prime.less,j,sep=":") are variables in x matrix
        # this is to get around situations where there is only interactions with the E variable
        j.prime.less.interaction <- intersect(paste(j.prime.less,j, sep = ":"), colnames(x))

        # need to get the main effects that are in j.prime.greater.interaction
        j.prime.less <- gsub("\\:(.*)", "", j.prime.less.interaction)


        # the if conditions in term1 and term2 are to check if there are
        # any variables greater or less than j
        # lapply is faster than mclapply here
        term_1 <- if (length(j.prime.less.interaction) != 0 ) {
          x[,j.prime.less.interaction] %*%
            (gamma_hat_next[j.prime.less.interaction,, drop = F] *
               beta_hat_next[j.prime.less, , drop = F])} else 0

        # term_1_list <- replace(term_1_list, not_converged, term_1_list_not_converged)

        # j' greater than j
        j.prime.greater <- index[which(index[,"index"] >
                                         index[which(index$main.effect.names == j),2]),
                                 "main.effect.names"]

        # need to make sure j.prime.greater is a variable in x matrix
        # this is to get around situations where there is only interactions with the E variable
        j.prime.greater.interaction <- intersect(paste(j,j.prime.greater,sep = ":"), colnames(x))

        # need to get the main effects that are in j.prime.greater.interaction
        j.prime.greater <- if (all(gsub("\\:(.*)", "", j.prime.greater.interaction) == j))
          gsub("(.*)\\:", "", j.prime.greater.interaction) else gsub("\\:(.*)", "", j.prime.greater.interaction)

        term_2 <- if (length(j.prime.greater) != 0) {
          x[,j.prime.greater.interaction] %*%
            (gamma_hat_next[j.prime.greater.interaction,, drop = F] *
               beta_hat_next[j.prime.greater,,drop = F]) } else 0

        x_tilde_2 <- x[,j, drop = F] +  term_1 + term_2

        # glmnet is giving weired results for this... and is slower than using my
        # soft function. use this. non-parallel version is faster
        # the result of this should give 1 beta for each tuningn parameter
        # This calculates for all tuning parameters
        beta_hat_next_j <- soft(x = x_tilde_2,
                                y = y_tilde_2,
                                weight = adaptive.weights[j,,drop=F],
                                lambda = lambda_beta)

        beta_hat_next[j,] <- beta_hat_next_j

      }

      Q[m + 1, 1] <- Q_theta(beta = beta_hat_next,
                             gamma = gamma_hat_next,
                             lambda.beta = lambda_beta,
                             lambda.gamma = lambda_gamma,
                             weights = adaptive.weights,
                             x = x, y = y,
                             main.effect.names = main.effect.names,
                             interaction.names = interaction.names)

      delta <- abs(Q[m,1] - Q[m + 1, 1])/abs(Q[m,1])
      converged <- as.numeric(delta <= threshold)
      if (verbose) print(paste("Iteration:", m))
      if (verbose) print(converged)

      m <- m + 1

      beta_hat_previous <- beta_hat_next

      # adaptive weight for each tuning parameter. currently this is the
      # same for iterations, but I am coding it here
      # for flexibility in case we want to change the weights at each iteration

      adaptive.weights <- update_weights(betas = beta_hat_previous,
                                         gammas = gamma_hat_next,
                                         main.effect.names = main.effect.names,
                                         interaction.names = interaction.names)

    }

    Betas_and_Alphas <- convert2(beta_hat_next, gamma_hat_next,
                                 main.effect.names = main.effect.names,
                                 interaction.names = interaction.names)

    dfbeta <- length(nonzero(Betas_and_Alphas[main.effect.names,]))
    dfalpha <- length(nonzero(Betas_and_Alphas[interaction.names,]))
    deviance <- crossprod(y - x %*% Betas_and_Alphas)
    devRatio <- 1-deviance/nulldev
    outPrint[LAMBDA,] <- c(if (dfbeta==0) 0 else dfbeta,
                           if (dfalpha==0) 0 else dfalpha,
                           deviance,
                           devRatio,
                           lambda_beta, lambda_gamma)

    betaMat[,lambdaIndex] <- beta_hat_next
    gammaMat[,lambdaIndex] <- gamma_hat_next

    # for a fixed lambda.beta, keep using the previous solution for each lambda.gamma
    # for the next fixed lambda.beta restart the calculation using the initial value
    # from the uni_fun function
    betaWarmStart <- if (lambdaIndex %ni% switchWarmStart$X1) beta_hat_next

    gammaWarmStart <- if (lambdaIndex %ni% switchWarmStart$X1) gamma_hat_next

    # uni_start <- rbind(beta_hat_next, gamma_hat_next)
    uni_start <- if (lambdaIndex %ni% switchWarmStart$X1) {
      rbind(betaWarmStart, gammaWarmStart) } else {
        uni_start_iteration1
      }

    # need to update weights also!
    adaptive.weights <- if (lambdaIndex %ni% switchWarmStart$X1) {
      update_weights(betaWarmStart, gammaWarmStart,
                     main.effect.names, interaction.names) } else {
                       adaptive.weights.start
                     }

    devianceDiff <- outPrint[lambdaIndex,"deviance"] - outPrint[lambdaIndex-1,"deviance"]
    #coefficientMat[,LAMBDA] <- Betas_and_Alphas

    if (length(devianceDiff)!=0 && !is.na(devianceDiff) && devRatio>1e-3) {
      # if (devianceDiff < 1e-50 | outPrint[LAMBDA,"percentDev"] > 0.999) break }
      if (outPrint[LAMBDA,"percentDev"] > 0.999) break }

    pb$tick()
  }

  # Rprof()
  # proftable(tmp)
  # summaryRprof(tmp)
  # b <- Sys.time()

  # outPrint[complete.cases(outPrint),]
  # coefficientMat
  # b-a

  beta_hat_next_list <- lapply(seq_len(ncol(betaMat)),
                               function(i) betaMat[,i,drop=F])
  # convert to original scale
  betas_original_scale_list <- lapply(beta_hat_next_list, function(i) i / sx[main.effect.names])

  gamma_hat_next_list <- lapply(seq_len(ncol(gammaMat)),
                                function(i) gammaMat[,i,drop=F])

  gammas_original_scale_list <- lapply(gamma_hat_next_list, function(i) i / sx[interaction.names])

  # convert gammas to alphas
  betas_alphas_original_scale <- mapply(convert2,
                                        beta = betas_original_scale_list,
                                        gamma = gammas_original_scale_list,
                                        MoreArgs = list(main.effect.names = main.effect.names,
                                                        interaction.names = interaction.names))

  dimnames(betas_alphas_original_scale) <- list(c(main.effect.names, interaction.names),
                                                paste0("s",1:nlambda))

  betas_original_scale <- betas_alphas_original_scale[main.effect.names, , drop = F]
  alphas_original_scale <- betas_alphas_original_scale[interaction.names, , drop = F]

  b0 <- vector(length = nlambda)
  for (lam in seq_len(nlambda)) {
    b0[lam] <- by - sum(betas_original_scale[,lam,drop = F] * bx[main.effect.names]) -
      sum(alphas_original_scale[,lam,drop=F] * bx[interaction.names])
  }
  names(b0) <- paste0("s",1:nlambda)

  gamma_final <- as(matrix(unlist(gammas_original_scale_list, use.names = F),
                           ncol = nlambda,
                           byrow = T,
                           dimnames = list(interaction.names, paste0("s",1:nlambda))),
                    "dgCMatrix")
  beta_final <- as(betas_original_scale,"dgCMatrix")
  alpha_final <- as(alphas_original_scale,"dgCMatrix")

  lambda.beta <- unlist(lambda_beta_list)
  names(lambda.beta) <- paste0("s",1:nlambda,".beta")
  lambda.gamma <- unlist(lambda_gamma_list)
  names(lambda.gamma) <- paste0("s",1:nlambda, ".gamma")

  out <- list(b0 = b0,
              beta = beta_final,
              alpha = alpha_final,
              gamma = gamma_final,
              lambda.beta = lambda.beta,
              lambda.gamma = lambda.gamma,
              tuning.parameters = tuning_params_mat,
              dfbeta = outPrint[,"dfBeta", drop = F],
              dfalpha = outPrint[,"dfAlpha", drop = F],
              dev.ratio = outPrint[,"percentDev", drop = F],
              deviance = outPrint[,"deviance", drop = F],
              converged = converged, x = x, y = y, bx = bx, by = by, sx = sx,
              center = center, normalize = normalize,
              nlambda.gamma = nlambda.gamma,
              nlambda.beta = nlambda.beta,
              nlambda = nlambda,
              interaction.names = interaction.names,
              main.effect.names = main.effect.names)
  class(out) <- "lspath"
  return(out)

}




