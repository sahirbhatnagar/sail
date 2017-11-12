#' Gaussian Response fitting function with warm starts
#'
lspath <- function(x, y, e, df,
                   group.penalty,
                   lambda.beta,
                   lambda.gamma,
                   weights,
                   lambda.factor,
                   nlambda.gamma,
                   nlambda.beta,
                   nlambda,
                   thresh,
                   maxit,
                   initialization.type,
                   center,
                   normalize,
                   verbose,
                   cores) {

  # rm(list = ls())
  # devtools::load_all()
  # options(scipen = 999, digits = 4)
  # DT <- gendata(n = 300, p = 50, df = 5)
  # x = DT$x; y = DT$y; e = DT$e
  # lambda.beta = NULL ; lambda.gamma = NULL
  # thresh = 1e-7 ; maxit = 500 ; initialization.type = "ridge";
  # nlambda.gamma = 10; nlambda.beta = 10;
  # nlambda = 100 ; lambda.factor = 0.001
  # cores = 1;
  # df = DT$df;
  # center=TRUE; normalize=TRUE; verbose = TRUE
# ==============================================================
# print(lambda.factor)
  # browser()
  y <- drop(y)
  e <- drop(e)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])

  # group membership
  group <- rep(seq_len(nvars), each = df)

  # Expand X's
  Phi_j <- do.call(cbind, lapply(seq_len(nvars), function(j) splines::bs(x[,j], df = df)))
  # Phi_j <- do.call(cbind, lapply(seq_len(nvars), function(j) splines::ns(x[,j], df = df)))

  main_effect_names <- paste(paste0("X", group), rep(seq_len(df), times = nvars), sep = "_")
  dimnames(Phi_j)[[2]] <- main_effect_names

  # the names of the main effects in list form.. Each element contains the main effect names
  # for an X
  list_group_main <- split(main_effect_names, group)

  # X_E x Phi_j
  XE_Phi_j <- e * Phi_j
  interaction_names <- paste(main_effect_names, "X_E", sep = ":")
  dimnames(XE_Phi_j)[[2]] <- interaction_names

  list_group_inter <- split(interaction_names, group)

  design <- cbind(Phi_j, "X_E" = e, XE_Phi_j)
  # design[1:5,1:10]
  # colnames(design)
  obj <- standardize(x = design, y = y, center = center, normalize = normalize)
  x <- obj$x
  y <- obj$y
  bx <- obj$bx
  by <- obj$by
  sx <- obj$sx

  nulldev <- as.numeric(crossprod(y))

  adaptive.weights <- ridge_weights(x = x,
                                    y = y,
                                    group = group,
                                    main.effect.names.list = list_group_main,
                                    interaction.names.list = list_group_inter,
                                    include.intercept = F)

  # used in the warm start strategy because we ONLY want to use warm starts
  # for each lambda_gamma for a fixed lambda_beta. For the next lambda_beta
  # we restart using the initial values
  adaptive.weights.start <- adaptive.weights
  # initialization
  betas_and_alphas <- uni_fun(x = x, y = y,
                              include.intercept = F,
                              type = initialization.type)
# browser()
  if (is.null(lambda.gamma) & is.null(lambda.beta)) {

    # the sequence needs to have beta fixed first and then iterate over
    # lambda_gamma. Ive tried it the other oway arpund and the solutions are
    # too sparse
    lamb <- if (nlambda.beta == 1) rev(lambda_sequence(x, y, nlambda = nlambda.gamma,
                                lambda.factor = lambda.factor)) else rev(lambda_sequence(x, y, nlambda = nlambda.beta,
                                                                                         lambda.factor = lambda.factor))
    lambda.beta <- rep(lamb, each = nlambda.beta)
    lambda.gamma <- rep(lamb, nlambda.beta)

    # lambda.gamma <- rep(lamb, each = nlambda.beta)
    # lambda.beta <- rep(lamb, nlambda.beta)

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
    # the cv.sail function which uses the same lambda sequences for each fold...
    lambda_gamma_list <- lapply(seq_len(length(lambda.gamma)),
                                function(i) lambda.gamma[i])

    lambda_beta_list <- lapply(seq_len(length(unlist(lambda.beta))),
                               function(i) unlist(lambda.beta)[i])
  }

  tuning_params_mat <- matrix(c(lambda_gamma_list, lambda_beta_list),
                              nrow = 2, ncol = nlambda, byrow = T)
  dimnames(tuning_params_mat)[[1]] <- c("lambda.gamma","lambda.beta")
  dimnames(tuning_params_mat)[[2]] <- paste0("s",seq_len(nlambda))


  print(tuning_params_mat)

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

  # this converts the alphas to gammas
  # uni_start <- convert(betas.and.alphas = betas_and_alphas,
  #                      # main.effect.names = c(main_effect_names,"X_E"),
  #                      # interaction.names = interaction_names,
  #                      group = group,
  #                      main.effect.names.list = list_group_main,
  #                      interaction.names.list = list_group_inter)

  uni_start <- matrix(c(betas_and_alphas[c(main_effect_names,"X_E"),], rep(0, length(unique(group)))))

  # used in the warm start strategy because we ONLY want to use warm starts
  # for each lambda_gamma for a fixed lambda_beta. For the next lambda_beta
  # we restart using the initial values
  uni_start_iteration1 <- uni_start

  # for all the x_tilde in zero_x_tilde, return the following matrix with 0 for each coefficient
  # this is like a place holder.
  coef_zero_gamma_matrix <- matrix(data = 0,
                                   nrow = length(unique(group)),
                                   ncol = 1,
                                   dimnames = list(paste0("X", unique(group))))

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
  # pb <- progress::progress_bar$new(
  #   format = "  fitting over all pairs of tuning parameters [:bar] :percent eta: :eta",
  #   total = 100, clear = FALSE, width= 90)
  # pb$tick(0)

  for (LAMBDA in lambdaNames) {
    # (LAMBDA <- lambdaNames[1])
    lambdaIndex <- which(LAMBDA==lambdaNames)
    lambda_beta <- tuning_params_mat["lambda.beta",LAMBDA][[1]]
    lambda_gamma <- tuning_params_mat["lambda.gamma",LAMBDA][[1]]

    if (verbose) {
      message(paste("Index:",LAMBDA, ", lambda_beta:", lambda_beta, ", lambda_gamma:", lambda_gamma))
    }

    beta_hat_previous <- betas_and_alphas[c(main_effect_names,"X_E"), , drop = F]
    gamma_hat_previous <- as.matrix(rep(0, length(unique(group))))
    dimnames(gamma_hat_previous)[[1]] <- paste0("X",unique(group))

    # gamma_hat_previous <- uni_start[interaction_names, , drop = F]

    # store likelihood values at each iteration in a matrix Q
    # rows: iteration number
    Q <- matrix(nrow = maxit + 1, ncol = 1)

    # store the value of the likelihood at the 0th iteration
    Q[1,1] <- Q_theta(x = x, y = y,
                      beta = beta_hat_previous,
                      gamma = gamma_hat_previous,
                      lambda.beta = lambda_beta,
                      lambda.gamma = lambda_gamma,
                      weights = adaptive.weights,
                      main.effect.names = list_group_main,
                      interaction.names = list_group_inter,
                      group = group)
    message(Q[1,1])

    #iteration counter
    m <- 1

    # to enter while loop
    converged <- FALSE

    while (!converged && m < maxit){

      Theta_init <- c(drop(beta_hat_previous), drop(gamma_hat_previous))

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update gamma (interaction parameter)
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      # this is a nsubjects x lambda matrix for each tuning parameter stored in a list
      # each element of the list corresponds to a tuning parameter
      # need to keep y_tilde_list and x_tilde_list of length nlambda

      # all(c(main_effect_names,"X_E") == rownames(beta_hat_previous))
      y_tilde <- y - x[,c(main_effect_names,"X_E"), drop = F] %*% beta_hat_previous

      x_tilde <- xtilde(x = x,
                        group = group,
                        main.effect.names.list = list_group_main,
                        interaction.names.list = list_group_inter,
                        beta.main.effects = beta_hat_previous)

      # indices of the x_tilde matrices that have all 0 columns
      zero_x_tilde <- is.null(colnames(check_col_0(x_tilde)))

      # this will store the results but will be shorter than nlambda
      gamma_hat_next <- if (zero_x_tilde) coef_zero_gamma_matrix else
        as.matrix(coef(glmnet::glmnet(
          x = x_tilde,
          y = y_tilde,
          penalty.factor = adaptive.weights[paste0("X", unique(group), ":X_E"),,drop=F],
          lambda = c(.Machine$double.xmax,lambda_gamma),
          standardize = F, intercept = F))[-1,2,drop = F])

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update beta (main effect parameter) step 4 of algortihm in Choi et al
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      # this is actually theta_j and also includes Beta_E
      beta_hat_next <- beta_hat_previous

      # update beta_j, j = 1,..., p first
      for (j in unique(group)) {

        # j = 2
        # ==================

        # determine the main effects not in j
        j_prime_not_in_j <- as.vector(do.call(c, list_group_main[setdiff(unique(group),j)]))

        y_tilde_2 <- y -
          beta_hat_next["X_E",] * x[,"X_E"] -
          x[,j_prime_not_in_j, drop = F] %*% beta_hat_next[j_prime_not_in_j, , drop = F] -
          rowSums(do.call(cbind, lapply(setdiff(unique(group),j), function(kk) {
            beta_hat_next["X_E", ] *
              gamma_hat_next[kk, ] *
              (x[,list_group_inter[[kk]], drop = F] %*% beta_hat_next[list_group_main[[kk]], , drop = F])
          })))

        x_tilde_2 <- x[,list_group_main[[j]], drop = F] + gamma_hat_next[j,] * beta_hat_next["X_E",] * x[,list_group_inter[[j]]]

        # browser()

        beta_hat_next_j <- switch(group.penalty,
                                  gglasso = coef(gglasso::gglasso(x = x_tilde_2,
                                                                  y = y_tilde_2,
                                                                  group = rep(1, df),
                                                                  pf = as.vector(adaptive.weights[paste0("X",j),]),
                                                                  lambda = lambda_beta,
                                                                  intercept = F))[-1,],
                                  MCP = grpreg::grpreg(X = x_tilde_2,
                                                       y = y_tilde_2,
                                                       group = rep(1, df),
                                                       penalty = "grMCP",
                                                       family = "gaussian",
                                                       group.multiplier = as.vector(adaptive.weights[paste0("X",j),]),
                                                       lambda = lambda_beta,
                                                       intercept = T)$beta[-1,],
                                  SCAD = grpreg::grpreg(X = x_tilde_2,
                                                        y = y_tilde_2,
                                                        group = rep(1, df),
                                                        penalty = "grSCAD",
                                                        family = "gaussian",
                                                        group.multiplier = as.vector(adaptive.weights[paste0("X",j),]),
                                                        lambda = lambda_beta,
                                                        intercept = T)$beta[-1,])

        beta_hat_next[list_group_main[[j]],] <- beta_hat_next_j

      }

      # Update Beta_E
      R <- y - x[, as.vector(do.call(c, list_group_main))] %*% beta_hat_next[as.vector(do.call(c, list_group_main)),] -
        rowSums(do.call(cbind, lapply(unique(group), function(kk) {
          gamma_hat_next[kk, ] *
            (x[,list_group_main[[kk]], drop = F] %*% beta_hat_next[list_group_main[[kk]], , drop = F])
        })))


      x_tilde_E <- x[,"X_E", drop = F] +
        rowSums(do.call(cbind, lapply(unique(group), function(kk) {
          gamma_hat_next[kk, ] *
            (x[,list_group_inter[[kk]], drop = F] %*% beta_hat_next[list_group_main[[kk]], , drop = F])
        })))


      beta_hat_next["X_E", ] <- soft(x = x_tilde_E,
                                     y = R,
                                     weight = adaptive.weights["X_E", , drop = F],
                                     lambda = lambda_beta)

      Q[m + 1, 1] <- Q_theta(x = x, y = y,
                             beta = beta_hat_next,
                             gamma = gamma_hat_next,
                             lambda.beta = lambda_beta,
                             lambda.gamma = lambda_gamma,
                             weights = adaptive.weights,
                             main.effect.names = list_group_main,
                             interaction.names = list_group_inter,
                             group = group)


      # Theta_next <- c(drop(beta_hat_next), drop(gamma_hat_next))
      converged <- abs(Q[m,1] - Q[m + 1, 1])/abs(Q[m,1]) < thresh
      converged <- if (is.na(converged)) FALSE else converged
      if (verbose) message(sprintf("Iteration: %f, Converged: %f, Crossprod: %f", m, converged, abs(Q[m,1] - Q[m + 1, 1])/abs(Q[m,1])))

      beta_hat_previous <- beta_hat_next
      gamma_hat_previous <- gamma_hat_next

      m <- m + 1
      # adaptive weight for each tuning parameter. currently this is the
      # same for iterations, but I am coding it here
      # for flexibility in case we want to change the weights at each iteration

      # adaptive.weights <- update_weights(betas = beta_hat_previous,
      #                                    gammas = gamma_hat_next,
      #                                    main.effect.names = main.effect.names,
      #                                    interaction.names = interaction.names)

    }

    Betas_and_Alphas <- convert2(beta = beta_hat_next,
                                 gamma = gamma_hat_next,
                                 main.effect.names = list_group_main,
                                 interaction.names = list_group_inter,
                                 group = group)

    dfbeta <- length(nonzero(Betas_and_Alphas[c(main_effect_names,"X_E"),]))
    dfalpha <- length(nonzero(Betas_and_Alphas[interaction_names,]))
    deviance <- crossprod(y - x %*% Betas_and_Alphas)
    devRatio <- 1 - deviance/nulldev
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
      update_weights(betas = betaWarmStart,
                     # gammas = gammaWarmStart,
                     alphas = Betas_and_Alphas[interaction_names,,drop=F],
                     main.effect.names = list_group_main,
                     interaction.names = list_group_inter,
                     group = group) } else {
                       adaptive.weights.start
                     }

    devianceDiff <- outPrint[lambdaIndex,"deviance"] - outPrint[lambdaIndex-1,"deviance"]
    #coefficientMat[,LAMBDA] <- Betas_and_Alphas

    if (length(devianceDiff)!=0 && !is.na(devianceDiff) && devRatio>1e-3) {
      # if (devianceDiff < 1e-50 | outPrint[LAMBDA,"percentDev"] > 0.999) break }
      if (outPrint[LAMBDA,"percentDev"] > 0.999) break }

    # pb$tick()
  }

  # Rprof()
  # proftable(tmp)
  # summaryRprof(tmp)
  # b <- Sys.time()

  # outPrint[complete.cases(outPrint),]
  # coefficientMat
  # b-a

  # browser()

  beta_hat_next_list <- lapply(seq_len(ncol(betaMat)),
                               function(i) betaMat[,i,drop=F])
  # convert to original scale
  betas_original_scale_list <- if (normalize) lapply(beta_hat_next_list, function(i) i / sx[c(main_effect_names, "X_E")]) else lapply(beta_hat_next_list, function(i) i )

  gamma_hat_next_list <- lapply(seq_len(ncol(gammaMat)),
                                function(i) gammaMat[,i,drop=F])

  gammas_original_scale_list <- if (normalize) lapply(gamma_hat_next_list, function(i) i / sx[interaction_names]) else lapply(gamma_hat_next_list, function(i) i )

  # convert gammas to alphas
  betas_alphas_original_scale <- mapply(convert2,
                                        beta = betas_original_scale_list,
                                        gamma = gammas_original_scale_list,
                                        MoreArgs = list(main.effect.names = list_group_main,
                                                        interaction.names = list_group_inter,
                                                        group = group))

  dimnames(betas_alphas_original_scale) <- list(c(main_effect_names,"X_E",interaction_names),
                                                paste0("s",1:nlambda))

  betas_original_scale <- betas_alphas_original_scale[c(main_effect_names,"X_E"), , drop = F]
  alphas_original_scale <- betas_alphas_original_scale[interaction_names, , drop = F]

  b0 <- vector(length = nlambda)
  for (lam in seq_len(nlambda)) {
    b0[lam] <- by - sum(betas_original_scale[,lam,drop = F] * bx[c(main_effect_names,"X_E")]) -
      sum(alphas_original_scale[,lam,drop=F] * bx[interaction_names])
  }
  names(b0) <- paste0("s", 1:nlambda)

  gamma_final <- as(matrix(unlist(gammas_original_scale_list, use.names = F),
                           ncol = nlambda,
                           byrow = T,
                           dimnames = list(interaction_names, paste0("s",1:nlambda))),
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
              group = group,
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
              design = design,
              df = df,
              interaction.names = interaction_names,
              main.effect.names = c(main_effect_names,"X_E"))
  class(out) <- "lspath"
  return(out)

}



