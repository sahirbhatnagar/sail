#' Gaussian Response fitting function with warm starts
#'
lspath <- function(x,
                   y,
                   e,
                   df,
                   degree,
                   group.penalty,
                   weights,
                   nlambda,
                   thresh,
                   maxit,
                   verbose,
                   alpha,
                   nobs,
                   nvars,
                   jd,
                   vp,
                   we,
                   wj,
                   wje,
                   flmin,
                   vnames,
                   ne, #dfmax
                   ulam) {

  # Basis Expansion and Design Matrix ---------------------------------------

  # Rprof(tmp <- tempfile())

  # group membership
  group <- rep(seq_len(nvars), each = df)

  # Expand X's
  Phi_j_list <- lapply(seq_len(nvars), function(j) splines::bs(x[,j], df = df, degree = degree))
  Phi_j <- do.call(cbind, Phi_j_list)
  main_effect_names <- paste(rep(vnames, each = df), rep(seq_len(df), times = nvars), sep = "_")
  dimnames(Phi_j)[[2]] <- main_effect_names

  # X_E x Phi_j
  XE_Phi_j_list <- lapply(Phi_j_list, function(i) e * i)
  XE_Phi_j <- do.call(cbind, XE_Phi_j_list)
  interaction_names <- paste(main_effect_names, "X_E", sep = ":")
  dimnames(XE_Phi_j)[[2]] <- interaction_names

  design <- cbind(Phi_j, "X_E" = e, XE_Phi_j)

  nulldev <- as.numeric(crossprod(y))

  # Initialize -------------------------------------------------------------

  # the initial values here dont matter, since at Lambda_max everything is 0
  b0 <- mean(y)
  betaE <- 0
  theta <- split(setNames(rep(0, length(main_effect_names)), main_effect_names), group)
  gamma <- rep(0, nvars)
  theta_next <- theta
  R.star <- y - b0
  # update this at the end once betaE and theta are updated. x_tilde is used for gamma update
  x_tilde <- matrix(0, nrow = nobs, ncol = nvars)
  add_back <- rep(0, nobs)

  Theta_init <- c(b0, betaE, do.call(c,theta), gamma)

  # Lambda Sequence ---------------------------------------------------------

  if (ulam == 0) {

    # R1 <- R2 <- y - b0 # this is used as the starting residual for Gamma and Theta update
    term1 <- (1 / we) * (crossprod(e, R.star))
    term2 <- (1 / wj) * sapply(Phi_j_list, function(i) l2norm(crossprod(i, R.star)))
    lambda_max <- (1 / (nobs * (1 - alpha))) * max(term1, max(term2))
    lambdas <- rev(exp(seq(log(flmin * lambda_max), log(lambda_max), length.out = nlambda)))
    lambdaNames <- paste0("s", seq_along(lambdas))

  } else {

    # not sure what to do yet, need to think about cv.sail and supplying the same lambda.sequence
    # or when using adaptive lasso?

  }

  # for all the x_tilde in zero_x_tilde, return the following matrix with 0 for each coefficient
  # this is like a place holder.
  coef_zero_gamma_matrix <- matrix(data = 0, nrow = nvars, ncol = 1,
                                   dimnames = list(vnames))

  # Objects to store results ------------------------------------------------

  # matrix to store results of betas and alphas on standardized scale
  coefficientMat <- matrix(nrow = length(c(main_effect_names,"X_E",interaction_names)),
                           ncol = nlambda,
                           dimnames = list(c(main_effect_names,"X_E", interaction_names),
                                           lambdaNames))

  betaMat <- matrix(nrow = length(c("b0",main_effect_names,"X_E")), ncol = nlambda,
                    dimnames = list(c("b0",main_effect_names,"X_E"),
                                    lambdaNames))

  gammaMat <- matrix(nrow = nvars, ncol = nlambda,
                     dimnames = list(c(paste0(vnames,"X_E")),
                                     lambdaNames))

  outPrint <- matrix(NA, nrow = nlambda, ncol = 5,
                     dimnames = list(lambdaNames,
                                     c("dfBeta","dfAlpha","deviance",
                                       "percentDev",
                                       "lambda")))

  # trying to implement this one at a time, so that
  # pb <- progress::progress_bar$new(
  #   format = "  fitting over all pairs of tuning parameters [:bar] :percent eta: :eta",
  #   total = 100, clear = FALSE, width= 90)
  # pb$tick(0)

  # Lambda Loop Start -------------------------------------------------------

  for (LAMBDA in lambdas) {

    # LAMBDA = lambdas[20]
    #======================

    lambdaIndex <- which(LAMBDA==lambdas)

    if (verbose) {
      message(sprintf("Index: %g, lambda: %0.2f", lambdaIndex, LAMBDA))
    }

    # store likelihood values at each iteration in a matrix Q
    # rows: iteration number
    Q <- vector("numeric", length = maxit + 1)

    # store the value of the likelihood at the 0th iteration
    Q[1] <- (1 / (2 * nobs)) * crossprod(R.star)

    #iteration counter
    m <- 1

    # to enter while loop
    converged <- FALSE

    # un-comment if we dont want warm starts
    # b0 <- mean(y)
    # theta <- split(setNames(rep(0, length(main_effect_names)), main_effect_names), group)
    # betaE <- 0
    # gamma <- rep(0, nvars)
    # R.star <- y - b0
    # b0_next <- b0 ;
    # theta_next <- theta

    # While loop for convergence at a given Lambda value ----------------------

# browser()
    while (!converged && m < maxit){

      # Theta_init <- c(drop(beta_hat_previous), drop(gamma_hat_previous))

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update gamma (interaction parameter)
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      # move this down to the betaE update, because we need this same computation for it
      # after updating the theta. on the first pass this is 0 anyways
      # R1 is initiated above

      # x_tilde <- betaE * do.call(cbind,
      #                            lapply(seq_along(XE_Phi_j_list),
      #                                   function(i) XE_Phi_j_list[[i]] %*% theta[[i]]))


      # add_back <- rowSums(sweep(x_tilde, 2, gamma, FUN = "*"))

      # gamma_Phi_tilde_theta_sum <- rowSums(sweep(Phi_tilde_theta, 2, gamma, FUN = "*"))
      # R1 <- y - b0 - betaE * e - rowSums(Phi_j_theta_j)
      # message("Update gamma")

      R <- R.star + add_back

      # indices of the x_tilde matrices that have all 0 columns
      zero_x_tilde <- dim(check_col_0(x_tilde))[2]
      # message(zero_x_tilde)
      gamma_next <- if (zero_x_tilde == 0) drop(coef_zero_gamma_matrix) else {
        coef(glmnet::glmnet(
          x = x_tilde,
          y = R,
          # thresh = 1e-10,
          penalty.factor = wje,
          lambda = c(.Machine$double.xmax, LAMBDA * alpha),
          standardize = F, intercept = F))[-1,2]
      }

      # Delta <- rowSums(
      #   do.call(cbind,
      #           lapply(seq_along(Phi_j_list),
      #                  function(i) (
      #                  ( ( gamma[i] - gamma_next[i] ) * betaE * XE_Phi_j_list[[i]] ) %*% theta[[i]]
      #                  )
      #           )
      #   )
      # )

      Delta <- rowSums(sweep(x_tilde, 2, (gamma - gamma_next), FUN = "*"))


      R.star <- R.star + Delta

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update theta (main effect parameters)
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      # message("Update theta")


      x_tilde_2 <- lapply(seq_along(Phi_j_list),
                          function(i) Phi_j_list[[i]] + gamma_next[i] * betaE * XE_Phi_j_list[[i]])

      # converged_theta <- FALSE
      # k <- 1
      # while (!converged_theta && k < maxit){

        for (j in seq_len(nvars)) {

          #   delta <- if (j == 1) 0 else {
          #     x_tilde_2[[j]] %*% theta_next[[j]] - x_tilde_2[[(j-1)]] %*% theta_next[[(j-1)]]
          #   }
          #
          #   R2 <- R2 + delta

          R <- R.star + x_tilde_2[[j]] %*% theta_next[[j]]

          # R2 <- y - b0_next - betaE * e - rowSums(do.call(cbind,
          #                                         lapply(seq_along(x_tilde_2),
          #                                                function(i) x_tilde_2[[i]] %*% theta_next[[i]]))) +
          #   x_tilde_2[[j]] %*% theta_next[[j]]

          theta_next_j <- switch(group.penalty,
                                 gglasso = coef(gglasso::gglasso(x = x_tilde_2[[j]],
                                                                 y = R,
                                                                 # eps = 1e-10,
                                                                 group = rep(1, df),
                                                                 pf = wj[j],
                                                                 lambda = LAMBDA * (1 - alpha),
                                                                 intercept = F))[-1,],
                                 MCP = grpreg::grpreg(X = x_tilde_2[[j]],
                                                      y = R,
                                                      group = rep(1, df),
                                                      penalty = "grMCP",
                                                      family = "gaussian",
                                                      group.multiplier = as.vector(wj[j]),
                                                      lambda = LAMBDA * (1 - alpha),
                                                      intercept = T)$beta[-1,],
                                 SCAD = grpreg::grpreg(X = x_tilde_2[[j]],
                                                       y = R,
                                                       group = rep(1, df),
                                                       penalty = "grSCAD",
                                                       family = "gaussian",
                                                       group.multiplier = as.vector(wj[j]),
                                                       lambda = LAMBDA * (1 - alpha),
                                                       intercept = T)$beta[-1,])


          Delta <- x_tilde_2[[j]] %*% (theta_next[[j]] - theta_next_j)

          theta_next[[j]] <- theta_next_j

          R.star <- R.star + Delta

        }

      #   criterion_theta <- l2norm(do.call(c,theta_next) - do.call(c,theta)) ^ 2
      #   converged_theta <- criterion_theta < 1e-2
      #   converged_theta <- if (is.na(converged_theta)) FALSE else converged_theta
      #   # if (verbose) message(sprintf("Iteration: %f, l2norm((theta_next - theta)^2: %f", k,
      #   #                              criterion_theta))
      #   k <- k + 1
      #   theta <- theta_next
      #
      # }

      # used to check convergence
      theta_next_vec <- do.call(c, theta_next)

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update betaE
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      # message("Update betaE")

      # this will be used for R1 (the residual for gamma update) also and R4
      # Phi_j_theta_j <- do.call(cbind,
      #                          lapply(seq_along(Phi_j_list),
      #                                 function(i) Phi_j_list[[i]] %*% theta_next[[i]]))

      # R3 <- y - b0_next - rowSums(Phi_j_theta_j) #- rowSums(sweep(Phi_j_theta_j, 2, gamma_next, FUN = "*"))

      # this can be used for betaE, b0 and gamma update!
      Phi_tilde_theta <- do.call(cbind,
                                 lapply(seq_along(XE_Phi_j_list),
                                   function(i) XE_Phi_j_list[[i]] %*% theta_next[[i]]))

      gamma_Phi_tilde_theta_sum <- rowSums(sweep(Phi_tilde_theta, 2, gamma_next, FUN = "*"))
      x_tilde_E <- e + gamma_Phi_tilde_theta_sum

      R <- R.star + betaE * x_tilde_E

      betaE_next <- SoftThreshold(x = (1 / (nobs * we)) * sum(x_tilde_E * R),
                                  lambda = LAMBDA * (1- alpha))

      Delta <- (betaE - betaE_next) * x_tilde_E

      R.star <- R.star + Delta

      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update beta0
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      # message(sprintf("Update beta0, lambda = %0.2f", LAMBDA))

      # this is the linear predictor without the intercept, used for b0 update and objective function value
      # lp <- betaE_next * e - rowSums(Phi_j_theta_j) - betaE_next * gamma_Phi_tilde_theta_sum
      #
      # R4 <- y - lp

      R <- R.star + b0
      b0_next <- mean(R)

      # used for gamma update
      x_tilde <- betaE_next * Phi_tilde_theta
      add_back <- rowSums(sweep(x_tilde, 2, gamma_next, FUN = "*"))
      # message(sprintf("betaE = %0.2f", betaE_next))

      Delta <- (b0 - b0_next)

      R.star <- R.star + Delta

      Q[m+1] <- Q_theta(R = R.star, nobs = nobs, lambda = LAMBDA, alpha = alpha,
                        we = we, wj = wj, wje = wje, betaE = betaE_next,
                        theta_list = theta_next, gamma = gamma_next)

      Theta_next <- c(b0_next, betaE_next, theta_next_vec, gamma_next)

      criterion <- abs(Q[m] - Q[m + 1])/abs(Q[m])
      # criterion <- l2norm(Theta_next - Theta_init)^2
      converged <- criterion < thresh
      converged <- if (is.na(converged)) FALSE else converged
      if (verbose) message(sprintf("Iteration: %f, Converged: %f, Crossprod: %f", m, converged,
                                   criterion))

      b0 <- b0_next
      betaE <- betaE_next
      theta <- theta_next
      gamma <- gamma_next
      Theta_init <- Theta_next



      m <- m + 1
      # adaptive weight for each tuning parameter. currently this is the
      # same for iterations, but I am coding it here
      # for flexibility in case we want to change the weights at each iteration

      # adaptive.weights <- update_weights(betas = beta_hat_previous,
      #                                    gammas = gamma_next,
      #                                    main.effect.names = main.effect.names,
      #                                    interaction.names = interaction.names)

    }

    message(paste0(capture.output(Theta_next), collapse = "\n"))

    # Betas_and_Alphas <- convert2(beta = beta_hat_next,
    #                              gamma = gamma_next,
    #                              main.effect.names = list_group_main,
    #                              interaction.names = list_group_inter,
    #                              group = group)
    #
    # dfbeta <- length(nonzero(Betas_and_Alphas[c(main_effect_names,"X_E"),]))
    # dfalpha <- length(nonzero(Betas_and_Alphas[interaction_names,]))
    # deviance <- crossprod(y - x %*% Betas_and_Alphas)
    # devRatio <- 1 - deviance/nulldev
    # outPrint[LAMBDA,] <- c(if (dfbeta==0) 0 else dfbeta,
    #                        if (dfalpha==0) 0 else dfalpha,
    #                        deviance,
    #                        devRatio,
    #                        lambda_beta, lambda_gamma)

    betaMat[,lambdaIndex] <- c(b0_next,theta_next_vec, betaE_next)
    gammaMat[,lambdaIndex] <- gamma_next

    # par(mfrow=c(2,1))
    # matplot(t(betaMat), type = "l")
    # matplot(t(gammaMat), type = "l")

    # for a fixed lambda.beta, keep using the previous solution for each lambda.gamma
    # for the next fixed lambda.beta restart the calculation using the initial value
    # from the uni_fun function
    # betaWarmStart <- if (lambdaIndex %ni% switchWarmStart$X1) beta_hat_next
    # trying without warm start

    # uni_start <- rbind(beta_hat_next, gamma_next)
    # uni_start <- if (lambdaIndex %ni% switchWarmStart$X1) {
    #   rbind(betaWarmStart, gammaWarmStart) } else {
    #     uni_start_iteration1
    #   }

    # uni_start <- uni_start_iteration1

    # need to update weights also!
    # adaptive.weights <- if (lambdaIndex %ni% switchWarmStart$X1) {
    #   update_weights(betas = betaWarmStart,
    #                  # gammas = gammaWarmStart,
    #                  alphas = Betas_and_Alphas[interaction_names,,drop=F],
    #                  main.effect.names = list_group_main,
    #                  interaction.names = list_group_inter,
    #                  group = group) } else {
    #                    adaptive.weights.start
    #                  }

    # adaptive.weights <- adaptive.weights.start

    # devianceDiff <- outPrint[lambdaIndex,"deviance"] - outPrint[lambdaIndex-1,"deviance"]
    #coefficientMat[,LAMBDA] <- Betas_and_Alphas

    # if (length(devianceDiff)!=0 && !is.na(devianceDiff) && devRatio>1e-3) {
      # if (devianceDiff < 1e-50 | outPrint[LAMBDA,"percentDev"] > 0.999) break }
      # if (outPrint[LAMBDA,"percentDev"] > 0.999) break }

    # pb$tick()
  }

  return(betaMat)

  browser()
  #
  #
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

  gamma_next_list <- lapply(seq_len(ncol(gammaMat)),
                                function(i) gammaMat[,i,drop=F])

  gammas_original_scale_list <- if (normalize) lapply(gamma_next_list, function(i) i / sx[interaction_names]) else lapply(gamma_next_list, function(i) i )

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



