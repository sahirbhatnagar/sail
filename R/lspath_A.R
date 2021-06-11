######################################
# R Source code file for least squares strong hierarchy
# this is where most of the work is being done
# not exported
# Author: Sahir Bhatnagar
# Created: 2016
# Updated: April 6, 2018
#####################################

lspathA <- function(x,
                   y,
                   e,
                   basis,
                   center.x,
                   center.e,
                   expand,
                   group,
                   group.penalty,
                   weights,
                   nlambda,
                   thresh,
                   fdev,
                   maxit,
                   verbose,
                   alpha,
                   nobs,
                   nvars,
                   vp, # penalty.factor
                   we, # we,wj,wje are subsets of vp
                   wj,
                   wje,
                   flmin,
                   vnames,
                   ne, # dfmax
                   ulam) {

  # Basis Expansion and Design Matrix ---------------------------------------

  expansion <- design_sail(
    x = x, e = e,expand = expand,weights=weights, group = group, basis = basis, nvars = nvars,
    vnames = vnames, center.x = center.x, center.e = center.e
  )

  y=y-weighted.mean(y,weights)

  # y <- drop(scale(y, center = TRUE, scale = FALSE))
  Phi_j_list <- expansion$Phi_j_list
  Phi_j <- expansion$Phi_j
  XE_Phi_j_list <- expansion$XE_Phi_j_list
  XE_Phi_j <- expansion$XE_Phi_j
  main_effect_names <- expansion$main_effect_names
  interaction_names <- expansion$interaction_names
  ncols <- expansion$ncols
  e <- expansion$E
  # group_list <- split(group, group)


  if (expand) {
    group <- rep(seq_len(nvars), each = ncols)
  }

  # this is used for the predict function
  design <- expansion$design

  nulldev <- as.numeric(crossprod(sqrt(weights)*y))

  # Initialize -------------------------------------------------------------
  # the initial values here dont matter, since at Lambda_max everything is 0
  # b0 <- mean(y)
  betaE <- 0
  theta <- split(stats::setNames(rep(0, length(main_effect_names)), main_effect_names), group)
  gamma <- rep(0, nvars)
  theta_next <- theta
  R.star <- y
  # R.star <- y - b0

  # update this at the end once betaE and theta are updated. x_tilde is used for gamma update
  x_tilde <- matrix(0, nrow = nobs, ncol = nvars)
  add_back <- rep(0, nobs)

  Theta_init <- c(betaE, do.call(c, theta), gamma)

  # Lambda Sequence ---------------------------------------------------------

  if (is.null(ulam)) {
    # R1 <- R2 <- y - b0 # this is used as the starting residual for Gamma and Theta update
    term1 <- (1 / we) * (crossprod(weights*e, R.star))
    term2 <- (1 / wj) * sapply(Phi_j_list, function(i) l2norm(crossprod(weights*i, R.star)))
    ###   r. star????
    lambda_max <- (1 / (nobs * (1 - alpha))) * max(term1[term1 != Inf], max(term2[term2 != Inf]))
    lambdas <- rev(exp(seq(log(flmin * lambda_max), log(lambda_max), length.out = nlambda)))
    lambdaNames <- paste0("s", seq_along(lambdas))
  } else {
    lambdas <- ulam
    lambdaNames <- paste0("s", seq_along(lambdas))
    # not sure what to do yet, need to think about cv.sail and supplying the same lambda.sequence
    # or when using adaptive lasso?
  }

  # for all the x_tilde in zero_x_tilde, return the following matrix with 0 for each coefficient
  # this is like a place holder.

  coef_zero_gamma_matrix <- matrix(
    data = 0, nrow = nvars, ncol = 1,
    dimnames = ifelse(expand, list(vnames), list(paste0("V", seq(nvars))))
  )


  # Objects to store results ------------------------------------------------

  # originalintercept <- stats::setNames(rep(0, nlambda), lambdaNames)

  a0 <- stats::setNames(rep(0, nlambda), lambdaNames)

  environ <- stats::setNames(rep(0, nlambda), lambdaNames)

  betaMat <- matrix(
    nrow = length(main_effect_names), ncol = nlambda,
    dimnames = list(
      main_effect_names,
      lambdaNames
    )
  )

  if (expand) {
    gammaMat <- matrix(
      nrow = nvars, ncol = nlambda,
      dimnames = list(
        c(paste0(vnames, "E")),
        lambdaNames
      )
    )
  } else {
    gammaMat <- matrix(
      nrow = nvars, ncol = nlambda,
      dimnames = list(
        c(paste0("V", seq(nvars))),
        lambdaNames
      )
    )
  }

  alphaMat <- matrix(
    nrow = length(c(main_effect_names)),
    ncol = nlambda,
    dimnames = list(
      paste(main_effect_names, "E", sep = ":"),
      lambdaNames
    )
  )

  converged <- stats::setNames(rep(FALSE, nlambda), lambdaNames)

  outPrint <- matrix(NA,
    nrow = nlambda, ncol = 5,
    dimnames = list(
      lambdaNames,
      c(
        "dfBeta", "dfAlpha", "dfEnviron", "deviance",
        "percentDev"
      )
    )
  )

  active <- vector("list", length = nlambda)
  # browser()

  # Lambda Loop Start -------------------------------------------------------

  lambdas[1] <- .Machine$double.xmax
  for (LAMBDA in lambdas) {
    lambdaIndex <- which(LAMBDA == lambdas)

    if (verbose >= 1) {
      message(sprintf("Index: %g, lambda: %0.4f", lambdaIndex, if (lambdaIndex==1) lambda_max else LAMBDA))
    }

    # store likelihood values at each iteration in a matrix Q
    # rows: iteration number
    Q <- vector("numeric", length = maxit + 1)

    # store the value of the likelihood at the 0th iteration
    Q[1] <- (1 / (2 * nobs)) * crossprod(sqrt(weights)*R.star)

    # iteration counter
    m <- 1


    # un-comment if we dont want warm starts for not converged lambdas

    # if (lambdaIndex > 1) {
    #   if (!converged[lambdaIndex - 1]) {
    #     b0 <- mean(y)
    #     theta <- split(setNames(rep(0, length(main_effect_names)), main_effect_names), group)
    #     betaE <- 0
    #     gamma <- rep(0, nvars)
    #     # R.star <- y - b0
    #     b0_next <- b0 ;
    #     theta_next <- theta
    #   }
    # }

    # While loop for convergence at a given Lambda value ----------------------

    while (!converged[lambdaIndex] && m < maxit) {

      # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update gamma (interaction parameter)
      # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      R <- R.star + add_back
      # indices of the x_tilde matrices that have all 0 columns
      zero_x_tilde <- dim(check_col_0(x_tilde))[2]

      gamma_next <- if (zero_x_tilde == 0) {
        drop(coef_zero_gamma_matrix)
      } else {
        coef(glmnet::glmnet(
          x = x_tilde,
          y = R,
          thresh = 1e-8,
          weights = weights,
          penalty.factor = wje,
          lambda = c(.Machine$double.xmax, LAMBDA * alpha),
          standardize = F, intercept = F
        ))[-1, 2]
      }

      Delta <- rowSums(sweep(x_tilde, 2, (gamma - gamma_next), FUN = "*"))

      R.star <- R.star + Delta

      # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update theta (main effect parameters)
      # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      x_tilde_2 <- lapply(
        seq_along(Phi_j_list),
        function(i) Phi_j_list[[i]]
      )


      ###  Note:  This is only use in the DTRs!!!!!!!!
      x_tilde_2=matrix(unlist(x_tilde_2), ncol=length(x_tilde_2))
      add_back <- rowSums(sweep(x_tilde_2, 2, unlist(theta), FUN = "*"))
      R <- R.star + add_back
      theta_next <- coef(glmnet::glmnet(
          x = x_tilde_2,
          y = R,
          thresh = 1e-8,

          weights = weights,

          penalty.factor = wj,
          lambda = c(.Machine$double.xmax, LAMBDA *(1- alpha)),
          standardize = F, intercept = F
        ))[-1, 2]

      Delta <- rowSums(sweep(x_tilde_2, 2, (unlist(theta) - unlist(theta_next)), FUN = "*"))

      R.star <- R.star + Delta

      # converged_theta <- FALSE
      # k <- 1
      # while (!converged_theta && k < maxit){
      # browser()

      # if (any(wj == 0)) {
      #   for (j in seq_len(nvars)) {
      #     R <- R.star + x_tilde_2[[j]] %*% theta_next[[j]]
      #     if (wj[j] != 0) {
      #       theta_next_j <- switch(group.penalty,
      #         gglasso = coef(gglasso::gglasso(
      #           x = x_tilde_2[[j]],
      #           y = R,
      #           weight=diag(weights),
      #           # eps = 1e-12,
      #           maxit = 100000,
      #           group = if (expand) rep(1, ncols) else rep(1, ncols[j]),
      #           pf = wj[j],
      #           lambda = LAMBDA * (1 - alpha),
      #           intercept = F
      #         ))[-1, ],
      #         grMCP = grpreg::grpreg(
      #           X = x_tilde_2[[j]],
      #           y = R,
      #           group = if (expand) rep(1, ncols) else rep(1, ncols[j]),
      #           penalty = "grMCP",
      #           family = "gaussian",
      #           group.multiplier = as.vector(wj[j]),
      #           lambda = LAMBDA * (1 - alpha),
      #           intercept = T
      #         )$beta[-1, ],
      #         grSCAD = grpreg::grpreg(
      #           X = x_tilde_2[[j]],
      #           y = R,
      #           group = if (expand) rep(1, ncols) else rep(1, ncols[j]),
      #           penalty = "grSCAD",
      #           family = "gaussian",
      #           group.multiplier = as.vector(wj[j]),
      #           lambda = LAMBDA * (1 - alpha),
      #           intercept = T
      #         )$beta[-1, ]
      #       )
      #     } else {
      #       theta_next_j <- stats::lm.fit(x_tilde_2[[j]], R, w=weights)$coef
      #     }

      #     Delta <- x_tilde_2[[j]] %*% (theta_next[[j]] - theta_next_j)
      #
      #     theta_next[[j]] <- theta_next_j
      #
      #     R.star <- R.star + Delta
      #   }
      # } else {
      #   for (j in seq_len(nvars)) {
      #     R <- R.star + x_tilde_2[[j]] %*% theta_next[[j]]
      #     theta_next_j <- switch(group.penalty,
      #       gglasso = coef(gglasso::gglasso(
      #         x = x_tilde_2[[j]],
      #         y = R,
      #
      #         # eps = 1e-12,
      #         weight=diag(weights),
      #         group = if (expand) rep(1, ncols) else rep(1, ncols[j]),
      #         pf = wj[j],
      #         lambda = LAMBDA * (1 - alpha),
      #         intercept = F
      #       ))[-1, ],
      #       grMCP = grpreg::grpreg(
      #         X = x_tilde_2[[j]],
      #         y = R,
      #         group = if (expand) rep(1, ncols) else rep(1, ncols[j]),
      #         penalty = "gel",
      #         family = "gaussian",
      #         group.multiplier = as.vector(wj[j]),
      #         lambda = LAMBDA * (1 - alpha),
      #         intercept = T
      #       )$beta[-1, ],
      #       grSCAD = grpreg::grpreg(
      #         X = x_tilde_2[[j]],
      #         y = R,
      #         group = if (expand) rep(1, ncols) else rep(1, ncols[j]),
      #         penalty = "grSCAD",
      #         family = "gaussian",
      #         group.multiplier = as.vector(wj[j]),
      #         lambda = LAMBDA * (1 - alpha),
      #         intercept = T
      #       )$beta[-1, ]
      #     )
      #
      #     Delta <- x_tilde_2[[j]] %*% (theta_next[[j]] - theta_next_j)
      #
      #     theta_next[[j]] <- theta_next_j
      #
      #     R.star <- R.star + Delta
      #   }
      # }



      # used to check convergence
      # theta_next_vec <- do.call(c, theta_next)   not use gglasso
      theta_next_vec <- theta_next

      # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update betaE
      # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      # this can be used for betaE, b0 and gamma update!
      Phi_tilde_theta <- do.call(
        cbind,
        lapply(
          seq_along(XE_Phi_j_list),
          function(i) XE_Phi_j_list[[i]]
        )
      )

      gamma_Phi_tilde_theta_sum <- rowSums(sweep(Phi_tilde_theta, 2, gamma_next, FUN = "*"))

      x_tilde_E <- e + gamma_Phi_tilde_theta_sum

      R <- R.star + betaE * x_tilde_E

      betaE_next =
        coef(glmnet::glmnet(
          x = cbind(0,x_tilde_E),
          y = R,

          thresh = 1e-8,
          weights = weights,

          penalty.factor = c(1,we),
          lambda = c(.Machine$double.xmax, LAMBDA *(1- alpha)),
          standardize = F, intercept = F
        ))[c(-1,-2), 2]

      Delta <- (betaE - betaE_next) * x_tilde_E

      R.star <- R.star + Delta

      # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # update beta0
      # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      # R <- R.star + b0
      # b0_next <- mean(weights*R)

      # used for gamma update
      x_tilde <- betaE_next * Phi_tilde_theta
      add_back <- rowSums(sweep(x_tilde, 2, gamma_next, FUN = "*"))

      # Delta <- (b0 - b0_next)
      #
      # R.star <- R.star + Delta

      Q[m + 1] <- Q_theta(
        R = R.star, nobs = nobs,  lambda = LAMBDA, alpha = alpha,weights = weights,
        we = we, wj = wj, wje = wje, betaE = betaE_next,
        theta_list = theta_next, gamma = gamma_next
      )

      Theta_next <- c( betaE_next, theta_next_vec, gamma_next)

      criterion <- abs(Q[m] - Q[m + 1]) / abs(Q[m])
      # criterion <- l2norm(Theta_next - Theta_init)
      converged[lambdaIndex] <- criterion < thresh
      converged[lambdaIndex] <- if (is.na(converged[lambdaIndex])) FALSE else converged[lambdaIndex]
      if (verbose >= 2) {
        message(sprintf(
          "Iteration: %f, Criterion: %f", m, criterion
        ))
      }
      R <- R.star + add_back
      # b0 <- b0_next
      betaE <- betaE_next
      theta <- theta_next
      gamma <- gamma_next
      Theta_init <- Theta_next

      m <- m + 1
    }

    # Store Results -----------------------------------------------------------
    # originalintercept[lambdaIndex]=b0_next
    environ[lambdaIndex] <- betaE_next
    betaMat[, lambdaIndex] <- theta_next_vec
    gammaMat[, lambdaIndex] <- gamma_next
    alphaMat[, lambdaIndex] <- do.call(c, lapply(seq_along(theta_next), function(i) betaE_next * gamma_next[i] ))
    # a0[lambdaIndex] <- b0_next

        # -crossprod(as.vector(expansion$mPhi_j),as.vector(do.call(cbind,theta_next))) -
        # expansion$mE * betaE_next -
        # crossprod(as.vector(expansion$mXE_Phi_j),alphaMat[,lambdaIndex])

    active[[lambdaIndex]] <- c(
      unique(gsub("\\_\\d*", "", names(which(abs(betaMat[, lambdaIndex]) > 0)))),
      unique(gsub("\\_\\d*", "", names(which(abs(alphaMat[, lambdaIndex]) > 0)))),
      if (abs(environ[lambdaIndex]) > 0) "E"
    )

    deviance <- crossprod(sqrt(weights)*R.star)
    devRatio <- 1 - deviance / nulldev
    dfbeta <- sum(abs(betaMat[, lambdaIndex]) > 0) / ifelse(expand, ncols, 1)
    dfalpha <- sum(abs(alphaMat[, lambdaIndex]) > 0) / ifelse(expand, ncols, 1)
    dfenviron <- sum(abs(environ[lambdaIndex]) > 0)


    outPrint[lambdaIndex, ] <- c(
      if (dfbeta == 0) 0 else dfbeta,
      if (dfalpha == 0) 0 else dfalpha,
      if (dfenviron == 0) 0 else dfenviron,
      deviance, devRatio
    )


    # dfmax
    if (sum(outPrint[lambdaIndex, c("dfBeta", "dfAlpha", "dfEnviron")]) > ne) break

    # dev.off()
    # par(mfrow=c(3,1), mai = c(0.2,0.2,0.2,0.2))
    # matplot(t(betaMat), type = "l")
    # matplot(t(gammaMat), type = "l")
    # matplot(t(alphaMat), type = "l")
    # browser()
    # devianceDiff <- outPrint[lambdaIndex,"deviance"] - outPrint[lambdaIndex-1,"deviance"]
    devianceDiff <- (outPrint[lambdaIndex, "percentDev"] - outPrint[lambdaIndex - 1, "percentDev"]) /
      outPrint[lambdaIndex - 1, "percentDev"]
    if (length(devianceDiff) != 0 && !is.na(devianceDiff) && devRatio > 1e-3) {
      if (devianceDiff < fdev | outPrint[lambdaIndex, "percentDev"] > 0.999) break
    }
    # if (outPrint[LAMBDA,"percentDev"] > 0.999) break #}
  }

  beta_final <- methods::as(betaMat, "dgCMatrix")
  alpha_final <- methods::as(alphaMat, "dgCMatrix")
  gamma_final <- methods::as(gammaMat, "dgCMatrix") # used for KKT check


  # browser()

  if (all(!converged)) warning("The algorithm did not converge for all values of lambda.\n
                               Try changing the value of alpha and the convergence threshold.")

  lambdas[1] <- lambda_max

  out <- list(
    # b0=originalintercept[converged],
    a0 = a0[converged],
    beta = beta_final[, converged, drop = FALSE],
    alpha = alpha_final[, converged, drop = FALSE],
    gamma = gamma_final[, converged, drop = FALSE],
    bE = environ[converged],
    active = active[converged],
    lambda = lambdas[converged],
    lambda2 = alpha,
    dfbeta = outPrint[converged, "dfBeta"],
    dfalpha = outPrint[converged, "dfAlpha"],
    dfenviron = outPrint[converged, "dfEnviron"],
    dev.ratio = outPrint[converged, "percentDev"],
    converged = converged,
    nlambda = sum(converged),
    design = design,
    # we = we,
    # wj = wj,
    # wje = wje,
    # Phi_j_list = Phi_j_list,
    # XE_Phi_j_list = XE_Phi_j_list,
    # Phi_j = Phi_j,
    # XE_Phi_j = XE_Phi_j,
    # x = x,
    # e = e,
    # y = y,
    nobs = nobs,
    nvars = nvars,
    vnames = vnames,
    ncols = ncols,
    center.x = center.x,
    center.e = center.e,
    basis = basis,
    expand = expand,
    group = group,
    interaction.names = interaction_names,
    main.effect.names = main_effect_names
  )
  class(out) <- "lspath"
  return(out)
}
