# Downloaded from http://arantxa.ii.uam.es/%7edhernan/RMGPC/ with small modifications

#########################################################################################################
# Module: epORGPC.R
# Date  : January 2011
# Author: Daniel Hernandez Lobato and José Miguel Hernández Lobato
#
# Implements a robust multi-class Gaussian Process classification
# using parallel EP updates for the likelihood and the prior. The Gaussian Process has a mechanism to detect points
# that lie at the wrong side of the decision boundary and are hence incorrectly classified. Basically, we
# have a latent variable for each instance that indicates werther of not that instance should not be considered
# for prediction.
#
#########################################################################################################
#
# EXPORTED FUNCTIONS: epORGPC
#
#########################################################################################################
#
# The main function that should be used in this module is "epORGPC". You have to call it with the arguments:
#
# 	X  -> Design matrix for the classification problem.
#       Y  -> Target vector for the classification problem. (must be an R factor)
#       kernel -> Either linear or Gaussian
#       a -> First Parameter of the beta prior
#       b -> Second Parameter of the beta prior
#       noise -> variance of the latent noise
#       sigma -> parameter for the gaussian kernel
#       intial_a -> initial EP approximation (NULL if not used)
#
# "epORGPC" returns the approximate distribution for the posterior as a list with several components.
#  Another posibility is to call optimizeEvidenceGaussianKernelNoNoise(X, Y) which will optimize the kernel
#  hyper-paramters and return the corresponding EP approximation.
#

#' @export
epORGPC <- function(X, Y, kernel = "linear", a = 1, b = 9, noise = rep(1, length(levels(Y))), sigma = rep(1, length(levels(Y))), initial_a = NULL, tol = 1e-3, verbose=F) {

	# We extend the data with an extra component and compute the number of classes, features and samples

	X <- t(cbind(X, rep(1, nrow(X))))

	nK <- length(levels(Y))
	d <- nrow(X)
	n <- ncol(X)

	if (length(Y) != n) {
	  stop("Y must have the same length as X")
	}

	# We initialize the approximate terms. Note that we are only approximating the likelihood.
	# The Gaussian Process Prior is approximated exactly. Each term of the likelihood is defined as
	# s_ik exp(-0.5 * vfHat[ i, ky_i ]^-1 (fiy_i - mfHat[ i, ky_i ])^2) exp(-0.5 * vfHat[ i, k ]^-1 (fik - mfHat[ i, k ])^2)
	# where y_i is the label of the data instance and k != y_i and fik is the k-th Gaussian Process evaluated on x_i.

	t0Hat <- list(mfHatk = matrix(0, n, nK), mfHatkyi = matrix(0, n, nK), vfHatk = matrix(0, n, nK), vfHatkyi = matrix(0, n, nK),
		sfHat = matrix(0, n, nK), pfHat = matrix(0, n, nK))

	t1Hat <- list(aHat = rep(1, n), bHat = rep(1, n), pHat = rep(logit(a / (a + b)), n), sHat = rep(0, n))

	C <- list()
	kernelFun <- list()

	# We compute the kernel matrix

	if (kernel == "linear") {
		for (k in 1 : nK) {
			C[[ k ]] <- linearKernel(X, noise[ k ])
			kernelFun[[ k ]] <- function(x, noise, sigma) linearKernel(x, noise = noise)
		}
	} else {
		for (k in 1 : nK) {
			C[[ k ]] <- gaussianKernel(X, sigma[ k ], noise[ k ])
			kernelFun[[ k ]] <- function(x, noise, sigma) gaussianKernel(x, sigma = sigma, noise = noise)
		}
	}

	Cinv <- lapply(C, solve)

	# We initialize the posterior approximation to the Gaussian priors or to the provided approximation

	if (! is.null(initial_a)) {
		a <- initial_a
		a$C <- C
		a$Cinv <- Cinv
		a$kernelFun <- kernelFun
		a$noise <- noise
		a$sigma <- sigma
		a$p <- a$p

		M <- list()

		for (k in 1 : nK) {

			lambda <- a$t0Hat$vfHatk[ , k ]
                	lambda[ Y == levels(Y)[ k ] ] <- apply(a$t0Hat$vfHatkyi[ Y == levels(Y)[ k ], ], 1, sum)

			upsilon <- a$t0Hat$mfHatk[ , k ]
			upsilon[ Y == levels(Y)[ k ] ]  <- apply(a$t0Hat$mfHatkyi[ Y == levels(Y)[ k ], ], 1, sum)

			M[[ k ]] <- solve(Cinv[[ k ]] + diag(lambda))

			a$vfik[ , k ] <- diag(M[[ k ]])^-1
			a$mfik[ , k ] <- (M[[ k ]] %*% upsilon) * a$vfik[ , k ]
		}
	} else {

		a <- list(t0Hat = t0Hat, t1Hat = t1Hat, indexNegative = NULL, C = C, mfik = matrix(0, n, nK), vfik = matrix(0, n, nK),
			Cinv = Cinv, kernelFun = kernelFun, p = rep(logit(a / (a + b)), n), O = (1 / nK)^(1 / (nK - 1)), M = list(),
			aPrior = a, bPrior = b, a = a, b = b, kernel = kernel)

		for (k in 1 : nK)
			a$vfik[ , k ] <- diag(C[[ k ]])^-1
	}

	# Main loop of EP

	j <- 1
	damping <- 0.5
	convergence <- FALSE

	while(!convergence && j < 1000) {

		aOld <- a

		# We process the approximate terms for the prior for the latent variables in parallel.

		for (i in 1 : n) {

			# We compute the old posterior

			pfiOldLogit <- a$p[ i ] - a$t1Hat$pHat[ i ]
			pfiOld <- sigmoid(pfiOldLogit)
			aBetaOld <- a$a - a$t1Hat$aHat[ i ] + 1
			bBetaOld <- a$b - a$t1Hat$bHat[ i ] + 1

			if (aBetaOld > 0 && bBetaOld > 0) {

				# We compute the new posterior

				epsOld <- aBetaOld / (aBetaOld + bBetaOld)
				Zi <- pfiOld * epsOld + (1 - pfiOld) * (1 - epsOld)
				pfiNew <- a$p[ i ] - a$t1Hat$pHat[ i ] + log(aBetaOld) - log(bBetaOld)

				# We compute the mean and the variance of the new posterior for the prior

				eE <- 1 / Zi * ((1 - pfiOld) * (1 - epsOld) * aBetaOld / (aBetaOld + bBetaOld + 1) + pfiOld * epsOld * (aBetaOld + 1) / (aBetaOld + bBetaOld + 1))
				eE2 <- 1 / Zi * ((1 - pfiOld) * (1 - epsOld) * aBetaOld  * (aBetaOld + 1) / (aBetaOld + bBetaOld + 1)  / (aBetaOld + bBetaOld + 2) +
					pfiOld * epsOld * (aBetaOld + 1) * (aBetaOld + 2) / (aBetaOld + bBetaOld + 1) / (aBetaOld + bBetaOld + 2))

				aNew <- (eE - eE2) / (eE2 - eE^2) * eE
				bNew <- (eE - eE2) / (eE2 - eE^2) * (1 - eE)

				# We compute the new approximate term

				pHatNew <- pfiNew - a$p[ i ] + a$t1Hat$pHat[ i ]
				aHatNew <- aNew - aBetaOld + 1
				bHatNew <- bNew - bBetaOld + 1

				# We update the approximate term

				a$t1Hat$pHat[ i ] <- damping * pHatNew + (1 - damping) * a$t1Hat$pHat[ i ]
				a$t1Hat$aHat[ i ] <- damping * aHatNew + (1 - damping) * a$t1Hat$aHat[ i ]
				a$t1Hat$bHat[ i ] <- damping * bHatNew + (1 - damping) * a$t1Hat$bHat[ i ]

				# We update the posterior approximation

				a$p[ i ] <- pfiOldLogit + a$t1Hat$pHat[ i ]
				a$a <- aBetaOld + a$t1Hat$aHat[ i ] - 1
				a$b <- bBetaOld + a$t1Hat$bHat[ i ] - 1
			}
		}

		# We process the approximate terms for the likelihood in parallel.
		# First we remove each term to compute the corresponding posterior approximation of fik, Q^\ik.

		lev_Y <- levels(Y)
		for (i in 1 : n) {
		  cur_ind <- (lev_Y != Y[i])
			for (k in (1 : nK)[cur_ind]) {
				yi <- which(!cur_ind)

				vfikOld <- a$vfik[ i, k ] - a$t0Hat$vfHatk[ i, k ]
				vfikyiOld <- a$vfik[ i, yi ] - a$t0Hat$vfHatkyi[ i, k ]

				mfikOld <- a$mfik[ i, k ] - a$t0Hat$mfHatk[ i, k ]
				mfikyiOld <- a$mfik[ i, yi ] - a$t0Hat$mfHatkyi[ i, k ]

				pfiOld <- sigmoid(a$p[ i ] - a$t0Hat$pfHat[ i, k ])

				# We compute the new posterior

				if (length(vfikyiOld^-1 + vfikOld^-1) == 0) {
				  print("ERROR")
				}
				if (vfikyiOld^-1 + vfikOld^-1 > 0) {

					u_ik <- (mfikyiOld * vfikyiOld^-1 - mfikOld * vfikOld^-1) / sqrt(vfikyiOld^-1 + vfikOld^-1)
					Z_ik <- (1 - pfiOld) * pnorm(u_ik) + pfiOld * a$O
					a_ik <- (1 - pfiOld) * dnorm(u_ik) / Z_ik * 1 / sqrt(vfikOld^-1 + vfikyiOld^-1)

	                                meanfikyiNew <- (mfikyiOld + a_ik ) / vfikyiOld
	                                meanfikNew <- (mfikOld - a_ik) / vfikOld

					# We compute the new approximate term

					vfHatkNew <- ((vfikOld^-1 + vfikyiOld^-1) / (a_ik * (meanfikyiNew - meanfikNew)) - vfikOld^-1)^-1
					vfHatkyiNew <- ((vfikOld^-1 + vfikyiOld^-1) / (a_ik * (meanfikyiNew - meanfikNew)) - vfikyiOld^-1)^-1
					mfHatkNew <- meanfikNew * vfHatkNew - a_ik
					mfHatkyiNew <- meanfikyiNew * vfHatkyiNew + a_ik

					pfHatNew <- log(a$O) - pnorm(u_ik, log.p = TRUE)

					# We update the approximate term

					a$t0Hat$vfHatk[ i, k ] <- (damping * vfHatkNew + (1 - damping) * a$t0Hat$vfHatk[ i, k ])
					a$t0Hat$mfHatk[ i, k ] <- (damping * mfHatkNew + (1 - damping) * a$t0Hat$mfHatk[ i, k ])
					a$t0Hat$vfHatkyi[ i, k ] <- (damping * vfHatkyiNew + (1 - damping) * a$t0Hat$vfHatkyi[ i, k ])
					a$t0Hat$mfHatkyi[ i, k ] <- (damping * mfHatkyiNew + (1 - damping) * a$t0Hat$mfHatkyi[ i, k ])
					a$t0Hat$pfHat[ i, k ] <-damping * pfHatNew + (1 - damping) * a$t0Hat$pfHat[ i, k ]
				}
			}
		}

		# We update the posterior approximation, which is recomputed as the product of the approximate terms and the Gaussian Processes Priors


		for (k in 1 : nK) {

			lambda <- a$t0Hat$vfHatk[ , k ]
    	lambda[ Y == levels(Y)[ k ] ] <- apply(a$t0Hat$vfHatkyi[ Y == levels(Y)[ k ], ], 1, sum)

			upsilon <- a$t0Hat$mfHatk[ , k ]
			upsilon[ Y == levels(Y)[ k ] ]  <- apply(a$t0Hat$mfHatkyi[ Y == levels(Y)[ k ], ], 1, sum)

			a$M[[ k ]] <- solve(a$Cinv[[ k ]] + diag(lambda))

			a$vfik[ , k ] <- diag(a$M[[ k ]])^-1
			a$mfik[ , k ] <- (a$M[[ k ]] %*% upsilon) * a$vfik[ , k ]
		}

		a$p <- apply(a$t0Hat$pfHat, 1, sum) + a$t1Hat$pHat

		# Annealed damping scheme

#		damping <- damping * 0.99

		# We check for convergence

		convergence <- checkConvergence(a, aOld, tol, verbose)$flag

		j <- j + 1
	}

	# We compute the evidence and store useful information (data and labels)

	a$sigma <- sigma
	a$noise <- noise
	a <- computeEvidence(a, Y, X)
	a$X <- X; a$Y <- Y
	a$prior <- a$a / (a$a + a$b)

	# We return the current approximation

	a
}

##
# Checks convergence of the EP algorithm.
#
# Input:
# 	aOld -> The previous approximation.
# 	aNew -> The new approximation.
# Output:
# 	TRUE if the values in aOld are differ from those in aNew by less than a small constant.
#

checkConvergence <- function(aNew, aOld, tol, verbose=F) {

	convergence <- max(max(abs(aNew$mfik * aNew$vfik^-1 - aOld$mfik * aOld$vfik^-1)))
	convergence <- max(convergence, max(abs(aNew$vfik^-1 - aOld$vfik^-1)))
	convergence <- max(convergence, max(abs(aNew$p - aOld$p)))
	convergence <- max(convergence, max(abs(aNew$a - aOld$a)))
	convergence <- max(convergence, max(abs(aNew$b - aOld$b)))

	if (verbose){
	  print(convergence)
	}

	if (convergence < tol)
		list(flag = TRUE, value = convergence)
	else
		list(flag = FALSE, value = convergence)
}

##
# Function that computes the log evidence
#

computeEvidence <- function(a, Y, X) {

	# First we compute the constant that multiplies each approximate factor

	log_ctes <- NULL

	# First the terms corresponding to the prior

	for (i in 1 : ncol(X))  {

		# We compute the old posterior

		pfiOld <- sigmoid(a$p[ i ] - a$t1Hat$pHat[ i ])
		aBetaOld <- a$a - a$t1Hat$aHat[ i ] + 1
		bBetaOld <- a$b - a$t1Hat$bHat[ i ] + 1

		# We compute the normalization constnat

		epsOld <- aBetaOld / (aBetaOld + bBetaOld)
		Zi <- pfiOld * epsOld + (1 - pfiOld) * (1 - epsOld)

		a$t1Hat$sHat[ i ] <- lbeta(aBetaOld, bBetaOld) - lbeta(a$a, a$b) # Simplified!! (sum 0): + log(Zi) +
										 # log(sigmoid(a$p[ i ]) / pfiOld + sigmoid(-a$p[ i ]) / (1 - pfiOld))
		log_ctes <- c(log_ctes, a$t1Hat$sHat[ i ])
	}

	# Then, the terms corresponding to the likelihood

	for (i in 1 : ncol(X))
		for (k in (1 : length(levels(Y)))[ levels(Y) != Y[ i ] ]) {

			yi <- which(levels(Y) == Y[ i ])

			vfikOld <- a$vfik[ i, k ] - a$t0Hat$vfHatk[ i, k ]
			vfikyiOld <- a$vfik[ i, yi ] - a$t0Hat$vfHatkyi[ i, k ]

			mfikOld <- a$mfik[ i, k ] - a$t0Hat$mfHatk[ i, k ]
			mfikyiOld <- a$mfik[ i, yi ] - a$t0Hat$mfHatkyi[ i, k ]

			pfiOld <- sigmoid(a$p[ i ] - a$t0Hat$pfHat[ i, k ])

			if (vfikyiOld^-1 + vfikOld^-1 > 0) {

				u_ik <- (mfikyiOld * vfikyiOld^-1 - mfikOld * vfikOld^-1) / sqrt(vfikyiOld^-1 + vfikOld^-1)
				Z_ik <- (1 - pfiOld) * pnorm(u_ik) + pfiOld * a$O
				a_ik <- (1 - pfiOld) * dnorm(u_ik) / Z_ik * 1 / sqrt(vfikOld^-1 + vfikyiOld^-1)

				meanfikyiNew <- (mfikyiOld + a_ik ) / vfikyiOld
				meanfikNew <- (mfikOld - a_ik) / vfikOld

				# We compute the normalization contant to guarantee that the approximate term integrates the same

				a$t0Hat$sfHat[ i , k ] <- log(pnorm(u_ik) + a$O) + 0.5 * log(1 + a$t0Hat$vfHatk[ i, k ] / vfikOld) +
					0.5 * log(1 + a$t0Hat$vfHatkyi[ i, k ] / vfikyiOld) +
					-0.5 * (a$mfik[ i, k ]^2 / a$vfik[ i, k ] - mfikOld^2 / vfikOld - a$t0Hat$mfHatk[ i, k ]^2 / a$t0Hat$vfHatk[ i, k ]) +
				-0.5 * (a$mfik[ i, yi ]^2 / a$vfik[ i, yi ] - mfikyiOld^2 / vfikyiOld - a$t0Hat$mfHatkyi[ i, k ]^2 / a$t0Hat$vfHatkyi[ i, k ])

			}

			log_ctes <- c(log_ctes, a$t0Hat$sfHat[ i , k ])
		}

	# Then, we compute the constants that arises from the multiplication of the approximate terms

	B <- NULL
	log_dets <- NULL

	for (k in (1 : length(levels(Y)))) {

		lambda <- a$t0Hat$vfHatk[ , k ]
               	lambda[ Y == levels(Y)[ k ] ] <- apply(a$t0Hat$vfHatkyi[ Y == levels(Y)[ k ], ], 1, sum)

		upsilon <- a$t0Hat$mfHatk[ , k ]
		upsilon[ Y == levels(Y)[ k ] ]  <- apply(a$t0Hat$mfHatkyi[ Y == levels(Y)[ k ], ], 1, sum)

		d <- a$t0Hat$mfHatk[ , k ]^2 / a$t0Hat$vfHatk[ , k ]
		d[ Y == levels(Y)[ k ] ]  <- apply(a$t0Hat$mfHatkyi[ Y == levels(Y)[ k ], ]^2 / a$t0Hat$vfHatkyi[ Y == levels(Y)[ k ], ], 1,
			function(x) sum(x, na.rm = TRUE))

		B <- c(B, -0.5 * (sum(d) - sum(a$mfik[, k ] / a$vfik[ , k ] * upsilon)))

		log_dets <- c(log_dets, - 0.5 * determinant(matrix(lambda, ncol(a$C[[ k ]]), ncol(a$C[[ k ]])) * a$C[[ k ]] + diag(ncol(a$C[[ k ]])))$modulus[[ 1 ]])
	}

	# We compute the constant that arises from the summation of the Bernoulli terms

	log_p <- NULL

        for (i in 1 : ncol(X))  {

		prodPos <- 1
		prodNeg <- 1

		for (k in (1 : length(levels(Y)))[ levels(Y) != Y[ i ] ]) {
			prodPos <- prodPos * sigmoid(a$t0Hat$pfHat[ i , k ])
			prodNeg <- prodNeg * (1 - sigmoid(a$t0Hat$pfHat[ i , k ]))
		}

		log_p <- c(log_p, log(sigmoid(a$t1Hat$pHat[ i ]) * prodPos + (1 - sigmoid(a$t1Hat$pHat[ i ])) * prodNeg))
	}

	# We add the constants that airse from the beta part of the approximation

	log_beta <- lbeta(a$a, a$b) - lbeta(a$aPrior, a$bPrior)

	a$evidence <- sum(log_ctes) + sum(B) + sum(log_dets) + sum(log_p) + log_beta

	# We compute the gradients of the log likelihood with respect to the parameters of the kernel (i.e. the latent noise)

	gradientNoise <- NULL
	gradientSigma <- NULL

	for (k in (1 : length(levels(Y)))) {

		G <- dSigmaGaussianKernel(X, a$sigma[ k ])

		lambda <- a$t0Hat$vfHatk[ , k ]
               	lambda[ Y == levels(Y)[ k ] ] <- apply(a$t0Hat$vfHatkyi[ Y == levels(Y)[ k ], ], 1, sum)

		upsilon <- a$t0Hat$mfHatk[ , k ]
		upsilon[ Y == levels(Y)[ k ] ]  <- apply(a$t0Hat$mfHatkyi[ Y == levels(Y)[ k ], ], 1, sum)

		U <- solve((matrix(lambda, ncol(a$C[[ k ]]), ncol(a$C[[ k ]])) * a$C[[ k ]]) + diag(ncol(a$C[[ k ]])))

		gradientNoise <- c(gradientNoise, - 0.5 * sum(diag(U %*% (matrix(lambda, ncol(a$C[[ k ]]), ncol(a$C[[ k ]])) * 2 * a$noise[ k ] * diag(ncol(a$C[[ k ]]))))) +
				0.5 * (upsilon %*% t(U) %*% (2 * a$noise[ k ] * diag(ncol(a$C[[ k ]]))) %*% U %*% upsilon)[ 1, 1 ])

		if (a$kernel != "linear")
			gradientSigma <- c(gradientSigma, -0.5 * sum(diag(U %*% (matrix(lambda, ncol(a$C[[ k ]]), ncol(a$C[[ k ]])) * G))) +
				0.5 * (upsilon %*% t(U) %*% G %*% U %*% upsilon)[ 1, 1])
	}

	a$gNoise <- gradientNoise
	a$gSigma <- gradientSigma
	a
}

##
# Function that computes the prediccions for new data instances
#
#' @export
predictORMCGP <- function(a, X) {

	n <- nrow(X)
	X <- t(cbind(X, rep(1, nrow(X))))
	nK <- length(levels(a$Y))
	# We compute the extended kernel data

	C <- list()
	K <- list()
	c <- list()

	for (k in 1 : nK) {

		if (is.null(a$sigma))
			C[[ k ]] <- a$kernelFun[[ k ]](cbind(a$X, X), noise = a$noise[ k ], sigma = NULL)
		else
			C[[ k ]] <- a$kernelFun[[ k ]](cbind(a$X, X), noise = a$noise[ k ], sigma = a$sigma[ k ])

		c[[ k ]] <- diag(C[[ k ]])[ (1 : n) + ncol(a$X) ]
		K[[ k ]] <- C[[ k ]][ 1 : ncol(a$X) ,(1 : n) + ncol(a$X) ]
	}

	# We compute the means and the variances of each Gaussian Process for each data instance

	means <- matrix(0, n, length(levels(a$Y)))
	variances <- matrix(0, n, length(levels(a$Y)))

	for (k in 1 : nK) {

		upsilon <- a$t0Hat$mfHatk[ , k ]
		upsilon[ a$Y == levels(a$Y)[ k ] ]  <- apply(a$t0Hat$mfHatkyi[ a$Y == levels(a$Y)[ k ], ], 1, sum)

		lambda <- a$t0Hat$vfHatk[ , k ]
               	lambda[ a$Y == levels(a$Y)[ k ] ] <- apply(a$t0Hat$vfHatkyi[ a$Y == levels(a$Y)[ k ], ], 1, sum)

		means[ , k ] <- t(K[[ k ]]) %*% a$Cinv[[ k ]] %*% a$M[[ k ]] %*% upsilon
		variances[ , k ] <- c[[ k ]] - rowSums((t(K[[ k ]]) %*% (a$Cinv[[ k ]] - a$Cinv[[ k ]] %*% a$M[[ k ]] %*% a$Cinv[[ k ]])) * t(K[[ k ]]))
	}

	# We compute the label of each instance using quadrature

	labels <- rep(a$Y[ 1 ], n)
	prob <- matrix(0, n, length(levels(a$Y)))
	colnames(prob) <- levels(a$y)

	for (i in 1 : n)  {

		maxProb <- 0

		for (k in 1 : length(levels(a$Y))) {

			pTmp <- integrate(function(u)
					apply(as.matrix(u), 1, function(u_tmp) {
					prod(pnorm((u_tmp - means[ i, -k ]) / sqrt(variances[ i, -k ]))) *
					dnorm(u_tmp, mean = means[ i, k ], sd = sqrt(variances[ i, k ]))})
				, -Inf, +Inf)$value

			prob[ i , k ] <- pTmp

			if (maxProb < pTmp) {
				labels[ i ] <- levels(a$Y)[ k ]
				maxProb <- pTmp
			}

		}
	}

	list(labels = labels, prob = prob, maxprob = apply(prob, 1, max))
}

##
# Logit and sigmoid function
#

sigmoid <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))

##
# Functions that computes a Gaussian Kernel very fast
#
# Arguments:
#		x <- matrix with the data. Each column is an instance
#		noise <- variance of the latent noise

gaussianKernel <- function(x, sigma, noise) {

	C <- matrix(0, ncol(x), ncol(x))

	for (i in 1 : nrow(x))
		C <- C + (matrix(rep(x[ i,  ], ncol(x)), ncol(x), ncol(x)) - t(matrix(rep(x[ i,  ], ncol(x)), ncol(x), ncol(x))))^2

	C <- C * -0.5 * sigma^2
	C <- exp(C)

	C <- C + diag(ncol(x)) * noise^2 + diag(ncol(x)) * 1e-3

	C
}

##
# Functions that computes the derivative of a Gaussian Kernel very fast
#
# Arguments:
#		x <- matrix with the data. Each column is an instance
#		noise <- variance of the latent noise

dSigmaGaussianKernel <- function(x, sigma) {

	C <- matrix(0, ncol(x), ncol(x))

	for (i in 1 : nrow(x))
		C <- C + (matrix(rep(x[ i,  ], ncol(x)), ncol(x), ncol(x)) - t(matrix(rep(x[ i,  ], ncol(x)), ncol(x), ncol(x))))^2

	C <- exp(C * - 0.5 * sigma^2) * -0.5 * C * 2 * sigma
	C
}

##
# Functions that computes a linear Kernel very fast
#
# Arguments:
#		x <- matrix with the data. Each column is an instance
#		noise <- variance of the latent noise

linearKernel <- function(x, noise) {
	t(x) %*% x + diag(ncol(x)) * noise^2 + diag(ncol(x)) * 1e-3
}


##
# Function that optimizes the evidence
#
# Arguments:
#
# Returns:
#	The estimated lower bound without constants

optimizeEvidenceGaussianKernelNoNoise <- function(X, Y, b = 9, verbose=F) {

	tol <- 1e-3
	a <- epORGPC(X, Y, kernel = "gaussian", noise = rep(0, length(levels(Y))), sigma = rep(1, length(levels(Y))), b = b, tol = tol)
	converged <- FALSE
	eps <- 1e-3

	while (! converged) {

		# We approximate the Snd derivative

		sigma <- a$sigma + eps * sum(a$gSigma)

		if (verbose) {
		  cat("Training with sigma: ", sigma, "\n")
		}

		anew <- epORGPC(X, Y, kernel = "gaussian", noise = rep(0, length(levels(Y))), sigma = sigma, initial_a = a, b = b, tol = tol)

		if (is.na(anew$evidence))
			return(a)

		if (abs(anew$evidence - a$evidence) < 1e-4)
			converged <- TRUE

		if (anew$evidence < a$evidence)
			eps <- eps * 0.5
		else
			eps <- eps * 1.1

		if (verbose) {
		  cat("New evidence:", anew$evidence, "Sigma:", anew$sigma, "Change:", anew$evidence - a$evidence, "Eps:", eps, "\n")
		}

		if (anew$evidence > a$evidence) {
		  a <- anew
		}
	}

	a
}

##
# Function that optimizes the evidence
#
# Arguments:
#
# Returns:
#	The estimated lower bound without constants

optimizeEvidenceGaussianKernel <- function(X, Y, b = 9, tol=1e-3, verbose=F, max.iter.num=100) {
	a <- epORGPC(X, Y, kernel = "gaussian", noise = rep(1, length(levels(Y))), sigma = rep(1, length(levels(Y))), b = b)

	converged <- FALSE
	eps <- 1e-3

	for (step in 1:max.iter.num) {

		sigma <- a$sigma + eps * sum(a$gSigma)
		noise <- a$noise + eps * a$gNoise

		if (verbose) {
		  cat("Training with sigma: ", sigma, "noise: ", noise, "\n")
		}

		anew <- epORGPC(X, Y, kernel = "gaussian", noise = noise, sigma = sigma, initial_a = a, b = b)

		if (is.na(anew$evidence))
			return(a)

		if (abs(anew$evidence - a$evidence) < tol)
			converged <- TRUE

		if (anew$evidence < a$evidence) {
		  eps <- eps * 0.5
		}
		else {
		  eps <- eps * 1.1
		}

		if (verbose) {
		  cat("New evidence:", anew$evidence, "Sigma:", anew$sigma, "Noise:", anew$noise, "Change:", anew$evidence - a$evidence,
		      "Eps:", eps, "\n")
		}

		if (converged)
		  break

		if (anew$evidence > a$evidence) {
		  a <- anew
		}
	}

	return(a)
}

##
# Function that optimizes the evidence
#
# Arguments:
#
# Returns:
#	The estimated lower bound without constants

optimizeEvidenceLinearKernel <- function(X, Y, b = 9) {

	tol <- 1e-3
	a <- epORGPC(X, Y, kernel = "linear", noise = rep(3, length(levels(Y))), b = b, tol = tol)
	converged <- FALSE
	eps <- 1e-3

	while (! converged) {

		noise <- a$noise + eps * a$gNoise

		cat("Training with noise: ", noise, "\n")

		anew <- epORGPC(X, Y, kernel = "linear", noise = noise, initial_a = a, b = b, tol = tol)

		if (is.na(anew$evidence))
			return(a)

		if (abs(anew$evidence - a$evidence) < 1e-4)
			converged <- TRUE

		if (anew$evidence < a$evidence)
			eps <- eps * 0.5
		else
			eps <- eps * 1.1

		cat("New evidence:", anew$evidence, "Noise:", anew$noise, "Change:", anew$evidence - a$evidence, "Eps:", eps, "\n")

		a <- anew
	}

	a
}
