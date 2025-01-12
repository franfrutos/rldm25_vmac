# Función para generar valores iniciales
generate_initial_values <- function(chain_id, J, model = "1", chains = 4) {
  out <- vector("list", chains) 
  for (i in 1:chains) {
    if (model == "1") {
      out[[i]] <- list(
        mu_intercept = rnorm(1, 0, 1),        # Prior ~ N(0, 1)
        mu_beta = rnorm(1, 0, 1),            # Prior ~ N(0, 1)
        mu_beta2 = rnorm(1, 0, 1),            # Prior ~ N(0, 1)
        alpha_0 = rnorm(1, 0, 1),            # Prior ~ N(0, 1)
        mu_lr = rnorm(1, 0, 1),              # Prior ~ N(0, 1)
        sigma = abs(rnorm(5, 0, 0.5)),       # Prior ~ N(0, 1) (truncated)
        z = matrix(rnorm(5 * J, 0, 1), nrow = 5, ncol = J) # Prior ~ N(0, 1)
      )
    } else {
      out[[i]] <- list(
        mu_intercept = rnorm(1, 0, 1),        # Prior ~ N(0, 1)
        mu_beta = rnorm(1, 0, 1),            # Prior ~ N(0, 1)
        mu_beta2 = rnorm(1, 0, 1),            # Prior ~ N(0, 1)
        mu_lr = rnorm(1, 0, 1),              # Prior ~ N(0, 1)
        mu_gamma = rnorm(1, 0, 1),           # Prior ~ N(0, 1)
        alpha_0 = rnorm(1, 0, 1),            # Prior ~ N(0, 1)
        sigma = abs(rnorm(6, 0, 1)),       # Prior ~ N(0, 1) (truncated)
        z = matrix(rnorm(6 * J, 0, 1), nrow = 6, ncol = J) # Prior ~ N(0, 1)
      )
    }
  }
  
  return(out)
}

# Average PPC over participants 
get_PPC <- function(fit, group_size = 512) {
  preds <- fit$draws(variables = "y_rep", format = "matrix")
  tmp <- vector("list", length = 2)
  tmp[[1]] <- preds[,1:(ncol(preds)/2)]
  tmp[[2]] <- preds[,((ncol(preds)/2)+1):ncol(preds)]

  out <- data.frame()
  
  for (i in tmp) {
    N <- nrow(i)       # Todas las filas
    block_size <- group_size  # Número de columnas por submatriz
    
    # Dividir en submatrices por columnas
    num_blocks <- ceiling(ncol(i) / block_size)
    
    # Crear una matriz tridimensional
    sub <- array(i, dim = c(N, block_size, num_blocks))
    
    
    preds <- apply(sub, c(1,2), mean)
    
    yMean <- colMeans(preds)
    HDI <- apply(preds, 2, HDIofMCMC)
    out <- rbind(out, data.frame(
      Mean = yMean,
      Low = HDI[1,],
      High = HDI[2,]
    ))
  }
  
  rm(tmp)
  return(out)
}

# Functions to create epochs of blocks:
create_epochs <- function(blocks, epoch = 2) {
  vapply(blocks, function(x, e = epoch) {
    ceiling(x / e)
  }, FUN.VALUE = numeric(1))
}

HDIofMCMC = function(sampleVec, 
                     credMass = 0.89 ) {
  if ( class(sampleVec) == "mcmc.list" ) {
    sampleVec = as.matrix(sampleVec)
  }
  sortedPts = sort( sampleVec )
  ciIdxInc = floor( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}
