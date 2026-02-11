#' Generate Fractional Brownian Motion
#' 
#' Generate Fractional Brownian Motion using the FFT method with long-range dependence.
#' 
#' @param H Hurst parameter (0 < H < 1). H > 0.5 indicates positive long-range dependence.
#' @param n Number of grid points
#' @param T Total time length
#' @param seed Random seed for reproducibility
#' 
#' @return A vector of length n+1 containing the FBM process
#' 
#' @details
#' This function implements the fast Fourier transform (FFT) method for generating
#' fractional Brownian motion with specified Hurst parameter H.
#' 
#' @examples
#' # Generate FBM with H = 0.7 (positive long-range dependence)
#' fbm_process <- fbm(H = 0.7, n = 1000, T = 1, seed = 123)
#' plot(fbm_process, type = "l", main = "Fractional Brownian Motion (H=0.7)")
#' 
#' @export
fbm <- function(H, n, T, seed = NULL) {
  
  # Input validation
  if (missing(H) || missing(n) || missing(T)) {
    stop("H, n, and T must be specified")
  }
  
  if (H <= 0 || H >= 1) {
    stop("Hurst parameter H must be between 0 and 1")
  }
  
  if (n <= 0) {
    stop("Number of grid points n must be positive")
  }
  
  if (T <= 0) {
    stop("Time length T must be positive")
  }
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- as.integer(n)
  
  if (n < 2) {
    n <- 2
  }
  
  # Generate covariance sequence
  r <- numeric(n + 1)
  r[1] <- 0
  
  # Calculate covariance using the exact formula
  for (i in 1:n) {
    r[i + 1] <- 0.5 * ((i + 1)^(2 * H) - 2 * i^(2 * H) + (i - 1)^(2 * H))
  }
  
  # Symmetric extension for FFT
  r_sym <- c(r, r[seq(n, 2, by = -1)])
  
  # Calculate eigenvalues using FFT
  lambda <- Re(stats::fft(r_sym)) / (2 * n)
  
  # Ensure eigenvalues are non-negative
  lambda <- pmax(lambda, 0)
  
  # Generate complex Gaussian random variables
  real_part <- stats::rnorm(2 * n)
  imag_part <- stats::rnorm(2 * n)
  complex_vec <- complex(real = real_part, imaginary = imag_part)
  
  # FFT transform to generate the process
  W <- stats::fft(sqrt(lambda) * complex_vec)
  W <- n^(-H) * cumsum(Re(W[1:(n + 1)]))
  W <- T^H * W
  
  return(W)
}