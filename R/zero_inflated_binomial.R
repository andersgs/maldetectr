#' Zero-Inflated Binomial Log-Likelihood
#' 
#' \code{get_zib_likelihood} retuns the log-likelihood of observed counts of 
#' detection at \code{n_sites} tested \code{n_tests} times each, given a certain
#' probability of occupancy and conditional probability of detection given site is
#' occupied.
#' 
#' The Likehood function is described in equation (2) of Royle and Nichols (2003).
#' 
#' Royle, J. A., & Nichols, J. D. (2003). Estimating abundance from repeated 
#' presence–absence data or point counts. Ecology, 84(3), 777–790. 
#' http://doi.org/10.1890/0012-9658(2003)084[0777:EAFRPA]2.0.CO;2
#' 
#' @param params Numeric vector of length 2, where the first value is the probability of occupancy and the second is the probability of detection conditional on occupancy
#' @param data Data.frame with at least two columns: one named obs_counts with the observed counts at each site, and once called n_tests with the number of tests per sites
#' 
#' @return The log-likelihood of the data given the parameter values under a model of Zero-Inflated Binomial
#' 
#' @examples 
#' data <- maldetectr::sim_zib_data( n_sites = 10, n_samples = 5, p_oc = 0.2, p_dt = 0.8 )
#' maldetectr::get_zib_likelihood( params = c( 0.5, 0.5 ), data = data )
#' 
#' @export

get_zib_likelihood <- function( params, data ){
  p_occupied <- params[1]
  p_detection <- params[2]
  log_p_occupied <- log( p_occupied )
  p_no_detection <- 1 - p_detection
  p_not_occupied <- 1 - p_occupied
  res <- data %>%
    dplyr::mutate( L = ifelse( obs_count > 0, 
                               {dbinom( x = obs_count, 
                                        size = n_tests, 
                                        prob = p_detection ) * p_occupied },
                               {p_occupied * ( p_no_detection )^n_tests + 
                                   p_not_occupied } )
                   ) %>%
    dplyr::summarise( total_L = prod( L ) )
  return( log( res$total_L ) )
}

#' Estimate Probability of Occupancy and Conditional Probability of Detection under a Zero-Inflated Binomial model
#' 
#' Maximum likelihood estimates for probability of occupancy and probability of
#' detection coditional on site being occupied using the \code{get_zib_likelihood}
#' function.
#' 
#' @param data Data.frame with at least two columns: one named obs_counts with the observed counts at each site, and once called n_tests with the number of tests per sites
#' @param ci Value for confidence interval (between 0 and 1)
#' 
#' @return A data.frame with maximum likelihood estimates for probability of occupancy and probability of detection given site is occupied, plus confidence interval around the estimates.
#' 
#' @seealso \code{\link{get_zib_likelihood}}
#' 
#' @examples 
#' data <- maldetectr::sim_zib_data( n_sites = 10, n_samples = 5, p_oc = 0.2, p_dt = 0.8 )
#' maldetectr::get_zib_ml( data = data, ci = 0.95 )
#' 
#'@export 

get_zib_ml <- function( data, ci = 0.95 ) {
  z_score <- 1 - ( 1 - ci )/2
 res <-  optim( c( 0.5, 0.5 ), 
       fn = maldetectr::get_zib_likelihood, 
       data = data, 
       control = list( fnscale = -1 ), 
       hessian = T )
 stderr <- qnorm( z_score ) * diag( solve( -res$hessian ) )
 output = data.frame( ML = res$par )
 output <- output %>%
   dplyr::mutate( lower_CI = ML - stderr, 
                  upper_CI = ML + stderr )
 row.names( output ) <- c( "P_Occupancy", "P_Detection")
 return( output )
}

#' Simulate zero-inflated binomial data
#' 
#' @param n_sites Number of sites to simulate
#' @param n_samples Number of observations at each site
#' @param p_oc Probability of occupancy at site
#' @param p_dt Probability of detection at site conditional on site being occupied
#' @param seed Set the seed for reproduciability
#' 
#' @return A data.frame of simulated data with \code{n_sites} rows and three columns. An \code{occupied_sites} column with 0 if site is not occupied and 1 if occupied; an \code{n_tests} column, which indicates the number of times each site was sampled (all values in this column will be equal to \code{n_samples}, and will be used in calculating the likelihood); and an \code{obs_counts} column, which has the number of times a detection occurred (values will be 0 to a maximum of \code{n_smaples})
#' @examples
#' maldetectr::sim_zib_data( n_sites = 10, n_samples = 5, p_oc = 0.2, p_dt = 0.8 )
#' @export 
sim_zib_data <- function( n_sites, n_samples, p_oc, p_dt, seed = 42 ) {
  set.seed( seed )
  res <- data.frame( occupied_sites = rbinom( n = n_sites, 
                                              size = 1, 
                                              prob = p_oc ),
                     n_tests = n_samples )
  res_pos <- res %>%
    dplyr::filter( occupied_sites == 1 ) %>%
    dplyr::mutate( obs_count = rbinom( n = n(), 
                                    size = n_samples, prob = p_dt ) )
  res_neg <- res %>%
    dplyr::filter( occupied_sites == 0 ) %>%
    dplyr::mutate( obs_count = 0 )
  return( dplyr::bind_rows( res_pos, res_neg ) )
}
