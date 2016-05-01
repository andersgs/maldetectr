#' Zero-Inflated Binomial
#' 
#' 

detect_zib <- function( params, data ){
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

p_oc <- 0.2
p_detec <- 0.8
n_sites <- 50
n_samples <- 5
set.seed( 42 )

generate_test_data <- function( n_sites, n_samples, p_oc, p_dt ) {
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

test_data <- generate_test_data( n_sites, n_samples, p_oc, p_detec )

detect_zib( params = c( 0.2, 0.8 ), data = test_data )

optim( c( 0.5, 0.5 ), 
       fn = detect_zib, 
       data = test_data, 
       control = list( fnscale = -1 ), 
       hessian = T )
