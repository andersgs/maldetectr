library( maldetectr )
context( "Zero-Inflated Binomial" )

n_sites <- 30
n_samples <- 5
p_occupancy <- 0.2
p_detection <- 0.8
data <- sim_zib_data( n_sites = n_sites, 
                      n_samples = n_samples, 
                      p_oc = p_occupancy, 
                      p_dt = p_detection, 
                      seed = 42 )

test_that( "sim_zib_data works", {
  expect_equal( data$obs_count[1], 
                3 )
  expect_equal( data$obs_count[5], 
                5 )
  expect_equal( sum( data$occupied_sites ), 
                11 )
  expect_equal( data$n_tests[1], 
                5 )  
} )

test_that( "get_zib_likelihood works", {
  expect_equal( get_zib_likelihood( params = c( 0.5, 0.5 ),
                                    data = data ), 
                expected = -42.68476, 
                tolerance = 0.002 )
  expect_equal( get_zib_likelihood( params = c( 0.2, 0.8 ),
                                    data = data ), 
                expected = -35.20236, 
                tolerance = 0.002 )
  expect_equal( get_zib_likelihood( params = c( 0.8, 0.2 ),
                                    data = data ), 
                expected = -73.35574, 
                tolerance = 0.002 )
} )