Sys.setenv("R_TESTS" = "")
library(testthat)
library(dexter)
test_check("dexter")
# to do, replace skip_on_cran tests in: 
# test_ability (uses pisa, can easily just simulate a large dataset)
# test_interaction (find out what the special data condition was)
# test_oplike (need scr and dat, possible to randomly change the data I guess)
