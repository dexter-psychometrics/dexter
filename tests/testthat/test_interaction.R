context('test interaction')

library(dplyr)


expect_no_error = function(object, info=NULL) expect_error(object, regexp=NA, info=info)


# a dataset that generated an error before
test_that('fit_inter gives no error',{
  skip_on_cran()
  if(dir.exists('skip_on_cran'))
  {
    db = open_project('skip_on_cran/eva/670.db')
    
    expect_no_error({f=fit_inter(db,booklet_id=='1648')}, 
                    info='check that fit_inter calibration error in 670 does not occur')
    

    dbDisconnect(db)
  }
})