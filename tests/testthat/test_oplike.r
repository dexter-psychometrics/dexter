context('test oplike')


test_that('start_new_project_from_oplm',
{
  skip_on_cran()
  if(dir.exists('../skip_on_cran'))
  {
    # user incorrectly assumes booklet vector is (start,length) or just start
    # expect errors mentioning booklet_position
    expect_error(
    {
      db = start_new_project_from_oplm(dbname=':memory:',
                                      	scr_path='../skip_on_cran/duits/DUITS.scr',
                                      	dat_path='../skip_on_cran/duits/DUITS.DAT',
                                      	booklet_position=c(128,1),
                                  	    response_length = 1,
                                  	    responses_start = 131)
    }, 'booklet_position' )
    
    expect_error(
      {
        db = start_new_project_from_oplm(dbname=':memory:',
                                         scr_path='../skip_on_cran/duits/DUITS.scr',
                                         dat_path='../skip_on_cran/duits/DUITS.DAT',
                                         booklet_position = 128,
                                         response_length = 1,
                                         responses_start = 131)
      }, 'booklet_position' )
    
    # should give invalid booklets
    expect_error({
    db = start_new_project_from_oplm(dbname=':memory:',
                                     scr_path='../skip_on_cran/duits/DUITS.scr',
                                     dat_path='../skip_on_cran/duits/DUITS.DAT',
                                     booklet_position = c(12,14),
                                     responses_start = 19,
                                     use_discrim=TRUE)
    },'booklet id')
    expect_warning({
    db = start_new_project_from_oplm(dbname=':memory:',
                                     scr_path='../skip_on_cran/duits/DUITS.scr',
                                     dat_path='../skip_on_cran/duits/DUITS.DAT',
                                     booklet_position = c(3,4),
                                     responses_start = 19,
                                     use_discrim=TRUE)},
    'too short')
    
    close_project(db)
    # to do, add test with proper oplm file
    
  } 
})
