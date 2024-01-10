context('check latent correlations')

library(dplyr)
library(tidyr)

RcppArmadillo::armadillo_throttle_cores(1)

test_that('latent correlations work',{
  
  db = open_project('../verbAggression.db')

  lt = latent_cor(db,'behavior')
  x = get_responses(db, columns=c('person_id','item_id','item_score','behavior')) |>
    group_by(person_id,behavior) |>
    summarise(score=sum(item_score)) |>
    ungroup() |>
    pivot_wider(names_from=behavior, values_from=score)
    
  xcor = cor(select(x,-person_id))
  ut = upper.tri(lt$cor)
  
  expect_true(all(lt$cor[ut] > xcor[ut]), 'latent correlations should be larger than score correlations')
  expect_true(all(lt$cor[ut] > lt$hpd_l[ut]), 'hpd l lower than cor')
  expect_true(all(lt$cor[ut] < lt$hpd_h[ut]), 'hpd h higher than cor')
  
  skip_on_cran()
  # may fail in rare cases due to random nbr uncertainty
  expect_equal((lt$cor + 1.97 * lt$sd)[ut], lt$hpd_h[ut], tolerance=.02, info='upper hd ~ 1.97*sd+cor')
  
  #to do: incomplete designs
})

RcppArmadillo::armadillo_reset_cores()
