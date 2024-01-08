# dexter 1.3.2

## breaking changes

* A small change in the interface for the plausible_values and plausible_scores functions, the slightly different argument `parms_draw` replaces `use_draw`.

* The argument `standard_errors` in functions `ability` and `ability_tables` has been removed. Standard errors are now always returned.

* argument `prior.dist` in plausible values has become `prior_dist` to maintain a consistent naming scheme in dexter.

* Internal changes in the plausible values function will cause a slight increase in confidence intervals for population estimates when using multiple draws. We think this is more realistic. On the other hand, population estimates based on plausible values will now be more stable over different seeds or consecutive calls of `plausible_values`. Use of the same seed no longer guarantees the exact same outcomes across different computers.

## other changes

* more graphical options in `plot.prms`

# dexter 1.2.2

* changed license from GPL-3 to LGPL-3
* tia_tables optionally includes statistics for distractors for MC questions
* colors for distractor plots can be specified per response
* updated error messages
* bugfix: the confidence intervals in the fit plots for the enorm were drawn too wide for polytomous items in dexter <= 1.2.1. They are now correct.

# dexter 1.2.1

## breaking change

* an argument 'design' has been added to add_response_data. This obligates the user to explicitly specify the design rather than having it inferred from the data, as was the old behavior. The old behavior led to a lot of confusion, nevertheless it can be replicated by using `design=distinct(x,booklet_id,item_id)` as the design argument, but it is better practice to use design information from a separate source to prevent errors.  

## other changes

* adapted to comply with a change in the upcoming dplyr version. 


# dexter 1.1.5

* slight change in implementation of DIF. If some items are not administered in both groups, they are now removed after calibration instead of before, making the calibration more robust in incomplete designs.
* ability for EAP with normal prior is faster and SE's are more precise, argument nPV is removed
* deprecated objects in result of tia_tables have now been removed.


# dexter 1.1.4

* changed sd_score for items in tia_tables() for type='averaged' to the actual sd rather than a weighted average of sd over booklets (which was the case in previous versions) 
* tia_tables now returns a list with data.frame elements booklets and items with snake case column names rather than elements testStats and itemStats. testStas and itemStats are also still included for the time being to not break existing code.
* corrected a minor bug which caused progress bars to hang or not complete to 100%
* get_resp_matrix no longer includes empty rows for missing factor levels

# dexter 1.1.3

* new website, see dexter-psychometrics.github.io/dexter
* new function `latent_cor` which does what the name implies
* speed improvement in `probability_to_pass` 
* bugfix in `pobability_to_pass` which occurred when not all scores were achieved

# dexter 1.1.2

* update per request from cran, no changes in functionality

# dexter 1.1.1

* profiles/profile_tables can now take a data.frame of item parameters
* bugfix: mistakenly renormalized parameters when using a data.frame of polytomous item parameters 

# dexter 1.0.8

* added extra argument to `coef.prms` to retrieve variance-covariance matrix or posterior of
item parameters
* bugfix: x-axis scaled incorrectly when plotting parameters for Bayesian calibration
* update examples for updated cran check in 3.4

# dexter 1.0.7

* adapting to an update in tibble

# dexter 1.0.6

* specific workarounds for (cran) solaris

# dexter 1.0.5

* Almost sure this release finally fixes the solaris issues

# dexter 1.0.4

* added systemrequirements to description in the hope that it will solve solaris issues

# dexter 1.0.3

* fixed bug that caused an error message on some apple computers when using a matrix datasource

* removed unit tests for Solaris as we have no way to test on solaris, in effect this means dexter does not support solaris anymore

# dexter 1.0.2

* solved bug: accidental change of input data.frame in ability

* solved bug: unnecessary error message when input parameters are data.frame

# dexter 1.0.1

* weighted likelihood estimate option in `ability()` and `ability_tables()`

* added function `r_score_IM()` for simulation of item responses according to the interaction model

* added function p_score which provides conditional score distributions

* ability, plausible_values and plausible_scores accept a tidy data.frame of parameters

* it is now possible to use `get()` in predicates

* bugfix for mixed up plots in the interaction model when using a matrix dataSrc

* bugfix: reinstated check for unobserved zero score categories in fit_enorm (check was accidentally removed in 1.0)

* bugfix: uppercase letters in variable names in profile_plot work again



# dexter 1.0.0

## Breaking changes

* the interfaces for the standard setting functions have changed, these functions have new names and different arguments

* the interface for probability_to_pass has changed

* the name of the column 'sumScore', outputted by a number of functions, has changed to 'booklet_score'

## other changes

* dataSrc may now also be a wide format data matrix (in addition to a database or a long format data.frame)

* fit_enorm, ability and plausible_values have an extra argument 'merge_within_persons' which determines if different booklets for the same person are combined

* plausible_values can use a normal or a mixture prior

* the dependency on dbplyr has been removed. This again allows the use of data.frame columns in predicates (this caused error messages when using dbplyr 1.4+)

* generally faster due to the excellent Rcpp and RcppArmadillo packages


# dexter 0.8.5

* new function `design_info()` returns extensive information about incomplete test designs. Functions `design_as_network()` and `design_is_connected()` are deprecated.

* correction for a bug which caused NA's in plausible values for booklets with 1 respondent and nPV>1

# dexter 0.8.4

* Tiny changes for compatibility with dependencies

# dexter 0.8.3

* correction of a bug in 0.8.2 which caused coef.prms to fail with calibration method Bayes

# dexter 0.8.2

* new functions `information()` and `expected_score()` for tests and items in the ENORM.

* change in the algorithm for MLE estimation method in ability and ability_tables. The new algorithm is
  more precise so small changes in (the lower digits of) computed theta values may occur with respect 
  to previous versions of dexter.
  
* corrected a bug in the reporting of standard errors for Bayesian estimation of the ENORM. The
  SE's reported in previous versions of dexter were too high.
  
* corrected bug in plot for DIF (item labels were reversed) introduced in 0.8.0


# dexter 0.8.1

* significant speed increases in plausible_values, plausible_scores, ability and ability_tables

* added plot method for the extended nominal response model

* new function `add_response_data()` for importing normalized data

# dexter 0.8.0

* `profiles()` and `profile_tables()` functions added for analysis of individual deviation profiles
  (misfit) compared to the fitted model.

* bugfix for solaris

* You may want to take a look at the new sister package, dextergui


# dexter 0.7.0

* support for RSQLite release 2.1.0+ which broke dexter 0.6.0. 

* `probability_to_pass()` function added for equating to a reference test that is 
  connected to a target test.


# dexter 0.6.0

* `ability_tables()` function added to generate score transformation tables 
  for arbitrary subsets of items. `ability()` and `ability_tables()` can now also 
  take Jeffreys prior and variatons of the normal prior and can optionally 
  output standard errors.
  
* `plausible_scores()` function added to generate scores on arbitrary itemsets
  (or the entire bank). 
  
* Updates to documentation.


# dexter 0.5.1

* `plausible_values()` can now work with covariates. 

* `fit_enorm()` can optionally take fixed_parameters into account. 

* Significant speed increases in profile_plot, iTIA and iModels and minor speed increases in other functions.

* Further includes more user control over plotting parameters, updates to documentation, several bugfixes and some
  smaller updates.


# dexter 0.4.4

* included Rdpack in Imports at the request of CRAN, and fixed 4 minor bugs

# dexter 0.4.2

We go straight from version 0.1.7 to version 0.4.2 in a major new release.
The functionality envisaged for two other related packages, enorm and 
roger, has now been all incorporated into dexter.

* The database has been completely rewritten. Databases created with 0.1.7 can
  be imported but are no longer maintained.
  
* Many functions now support flexible subsetting of data through a predicate
  expression; furthermore, the data source can be either the dexter database
  (more usual), or a data frame or tibble. The latter can be useful in simple
  simulations or in some non-standard situations.
  
* Basic test and item analysis has been extended with two interactive tools
  based on shiny, and with an exploratory tool to search for DIF. 

* Version 0.4.2 estimates the Extended NOminal Response Model (ENORM) by
  either CML or a Gibbs sampler, computes person scores by either MLE, EAP,
  or plausible values, and provides a new test for the hypothesis of no true
  individual differences. In addition, it supports a method for standard
  settings known as 3DC. 

* Several vignettes have been added to explain the new features in detail.
