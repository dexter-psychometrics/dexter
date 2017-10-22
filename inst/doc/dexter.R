## ---- message=FALSE------------------------------------------------------
library(dexter)

## ------------------------------------------------------------------------
head(verbAggrRules)
db = start_new_project(verbAggrRules, "verbAggression.db", covariates=list(gender="<unknown>"))

## ------------------------------------------------------------------------
head(verbAggrData)
add_booklet(db, verbAggrData, "agg")

## ------------------------------------------------------------------------
get_booklets(db)
get_items(db)

## ------------------------------------------------------------------------
data("verbAggrProperties")
head(verbAggrProperties)
add_item_properties(db, verbAggrProperties)
get_item_properties(db)
get_person_properties(db)

## ------------------------------------------------------------------------
tt = tia_tables(db)
knitr::kable(tt$itemStats, digits=3)
knitr::kable(tt$testStats, digits=3)

## ------------------------------------------------------------------------
m = fit_inter(db, booklet_id=='agg')
print(m)

## ------------------------------------------------------------------------
plot(m, "S1DoScold", show.observed=TRUE)

## ------------------------------------------------------------------------
plot(m, 'S1DoCurse', summate=FALSE)

## ------------------------------------------------------------------------
mSit = fit_domains(db, item_property= "situation")
plot(mSit)

## ------------------------------------------------------------------------
profile_plot(db, item_property='mode', covariate='gender')

## ------------------------------------------------------------------------
parms = fit_enorm(db)

## ------------------------------------------------------------------------
parms_gibbs = fit_enorm(db, method='Bayes')

## ------------------------------------------------------------------------
pv = plausible_values(db, parms)
plot(density(pv$PV1))

## ------------------------------------------------------------------------
pv = merge(pv, get_persons(db))
boxplot(PV1~gender, data=pv)

## ---- eval=FALSE---------------------------------------------------------
#  library(dexter)
#  library(dplyr)
#  library(tidyr)
#  library(readr)
#  library(SAScii)
#  
#  # Fetching the data from the OECD site requires a certain ... dexterity
#  # for which we are indebted to a kind soul on stackexchange.
#  
#  # Download the data dictionary and read it in with SAScii
#  temp = tempfile()
#  download.file("https://www.oecd.org/pisa/pisaproducts/PISA2012_SAS_scored_cognitive_item.sas", temp)
#  dict_scored = parse.SAScii(sas_ri = temp)
#  unlink(temp)
#  
#  # Download the scored cognitive data
#  temp = tempfile()
#  download.file("https://www.oecd.org/pisa/pisaproducts/INT_COG12_S_DEC03.zip",temp, mode="wb")
#  unzip(temp, "INT_COG12_S_DEC03.txt")
#  scored = read_fwf(file = 'INT_COG12_S_DEC03.txt',
#                    col_positions = fwf_widths(dict_scored$width), progress = TRUE)
#  colnames(scored) = dict_scored$varname
#  unlink(temp)
#  
#  # Keep only the maths booklets and items
#  pisa12_M = scored %>%
#    filter(BOOKID %in% 1:13) %>%
#    select(CNT, BOOKID, starts_with('PM'))
#  
#  rm(scored)
#  
#  # Items missing by design are coded with 7, and actual nonresponse with 8
#  # Will replace both with NA for simplicity.
#  pisa12_M$BOOKID = paste0('B',pisa12_M$BOOKID)
#  pisa12_M[pisa12_M==7] = NA
#  pisa12_M[pisa12_M==8] = NA
#  
#  # remove columns that contain only NA values
#  pisa12_M = select_if(pisa12_M, function(x) !all(is.na(x)))
#  
#  mathItems = grep("PM",names(pisa12_M), value=TRUE)
#  
#  # prepare the scoring rules
#  # unique combinations of items and responses
#  # code NA as 0 and ontherwise make score equal to response (since we have scored data)
#  rules = gather(pisa12_M, key='item_id', value='response', starts_with('PM')) %>%
#    distinct(item_id, response) %>%
#    mutate(item_score = ifelse(is.na(response), 0, response))
#  
#  
#  # create the new project
#  db = start_new_project(rules, "pisa_math.db", covariates=list(cnt='<unknown country>'))
#  
#  # add all booklets one by one, deleting columns that may be all NA
#  pisa12_M %>%
#    group_by(BOOKID) %>%
#    do({
#      rsp = select_if(., function(x) !all(is.na(x)))
#      booklet_id = .$BOOKID[1]
#      add_booklet(db, rsp, booklet_id)
#      # return an empty data.frame since we're using the do construct only for its side-effects
#      data.frame()
#    })
#  
#  rm(pisa12_M)
#  
#  
#  # add some item properties -- we have supplied a data set to
#  # make things a bit easier
#  items = get_items(db)
#  
#  domain = merge(items, PISA_item_class, by.x="item_id", by.y="ItemCode") %>%
#    select(item_id, Content) %>%
#    mutate(isSaS = ifelse(Content=="Space and shape",'SaS','notSas'))
#  
#  add_item_properties(db, domain)
#  
#  # an overview of person and item properties in the project
#  get_item_properties(db)
#  get_person_properties(db)
#  
#  # Fit the interaction model for booklet 1, all countries
#  m <- fit_inter(db,booklet_id=='B1')
#  plot(m)
#  
#  # Analyse by domain
#  md = fit_domains(db, 'content', booklet_id=='B1')
#  plot(md, nr=2, nc=2)
#  
#  # Compare three countries on 'Space and shape' vs NOT 'Space and shape' in booklet 1
#  profile_plot(db, item_property = "isSaS", covariate="cnt", predicate=(cnt %in% c("DEU","FRA","ITA") & booklet_id=='B1'))
#  

## ---- show=FALSE---------------------------------------------------------
dbDisconnect(db)

