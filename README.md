<!-- README.md is generated from README.Rmd. Please edit that file -->
Dexter
======

Dexter is an R package for psychometric analysis of data from educational and psychological tests. Dexter typically works with project database files saved on disk.

Installation
------------

``` r
install.packages('dexter')
```

If you encounter a bug, please post a minimal reproducible example on [github](https://github.com/jessekps/dexter/issues). We post news and examples on a [blog](http://dexterities.netlify.com), it's also the place for general questions.

Example
-------

``` r
library(dexter)
# start a project and fill it with data
# verbAggrRules and verbAggrData are example datasets provided with dexter

db = start_new_project(verbAggrRules, "verbAggression.db")
add_booklet(db, verbAggrData, booklet_id = "verb_agg")

# Classical test theory
tia = tia_tables(db)

tia$testStats
```

| booklet\_id |  nItems|  alpha|  meanP|  meanRit|  meanRir|  maxTestScore|    N|
|:------------|-------:|------:|------:|--------:|--------:|-------------:|----:|
| verb\_agg   |      24|  0.888|  0.339|    0.527|    0.468|            48|  316|

``` r
head(tia$itemStats)
```

| booklet\_id | item\_id    |  meanScore|  sdScore|  maxScore|  pvalue|    rit|    rir|    n|
|:------------|:------------|----------:|--------:|---------:|-------:|------:|------:|----:|
| verb\_agg   | S1DoCurse   |      1.082|    0.808|         2|   0.541|  0.582|  0.519|  316|
| verb\_agg   | S1DoScold   |      0.832|    0.817|         2|   0.416|  0.651|  0.596|  316|
| verb\_agg   | S1DoShout   |      0.468|    0.710|         2|   0.234|  0.520|  0.460|  316|
| verb\_agg   | S1WantCurse |      1.123|    0.828|         2|   0.562|  0.537|  0.468|  316|
| verb\_agg   | S1WantScold |      0.930|    0.852|         2|   0.465|  0.593|  0.528|  316|
| verb\_agg   | S1WantShout |      0.712|    0.778|         2|   0.356|  0.529|  0.464|  316|

``` r
# IRT, extended nominal response model
f = fit_enorm(db)

head(coef(f))
```

| item\_id  |  item\_score|        beta|      SE\_b|
|:----------|------------:|-----------:|----------:|
| S1DoCurse |            1|  -1.3422140|  0.1541565|
| S1DoCurse |            2|  -0.6375015|  0.1418423|
| S1DoScold |            1|  -0.6702036|  0.1429057|
| S1DoScold |            2|  -0.2589855|  0.1579467|
| S1DoShout |            1|   0.3254326|  0.1480166|
| S1DoShout |            2|   0.3687574|  0.2099654|

``` r
# ability estimates per person
abl = ability(db, parms = f)
head(abl)
```

| booklet\_id | person\_id |  sumScore|       theta|
|:------------|:-----------|---------:|-----------:|
| verb\_agg   | dxP1       |        13|  -1.0218323|
| verb\_agg   | dxP10      |         9|  -1.4832759|
| verb\_agg   | dxP100     |        14|  -0.9204799|
| verb\_agg   | dxP101     |         0|        -Inf|
| verb\_agg   | dxP102     |        12|  -1.1276084|
| verb\_agg   | dxP103     |         8|  -1.6191967|

``` r
# ability estimates without item S1DoScold
abl2 = ability(db, parms = f, item_id != "S1DoScold")

# plausible values

pv = plausible_values(db, parms = f, nPV = 5)
head(pv)
```

| booklet\_id | person\_id |  sumScore|         PV1|         PV2|         PV3|         PV4|         PV5|
|:------------|:-----------|---------:|-----------:|-----------:|-----------:|-----------:|-----------:|
| verb\_agg   | dxP1       |        13|  -0.8830828|  -0.7066814|  -0.7703465|  -1.0441590|  -0.8918823|
| verb\_agg   | dxP10      |         9|  -1.5388849|  -1.0014150|  -1.7249607|  -1.3372498|  -1.5627987|
| verb\_agg   | dxP100     |        14|  -0.6678893|  -1.4039775|  -0.6161149|  -0.9249503|  -0.7607182|
| verb\_agg   | dxP101     |         0|  -2.3034185|  -3.0623911|  -2.3793639|  -3.1658178|  -3.3999772|
| verb\_agg   | dxP102     |        12|  -1.0434515|  -1.5317991|  -1.1518701|  -1.0103475|  -1.0057779|
| verb\_agg   | dxP103     |         8|  -2.0467553|  -1.5357100|  -2.0598186|  -1.2366479|  -1.7367032|

Contributing
------------

Contributions are welcome but please check with us first about what you would like to contribute.
