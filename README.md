<!-- README.md is generated from README.Rmd. Please edit that file -->
Dexter
======

Dexter is an R package for and psychometric analysis of data from educational and psychological tests. Dexter typically works with project database files saved on disk.

Installation
------------

``` r
install.packages('dexter')
```

If you encounter a bug, please post a minimal reproducible example on [github](https://github.com/jessekps/dexter/issues). We post news and examples on a [blog](dexterities.netlify.com), it's also the place for general questions.

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

| booklet\_id | person\_id |  sumScore|         PV1|         PV2|        PV3|        PV4|         PV5|
|:------------|:-----------|---------:|-----------:|-----------:|----------:|----------:|-----------:|
| verb\_agg   | dxP1       |        13|  -1.0214890|  -0.9521324|  -1.252551|  -1.669846|  -0.9245417|
| verb\_agg   | dxP10      |         9|  -1.2626779|  -1.2457954|  -1.398498|  -1.378339|  -1.5915501|
| verb\_agg   | dxP100     |        14|  -1.3186758|  -0.7754019|  -1.266446|  -1.222720|  -0.6871353|
| verb\_agg   | dxP101     |         0|  -3.2277562|  -2.9269698|  -3.819951|  -3.402195|  -3.4895825|
| verb\_agg   | dxP102     |        12|  -1.2779256|  -0.7690105|  -1.132858|  -1.274146|  -0.7648520|
| verb\_agg   | dxP103     |         8|  -0.7125214|  -1.2058239|  -1.366199|  -1.396237|  -1.2682042|

Contributing
------------

Contributions are welcome but please check with us first about what you want to contribute to prevent disappointment. Please also note that dexter uses CML estimation and we will not include models that require MML. There are already other R packages for that, most notably [mirt](https://cran.r-project.org/package=mirt).
