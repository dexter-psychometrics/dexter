<!-- README.md is generated from README.Rmd. Please edit that file -->

Dexter
======

Dexter is an R package for psychometric analysis of data from
educational and psychological tests. Dexter typically works with project
database files saved on disk.

Installation
------------

``` r
install.packages('dexter')
```

If you encounter a bug, please post a minimal reproducible example on
[github](https://github.com/dexter-psychometrics/dexter/issues). We post
news and examples on a [website and
blog](https://dexter-psychometrics.github.io/dexter).

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
| verb\_agg   | S1DoCurse   |      1.082|    0.807|         2|   0.541|  0.582|  0.519|  316|
| verb\_agg   | S1DoScold   |      0.832|    0.815|         2|   0.416|  0.651|  0.596|  316|
| verb\_agg   | S1DoShout   |      0.468|    0.709|         2|   0.234|  0.520|  0.460|  316|
| verb\_agg   | S1WantCurse |      1.123|    0.827|         2|   0.562|  0.537|  0.468|  316|
| verb\_agg   | S1WantScold |      0.930|    0.850|         2|   0.465|  0.593|  0.528|  316|
| verb\_agg   | S1WantShout |      0.712|    0.777|         2|   0.356|  0.529|  0.464|  316|

``` r
# IRT, extended nominal response model
f = fit_enorm(db)

head(coef(f))
```

| item\_id  |  item\_score|        beta|   SE\_beta|
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

| booklet\_id | person\_id  |  booklet\_score|       theta|
|:------------|:------------|---------------:|-----------:|
| verb\_agg   | dx\_0000001 |              13|  -1.0238738|
| verb\_agg   | dx\_0000002 |              28|   0.3124831|
| verb\_agg   | dx\_0000003 |               4|  -2.3748882|
| verb\_agg   | dx\_0000004 |              19|  -0.4630604|
| verb\_agg   | dx\_0000005 |               7|  -1.7721275|
| verb\_agg   | dx\_0000006 |              25|   0.0512826|

``` r
# ability estimates without item S1DoScold
abl2 = ability(db, parms = f, item_id != "S1DoScold")

# plausible values

pv = plausible_values(db, parms = f, nPV = 5)
head(pv)
```

| booklet\_id | person\_id  |  booklet\_score|         PV1|         PV2|         PV3|         PV4|         PV5|
|:------------|:------------|---------------:|-----------:|-----------:|-----------:|-----------:|-----------:|
| verb\_agg   | dx\_0000001 |              13|  -1.0616680|  -0.9219415|  -0.8161716|  -0.8717352|  -0.6656054|
| verb\_agg   | dx\_0000002 |              28|  -0.1858848|  -0.0366271|   0.3216584|   0.0263411|   0.3838943|
| verb\_agg   | dx\_0000003 |               4|  -1.9804484|  -2.4675297|  -2.5928338|  -2.1755638|  -1.7619598|
| verb\_agg   | dx\_0000004 |              19|  -0.2550607|  -0.4989645|  -0.4322568|  -0.2495976|  -0.3305588|
| verb\_agg   | dx\_0000005 |               7|  -2.4189426|  -1.5473985|  -1.2680878|  -1.3221576|  -1.3139852|
| verb\_agg   | dx\_0000006 |              25|   0.3409678|  -0.2472593|  -0.0050744|  -0.3793325|  -0.2298762|
