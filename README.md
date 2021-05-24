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
blog](https://dexter-psychometrics.github.io/dexter/).

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

tia$booklets
```

<table>
<colgroup>
<col style="width: 13%" />
<col style="width: 9%" />
<col style="width: 7%" />
<col style="width: 14%" />
<col style="width: 10%" />
<col style="width: 10%" />
<col style="width: 21%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">booklet_id</th>
<th style="text-align: right;">n_items</th>
<th style="text-align: right;">alpha</th>
<th style="text-align: right;">mean_pvalue</th>
<th style="text-align: right;">mean_rit</th>
<th style="text-align: right;">mean_rir</th>
<th style="text-align: right;">max_booklet_score</th>
<th style="text-align: right;">n_persons</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">verb_agg</td>
<td style="text-align: right;">24</td>
<td style="text-align: right;">0.888</td>
<td style="text-align: right;">0.339</td>
<td style="text-align: right;">0.527</td>
<td style="text-align: right;">0.468</td>
<td style="text-align: right;">48</td>
<td style="text-align: right;">316</td>
</tr>
</tbody>
</table>

``` r
head(tia$items)
```

| booklet\_id | item\_id    |  mean\_score|  sd\_score|  max\_score|  pvalue|    rit|    rir|  n\_persons|
|:------------|:------------|------------:|----------:|-----------:|-------:|------:|------:|-----------:|
| verb\_agg   | S1DoCurse   |        1.082|      0.807|           2|   0.541|  0.582|  0.519|         316|
| verb\_agg   | S1DoScold   |        0.832|      0.815|           2|   0.416|  0.651|  0.596|         316|
| verb\_agg   | S1DoShout   |        0.468|      0.709|           2|   0.234|  0.520|  0.460|         316|
| verb\_agg   | S1WantCurse |        1.123|      0.827|           2|   0.562|  0.537|  0.468|         316|
| verb\_agg   | S1WantScold |        0.930|      0.850|           2|   0.465|  0.593|  0.528|         316|
| verb\_agg   | S1WantShout |        0.712|      0.777|           2|   0.356|  0.529|  0.464|         316|

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
| verb\_agg   | dx\_0000001 |              13|  -1.1129358|  -0.3718473|  -0.6251660|  -1.0861430|  -1.4627665|
| verb\_agg   | dx\_0000002 |              28|  -0.0643628|   0.1312981|  -0.6255087|  -0.1063878|   0.3527158|
| verb\_agg   | dx\_0000003 |               4|  -1.8603752|  -2.5817701|  -1.6376942|  -1.7172643|  -1.8129040|
| verb\_agg   | dx\_0000004 |              19|  -1.0663677|  -0.6810231|  -1.1512542|  -0.9481833|  -0.8004328|
| verb\_agg   | dx\_0000005 |               7|  -2.0083386|  -1.5279208|  -1.3307093|  -1.5430697|  -1.3224421|
| verb\_agg   | dx\_0000006 |              25|   0.2020330|   0.2050144|   0.2221347|  -0.4249511|   0.0809038|
