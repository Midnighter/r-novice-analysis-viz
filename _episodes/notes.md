Since our data is in a data frame however, the better option is to use the `dplyr` package for this sort of operations. `dplyr` has many similarities with Python pandas and we could for example


```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
group_by(growth, strain) %>%
    summarize(max_od=max(od))
```

```
## Error in group_by(growth, strain): object 'growth' not found
```
