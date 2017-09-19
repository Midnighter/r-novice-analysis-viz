---
title: "Reading data to a data frame and performing basic operations"
teaching: 30
exercises: 0
questions:
- "What is a data frame?"
- "How do I access subsets of a data frame?"
- "How do I calculate simple statistics like mean and median?"
- "How can I access documentation?"
- "How do I install a package?"
- "How do I reshape data?"
objectives:
- "Select individual values and subsections from data."
- "Perform operations on a data frame of data."
- "Be able to install a package from CRAN."
- "Be aware of useful packages for working with data frames."
keypoints:
- "Use `read.table` and `write.table` to import / export data."
- "The function `str` describes the data frame."
- "Use `object[x, y]` to select a single element from a data frame."
- "Use `from:to` to specify a sequence that includes the indices from `from` to `to`."
- "All the indexing and slicing that works on data frames also works on vectors."
- "Use `#` to add comments to programs."
- "Use `mean`, `max`, `min` and `sd` to calculate simple statistics."
- "Use split-apply to calculate statistics across the groups in a data frame."
- "Use dplyr/tidyr in R for manipulating data frames"
---




## Read data into R <!-- 7 -->

Now that we know how to assign things to variables and use functions, let's read some yeast OD growth data into R using `read.table` and briefly examine the dataset.
 


~~~
growth <- read.table(file = "data/yeast-growth.csv", header = TRUE, sep = ",")
~~~
{: .r}

> ## Loading Data Without Headers
>
> What happens if you put `header = FALSE`? The default value is `header = TRUE`?. What do you expect will happen if you leave the default value? Before you run any code, think about what will happen to the first few rows of your data frame, and its overall size. Then run the following code and see if your expectations agree:
>
> 
> ~~~
> head(read.table(file = "data/yeast-growth.csv", header = FALSE, sep = ","))
> ~~~
> {: .r}
> 
> 
> 
> ~~~
>     V1        V2    V3                  V4            V5
> 1 well timepoint    od concentration_level concentration
> 2    a         1 0.017                 low          0.01
> 3    b         1 0.017                 low          0.03
> 4    c         1 0.018              medium             1
> 5    d         1 0.017              medium             3
> 6    e         1 0.017              medium            30
> ~~~
> {: .output}
{: .challenge}


> ## Where is that file? Or, what is my working directory?
>
> R is always running inside a specific directory, the *working
> directory*. Paths can be given relative to that directory so with
> `data/yeast-growth.csv` we mean 'the file `yeast-growth.csv` in the
> `data` directory that is right at the working directory. Set the
> working directory using RStudio `Session` > `Set Working Directory..` or `setwd()`
{: .callout}


> ## Reading Files from Different 'Locales'
>
> Depending on what countrys standard your computer is set to (the 'locale'), software such as Excel will use different characters to separate fields. E.g., the default for a computer with UK defaults will be to use ';' to separate fields and ',' to separate thousands. Try finding the right arguments to `read.table` to get something sensible out of `data/example-us.txt` and `data/example-dk.txt`.
>
> > ## Solution
> > 
> > ~~~
> > read.table('data/example-uk.txt', sep=',', header=TRUE)
> > ~~~
> > {: .r}
> > 
> > 
> > 
> > ~~~
> >          name height age income
> > Ian       183     27  12  1e-02
> > Peter     162     28  11  5e-01
> > Bernhard  173     30  10  0e+00
> > Steven    163     32  12  5e+02
> > ~~~
> > {: .output}
> > 
> > 
> > 
> > ~~~
> > read.table('data/example-dk.txt', sep=';', dec=',', header=TRUE)
> > ~~~
> > {: .r}
> > 
> > 
> > 
> > ~~~
> >       name height age   income
> > 1      Ian    183  27 12000.01
> > 2    Peter    162  28 11100.50
> > 3 Bernhard    173  30 11000.00
> > 4   Steven    163  32 12500.00
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}

Reading excel files is not natively supported in R so we need to use a special package for that, `readxl` is recommended.


~~~
library(readxl)
read_xlsx('data/example-dk.xlsx')
~~~
{: .r}



~~~
# A tibble: 4 x 4
      name height   age income
     <chr>  <dbl> <dbl>  <dbl>
1      Ian    183    27  12000
2    Peter    162    28  11100
3 Bernhard    173    30  11000
4   Steven    163    32  12500
~~~
{: .output}

## Working with data frames <!-- 10 -->

Now that our data is loaded in memory, we can start doing things with it.
First, let's ask what type of thing `growth` is:


~~~
head(growth)
~~~
{: .r}



~~~
  well timepoint    od concentration_level concentration
1    a         1 0.017                 low         1e-02
2    b         1 0.017                 low         3e-02
3    c         1 0.018              medium         1e+00
4    d         1 0.017              medium         3e+00
5    e         1 0.017              medium         3e+01
6    f         1 0.016                high         1e+02
~~~
{: .output}



~~~
str(growth) # what data types are the different columns?
~~~
{: .r}



~~~
'data.frame':	455 obs. of  5 variables:
 $ well               : Factor w/ 7 levels "a","b","c","d",..: 1 2 3 4 5 6 7 1 2 3 ...
 $ timepoint          : int  1 1 1 1 1 1 1 2 2 2 ...
 $ od                 : num  0.017 0.017 0.018 0.017 0.017 0.016 0.015 0.015 0.018 0.021 ...
 $ concentration_level: Factor w/ 3 levels "high","low","medium": 2 2 3 3 3 1 1 2 2 3 ...
 $ concentration      : num  1e-02 3e-02 1e+00 3e+00 3e+01 1e+02 3e+02 1e-02 3e-02 1e+00 ...
~~~
{: .output}



~~~
class(growth)
~~~
{: .r}



~~~
[1] "data.frame"
~~~
{: .output}

The output tells us that is a data frame. Think of this structure as a spreadsheet in MS Excel that many of us are familiar with.
Data frames are very useful for storing data and you will find them elsewhere when programming in R. A typical data frame of experimental data contains individual observations in rows and variables in columns.

We can see the shape, or [dimensions]({{ site.github.url }}/reference/#dimensions), of the data frame with the function `dim`:


~~~
dim(growth)
~~~
{: .r}



~~~
[1] 455   5
~~~
{: .output}

This tells us that our data frame, `growth`, has 455 rows and 5 columns.

If we want to get a single value from the data frame, we can provide an [index]({{ site.github.url }}/reference/#index) in square brackets, just as we do in math:


~~~
 # first value in dat
growth[1, 1]
~~~
{: .r}



~~~
[1] a
Levels: a b c d e f g
~~~
{: .output}



~~~
 # middle value in dat
growth[30, 2]
~~~
{: .r}



~~~
[1] 5
~~~
{: .output}

An index like `[30, 2]` selects a single element of a data frame, but we can select whole sections as well.
For example, we can select the first ten days (columns) of values for the first four patients (rows) like this:


~~~
growth[1:4, 1:2]
~~~
{: .r}



~~~
  well timepoint
1    a         1
2    b         1
3    c         1
4    d         1
~~~
{: .output}

We can use the function `c`, which stands for **c**oncatenate, to select non-contiguous values:


~~~
growth[c(3, 8, 37, 56), c(1, 3)]
~~~
{: .r}



~~~
   well    od
3     c 0.018
8     a 0.015
37    b 0.020
56    g 0.024
~~~
{: .output}

We also don't have to provide a subset for either the rows or the columns.
If we don't include a subset for the rows, R returns all the rows; if we don't include a subset for the columns, R returns all the columns.
If we don't provide a subset for either rows or columns, e.g. `growth[, ]`, R returns the full data frame.


~~~
growth[5, ]
~~~
{: .r}



~~~
  well timepoint    od concentration_level concentration
5    e         1 0.017              medium            30
~~~
{: .output}

> ## Addressing Columns by Name (the better way)
>
> Columns can also be addressed by name, with either the `$` operator (ie. `growth$medium`) or square brackets (ie. `growth[,"medium"]`).
> You can learn more about subsetting by column name in this supplementary [lesson]({{ site.github.url }}/10-supp-addressing-data/).
{: .callout}

Particularly useful is also to user other vectors as filters and only return the rows that evaluate to `TRUE`. Here, `growth$well == "a"` gives a vector with `TRUE` or `FALSE` for every element in `growth$well` that is equal to `"a"`.


~~~
head(growth[growth$well == "e",])
~~~
{: .r}



~~~
   well timepoint    od concentration_level concentration
5     e         1 0.017              medium            30
12    e         2 0.019              medium            30
19    e         3 0.016              medium            30
26    e         4 0.022              medium            30
33    e         5 0.022              medium            30
40    e         6 0.022              medium            30
~~~
{: .output}

Now let's perform some common mathematical operations to learn about our growth curves.


~~~
max(growth[growth$well == "e", "od"])
~~~
{: .r}



~~~
[1] 0.065
~~~
{: .output}

> ## Forcing Conversion
>
> Note that R may return an error when you attempt to perform similar calculations on 
> subsetted *rows* of data frames. This is because some functions in R automatically convert 
> the object type to a numeric vector, while others do not (e.g. `max(growth[1, ])` works as 
> expected, while `mean(growth[1, ])` returns an error). You can fix this by including an 
> explicit call to `as.numeric()`, e.g. `mean(as.numeric(growth[1, ]))` (but mostly this is not what you want to do). By contrast,  calculations on subsetted *columns* always work as expected, since columns of data frames  are already defined as vectors.
{: .callout}

Particularly useful is also to user other vectors as filters and only return the rows that evaluate to `TRUE`. Here, `growth$well` 
R also has functions for other common calculations, e.g. finding the minimum, mean, and standard deviation of the data:


~~~
min(growth[growth$well == "e", "od"])
~~~
{: .r}



~~~
[1] 0.016
~~~
{: .output}



~~~
mean(growth[growth$well == "e", "od"])
~~~
{: .r}



~~~
[1] 0.04347692
~~~
{: .output}



~~~
sd(growth[growth$well == "e", "od"])
~~~
{: .r}



~~~
[1] 0.01259974
~~~
{: .output}

We may want to compare the different wells and for that we can use the split-apply approach which is very common in R. A common approach is to first split the data:


~~~
splitData <- split(growth$od, growth$well)
 ## and then apply a function
sapply(splitData, max)
~~~
{: .r}



~~~
    a     b     c     d     e     f     g 
0.237 0.233 0.231 0.221 0.065 0.040 0.034 
~~~
{: .output}



~~~
 ## or in one go
tapply(growth$od, growth$well, max)
~~~
{: .r}



~~~
    a     b     c     d     e     f     g 
0.237 0.233 0.231 0.221 0.065 0.040 0.034 
~~~
{: .output}

<img src="../fig/split-apply.svg" alt="the split apply approach, divide data to chunks, then run a given function on each ot the chunk sepearately" />

There are many more `apply` style functions among which `lapply` for applying functions to elements of lists, `apply` for applying functions to rows or columns of a matrix.

## Gentle introduction to dplyr and tidyr <!-- 10 -->
Two great packages for doing much more advanced things with data frame are `dplyr` and `tidyr` which together overlaps a lot with Python pandas but it is not practical to use R completely without these so let's cover the basics.

As an example, let's take some mildly messy data and make it easier to use in R. The exact process for cleaning up data is of course entirely dependent on the problem with your data but this example demonstrates the basic process.

Original data looks like this


~~~
messy <- read.csv("data/yeast-growth-messy.csv")
messy[, 1:10]
~~~
{: .r}



~~~
     V1     V2    V3    V4    V5    V6    V7    V8    V9   V10
1 Test1    low 1e-02 0.017 0.015 0.016 0.018 0.022 0.021 0.025
2 Test1    low 3e-02 0.017 0.018 0.015 0.019 0.021 0.020 0.024
3 Test1 medium 1e+00 0.018 0.021 0.016 0.020 0.024 0.023 0.025
4 Test1 medium 3e+00 0.017 0.015 0.017 0.018 0.022 0.022 0.024
5 Test1 medium 3e+01 0.017 0.019 0.016 0.022 0.022 0.022 0.024
6 Test1   high 1e+02 0.016 0.015 0.015 0.018 0.021 0.020 0.023
7 Test1   high 3e+02 0.015 0.016 0.015 0.018 0.021 0.019 0.021
~~~
{: .output}

Problems
- OD measurements for each replicate are on different columns
- There are no sensical header so we would have to work with "messy$V1" etc

Let's fix this using dplyr and tidyr. 


~~~
library(dplyr)
library(tidyr)
~~~
{: .r}

Both these packages make use of the `%>%` operator which allows you to chain functions with each other, e.g. instead of.


~~~
mean(rnorm(10))
~~~
{: .r}

we can write


~~~
rnorm(10) %>%
    mean()
~~~
{: .r}

For this simple example, the first might look easier but with many calls it really helps readability.


~~~
tidy <- as.tbl(messy) %>%
    dplyr::mutate(timepoint=1:7)  %>%
    tidyr::gather(well, od, -V1, -V2, -V3, -timepoint) %>%
    dplyr::rename(concentration_level=V2) %>%
    dplyr::rename(concentration=V3) %>%
    dplyr::select(-V1)
tidy
~~~
{: .r}



~~~
# A tibble: 455 x 5
   concentration_level concentration timepoint  well    od
                <fctr>         <dbl>     <int> <chr> <dbl>
 1                 low         1e-02         1    V4 0.017
 2                 low         3e-02         2    V4 0.017
 3              medium         1e+00         3    V4 0.018
 4              medium         3e+00         4    V4 0.017
 5              medium         3e+01         5    V4 0.017
 6                high         1e+02         6    V4 0.016
 7                high         3e+02         7    V4 0.015
 8                 low         1e-02         1    V5 0.015
 9                 low         3e-02         2    V5 0.018
10              medium         1e+00         3    V5 0.021
# ... with 445 more rows
~~~
{: .output}

How could we reverse the process if we for some reason wanted a wide format again?


~~~
unite(tidy, key, concentration_level, concentration, timepoint) %>%
    spread(key, od)
~~~
{: .r}



~~~
# A tibble: 65 x 8
    well high_100_6 high_300_7 low_0.01_1 low_0.03_2 medium_1_3 medium_3_4
 * <chr>      <dbl>      <dbl>      <dbl>      <dbl>      <dbl>      <dbl>
 1   V10      0.023      0.021      0.025      0.024      0.025      0.024
 2   V11      0.024      0.024      0.028      0.027      0.027      0.026
 3   V12      0.023      0.022      0.026      0.025      0.025      0.025
 4   V13      0.025      0.024      0.028      0.028      0.028      0.028
 5   V14      0.027      0.025      0.030      0.031      0.029      0.030
 6   V15      0.029      0.027      0.032      0.032      0.031      0.032
 7   V16      0.029      0.027      0.033      0.033      0.035      0.034
 8   V17      0.032      0.028      0.037      0.037      0.038      0.036
 9   V18      0.033      0.028      0.041      0.038      0.038      0.038
10   V19      0.035      0.031      0.042      0.042      0.041      0.041
# ... with 55 more rows, and 1 more variables: medium_30_5 <dbl>
~~~
{: .output}
