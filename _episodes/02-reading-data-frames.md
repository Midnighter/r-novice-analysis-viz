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
objectives:
- "Select individual values and subsections from data."
- "Perform operations on a data frame of data."
- "Be able to install a package from CRAN."
keypoints:
- "Use `read.table` and `write.table` to import / export data."
- "The function `str` describes the data frame."
- "Use `object[x, y]` to select a single element from a data frame."
- "Use `from:to` to specify a sequence that includes the indices from `from` to `to`."
- "All the indexing and slicing that works on data frames also works on vectors."
- "Use `#` to add comments to programs."
- "Use `mean`, `max`, `min` and `sd` to calculate simple statistics."
- "Use split-apply to calculate statistics across the groups in a data frame."
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
## Gentle introduction to dplyr and tidyr
Two great packages for doing much more advanced things with data frame are `dplyr` and `tidyr` which together overlaps a lot with Python pandas but it is not practical to use R completely without these so let's cover the basics.
