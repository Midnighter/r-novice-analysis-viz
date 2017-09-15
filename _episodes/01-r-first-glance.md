---
title: "A brief introduction to R"
teaching: 40
exercises: 0
questions:
- "What is R and when should I use it?"
- "How do I use the RStudio graphical user interface?"
- "How do I assign variables?"
objectives:
- "Get introduced to R and RStudio."
- "Learn basic R syntax."
- "Assign values to variables."
- "Be able to install a package from CRAN."
keypoints:
- "R is a strong statistical computing environment"
- "Thousands of packages for R"
- "Use `variable <- value` to assign a value to a variable in order to
  record it in memory."
- "Objects are created on demand whenever a value is assigned to them."
---




## What is R? <!-- 2 -->

R is a statistical computing environment that has become a widely adopted standard for data analysis various fields - notably bioinformatics.

- R is a dialect of the old (commercial) S language
- R is relatively easy to learn to use - tons of built-in functions
  that are just right for data analysis, focused on interactive use
  and often very high-level
- There are thousands of package implementing algorithms, procedures,
  plots - for life science in particular
- R is also a programming language

Consider using R to get access to packages that implement solution to a given problems like

- differential gene expression analysis, gene set enrichment analysis,
- metabolomics data analysis,
- gene clustering,
- access to online databases on gene expression / annotation etc,
- qRT-PCR,
- multivariate statistics, mixed linear models,
- and much much more..

In this one-day course we will learn

- How to use RStudio, the most diverse and mature integrated
  development environment for R
- R's basic syntax
- Re..
- Visualising data using ggplot2

## RStudio - a graphical interface to R <!-- 1 -->

While R can be used directly in the shell, it is much nicer with a graphical interface. RStudio is by far the best one - let's learn how to use it.

* Point to the different panels: Console, Scripts, Environments, Plots.
* Code and workflow are more reproducible if we can document everything that we do.
* Our end goal is not just to "do stuff" but to do it in a way that anyone can
  easily and exactly replicate our workflow and results.
* The best way to achieve this is to write scripts. RStudio provides an
  environment that allows you to do that.

## Interacting with R <!-- 5 -->

There are two main ways of interacting with R: using the console or by using
script files (plain text files that contain your code).

The console window (in RStudio, the bottom left panel) is the place where R is
waiting for you to tell it what to do, and where it will show the results of a
command.  You can type commands directly into the console, but they will be
forgotten when you close the session. It is better to enter the commands in the
script editor, and save the script. This way, you have a complete record of what
you did, you can easily show others how you did it and you can do it again later
on if needed. You can copy-paste into the R console, but the Rstudio script
editor allows you to 'send' the current line or the currently selected text to
the R console using the `Ctrl-Enter` shortcut.

At some point in your analysis you may want to check the content of variable or
the structure of an object, without necessarily keep a record of it in your
script. You can type these commands directly in the console. RStudio provides
the `Ctrl-1` and `Ctrl-2` shortcuts allow you to jump between the script and the
console windows.

If R is ready to accept commands, the R console shows a `>` prompt. If it
receives a command (by typing, copy-pasting or sent from the script editor using
`Ctrl-Enter`), R will try to execute it, and when ready, show the results and
come back with a new `>`-prompt to wait for new commands.

If R is still waiting for you to enter more data because it isn't complete yet,
the console will show a `+` prompt. It means that you haven't finished entering
a complete command. This is because you have not 'closed' a parenthesis or
quotation. If you're in RStudio and this happens, click inside the console
window and press `Esc`; this should help you out of trouble.

## Basic R syntax and loading a package <!-- 7 -->

Just like with python, we can perform simple operations using the R console and assign the output to variables


~~~
1 + 1
~~~
{: .r}



~~~
[1] 2
~~~
{: .output}

`<-` is the assignment operator. It assigns values on the right to objects on
the left. So, after executing `x <- 3`, the value of `x` is `3`. The arrow can
be read as 3 **goes into** `x`.  You can also use `=` for assignments but not in
all contexts so it is good practice to use `<-` for assignments. `=` should only
be used to specify the values of arguments in functions, see below.

> ## Save some keystrokes..
> In RStudio, typing `Alt + -` (push `Alt`, the key next to your space bar at the
> same time as the `-` key) will write ` <- ` in a single keystroke.
{: .callout}


~~~
foo <- 1 + 1
foo + 1
~~~
{: .r}



~~~
[1] 3
~~~
{: .output}

> ## Commenting
>
> We can add comments to our code using the `#` character. It is
> useful to document our code in this way so that others (and us the
> next time we read it) have an easier time following what the code is
> doing.
{: .callout}

We can create some random numbers from the normal distribution and calculate the mean by using two built-in functions


~~~
randomNumbers <- rnorm(10)
mean(randomNumbers)
~~~
{: .r}



~~~
[1] 0.2680998
~~~
{: .output}

> > ## Variable Naming Conventions
>
> Historically, R programmers have used a variety of conventions for naming variables. The `.` character
> in R can be a valid part of a variable name; thus the above assignment could have easily been `random.numbers <- rnorm(10)`.
> This is often confusing to R newcomers who have programmed in languages where `.` has a more significant meaning (like in Python).
> There are many 'standards' in use e.g. `random.numbers`, `random_numbers` or `randomNumbers`. Choose what you prefer, but, just as with british or american spelling, the rule is to be consistent!
{: .callout}


> ## All those built-ins..
>
> In contrast to Python, R has many functions available without
> loading any extra packages (2383 to be exact). These are mostly
> functions and becoming proficient with R is a lot about learning how
> these work.
{: .callout}

While plain R already comes with many functions, one almost always wants to make use of code contributed as packages. The concept is very similar to Python but has some important differences. Say we want to use the `readxl` package which allows us to read Excel sheets directly to R. If you by any chance do not have the `readxl` package on your system you have to install it first. This can be done via the menues in RStudio, or by typing:


~~~
install.packages("readxl")
~~~
{: .r}

After that, we can load the package by:


~~~
library(readxl) # equivalent in Python would have been 'from readxl import *'
~~~
{: .r}

To see what became available we can


~~~
ls('package:readxl') # or read the documentation
help(package='readxl')
~~~
{: .r}

Sometimes you just want to use one function and do not want to load all functions. Then you can, without loading the package call the functions using the `::` syntax.


~~~
readxl::readxl_example()
~~~
{: .r}



~~~
 [1] "clippy.xls"    "clippy.xlsx"   "datasets.xls"  "datasets.xlsx"
 [5] "deaths.xls"    "deaths.xlsx"   "geometry.xls"  "geometry.xlsx"
 [9] "type-me.xls"   "type-me.xlsx" 
~~~
{: .output}
## Reading function documentation <!-- 1 -->

All R functions are documented and you can read about them using RStudio documentation pane, or typing `?object`, eg


~~~
?mean
~~~
{: .r}
## R data types <!-- 5 -->
R has three main data types that we need to know about, the two main ones are `numeric` which is both integers and floats, `character`, and `factor` which is like integer but with a character label attached to number (or `level`). 


~~~
is.numeric(1)
~~~
{: .r}



~~~
[1] TRUE
~~~
{: .output}



~~~
is.character("foo")
~~~
{: .r}



~~~
[1] TRUE
~~~
{: .output}



~~~
is.factor(factor(c("a", "b", "c", "c")))
~~~
{: .r}



~~~
[1] TRUE
~~~
{: .output}

There are numerous data structures and classes but mainly three we need to care about.

**Vectors** Scalars act like vectors of length 1.

~~~
foo <- c(1, 2, 3)
foo[2]
~~~
{: .r}



~~~
[1] 2
~~~
{: .output}



~~~
bar <- 1
bar[1]
~~~
{: .r}



~~~
[1] 1
~~~
{: .output}

> ## Indexing starts at 1!
>
> Unlike python (and indeed most other programming languages) indexing starts at 1. `bar[0]` is however still syntactic correct, it just selects the null-vector from vector `bar` thus always a vector of length 0.
{: .callout}

**Matrices** work like one would expect, values in two dimensions. 


~~~
(foo <- matrix(c(1, 2, 3, 4, 5, 6), nrow=2))
~~~
{: .r}



~~~
     [,1] [,2] [,3]
[1,]    1    3    5
[2,]    2    4    6
~~~
{: .output}

**Lists** are widely used and work like vectors except that each element can be any data structure - so you can have lists of matrices, lists of lists etc.


~~~
(foo <- list(bar=c(1, 2, 3), baz=matrix(c(1, 2, 3, 4, 5, 6), nrow=2)))
~~~
{: .r}



~~~
$bar
[1] 1 2 3

$baz
     [,1] [,2] [,3]
[1,]    1    3    5
[2,]    2    4    6
~~~
{: .output}



~~~
foo[[2]]
~~~
{: .r}



~~~
     [,1] [,2] [,3]
[1,]    1    3    5
[2,]    2    4    6
~~~
{: .output}



~~~
foo[2]
~~~
{: .r}



~~~
$baz
     [,1] [,2] [,3]
[1,]    1    3    5
[2,]    2    4    6
~~~
{: .output}



~~~
foo$baz
~~~
{: .r}



~~~
     [,1] [,2] [,3]
[1,]    1    3    5
[2,]    2    4    6
~~~
{: .output}

Finally, **data frames** but they are so important they deserve a small section on their own.

For a rough comparison: 

- R vectors are like lists (`[]`) in Python 
- R lists are like dicts (`{}`) in Python
- R data frames are very much like `pandas.DataFrame` in Python

    
