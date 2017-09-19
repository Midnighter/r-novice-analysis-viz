---
title: "Getting help"
teaching: 15
exercises: 0
questions:
- "How do I use R to do X?"
- "What to do with error X?"
objectives:
- "Learn how to efficiently search for solutions."
- "Learn how to ask questions and help the helpers."
keypoints:
- "Don't panic, try to diagnose your problem."
- "Others may have had the same problem and poster answer online."
- "Good questions are more likely to get good answers."
---



## Find your package

CRAN has >10,000 packages and Bioconductor adds 1383. For most of the analysis you may want to do, someone will have already made a package to help you with this. How do you best find the right package? In addition to google and pubmed, you may want to consider 

- browsing the [bioconductor task views](https://bioconductor.org/packages/3.5/BiocViews.html)
- browsing the [CRAN task views](https://cran.r-project.org/web/views/)
- use the `findFn` function from package `sos` (might not be pretty but can be useful)


~~~
install.packages('sos')
findFn('phylogenetic tree')
~~~
{: .r}
## Finding the right forum

To get introduced to a package we already discussed the vignette's and the function documentation, note in particular the examples. When you get stuck, which will undoubtedly happen at some point, try googling for your question, more often than not, someone has posted an answer on stackoverflow. If no luck, ask for help. The vast majority of package authors / maintainers appreciate hearing from users as it is a proof of the value of their efforts. That said, the work that is being put in is almost always unpaid; make support requests sparingly and, most of all, efficiently.

The 10,000 R packages are written by some 5000? authors. They do not have a coherent forum where you can ask for help so efficiently asking for help may require a (small) bit of investigation. I recommend the following search order

1. Does the package have a webpage? Search for instructions on how to
   ask questions.
2. Is the package on github and you think you found a bug, e.g. and
   exception where there probably shouldn't have been one, *raise an
   issue*
3. Is it a bioconductor package? Make post on
   https://support.bioconductor.org and cc the package maintainer (see
   `DESCRIPTION` file)
4. Make a post to https://stackoverflow.com/ or
   https://www.biostars.org/ and cc the package maintainer.
5. Send an e-mail to
   [R-help](https://stat.ethz.ch/mailman/listinfo/r-help) and cc the
   package maintainer.
6. Last option: just email the package maintainer.

It may seem attractive to immediately asking the package maintainer, but by doing so you also ensure that that person may have to answer your question multiple times as the next person can't google their way to the answer. 

## Good questions are more likely to get good answers

Stackoverflow has a [good article](https://stackoverflow.com/help/how-to-ask) on how to write your question so that you *help the people you are asking for to help you*. It boils down to

- Formulate your question so that it is concise, precise and easy to understand. 
- Has all the details necessary to reproduce your problem -
  particularly important if you are reporting something you think is a
  bug.
- Goes without saying but, be nice and friendly, you're talking to humans now :)

A good question could be formulated like this:

```
*Posted to stackoverflow and/or R-help mailing*

## All package installs on Windows fail with Warning: unable to move temporary installation

All package installations I try fail e.g.:

> install.packages('MASS')
> Warning: unable to move temporary installation 'C:\Program Files\R\R-3.4.1\library\file6cae3bcf\MASS' to 'C:\Program Files\R\R-3.4.1\library\MASS'

I have read-write permissions to target directory but R appears unable to write to that directory.

> sessionInfo()
R version 3.4.1 (2017-04-21)
Platform: Windows 7
Running under: Windows 7

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] sos_2.0-0             brew_1.0-6            ape_4.1
 [4] hexbin_1.27.1         pathview_1.16.5       org.Hs.eg.db_3.4.1

loaded via a namespace (and not attached):
 [1] KEGGgraph_1.38.1  Rcpp_0.12.12      locfit_1.5-9.1    lattice_0.20-35
> 

Thanks in advance for help!
```

> ## What is good about the question above?
>
> Discuss! Anything to improve?
{: .callout!}

> ## Study other peoples questions and discuss
> 
> Look at these questions
> 
> https://support.bioconductor.org/p/100518/
> 
> https://support.bioconductor.org/p/100325/
> 
> https://support.bioconductor.org/p/100541/
> 
> https://support.bioconductor.org/p/100300/
> 
> Is the question clear? Well described? Does it immediately require a question in return or can you answer directly? Discuss with your neighbor!
{: .challenge}
