---
title: "Using ggplot2 to visualize your data"
teaching: 90
exercises: 0
questions:
- "What is ggplot2 and why should I use it?"
- "What are aes and geom?"
- "How do I make a basic plot?"
- "What geoms are available?"
- "How can I best visualize groups of data?"
objectives:
- "Get introduced to ggplot2 syntax and philosophy."
- "Learn how to use most common geoms."
- "Learn how to use facets."
- "Learn to change common defaults."
keypoints:
- "ggplot2 is highlevel but so flexible that it replaces most needs for R's base graphics"
- "Generating intuitive, publication-ready graphs in ggplot2 is easy once you get the hang of it"
---




## The grammar of graphics

While plain R comes with its own plotting capabilities, these functions are quiet cumbersome to use as you typically have to write code for every little change you want and the code you write mostly is not very re-usable. Based on the [Grammar of graphics](https://www.amazon.com/exec/obidos/ASIN/0387245448/7210-20), the ggplot2 package implements the idea that graphs can be expressed quite generally using the right way to describe them. For any given plot we have two basic concepts

- **Aesthetic:** What we map data to, like the x or y axis, or higher values being a darker color or bigger circles.
- **Geometries:** How we draw the aesthetic, e.g. take the aesthetic height and create a rectangle with a the given height.

ggplot2 provides a large set of geometries and the means to map aesthetics to these along with capability to arranging plots nicely.
## Input data
Your data must be in a data frame to be really useful with ggplot2. Ideally, the data should also be fairly normalized, i.e. each column should have all the values that go on each asthetic, not spread over multiple columns (you may want to do this using Python pandas) e.g.

| asdf  | asfs   |
|-------|--------|
| asdfd | sdfdsf |

## A first plot
`ggplot2` works well with data frames, particularly when formatted in the 'long' format that our growth data is already in. Plots are initialized with the `ggplot()` function and then we add layers to it to represent the data. Let's first make a simple scatter plot.


~~~
ggplot(growth, aes(x=timepoint, y=od)) +
    geom_point()
~~~
{: .r}



~~~
Error in ggplot(growth, aes(x = timepoint, y = od)): could not find function "ggplot"
~~~
{: .error}

Let's add another layer, a line this time.


~~~
ggplot(growth, aes(x=timepoint, y=od)) +
    geom_point() +
    geom_line()
~~~
{: .r}



~~~
Error in ggplot(growth, aes(x = timepoint, y = od)): could not find function "ggplot"
~~~
{: .error}

Oops, that looks funny. Why? Because we haven't informed ggplot about the strains that each should make up a trajectory in our plot. We can do that by simply adding strain as another aesthetic. 


~~~
ggplot(growth, aes(x=timepoint, y=od, color=strain)) +
    geom_point() +
    geom_line()
~~~
{: .r}



~~~
Error in ggplot(growth, aes(x = timepoint, y = od, color = strain)): could not find function "ggplot"
~~~
{: .error}

Plotting each line in separat facet would have been another option


~~~
ggplot(growth, aes(x=timepoint, y=od)) +
    geom_point() +
    geom_line() +
    facet_wrap(~strain)
~~~
{: .r}



~~~
Error in ggplot(growth, aes(x = timepoint, y = od)): could not find function "ggplot"
~~~
{: .error}

`ggplot2` can present data in a large number of ways, explore the
[online documentation](http://docs.ggplot2.org) or the
[R graph gallery](http://www.r-graph-gallery.com/portfolio/ggplot2-package/)
for inspiration.
