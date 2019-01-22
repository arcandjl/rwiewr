--- 
title: "Real World Impact Evaluation with R"
author:
- Department of International Economics
- The Graduate Institute | Geneva
- jean-louis.arcand@graduateinstitute.ch
date: "2019-01-22"
output:
  html_document:
    df_print: paged
description: Real world impact evaluation with R.
documentclass: book
github-repo: arcandjl/rwIEwr
link-citations: yes
bibliography: my_refs.bib
site: bookdown::bookdown_site
subtitle: Jean-Louis Arcand
biblio-style: apalike
url: null
---

# Welcome to Real World Impact Evaluation with R {-}

<!-- ```{r echo=FALSE, out.width="33%"} -->
<!-- knitr::include_graphics("figures/f0_web.png") -->
<!-- ``` -->

This is the [online version](https://csgillespie.github.io/efficientR/) of the O'Reilly book: [Efficient R programming](http://shop.oreilly.com/product/0636920047995.do). Pull requests and general comments are welcome.

## Author {-}

[Jean-Louis Arcand](http://www.graduateinstitute.ch) is Professor and Chair of the Department of International Economics at the Graduate Institute in Geneva.

# Preface {-}

Statistically assessing the causal impact of development policies has now become an extremely large industry. While there will always be five available identification strategies, doing things right in a policy-relevant manner is neither obvious, nor easy.[^1]  And the tendency to evaluate stuff simply because it is possible to do so rigorously --instead of trying to evaluate what deserves to be evaluated from the development perspective, has sometimes led to what I\ call the ``tail wagging the dog'' syndrome.
 There is also the minor issue of scientific progress and actually accumulating a useful body of knowledge.
 
[^1]: In case you're wondering, the five identification strategies are: (i) wishing the problem away by assuming selection on observables and applying either OLS or some type of  matching estimator (I will have little if anything to say about this in what follows), (ii) randomization, (iii) instrumental variables, (iv) regression discontinuity design and (v) some sort of covariance transformation that you can implement because you are lucky enough to have panel data.  Convex combinations of these are, of course, very common as well.

The literature on impact evaluation is a subset of econometrics, sometimes with a vocabulary of its own.  As such, econometric methods that you have learnt will figure prominently in what follows.  There is no textbook for this course.  However, there are several great surveys of impact evaluation methods, by masters in the field.  Two of my favorites, which adopt diametrically opposite philosophical stances, are 
@Imbens2009 and @Heckman2007.[^2] 

[^2]: The first views the world through a-theoretical Rubin lenses, the second embraces economic theory.

A very nice survey of more recent approaches, including synthetic cohorts and machine learning methods is given by @Athey2017.  

For my own (cynical) view on all that follows, see @Arcand2014a. But mere cynicism is not going to stop me from enjoying teaching this class.  

There is a plethora of websites dedicated to impact evaluation-related topics that you should consult on a regular basis.  A non-random sample of these includes:

- https://blogs.worldbank.org/impactevaluations/
- http://www.3ieimpact.org/
- https://www.povertyactionlab.org/
- http://www.poverty-action.org/
- http://www.worldbank.org/en/research/dime

Plus my favorite statistics blog: http://andrewgelman.com/
