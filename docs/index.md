--- 
title: "Real World Impact Evaluation with R"
subtitle: 'Jean-Louis Arcand'  
author: 
- 'Department of International Economics'
- 'The Graduate Institute | Geneva'
- 'jean-louis.arcand@graduateinstitute.ch'
date: "2019-01-22"
site: bookdown::bookdown_site
documentclass: book
bibliography: [my_refs.bib]
biblio-style: apalike
link-citations: yes
github-repo: arcandjl/rwIEwr
url: 
description: "Real world impact evaluation with R."
---

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
