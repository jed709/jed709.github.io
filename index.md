---
title: Home
layout: home
---

===============
Introduction to Bayesian Parameter Estimation
===============

This tutorial uses the following _R_ packages:

```R
library(tidyverse)
library(BayesFactor)
library(brms)
library(tidybayes)
```
---
Why use parameter estimation?
---

Bayes Factors are an intuitive and easily interpretable alternative to null hypothesis significance testing. However,
these ratios do not necessarily provide a complete picture of our beliefs. To demonstrate this, let's work with an example
using the `ToothGrowth` dataset.

First, we assign the built-in dataset to a variable so we can modify it and recover the original if necesssary:

```R
d <- ToothGrowth 
```

```R
mt.2 <- brm(len~supp*dose,
            data = TG_noMid,
            chains = 4,
            cores = 4,
            prior = c(prior(normal(10, 10), class = 'Intercept'),
                      prior(normal(0,20), class = 'b'),
                      prior(normal(0,10), class = 'sigma')),
            sample_prior = 'yes',
            file = 'toothmod_2')
```

![ graph ](assets/images/e4-faceted-c.png)
