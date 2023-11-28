---
title: Home
layout: home
---

===============
Practical Introduction to Bayesian Parameter Estimation
===============

This tutorial uses the following _R_ packages:

```R
library(tidyverse)
library(BayesFactor)
library(brms)
library(tidybayes)
```
---
Part I: Why use parameter estimation?
---

Bayes Factors are an intuitive and easily interpretable alternative to null hypothesis significance testing. However,
these ratios do not necessarily provide a complete picture of our beliefs. To demonstrate this, let's work through an example
using the built-in `sleep` dataset.

First, we assign the built-in dataset to a variable so we can modify it and recover the original if necesssary:

```R
d <- sleep
```

Next, let's conduct a Bayes Factor t-test using the `ttestBF` function from the `BayesFactor` package. For this example, we'll compare the extra hours
of sleep that each condition receieved. Note that this paired data, despite the fact that the condition variable is named `group`.

```R
ttestBF(d[d$group == 1,]$extra,
        d[d$group == 2,]$extra,
        paired = TRUE,
        data = d)
```

Output:

```R
Bayes factor analysis
--------------
[1] Alt., r=0.707 : 17.25888 ±0%

Against denominator:
  Null, mu = 0 
---
Bayes factor type: BFoneSample, JZS
```

According to our trusty rules of thumb, the evidence here is pretty solid. But where exactly do we derive this from? 
It might help to know what the posterior actually looks like. Let's repeat this test, but this time we'll generate a posterior
distribution and save it as a variable. For the purposes of this example, we'll draw 10000 samples from the posterior. Because 
sampling from the posterior is an iterative process, the samples drawn will vary slightly everytime this test is repeated, even if
we are working with the exact same dataset. To ensure we all get the exact same values, I include a call to `set.seed()` before
running the test.

```R
set.seed(999)

sleep.post <- ttestBF(d[d$group == 1,]$extra,
                     d[d$group == 2,]$extra,
                     paired = TRUE,
                     posterior = TRUE,
                     iterations = 1e5,
                     data = d)
```

We've successfully saved our posterior in the `sleep.post` variable. Let's take a look at it:

```R
head(supp.post)
```

Output:

```R
Markov Chain Monte Carlo (MCMC) output:
Start = 1 
End = 7 
Thinning interval = 1 
            mu      sig2      delta         g
[1,] -1.630122 2.9337786 -0.9517138 0.1575195
[2,] -1.206233 1.5606834 -0.9655476 0.5256748
[3,] -1.861915 1.8127548 -1.3828983 4.7538352
[4,] -1.474748 0.9622143 -1.5034260 0.6653339
[5,] -1.049121 1.3211610 -0.9127415 1.7021058
[6,] -1.265170 1.8386865 -0.9330288 5.1908654
[7,] -1.277349 2.6216533 -0.7889002 0.8153779
```

There's a lot going on there. For the purposes of this tutorial, we're going to ignore most of it. The parameter
of interest to us is `mu`. This parameter reflects posterior samples - 10000 of them, to be precise - for the unstandardized difference in `extra` between conditions. That might sound somewhat absract - what does it _mean_?. To make it concrete, let's visualize the distribution of our posterior using `ggplot2`:

```R
sleep.post %>%
  as.data.frame() %>% # this makes the data easier to work with using tidyverse
  ggplot(aes(x = mu)) +
  geom_density()
```

You'll get a figure that looks something like this:

![ graph ](assets/images/sleep-dens.png)

Pretty neat, but what does this distribution tell us? 

Essentially, it tells us how our beliefs about `mu` are distributed: The higher the density of a given value, the higher the
posterior probability of that value. As we can see, it appears that our posterior density is highest for values of around 1.5 or so,
meaning values in this region are the most probable. This should be reflected if we calculate descriptive statistics for `mu`:

```R
mean(sleep.post[,'mu'])

median(sleep.post[,'mu'])
```

Output:

```R
[1] -1.408298

[1] -1.415092
```

That checks out, then: It seems that the most likely values for mu are around 1.4. 

But why leave it at that? One of the main advantages of estimating parameters is that your posterior quantifies your _distribution_
of beliefs about the data. Accordingly, a single estimate for `mu` does not necessarily represent this distribution adequately. Looking
at our figure, we can see our posterior includes values ranging from below -4 to upwards of 1. That's a lot of uncertainty about the value of `mu`! 
Of course, not all of these extreme values are meaningful. If we look at the tail ends of the distribution, we can see that posterior density - 
and therefore probability - is quite low. But what about values around -2, or around -0.5? The probability of these values is still fairly
high, and a single estimate derived from the posterior doesn't really capture this. 

Thus, computing a single estimate for `mu` means that we end up ignoring a lot of uncertainty and some fairly probable values. However, if we
use the full posterior, we might overemphasize extreme values that are very unlikely. So what is the solution? 

One common approach is to compute a _Credible Interval_, which essentially trims the most extreme values off of
the posterior distribution. This will provide a more complete representation of the posterior than a single estimate,
but it also won't emphasize extremely improbable values that exist on either end of the range. It will also allow us to make
inferences from the posterior, as I will explain below.

There are several possible approaches we could take to compute such an interval, but I'll focus on the Highest Posterior Density Interval (HPDI), which contains
the highest-probability values for a given confidence level. The confidence level is entirely arbitrary. Although 95% is a common cutoff, there
is no reason this cutoff must be used. We could use 99%, 50%, 89% (as McElreath likes to do), or any other value that we believe is justified. To keep things
simple and straightforward, I'll stick with 95%; however, keep in mind that you don't have to.

With that in mind, let's go back into _R_ and do something concrete. Now that we know what the HDPI is, let's compute this interval for `mu` 
using the `tidybayes` package:

```R
sleep.post %>%
  as.data.frame() %>%
  median_hdi(mu)
```

Output:

```R
# A tibble: 1 × 6
     mu .lower .upper .width .point .interval
  <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
1 -1.42  -2.29 -0.501   0.95 median hdi 
```

Our estimate for `mu` has not changed much. However, we've also computed lower and upper bounds for the HDPI: -2.29 to -0.50.
Let's visualize the HDPI, with the 95% confidence region shaded in blue:

![hdpi](https://github.com/jed709/jed709.github.io/assets/87210399/b849d620-f71b-4407-b229-cec9356e0c65)

As we can see, the HDPI cuts off all the extreme values that don't fall within our arbitrary confidence region. Thus, this interval
provides a good representation of the values in our posterior that we deem credible. Now that we have a Credible Interval, making inferences about
the parameter is easy: 95% of values in our HDPI are below zero, so we are 95% confident that `mu` is negative. To phrase this in terms
of the dataset we are working with, we are 95% certain the the difference between condition 1 and condition 2 falls between -2.29 and -0.50. 
We can say, then, that the effect of condition on sleep is _credible_. 

Accordingly, we can use parameter estimation to accomplish things that Bayes Factors cannot. Rather than a ratio of evidence,
we can compute a distribution of beliefs about a model that incorporates uncertainty. But we can only do so much with the `BayesFactor` package.
What if we had a more complicated design? Further, what if we wanted to incorporate prior knowledge into our models? With this package, the latter is possible,
but options for specifying priors are limited. To address these problems, we must turn to the `brms` package. 

---
Estimating Parameters Using _brms_
---

Although we will work with a more complex design in this section, the basics we learned about working with the posterior still apply here. For our
purposes, let's work with the `ToothGrowth` dataset.

```R
nd <- ToothGrowth
```

This dataset reflects an experiment wherein the dependent measure the length of odontoblasts (`len`) in guinea pigs. This doesn't feel intuitive, so we'll just pretend this is tooth length. The guinea pigs were administered vitamin C using one of two methods, either ascorbic acid (`VC`) or orange juice (`OJ`). Each guinea pig received vitamin C at one of three possible levels of `dose`, `0.5`, `1.0`, or `2.0`. Thus, the design of the experiment is 2 x 3, with both fixed effects manipulated between-subject; this is perfect for our purposes, because multilevel modelling is beyond the scope of this tutorial. 

So, how would we model this experiment? We can think about it in a pretty similar manner to how we would think about an ANOVA: We have two fixed effects we are interested in, and we are also probably interested in the interaction between them. `brms` uses the same formula syntax as the `lme4` package and most _R_ functions that can be used to fit linear models. So, we can translate our design into a linear formula, like so:

```R
len ~ supp * dose
```

Passed to `brms`, this syntax will model the effect of each predictor (and their interaction) on `len`. Before we do that, however, we should prepare our dataset:
We want dose to be a categorical predictor, so we should convert it into a factor like so:

```R
nd = nd %>%
  mutate(dose = factor(dose))
```

_Now_ we can fit a linear model of the data using `brms`. Here's how the code looks:

```R
m.1 <- brm(len~supp*dose,
           data = nd,
           chains = 4,
           cores = 4)
```

Simple as that! This might take a minute to compile and run.

Now let's call the summary function on the output:

```R
 Family: gaussian 
  Links: mu = identity; sigma = identity 
Formula: len ~ supp * dose 
   Data: nd (Number of observations: 60) 
  Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup draws = 4000

Population-Level Effects: 
             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept       13.23      1.16    10.98    15.52 1.00     2102     2638
suppVC          -5.22      1.63    -8.38    -1.96 1.00     1972     2676
dose1            9.49      1.66     6.21    12.72 1.00     1989     2189
dose2           12.83      1.64     9.57    16.09 1.00     2400     2817
suppVC:dose1    -0.71      2.36    -5.31     4.11 1.00     1919     2404
suppVC:dose2     5.30      2.31     0.83     9.79 1.00     2185     2825

Family Specific Parameters: 
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sigma     3.71      0.37     3.06     4.49 1.00     2993     2784

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

There's a lot of output there, but it is very similar to the output you would get from a typical linear model. Let's start with the population-level effects. Because our predictors are categorical,
the model intercept represents tooth length at our reference level, which is `supp = OJ` and `dose = 0.5`. The model slopes represent changes in tooth length
for other levels of our predictors _relative to our reference level_. So, a coefficient of 9.49 for `dose1` means that at levels of `dose = 1` and `supp = VC`, tooth length increased by an average of 9.49 relative to the reference level. See? Simple enough. We also get 95% Credible Intervals for each coefficient. Given that none of these intervals contain 0, it appears that our intercept and all of our slopes are credible. Finally, the model also gives us R-hat statistics and effective sample sizes for each estimate. These are indices of model convergence, which we will not worry about in this tutorial. 
