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

For this tutorial, I will assume that readers have some baseline knowledge of (1) Bayesian reasoning and statistics and (2) linear regression. This guide is intended to be practical and to show you that you - yes, you - can estimate parameters using Bayesian methods. I will briefly cover theoretical concepts where necessary, but my approach throughout this guide will be largely hands-on. 

With that out of the way, let's get started. 

---
Part I: Why use parameter estimation?
---

Assuming you are at least slightly familiar with Bayesian statistics, you are probably familiar with Bayes Factors. Bayes Factors are an intuitive and easily interpretable alternative to approaches that our statistics courses have taught us to distrust, like null hypothesis significance testing. So why would we want to use paramater estimation? The output of regression models is messy and more difficult to interpet. Can't we just stick with model comparison using Bayes Factors? The issue here is that Bayes Factors do not necessarily provide a complete, nuanced picture of our beliefs. Like everything else in life, statistics is not black and white. It is difficult - and probably not best practice - to condense our beliefs about something down to a single number. This is where Bayesian parameter estimation comes in: Using this approach, we can map out _distributions_ of beliefs that quantify our uncertainty about the data.

To demonstrate this, let's start with using Bayes Factors. We'll work through an example using the built-in `sleep` dataset.

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

According to our trusty rules of thumb, the evidence here is pretty solid: A difference between groups is favored over a null model by a factor of about 17. But where exactly do we derive this from? Does this Bayes Factor provide a full picture of our belief in the effect? Probably not.

To demonstrate this, it might help to know what our posterior for this test would actually look like. Let's repeat this test, but this time we'll generate a posterior
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
Because our Credible Interval does not contain zero, we can say that the effect of condition on sleep is _credible_. 

Let's return to the Bayes Factor we computed earlier. The inferences we derive from the Bayes Factor and parameter estimation are the same, in that we have good reason to believe there is a difference between conditions. However, the Bayes Factor did not provide a complete picture of the effect. With parameter estimation, then, we can accomplish things that Bayes Factors cannot. Rather than a ratio of evidence,
we can compute a distribution of beliefs about a model that incorporates uncertainty. 

But we can only do so much with the `BayesFactor` package.
What if we had a more complicated design? Further, what if we wanted to incorporate prior knowledge into our models? With this `BayesFactor`, the latter is possible, but options for specifying priors are limited. To address these problems, we must turn to the `brms` package. 

---
Part II: The Basics of Estimating Parameters Using _brms_
---

Although we will work with a more complex design in this section, the basics we learned about working with the posterior still apply here. For our
purposes, let's work with the `ToothGrowth` dataset.

```R
nd <- ToothGrowth
```

This dataset reflects an experiment wherein the dependent measure the length of odontoblasts (`len`) in guinea pigs. This doesn't feel intuitive, so we'll just pretend this is tooth length. The guinea pigs were administered vitamin C using one of two methods, either ascorbic acid (`VC`) or orange juice (`OJ`). Each guinea pig received vitamin C at one of three possible levels of `dose`, `0.5`, `1.0`, or `2.0`. Thus, the design of the experiment is 2 x 3, with both fixed effects manipulated between-subject; this is perfect for our purposes, because multilevel modelling is beyond the scope of this tutorial. 

So, how would we model this experiment? We can think about it in a pretty similar manner to how we would think about an ANOVA: We have two fixed effects we are interested in, and we are also probably interested in the interaction between them. `brms` uses the same formula syntax as the `lme4` package and most _R_ functions that can be used to fit linear models (e.g., `lm()`). So, it is fairly easy to translate our design into a linear formula, like so:

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

We can also directly access the full posterior for any model term. This can be accomplished using the `as_draws_df` function from the `brms` package. Let's try that and see what we get:

```R
m.1 %>%
  as_draws_df() %>%
  head()
```

Output:

```R
# A draws_df: 6 iterations, 1 chains, and 9 variables
  b_Intercept b_suppVC b_dose1 b_dose2 b_suppVC:dose1 b_suppVC:dose2 sigma lprior
1          12     -4.6    10.9      13           0.15            5.3   3.6   -5.8
2          14     -6.8     8.0      13           1.57            4.7   3.7   -5.8
3          15     -8.3     8.1      10           2.23            8.7   3.6   -5.8
4          14     -6.1     6.7      10           1.92            7.8   4.0   -5.8
5          15     -7.0     8.0      12           2.58            6.8   3.6   -5.8
6          15     -6.6     6.2      12           3.10            6.2   3.1   -5.8
# ... with 1 more variables
# ... hidden reserved variables {'.chain', '.iteration', '.draw'}
```

All of our model coefficients here. The values in this dataframe represent draws or samples from the posterior for each coefficient. This is pretty similar to what we did using `ttestBF`, except we have done it for a much more complex model. As an aside, we can see in our output from `as_draws_df` that `brms` uses weird naming conventions for model coefficients. If you're not familiar with these conventions, and if you don't want to generate a dataframe every time you need to know what the name of a model term is, the `get_variables` function from the `tidybayes` package is very handy:

```R
get_variables(m.1)
```

Output:

```R
 [1] "b_Intercept"    "b_suppVC"       "b_dose1"        "b_dose2"        "b_suppVC:dose1" "b_suppVC:dose2"
 [7] "sigma"          "lprior"         "lp__"           "accept_stat__"  "stepsize__"     "treedepth__"   
[13] "n_leapfrog__"   "divergent__"    "energy__"  
```

Easy!

Now, to get comfortable working with posteriors from complex models, let's look at the posterior for the `supp = OJ` and `dose = 0.5` condition (i.e., our intercept). Remember, to access the posterior for this condition, we will need to use the model intercept, i.e., the `b_Intercept` term. Let's glimpse at the posterior distribution for the intercept:

```R
m.1 %>%
  as_draws_df() %>%
  ggplot(aes(x = b_Intercept)) +
  geom_density() +
  theme_classic()
```

Output:

![intercept-post](https://github.com/jed709/jed709.github.io/assets/87210399/6e6769f0-1f50-481e-bbf0-95caa70eb894)

This is a similar exercise to what we did with our posterior from the sleep dataset. As we can see, the most probable values are around 13 or so, which matches up well with the model output. We can also use `as_draws_df` in conjunction with `median_hdi` from `tidybayes` to calculate a credible interval for any posterior distribution derived from our model:

```R
m.1 %>%
  as_draws_df() %>%
  median_hdi(b_Intercept)
```

Output:

```R
# A tibble: 1 × 6
  b_Intercept .lower .upper .width .point .interval
        <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
1        13.2   11.1   15.6   0.95 median hdi 
```

This might not seem that exciting, because our `brms` output already gave us a Credible Interval for this term, and it is quite similar to the HDI. However, something more exciting is that we can use this general method to conduct any number of arbitrary comparisons and evaluate whether or not they are meaningful. For example, what if we wanted to check whether tooth length in the `VC` condition was higher for `dose = 2` relative to `dose = 1`? This is very easy to do with `as_draws_df`. Here's a demonstration: 

```R
m.1 %>%
  as_draws_df() %>%
  
  # Remember, the model slopes reflect changes relative to the reference level.
  # To get the mean for a given condition, we need to add the slopes to the intercept
  
  mutate(contrast = (b_Intercept + b_dose2) - (b_Intercept + b_dose1)) %>%
  
  # Then we can compute the desired contrast between conditions and produce a
  # credible interval
  
  select(contrast) %>%
  median_hdi()
```

Output:

```R
# A tibble: 1 × 6
  contrast .lower .upper .width .point .interval
     <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
1     3.34  0.135   6.56   0.95 median hdi  
```

As we can see, our median posterior estimate is 3.34, with the HDI ranging from 0.14 to 6.56. In other words, tooth length was credibly higher for the `supp = OJ` and `dose = 2` condition relative to the `dose = 1` condition for the same method of delivery. This is what I find most intutitive about working with Bayesian models. When it comes to contrasts, the world is your oyster. You can generate a contrast for any relationship you are interested in - no need for special statistical tests to evaluate pairwise comparisons. 

The contrasts and condition means you compute can also be translated into informative visualizations. For example, what if we were interested in depicting the contrast between high and medium doses (`dose2` - `dose1`) and the contrast between high and low doses (`dose2` - `dose0.5`). We can make a pretty cool plot out of this using the `stat_halfeye` function from `tidybayes`. This function is imported from the `ggdist` package, which you'll want to install if you'll be doing a lot of this sort of visualization. For our purposes, however, the function `tidybayes` has will suffice. Here's how we could plot these two contrasts together:

```R
m.1 %>%
  as_draws_df() %>%
  mutate(HighMid = (b_Intercept + b_dose2) - (b_Intercept + b_dose1),
         
         # dose0.5 is our intercept, so the difference between dose2 and
         # dose0.5 is just equal to the slope of dose2
         
         HighLow = b_dose2) %>% 
  
  pivot_longer(cols = c(HighMid, HighLow), 
               names_to = 'contrast', 
               values_to = 'post') %>%
  ggplot(aes(x = post, 
             y = contrast)) +
  stat_halfeye() + 
  theme_classic()
```

Output:

![con_plot_1](https://github.com/jed709/jed709.github.io/assets/87210399/083e0249-17b7-4801-8c97-020f4822ce07)

Cool plot, right? It contains a lot of information. The point is your median posterior estimate. The thick line surrounding the point is the 50% HDI, and the thin line is the 95% HDI. On top of the interval, we get a visualization of the entire posterior distribution of the estimate. These contrasts might not be the most interesting to visualize, but it's a simple proof-of-concept that you can build off. 

Now that we understand the basics of fitting models with `brms` and working with the posteriors they generate, let's go through a more complex - but more practical - example, using the same dataset.

---
Part III: A Practical Example
---

You might have noticed that we fit our model of the `ToothGrowth` dataset without specifying any priors. This is not great. When we don't specify priors, `brms` uses default priors. What does this mean, exactly? To understand, we need to take a look at the priors the model used, which we can access using `prior_summary`:

```R
prior_summary(m.1)
```

Output:

```R
                 prior     class         coef group resp dpar nlpar lb ub       source
                (flat)         b                                               default
                (flat)         b        dose1                             (vectorized)
                (flat)         b        dose2                             (vectorized)
                (flat)         b       suppVC                             (vectorized)
                (flat)         b suppVC:dose1                             (vectorized)
                (flat)         b suppVC:dose2                             (vectorized)
 student_t(3, 19.2, 9) Intercept                                               default
    student_t(3, 0, 9)     sigma                                     0         default
```

We can see that by default, `brms` puts "flat" priors on our slopes. These are the priors we will focus on, for now. But what would "flat" priors look like? Something like this:

![unif](https://github.com/jed709/jed709.github.io/assets/87210399/eb1f5f3e-4d2c-47dd-be99-be6674d95f3a)

I've added tails so the "curve" is visible. The values on each axis aren't important. What is important is that all values within the distribution (within the range of -5 to 5, for the purposes of my visualization) have equal probability. So, this is the prior knowledge that `brms` puts on our slopes if we don't specify anything. This is about as uninformative as a prior can get: We are telling the model that any given value of our slopes is just as likely as any other given value. This is also nearly equivalent to the approach taken by Frequentist methods.

It is good not to use priors that are _too_ informative; this could tip the scale one way or the other, which we want to avoid. In practice, however, there is no reason that priors need to be _so_ uninformative as to afford equal probablity to all possible values. What reasonable constraints might we apply to the data?

In this example, our dependent variable is tooth length. The default prior used by `brms` will consider values below zero to be just as likely as values above zero. This makes no sense. A guinea pig could have a tooth length of zero, hypothetically, but lower values are not possible: Guinea pigs cannot have negative teeth. Accordingly, it would be very reasonable of us to tell the model not to consider values below zero. Similarly, the fact that we are working with guinea pigs tells us that we can place some reasonable constraints on how long their teeth can be. Guinea pigs are quite small, so their teeth are not going to be comparable in size to the tusks of an elephant. Using flat priors, we are essentially saying that tusk-sized teeth are just as likely as something more reasonable. 

With that in mind, let's fit a better model that incorporates some reasonable prior knowledge. Before we do that, however, let's modify the parameterization of the model so we can more easily specify priors for each condition. Because we are using categorical predictors, we can remove the model intercept and compute slopes that estimate the mean in each condition. Using this approach, our reference level will still be modelled, and the output will be much more intutive to interpret. Rather than thinking about the slopes as differences between conditions, we will instead simply get an estimate for each cell. This will also make contrasts easier to compute. We can remove the model intercept using the same syntax we would use in `lme4`. Here's the formula we'll use:

```R
len~supp:dose-1
```

Now that we have our formula, we can look at the parameterization of the model and decide how we will specify our priors. We can access this easily without fitting a model, using the very handy `get_prior` function from `brms`. If we were fitting a different type of model, we'd need to give it some additional information. But for our purposes, all we need is the formula and the dataframe we're using:

```R
get_prior(formula = len~supp:dose-1,
          data = nd)
```

Output:

```R
              prior class           coef group resp dpar nlpar lb ub       source
             (flat)     b                                                 default
             (flat)     b suppOJ:dose0.5                             (vectorized)
             (flat)     b   suppOJ:dose1                             (vectorized)
             (flat)     b   suppOJ:dose2                             (vectorized)
             (flat)     b suppVC:dose0.5                             (vectorized)
             (flat)     b   suppVC:dose1                             (vectorized)
             (flat)     b   suppVC:dose2                             (vectorized)
 student_t(3, 0, 9) sigma                                       0         default
```

Now we can see all the terms in our model. As we can see, we'll compute a slope corresponding to tooth length for every possible combination of our fixed effects. Pretty intuitive, right? This also means we don't have to specify priors on our slopes with respect to the _difference_ we expect relative to the reference level. Instead, we can specify priors that reflect our belief about what tooth length should be in each condition. 

Now, with this parameterization, we could give the model a specific prior for each condition. However, we do want to avoid tipping the scale, so we'll keep our priors pretty general and mostly uninformative. There are two pieces of prior knowledge we want to incorporate here:

1. Tooth length cannot be negative
2. Guinea pig teeth won't be particularly long

To do this, we'll have to come up with a belief about guinea pig tooth length. I'm no expert, but you can pretend I researched this and came up with the principled hypothesis that tooth length in any condition should fall between 0 and 16, with a mean of 8. In practice, your priors should be calibrated to principled beliefs. We're also going to set a hard lower bound on tooth length in all conditions so the model doesn't afford probability to values below zero. In `brms`-speak, all of that would look like this:

```R
prior(normal(8, 4), class = 'b', lb = 0)
```

This sets a normal prior - with a mean of 8 and a standard deviation of 4 - on all of our model slopes. If you wanted to set priors on individual model coefficients, you could do it like this:

```R
prior(normal(8, 4), coef = 'suppOJ:dose2')
```

But you would only want to do that if you had good reason to do so. We won't do that today.

Another thing we won't get into is addressing the priors for sigma. By default, `brms` uses a half-Cauchy prior, which is a Cauchy distribution truncated at zero. Sigma cannot be negative, so this is good. Half-Cauchy priors are uninformative and are recommended for sigma (e.g., Gelman, 2006). We don't have any principled reason to mess with this prior in this tutorial, so we'll leave it as is.

Now that we've gone over all of that, let's fit a new model, minus the intercept and inclusive of our principled priors on the model slopes. Here's how we could do that:

```R
m.2 <- brm(len~supp:dose-1,
           prior = prior(normal(8, 4), class = 'b', lb = 0),
           data = nd,
           chains = 4,
           cores = 4)
```

And our output:

```R
 Family: gaussian 
  Links: mu = identity; sigma = identity 
Formula: len ~ supp:dose - 1 
   Data: nd (Number of observations: 60) 
  Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup draws = 4000

Population-Level Effects: 
               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
suppOJ:dose0.5    12.75      1.21    10.37    15.12 1.00     5051     2749
suppVC:dose0.5     7.95      1.20     5.57    10.28 1.00     5438     2539
suppOJ:dose1      21.41      1.23    18.86    23.72 1.00     5609     2585
suppVC:dose1      16.00      1.18    13.57    18.25 1.00     5436     2978
suppOJ:dose2      24.51      1.21    22.07    26.78 1.00     4986     3050
suppVC:dose2      24.57      1.21    22.07    26.80 1.00     4728     2760

Family Specific Parameters: 
      Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sigma     3.88      0.42     3.16     4.78 1.00     3734     3284

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

Nice - much easier to interpret than our previous model! You'll also notice that tooth length in our reference level is a bit lower than before. This is because our priors pulled it a little bit closer to zero. Once again, we get a 95% Credible Interval for each estimate, but in this case, these are less meaningful. All they tell us is that tooth length in each condition is credibly different from zero, which we would hope it should be (unless we had toothless guinea pigs). 

But the real strengths of this model arise when we want to compute contrasts. The way we have paramterized the model makes this incredibly easy. For example, let's see if tooth length was greater in the `suppOJ:dose1` condition relative to the `suppVC:dose1` condition. Note that I wrap the column names in backticks. _R_ will interpret the `:` as an operator if you don't do this, so this is useful to know if you'll be working with models with crossed effects. Anyway, here's the code:

```R
m.2 %>%
  as_draws_df() %>%
  mutate(contrast = `b_suppOJ:dose1` - `b_suppVC:dose1`) %>%
  median_hdi(contrast)
```

Output:

```R
# A tibble: 1 × 6
  contrast .lower .upper .width .point .interval
     <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
1     5.40   1.96   8.53   0.95 median hdi 
```

Indeed, it's credibly higher for the former. Let's try something else - what if we wanted to know if average tooth length was higher for the `OJ` group overall relative to the `VC` group? We could check by averaging the lengths for each group and then subtracting one from the other, like this:

```R
m.2 %>%
  as_draws_df() %>%
  mutate(contrast = 
           (`b_suppOJ:dose0.5` + `b_suppOJ:dose1` + `b_suppOJ:dose2`)/3 - 
           (`b_suppVC:dose0.5` + `b_suppVC:dose1` + `b_suppVC:dose2`)/3) %>%
  median_hdi(contrast)
```

Output:

```R
# A tibble: 1 × 6
  contrast .lower .upper .width .point .interval
     <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
1     3.39   1.58   5.37   0.95 median hdi   
```

Credibly higher for the `OJ` group - analogous to a main effect of group that you would get in an ANOVA. We could do the same thing for levels of dose, if we chose to. It's also easy to check for interaction between the fixed effects, and one appears to be obvious from our model output. We already know that tooth length is higher for the `OJ` relative to the `VC` group at `dose = 1`, so let's compute the contrasts for the other levels of `dose`:

```R
m.2 %>%
  as_draws_df() %>%
  mutate(contrast = `b_suppOJ:dose0.5` - `b_suppVC:dose0.5`) %>%
  median_hdi(contrast)
```

Output: 

```R
# A tibble: 1 × 6
  contrast .lower .upper .width .point .interval
     <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
1     4.79   1.46   8.25   0.95 median hdi 
```

Same pattern for low doses. What about high doses?

```R
m.2 %>%
  as_draws_df() %>%
  mutate(contrast = `b_suppOJ:dose2` - `b_suppVC:dose2`) %>%
  median_hdi(contrast)
```

Output:

```R
# A tibble: 1 × 6
  contrast .lower .upper .width .point .interval
     <dbl>  <dbl>  <dbl>  <dbl> <chr>  <chr>    
1  -0.0765  -3.40   3.30   0.95 median hdi   
```

Just as I suspected - no credible difference at higher doses. Thus, it seems that `supp` and `dose` interact such that tooth length is higher for the `OJ` group relative to the `VC` group and low and medium doses, but length is similar at high doses. If this study wasn't so mundane, maybe this would be a very interesting finding. 

The reason I've gone through all these arbitrary contrasts is to show you that you can really do whatever you please with the posterior when estimating parameters. This is a verstaile approach and I believe working with posteriors in this manner is about as intuitive as statistics can be, once you get past the more intimidating parts of the code. As a final exercise, let's put together some visualizations of our more complex model. One option would be to plot the estimates for each condition by group. This should do a fairly good job of conveying the pattern of results we observed. Here's how we could do it:

```R
m.2 %>%
  as_draws_df() %>%
  pivot_longer(cols = `b_suppOJ:dose0.5`:`b_suppVC:dose2`,
               names_to = 'Dose',
               values_to = 'Length') %>%
  mutate(Supp = if_else(grepl('VC', Dose), 'VC', 'OJ')) %>%
  mutate(Dose = case_when(grepl('0.5', Dose) ~ 'Low',
                               grepl('1', Dose) ~ 'Mid',
                               grepl('2', Dose) ~ 'High')) %>%
  mutate(Dose = factor(Dose, levels = c('Low', 'Mid', 'High'))) %>%
  ggplot(aes(x = Length, y = Dose, fill = Supp)) +
  stat_halfeye(position = position_dodge(width = 0.9), 
               scale = 0.7, slab_alpha = 1,
               point_interval = median_hdi) +
  scale_fill_brewer(palette = 'Dark2') +
  theme_classic()
```

Which gives you this:

![complex-model-plot1](https://github.com/jed709/jed709.github.io/assets/87210399/f0aa1e76-195d-4f2a-aad8-b4e98f453e79)

Pretty cool, right? And it's informative, too! You can see the HDI, the full posterior, and the pattern of results is pretty easy to pick out. Of course, there are other ways we could visualize the data. For example, if we thought the `OJ` group would generally have much higher tooth length relative to the `VC` group, we might want to graphically highlight the contrast between levels of `supp` at different levels of `dose`. Hypothetical, of course, but here's how you could do that:

```R
m.2 %>%
  as_draws_df() %>%
  mutate('OJ Low - VC Low' = `b_suppOJ:dose0.5` - `b_suppVC:dose0.5`,
         'OJ Mid - VC Mid' = `b_suppOJ:dose1` - `b_suppVC:dose1`,
         'OJ High - VC High' = `b_suppOJ:dose2` - `b_suppVC:dose2`) %>%
  pivot_longer(cols = `OJ Low - VC Low`:`OJ High - VC High`, 
               names_to = 'Contrast', values_to = 'Length') %>%
  mutate(Contrast = factor(Contrast,
                           levels = c('OJ Low - VC Low',
                                      'OJ Mid - VC Mid',
                                      'OJ High - VC High'))) %>%
  ggplot(aes(x = Length, y = Contrast)) +
  stat_halfeye(slab_fill = 'firebrick4',
               slab_alpha = 0.8, 
               point_interval = median_hdi) + 
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_classic()
```

That will give you a plot like this:

![complex-model-con](https://github.com/jed709/jed709.github.io/assets/87210399/e36cfa10-9864-4a44-beab-0cc99a126764)

Here, I added a line at 0 so it's easy to tell if a contrast is credible or not. If this was the comparison we were interest it, a plot like this would do a very good job communicating it. 

Hopefully, this tutorial has conveyed that you can do a great deal with Bayesian parameter estimation. Unlike other approaches, you are not at all limited in the comparisons you can make. For even more complex designs, the general principles demonstrated here still apply. Hopefully you can transfer this knowledge to your own research - good luck!
