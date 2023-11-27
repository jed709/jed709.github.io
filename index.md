---
title: Home
layout: home
---

---
Introduction to Bayesian Parameter Estimation
---


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
