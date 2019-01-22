---
knit: "bookdown::preview_chapter"
---



# Setting up an evaluation

We describe our methods in this chapter.

## Statistical power, survey data and just doing it...

@Bloom1995
Optimal Design software: http://hlmsoft.net/od/

### Blocking

Consider a simple RCT regression ($i=1,...N$), with average treatment effect $\gamma$ and idiosyncractic treatment effects $\theta_{i}$, $E[\theta_{i}]=0$:
$$
Y_{i}=X_{i}\eta+D_{i}\gamma+D_{i}\theta_{i}+\epsilon_{i}
$$
Including covariates (*ex post*) or ``blocking" (*ex ante*), will reduce $\text{var}[U_{i}]$ and therefore improve precision:
$$
\text{var}[\widehat{\gamma}]=\frac{\text{var}[\epsilon_{i}]}{N\text{var}[D_{i}](1-R^{2}_{D|X})}=\frac{\text{var}[\epsilon_{i}]}{N\text{var}[D_{i}]}
$$

### Repeated observations

Add a time dimension $t$:
$$
Y_{it}=X_{it}\eta+D_{it}\gamma+D_{it}\theta_{i}+\lambda_{i}+\epsilon_{it}
$$
Gain in efficiency **as long as** $\theta_{i}\perp\lambda_{i}$ 
Standard panel data formula: 
$$
\text{var}[\widehat{\gamma}^{B}]-\text{var}[\widehat{\gamma}^{W}]=2N^{-1}\text{var}[\lambda_{i}]
$$
Disadvantage of within-subject designs: treatment effects may be time-dependent (history / learning):
$$
Y_{i2}=...+D_{i1}\gamma_{1}+D_{i2}\gamma_{2}+D_{i1}\theta_{i1}+D_{i2}\theta_{i2}...
$$

### Glossary and intuition


**Type I error** ($p-$ value or significance level): probability of falsely rejecting the null hypothesis

**Type II error**: probability of falsely not rejecting the null hypothesis = $1-$ power
 
**Power**: probability of correctly rejecting the null

**Effect size**: the magnitude of the treatment effect that you want to be able to detect

**Notation**: $Y_{ij}|X_{i}\sim N(\mu_{j},\sigma^{2}_{j}),j=0,1$; covariates $X_{i}$ have been "partialled out"

**Intuition**: $\text{H}_0:\mu_{0}=\mu_{1}\ +$  explicitly specify $\text{H}_1:\mu_{0}-\mu_{1}=\delta=$ **Minimum Detectable Effect Size (MDES)**

@Deaton97
@List2011
@Ranganathan2015

## Stuff to keep in mind:  Survey bias and Hawthorne effects

@Zwane2011
@Levitt2011
