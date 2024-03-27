Prentice definition (originally from Stacy 1962) of PDF for generalised gamma:

$\left \{ \frac{b}{\Gamma k} \right \}a^{-bk}x^{bk-1}exp{(-(\frac{x}{a})^b)}$ 

for $x > 0$

Wikipedia PDF definition = 
$\frac{p / a^{d}}{\Gamma (d / p)}x^{d - 1}e^{(-\frac{x}{a})^b}$

|              | Stacy & Mihram                                                                                                                         | Prentice                                                                      | Wikipedia                                                       |
| ------------ | -------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- | --------------------------------------------------------------- |
| PDF          | $\frac{\|p\|x^{p\nu-1}exp{(-(\frac{x}{a})^\nu)}}{a^{p\nu}\Gamma(\nu)}$ eqn(1)                                                          | $\left \{ \frac{b}{\Gamma k} \right \}a^{-bk}x^{bk-1}exp{(-(\frac{x}{a})^b)}$ | $\frac{p / a^{d}}{\Gamma (d / p)}x^{d - 1}e^{(-\frac{x}{a})^b}$ |
| *r*th moment | $E(X^r) = \begin{cases} a^r\frac{\Gamma((p\nu + r)/p)}{\Gamma(\nu)}, & r/p > -\nu\\ \infty,  & \textup{otherwise} \end{cases}$ eqn (3) |                                                                               |                                                                 |
| Mean         | $E(X) = \begin{cases} a\frac{\Gamma((p\nu + 1)/p)}{\Gamma(\nu)}, & r/p > -\nu\\ \infty,  & \textup{otherwise} \end{cases}$             |                                                                               | $a\frac{\Gamma((d + 1)/p)}{\Gamma(d/p)}$                        |
| Mode         |                                                                                                                                        |                                                                               | $a(\frac{d-1}{p})^{\frac{1}{p}}$ for $d > 1$, otherwise $0$     |
| Model        |                                                                                                                                        | $\begin{array}{c} y = log(x) \\ log(a) + b^{-1}w \end{array}$                 |                                                                 |

| Stacy to Prentice | Prentice to Wiki  |    Stacy to Wiki    |  Lawless to Wiki  |
| :---------------: | :---------------: | :-----------------: | :---------------: |
|      $p = b$      |     $b = p $      |       $p = p$       | $b = \frac{1}{p}$ |
|     $\nu = k$     | $k = \frac{d}{p}$ | $\nu = \frac{d}{p}$ | $k = \frac{d}{p}$ |
|         -         |     $bk = d$      |          -          |   $\alpha = a$    |


--- 
Lawless 1980: 

Starting with Prentice PDF defined as: 

$\left \{ \frac{\beta}{\Gamma k} \right \}\alpha^{-\beta k}t^{\beta k-1}exp{(-(\frac{t}{\alpha})^\beta)}, ~~~t > 0$

In log link space:

| Substitutions                                              | Lawless to Prentice |
| ---------------------------------------------------------- | :-----------------: |
| $Y = log(T);  ~~~~  e^{y} = t$                             |       $t = x$       |
| $u = log(\alpha); ~~~~~  e^u = \alpha$                     |     $\beta = b$     |
| $b = \frac{1}{\beta}; ~~~~~~~~~~~~~~~ \frac{1}{b} = \beta$ |    $\alpha = a$     |

$\frac{\frac{1}{b}}{\Gamma k}{(e^u)}^{-\frac{1}{b} k}{(e^y)}^{\frac{1}{b} k-1}exp{(-(\frac{e^{y}}{e^u})^\frac{1}{b})}$

$\frac1{b\Gamma k}{(e^u)}^{-\frac{k}{b}}{(e^y)}^{\frac{k}{b}-1}exp{(-(e^{y-u})^\frac{1}{b})}$

$\frac1{b\Gamma k}(e^{-\frac{uk}{b}}e^{\frac{yk}{b}-y}exp{\left(-e^{\frac{y-u}{b}} \right)})$

$\frac1{b\Gamma k}exp[\frac{yk}{b} -\frac{uk}{b} -y + \left(-e^{\frac{y-u}{b}} \right)]$

$\frac1{b\Gamma k}exp[k(\frac{y - u}{b}) -y + \left(-e^{\frac{y-u}{b}} \right)]$ I have an extra -y and I don't know why...



|     | Lawless 1980                                                                                              |
| --- | --------------------------------------------------------------------------------------------------------- |
| PDF | $\frac{1}{b\Gamma(k)}exp\left [ k\left(\frac{y - u}{b}\right) -e^{(y-u)/b}\right ], -\infty < y < \infty$ |


Original special cases from Prentice, where $b = \beta$ if referring to OG Prentice paper, but I am using Lawless 1980 terms: 

| Special cases | Conditions                    |
| ------------- | ----------------------------- |
| Exponential   | $\beta = k = 1$               |
| Weibull       | $k = 1$                       |
| Gamma         | $\beta = 1 = \sigma \sqrt(k) = \frac{\sigma}{Q}$ |
| ~ Lognormal   | $k -> \infty$                 |

$k = \frac{1}{Q^2} = Q^{-2}$

$Q = \frac{1}{\sqrt(k)} = k^{-1/2}$


Not sure what to do with this note or why I wrote this down: 
- With Prentice using Q parameterisation: $\left \{ \frac{b}{\Gamma \frac{1}{Q^{2}}} \right \}a^{\frac{-b}{Q^2}}x^{\frac{b}{Q^2}-1}exp{(-(\frac{x}{a})^b)}$



## Mean <-> Mu

I think that this is only an issue for the PDF? Not for `flexsurv::rgengamma`. 
I do not know about for anything else. 

From flexsurv documentation:

```
# Where mu is the mean
dgengamma(x, mu, sigma, Q=0)  =  dlnorm(x, mu, sigma)
dgengamma(x, mu, sigma, Q=1)  =  dweibull(x, shape=1/sigma, scale=exp(mu))
dgengamma(x, mu, sigma, Q=sigma)  =  dgamma(x, shape=1/sigma^2, rate=exp(-mu) / sigma^2)
```

Stacy & Mihram 1965 definition of the *r*th moment:

$E(X^r) = \begin{cases} 
a^r\Gamma((p\nu + r)/p), & r/p > -\nu\\ 
\infty,  & \textup{otherwise}
\end{cases}$

Note to self, for the mean, $r = 1$ ('moment ordinal in wikipedia table'). 

If we rewrite the equation from Stacy & Mihram 1965 using the notation of Lawless 1980 (and match what is the gengamma MS draft: (i.e., $p = \beta$, $\nu = k$)):

$mean = a \frac{\Gamma (\frac{k\beta + 1}{\beta})}{\Gamma k}$

And then use the parameterisations: 

-  Lawless eqn (2): $u = log(a)$ such that $a = e^u$
- And the Lawless parameterisation in eqn(3) and also using the Q parameterisation: $Q = k^{-1/2}$ also $k = Q^{-2}$ we can get:

| sigma                         | mu                        |
| ----------------------------- | ------------------------- |
| $\sigma = \frac{b}{\sqrt{k}}$ | $\mu = u + b\log(k)$      |
| $\sigma = bQ$                 | $u = \mu - b\log(Q^{-2})$ |

-----> **OOOOH** $\theta = u$ <------

A reminder that $\beta = \frac{1}{b} = \frac{Q}{\sigma}$, 


We get:

$\begin{align}
mean = \frac{\Gamma (\frac{k\beta + 1}{\beta})}{\Gamma k}\cdot e^{\mu - \frac{\sigma\log(Q^{-2})}{Q}} \nonumber \end{align}$
alternatively, to solve for $\mu$,
$\begin{align}
\log(mean) = \log(\frac{\Gamma (\frac{k\beta + 1}{\beta})}{\Gamma k}) + \log(e^{\mu - \frac{\sigma\log(Q^{-2})}{Q}}) \nonumber \\
\log(mean) = \log\Gamma \left(\frac{k\beta + 1}{\beta}\right) - \log(\Gamma k) + \mu - \frac{\sigma\log(Q^{-2})}{Q} \nonumber \\
\log(mean) - \log\Gamma \left(\frac{k\beta + 1}{\beta}\right) + \log(\Gamma k) + \frac{\sigma\log(Q^{-2})}{Q} = \mu  \nonumber \\
\end{align}$


```{r}
k <- q^(-2)
beta <- Q / sigma
# This works
mu <- log(mean) - lgamma((k * beta + 1) / beta) + lgamma(k) + (sigma * log(q^(-2)) /q)
mean <- (gamma((k * beta + 1) / beta) / gamma(k)) * exp(mu - (sigma * log(Q^(-2)) / Q))
```

---
Calculate PDF in R/cpp

```{r}
# In R
Q <- 0.2
sigma <- 0.9
mean <- 0.34
mu <- mu_from_mean(mean = mean, Q = Q, sigma = sigma)
x <- 2

# From flexsurv
y = log(x)
w = (y - mu) / sigma
abs_Q = abs(Q)
qi = 1 / (Q^2)
qw = Q * w

-log(sigma * x) + log(abs_Q) * (1 - 2 * qi) + qi * (qw - exp(qw)) - lgamma(qi)

# Rewrite in terms of Q
-log(sigma * x) + log(abs(Q)) * (1 - 2 * Q^-2) + (Q^-2) * (Q * w - exp(Q * w)) - lgamma(Q^-2)

# Simplify
- log(sigma * x) + log(abs(Q)^(1 - 2 * Q^-2)) + (Q^-2) * (Q * w - exp(Q * w)) - lgamma(Q^-2)
log(abs(Q)^(1 - 2 * Q^-2) / (sigma * x)) + (Q^-2) * (Q * w - exp(Q * w)) - lgamma(Q^-2)
log(abs(Q)^(1 - 2 * Q^-2) / (sigma * x * gamma(Q^-2))) + (Q^-2) * (Q * w - exp(Q * w))
log(abs(Q)^(1 - 2 * Q^-2) / (sigma * x * gamma(Q^-2)) * exp((Q^-2) * (Q * w - exp(Q * w))))

exp(log(abs(Q)^(1 - 2 * Q^-2) / (sigma * x * gamma(Q^-2)) * exp((Q^-2) * (Q * w - exp(Q * w)))))
(abs(Q)^(1 - 2 * Q^-2) / (sigma * x * gamma(Q^-2))) * exp((Q^-2) * (Q * w - exp(Q * w)))

((abs(Q)^(1) / (abs(Q)^ (2 * Q^-2))) / (sigma * x * gamma(Q^-2))) * exp((Q^-2) * (Q * w - exp(Q * w)))
# Use abs(Q) = (Q^2) ^ (1/2)
((abs(Q)^(1) / (((Q^2)^(1/2))^(2 * Q^-2))) / (sigma * x * gamma(Q^-2))) * exp((Q^-2) * (Q * w - exp(Q * w)))
# Simplify exponents
((abs(Q) / (Q^ (2 * (Q^-2)))) / (sigma * x * gamma(Q^-2))) * exp((Q^-2) * (Q * w - exp(Q * w)))
# Match gen gamma MS draft Table 5 PDF equation: 
((abs(Q) * (Q^ -(2 * (Q^-2)))) / (sigma * x * gamma(Q^-2))) * exp((Q^-2) * (Q * w - exp(Q * w)))

# Question for Jim, I don't clearly see where this equation comes from based on Stacy & Mihram, Prentice, and Lawless

```

$\log(\sigma * x) + (1 - 2Q)\log(\left |Q\right |) + Q(Qw - e^{Qw}) - \log\Gamma(Q)$

simplify before exponentiating: 

$\frac{\log(\sigma x\left (|Q\right |^{(1 - 2Q)}))}{\log\Gamma(Q)} + Q(Qw - e^{Qw})$

exponentiate:

$\frac{\sigma x\left (|Q\right |^{(1 - 2Q)})}{\Gamma(Q)} + e^{Q(Qw - e^{Qw})}$

---

From Prentice: transformation of k can give us $Q = \frac{1}{k^2}$, such that
$k = \frac{1}{Q^2} = Q^{-2}$


Model:

$y = \log(x)$

Can also be written as:

$y = \log(a) + \frac{w}{b}$

where density function for $w$ is:
$\left \{ \Gamma (k) \right \}^{-1}exp(kw)exp(-e^w)$

$\frac{1}{\left \{ \Gamma (k) \right \}}exp(kw)exp(-e^w)$

$\frac{1}{\left \{ \Gamma (k) \right \}}e^{kw}e^{-e^w}$


Variance function:





## Scale simulated values to CV

$cv = \frac{sd}{mean}$

**Gamma distribution:**  
$shape = k$  
$mean = k\theta$  
$variance = \sigma^2 = k\theta^2$  
$sd = \sqrt{k\theta^2} = \theta\sqrt{k}$  
$cv = \frac{\theta\sqrt{k}}{k\theta} = \frac{1}{\sqrt{k}}$  

sdmTMB gamma parameterisation:  
s1 = exp(ln_phi(m));        // shape  
$k = phi$
$cv = \frac{1}{\sqrt{phi}}$  
$phi = \frac{1}{cv^2}$

**Lognormal distribution:**

$mean = \exp\left(\ \mu +{\frac {\sigma ^{2}}{2}}\ \right)$

$variance = {\left[\ \exp(\sigma ^{2})-1\ \right]\exp \left(2\mu +\sigma ^{2}\right)\ }$


$cv = \frac{\sqrt{\left[\ \exp(\sigma ^{2})-1\ \right]\exp \left(2\mu +\sigma ^{2}\right)\ }}{\exp\left(\ \mu +{\frac {\sigma ^{2}}{2}}\ \right)}$  


```{r}
mu = 1
sigma2 = 0.5

.mean = exp(mu + sigma2 / 2)
.var = (exp(sigma2) - 1) * exp(2 * mu + sigma2)
.sd = sqrt((exp(sigma2) - 1) * exp(2 * mu + sigma2))

sd / mean

```