source(here::here('R', '00-utils.R'))
# ------------------------------------------------------------------------------

# Get sigma for desired_cv
f <- function(sigma, q, desired_cv) {
  print(sigma)
  set.seed(1)
  x <- rgengamma(3e6, mean = 1, sigma = sigma, Q = q)
  cv <- sd(x)/mean(x)
  (cv - desired_cv)^2
}

out <- optimize(f, c(0.0001, 3), tol = 0.00001, desired_cv = 2, q = 0.5)
out

x <- rgengamma(3e6, mean = 1, sigma = out$minimum, Q = 0.5)
sd(x)/mean(x)

hist(x, breaks = 200, xlim = c(0, 20))

# Check that for CV equals 2 and sigma = Q, that rgengamma matches gamma
f2 <- function(sigma, desired_cv) {
  print(sigma)
  set.seed(1)
  x <- rgengamma(3e6, mean = 1, sigma = sigma, Q = sigma)
  cv <- sd(x)/mean(x)
  (cv - desired_cv)^2
}

out <- optimize(f2, c(0.0001, 3), tol = 0.00001, desired_cv = 2)
out

x <- rgengamma(3e6, mean = 1, sigma = out$minimum, Q = out$minimum)
sd(x)/mean(x)

par(mfrow = c(1, 2))
hist(x, breaks = 200, xlim = c(0, 20), main = "rgengamma, CV = 2, sigma = Q")

# Get phi for gamma
# phi = 1 / cv^2
# phi * cv^2 = 1
# cv^2 = 1 / phi
# cv = sqrt(1 / phi)

cv <- 2
phi <- 1 / cv^2
#cv <- 1 / sqrt(phi)
x <- rgamma(3e6, shape = phi, scale = 1 / phi)

mean(x)
sd(x)/mean(x)

hist(x, breaks = 200, xlim = c(0, 20), main = "gamma, CV = 2")

# Lognormal
mconv <- function(m, sd) log(m) - 0.5 * log(1 + ((sd^2)/(m^2)))
sdconv <- function(m, sd) (log(1 + ((sd^2)/(m^2))))^0.5
mu <- 1
cv <- 2

mconv(mu, mu * cv)
sdconv(mu, mu * cv)

x <- rlnorm(3e6, mconv(mu, mu * cv), sdconv(mu, mu * cv))

x <- rlnorm(3e6, 1, 0.8)
sd(x)/mean(x)

hist(x, breaks = 200, xlim = c(0, 20))
plot(density(x), xlim = c(0, 50))

out <- optimize(f, c(0.0001, 3), tol = 0.00001, desired_cv = 2, q = 0.000001)

x <- rgengamma(3e6, mean = 1, sigma = out$minimum, Q = 0.000001)
mean(x)
sd(x)/mean(x)

# hist(x, breaks = 200, xlim = c(0, 20))

plot(density(x), xlim = c(0, 50))

f3 <- function(phi, desired_cv) {
  print(phi)
  set.seed(1)
  mu <- 1
  x <- rlnorm(3e6, meanlog = log(mu) - 0.5 * phi^2, sdlog = phi)
  cv <- sd(x)/mean(x)
  (cv - desired_cv)^2
}

out <- optimize(f3, c(0.0001, 3), tol = 0.00001, desired_cv = 2)

out

sdconv <- function(m, sd) (log(1 + ((sd^2)/(m^2))))^0.5
mu <- 9
cv <- 3
sdconv(mu, mu * cv) # to get input phi


# Tweedie
f4 <- function(phi, mu, desired_cv, p = 1.5) {
  print(phi)
  set.seed(1)
  x <- fishMod::rTweedie(1e6, mu = mu, phi = phi, p = p)
  cv <- sd(x)/mean(x)
  (cv - desired_cv)^2
}

x <- fishMod::rTweedie(1e6, mu = 1.45, phi = 1.2, p = 1.7)
mean(x)

out <- optimize(f4, c(0.0001, 3), tol = 0.00001, desired_cv = 2, mu = 1)

out

p <- 1.5
phi <- 3
mu <- 1
.variance <- phi * mu^p
.sd <- sqrt(phi * mu^p)

.cv <- 2
p <- 1.5
mu <- 1
# .cv <- sqrt(phi * mu^p) / mu
phi <- ((mu * .cv)^2)/(mu^p)

# binomial
x <- rbinom(3e6, p = 0.7)
sd(x) / mean(x)


get_phi <- function(cv, family, mu, p, Q) {
  # Numeric solution to get sigma for desired_cv
  f <- function(sigma, q, desired_cv = cv) {
    print(sigma)
    set.seed(1)
    x <- rgengamma(3e6, mean = mu, sigma = sigma, Q = Q)
    cv <- sd(x)/mean(x)
    (cv - desired_cv)^2
  }

  switch(family,
    "gamma" = 1 / cv^2,
    "lognormal" = log(1 + cv^2)^0.5, # simplification of `sdconv()`- independent of mu
    "tweedie" = ((mu * cv)^2) / (mu^p),
    "gengamma" = optimize(f, c(0.0001, 3), tol = 0.00001, desired_cv = cv, q = Q)$minimum
    )
}

mconv <- function(m, sd) log(m) - 0.5 * log(1 + ((sd^2)/(m^2)))
sdconv <- function(m, sd) (log(1 + ((sd^2)/(m^2))))^0.5

.mu <- 1
.Q <- 5
.cv <- 0.83
.p <- 1.5

.phi <- get_phi(cv = .cv, family = "lognormal")
# Lognormal
# works less well as .cv gets large
x <- rlnorm(10e6, log(.mu) - 0.5 * .phi^2, .phi)
sd(x) / mean(x)
.cv
.cv - (sd(x) / mean(x))

.phi <- get_phi(cv = .cv, family = "gamma")
x <- rgamma(3e6, shape = .phi, scale = 1 / .phi)
sd(x) / mean(x)
.cv

.phi <- get_phi(cv = .cv, family = "tweedie", mu = .mu, p = .p)
x <- fishMod::rTweedie(3e6, mu = .mu, phi = .phi, p = .p)
sd(x) / mean(x)
.cv

.phi <- get_phi(cv = .cv, family = "gengamma", mu = .mu, Q = .Q, p = .p)
x <- rgengamma(5e6, mean = .mu, sigma = .phi, Q = .Q)
sd(x) / mean(x)
.cv


# .mu <- 1
# .cv <- 1:10#c(seq(from = 0.1, to = 2, by = 0.2), seq(from = 3, to = 15, by = 1))
# n_iterations <- 5

# # Define a function to generate data for a given cv value
# generate_data <- function(cv_value) {
#   x <- map_df(rep(1, n_iterations), ~ {
#     .phi <- get_phi(cv = cv_value, family = "lognormal")
#     x <- rlnorm(10e6, log(.mu) - 0.5 * .phi^2, .phi)
#     sd_mean_ratio <- sd(x) / mean(x)
#     diff <- cv_value - sd_mean_ratio
#     tibble(cv = cv_value, sd_mean_ratio = sd_mean_ratio, diff = diff)
#   })
# }

# results_df <- lapply(.cv, generate_data) |> bind_rows()
# beepr::beep()
# results_df

# ggplot(data = results_df, aes(x = cv, y = diff)) +
#   geom_point()