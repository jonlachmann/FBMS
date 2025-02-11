# Title     : TODO
# Objective : TODO
# Created by: jonlachmann
# Created on: 2021-02-25

context("MJMCMC")

test_that("Testing MJMCMC algorithm", {
  # Dummy test likelihood function
  loglik.tester <- function (y, x, model, complex, params) {
    lmmod <- lm.fit(x[, model, drop = FALSE], y)
    n <- length(y)
    k <- length(lmmod$coefficients)
    rss <- sum(resid(lmmod)^2)
    aic <- n * log(rss / n) + 2 * k
    return(list(crit = aic, coefs = lmmod$coefficients))
  }

  data <- matrix(rnorm(600), 100)
  resm <- mjmcmc(data, loglik.tester)
  summary(resm, labels = c("a", "b", "c", "d", "e"))
  plot(resm)
  predm <- predict(resm, data[, -1, drop = FALSE])

  resg <- gmjmcmc(data, loglik.tester, NULL, c("p0", "exp_dbl"))
  summary(resg)
  plot(resg)
  prediction <- predict(resg, data[, -1, drop = FALSE])

  respm <- mjmcmc.parallel(2, 2, data, loglik.tester)
  summary(respm)
  plot(respm)
  pred_pm <- predict(respm, data[, -1, drop = FALSE])

  respg <- gmjmcmc.parallel(2, 2, NULL, data, loglik.tester, NULL, c("p0", "exp_dbl"))
  summary(respg)
  plot(respg)
  pred_pg <- predict(respg, data[, -1, drop = FALSE])

  fbms_result <- fbms(
   X1 ~ .,
   family = "gaussian",
   method = "gmjmcmc.parallel",
   data = data.frame(matrix(rnorm(600), 100)),
   transforms = c("sin","cos"),
   P = 10,
   runs = 1,
   cores = 1
  )

  # Dummy expect to run the test
  expect_true(is.list(fbms_result))
})

