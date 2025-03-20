nus <- seq(-.1,.1,.01)
p <- .001
mu <- mu[1]
sigma <- sigma[1]
nu <- seq(-1,1,.0001)
nu <- nu[nu != 0]
yy <- qggamma(p[1], mu[1], sigma[1], nu)
plot(nu, yy, type="l")

xi_inv <- (sigma*nu)^2
plot(nu, gamlss.dist::qGG(p, mu, sigma, nu), type="l")
lines(nu, qggamma(p, mu, sigma, nu), lty=2, col=2)

