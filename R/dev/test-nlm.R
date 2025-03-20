library(eggamma)

n <- 5e5
mu <- 1.25
sigma <- .25
nu <- 0
x <- gamlss.dist::rGG(n, mu, sigma, nu)
links <- c(mu="log", sigma="log", nu="identity")
exp12 <- \(x) c(exp(x[1:2]), x[3])
log12 <- \(x) c(log(x[1:2]), x[3])

opt <- nlm(\(p0) nllggamma(p0, x, 
                           w = 1, links=links, 
                           control = expansionControl(epsilon = c(.1,.5))),
           p = c(0,0,0), fscale = n, ndigit = 10, steptol = 1e-10, hessian = T,    
           stepmax = 1, print.level = 0)
opt$code
opt$estimate; 
rbind(exp12(opt$estimate), c(mu, sigma, nu))

opt$gradient; opt$hessian
