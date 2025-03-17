# eggamma
Implementation of the extended generalized gamma distribution.

# Extended generalized gamma distribution
The distribution implemented is due to Prentice (1974), although the latter considers a log-gamma (log generalized gamma, really) distribution. Note that the implementation proposed here uses the same parametrization as that of the gamlss package; see Section~19.4.3 of Rigby, Stasinopoulos, Heller, & De Bastiani (2019). The parameters are named mu, sigma and nu.

The functions provided make use of the Taylor series of (x/mu)^nu around nu=0 and Stirling's asymptotic formulas for lgamma(xi), digamma(xi), and trigamma(xi), where xi=1/(sigma*nu)^2. They are designed to improve numerical stability around nu=0. The method is described in a paper by myself (S. Perreault) to be submitted shortly. Stay tuned for more details (i.e., a link to the preprint).

# References
Prentice, R. L. (1974). A Log Gamma Model and Its Maximum Likelihood Estimation. Biometrika. 61(3):539--544.  
Stacy, E. W., & Mihram, G. A. (1965). Parameter estimation for a generalized gamma distribution. Technometrics, 7(3):349--358.  
Rigby, R. A., Stasinopoulos, M. D., Heller, G. Z., & De Bastiani, F. (2019). Distributions for modeling location, scale, and shape: Using GAMLSS in R. Chapman and Hall/CRC.

