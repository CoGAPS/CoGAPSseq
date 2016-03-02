# model for normal data known variance, unknown mean
set.seed(1)

# number of iterations
iter <- 5000

# allocate memory for chain
chain <- numeric(iter)

# generate data
D <- rnorm(1000, 
           # parameter of interest
           mean=15, 
           # known
           sd=4)

# priors
s.20 <- 1000
mu.0 <- 0

# data info
n <- length(D)
x.bar <- mean(D)
s.2 <- 4 ^ 2

lambda <- 1 / s.2
lambda.0 <- 1 / s.20
lambda.n <- lambda.0 + n * lambda

# full conditionals from here
# http://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf

for (s in 1:iter) {
    mu.n <- (x.bar * n * lambda + mu.0 * lambda.0) / lambda.n

    chain[s] <- rnorm(1, mu.n, sqrt(1 / lambda.n))
}

# converges immediately, but let's
# remove burnins and thin
chain <- chain[1001:5000]
chain <- chain[1:4000 %% 4 == 0]

# check convergence
ts.plot(chain)
mean(chain)
sd(chain)

# now let's generate some poisson data
D <- rpois(1000, 15)

# allocate memory for chain
chain <- numeric(iter)

# data info
n <- length(D)
x.bar <- mean(D)
s.2 <- 15

lambda <- 1 / s.2
lambda.0 <- 1 / s.20
lambda.n <- lambda.0 + n * lambda

# we keep two running values for mu.n
mu.n <- numeric(2)

# track acceptance ratios 
log.r.track <- numeric(s-1)

# generate first observation from prior, transform, and calculate mu_n
chain[1] <- rnorm(1, mu.0, sqrt(1 / lambda.0))
chain[1] <- chain[1] ^ 2
mu.n[1] <- (x.bar * n * lambda + mu.0 * lambda.0) / lambda.n

for (s in 2:iter) {
    # calculate posterior normal mean
    mu.n[2] <- (x.bar * n * lambda + mu.0 * lambda.0) / lambda.n

    # propose value
    z <- rnorm(1, mu.n, sqrt(1 / lambda.n))

    # transform proposed value (has to be nonnegative)
    z <- z ^ 2

    # acceptance probability
    log.r <- (sum(dpois(D, z, log=TRUE)) +
              dnorm(chain[s], mu.n[2], sqrt(1 / lambda.n), log=TRUE)) -
             (sum(dpois(D, chain[s-1], log=TRUE)) +
              dnorm(z, mu.n[1], sqrt(1 / lambda.n), log=TRUE))

    log.r.track[s-1] <- log.r

    if (log(runif(1)) < log.r) {
        chain[s] <- z
    } else {
        chain[s] <- chain[s-1]
    }

    # reset mu.n values
    mu.n[1] <- mu.n[2]
}

# remove burnins and thin
chain <- chain[1001:5000]
chain <- chain[1:4000 %% 4 == 0]

# check convergence
ts.plot(chain)
mean(chain)
sd(chain)
