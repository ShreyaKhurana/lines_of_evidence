library(raster)
library(rgdal)
library(MASS)

# Load genetic data
gen.data <- data.frame(readxl::read_excel("Data for alnus integrative model/Genetic_Data_For_Integrative_Model.xlsx"))
coord.g <- subset(gen.data, select=c(lat, long))

n.g <- 47
set.seed(2)
colnames(coord.g) <- c('lat', 'lon')

# Load sdm data
r <- raster("Data for alnus integrative model/FRU_lgmEM.asc.txt")
r.pts <- rasterToPoints(r, spatial = TRUE)
r.df <- data.frame(r.pts)
n.m <- 500
# Choose n.m rows
set.seed(3)
r.df <- r.df[sample(nrow(r.df), n.m), ]
coord.m <- subset(r.df, select=c(y,x))

colnames(coord.m) <- c('lat', 'lon')

# Load Pollen data
pol.data <- data.frame(readxl::read_excel('Data for alnus integrative model/Pollen_Data_For_Inegrative_Model_v2.xlsx'))
pol.data <- subset(pol.data, Pollen.Data != 0)
coord.p <- subset(pol.data, select=c(Latitude, Longitude))
n.p <- 13
colnames(coord.p) <- c('lat', 'lon')

# Union of all coordinates
n <- n.p + n.m + n.g
coords <- rbind(coord.m, coord.g, coord.p)

# M, G, P matrices
M <- diag(n)

G <- M
P <- M


# Initial values
tt <- 100000
tau <- 25
A <- 0.1
B <- 5

mu <- 0

# Arrays to hold samples
beta0 <- matrix(0, nrow = tt, ncol = 1)
beta1 <- matrix(0, nrow = tt, ncol = 1)
alpha0 <- matrix(0, nrow = tt, ncol = 1)
alpha1 <- matrix(0, nrow = tt, ncol = 1)
gamma0 <- matrix(0, nrow = tt, ncol = 1)
gamma1 <- matrix(0, nrow = tt, ncol = 1)

sigma2.m <- matrix(0, nrow = tt, ncol = 1)
sigma2.g <- matrix(0, nrow = tt, ncol = 1)
sigma2.p <- matrix(0, nrow = tt, ncol = 1)
sigma.2 <- matrix(0, nrow = tt, ncol = 1)

phi <- matrix(0.25, nrow = tt, ncol = 1)

beta0[1] <- 0
alpha0[1] <- 0
gamma0[1] <- 0

beta1[1] <- 1
alpha1[1] <- 1
gamma1[1] <- 1


sigma2.m[1] <- 1
sigma2.g[1] <- 1
sigma2.p[1] <- 1
sigma.2[1] <- 1


x <- matrix(0, nrow = n, ncol = tt)

corr = 0.826
# Some helper functions
# Scaling the distance matrix
x.dist <- as.matrix(dist(coords))
minx = min(x.dist)
rangex = diff(range(x.dist))
x.dist <- apply(x.dist, MARGIN = 1, FUN = function(X) (X - minx)/rangex)
V <- function(phi) {return(exp(-x.dist/phi) + corr*diag(n))}
Sigma.inv <- function(sigma.2, V) {return(solve(sigma.2*V))}

# Simulate actual data here, fix slopes at 3 and intercepts at 1
beta0.a <- 1
beta1.a <- 3
alpha0.a <- 1
alpha1.a <- 3
gamma0.a <- 1
gamma1.a <- 3
phi.a <- 0.25
sigma2.a <- 6
sigma2g.a <- 5
sigma2m.a <- 4
sigma2p.a <- 8

x.a <- mvrnorm(mu = rep(0,n), Sigma = sigma2.a*exp(-x.dist/phi.a) )
x.m <- beta0.a + beta1.a*x.a + rnorm(nrow(M), mean = 0, sd = sqrt(sigma2m.a))
x.g <- alpha0.a + alpha1.a*x.a + rnorm(nrow(G), mean = 0, sd = sqrt(sigma2g.a))
x.p <- gamma0.a + gamma1.a*x.a + rnorm(nrow(P), mean = 0, sd = sqrt(sigma2p.a))

n.m <- length(x.m)
n.g <- length(x.g)
n.p <- length(x.p)

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

ptm <- proc.time()
for (i in 2:tt) {
  x[,i] <- x.a
  
  # Sample sigma2.m, sigma2.g, sigma2.p
  aa <- n.m/2 + a
  bb <- b + 0.5*sum((x.m - beta0[i-1] - beta1[i-1]*x[,i])^2)
  sigma2.m[i] <- 1/rgamma(1,aa,bb)
  
  aa <- n.g/2 + a
  bb <- b + 0.5*sum((x.g - alpha0[i-1] - alpha1[i-1]*x[,i])^2)
  sigma2.g[i] <- 1/rgamma(1,aa,bb)
  
  aa <- n.p/2 + a
  bb <- b + 0.5*sum((x.p - gamma0[i-1] - gamma1[i-1]*x[,i])^2)
  sigma2.p[i] <- 1/rgamma(1,aa,bb)
  
  # Sample beta1
  c.k <- 1/n.m * (t(x.m) %*% x[,i])
  s.k <- 1/n.m * (t(x[,i]) %*% t(M) %*% M %*% x[,i])
  eta.k <- n.m * s.k / sigma2.m[i] + 1/tau
  
  mn <- n.m/sigma2.m[i] * (c.k - beta0[i-1]*mean(x[,i])) + 1/tau
  
  beta1[i] <- rnorm(1, mean = mn / eta.k, sd = sqrt(1/eta.k))
  
  # Sample alpha1
  
  c.k <- 1/n.g * (t(x.g) %*% x[,i])
  s.k <- 1/n.g * (t(x[,i]) %*% t(G) %*% G %*% x[,i])
  eta.k <- n.g * s.k / sigma2.g[i] + 1/tau
  
  mn <- n.g/sigma2.g[i] * (c.k - alpha0[i-1]*mean(x[,i])) + 1/tau
  
  alpha1[i] <- rnorm(1, mean = mn / eta.k, sd = sqrt(1/eta.k))
  
  # Sample gamma1
  
  c.k <- 1/n.p * (t(x.p) %*% x[,i])
  s.k <- 1/n.p * (t(x[,i]) %*% t(P) %*% P %*% x[,i])
  eta.k <- n.p * s.k / sigma2.p[i] + 1/tau
  
  mn <- n.p/sigma2.p[i] * (c.k - gamma0[i-1]*mean(x[,i])) + 1/tau
  
  gamma1[i] <- rnorm(1, mean = mn / eta.k, sd = sqrt(1/eta.k))
  
  # Sample beta0
  rho.k <- n.m/sigma2.m[i] + 1/tau
  mn <- n.m/sigma2.m[i]*(mean(x.m) - beta1[i]*mean(x[,i]))
  beta0[i] <- rnorm(1, mean = mn/rho.k, sd = 1/sqrt(rho.k))
  
  # Sample alpha0
  rho.k <- n.g/sigma2.g[i] + 1/tau
  mn <- n.g/sigma2.g[i]*(mean(x.g) - alpha1[i]*mean(x[,i]))
  alpha0[i] <- rnorm(1, mean = mn/rho.k, sd = 1/sqrt(rho.k))
  
  # Sample gamma0
  rho.k <- n.p/sigma2.p[i] + 1/tau
  mn <- n.p/sigma2.p[i]*(mean(x.p) - gamma1[i]*mean(x[,i]))
  gamma0[i] <- rnorm(1, mean = mn/rho.k, sd = 1/sqrt(rho.k))
  

}

print(proc.time() - ptm)

# See the trace plots

par(mfrow = c(2,3))
plot(alpha0, type = "l")
plot(beta0,type="l")
plot(gamma0,type="l")
plot(alpha1,type="l")
plot(beta1,type="l")
plot(gamma1,type="l")

par(mfrow = c(2,2))
plot(sigma2.m,type="l")
plot(sigma2.g,type="l")
plot(sigma2.p,type="l")

mean(beta0[-(1:25000)])
mean(alpha0[-(1:25000)])
mean(gamma0[-(1:25000)])

mean(beta1[-(1:25000)])
mean(alpha1[-(1:25000)])
mean(gamma1[-(1:25000)])

mean(sigma2.m[-(1:25000)])
mean(sigma2.g[-(1:25000)])
mean(sigma2.p[-(1:25000)])
mean(sigma.2[-(1:1000)])

#par(mfrow = c(2,1))
#plot(mu,type="l")
plot(sigma.2,type="l")
#plot(phi,type="l")

# x.final <- apply(x[,25000:tt],1,mean)
#prob.map = pnorm(x.final)
# hist(x.final-x.a)
#hist(prob.map)
#cr <- colorRampPalette(c("black", "red"))
#colors <- rgb(cr(prob.map / max(prob.map)), maxColorValue = 255)

#plot(coords[,1], coords[,2], col=colors, xlab="lat", ylab="lon", main="Probability map of Fructiosa ")
alpha0 = alpha0[-(1:25000)]
alpha1 = alpha1[-(1:25000)]
beta0 = beta0[-(1:25000)]
beta1 = beta1[-(1:25000)]
gamma0 = gamma0[-(1:25000)]
gamma1 = gamma1[-(1:25000)]
sigma2.g = sigma2.g[-(1:25000)]
sigma2.m = sigma2.m[-(1:25000)]
sigma2.p = sigma2.p[-(1:25000)]


par(mfrow = c(2,3))
acf(alpha0)
acf(alpha1)
acf(beta0)
acf(beta1)
acf(gamma0)
acf(gamma1)


par(mfrow= c(2,2))
acf(sigma2.m)
acf(sigma2.g)
acf(sigma2.p)

plot(beta0,type="l")
plot(gamma0,type="l")
plot(alpha1,type="l")
plot(beta1,type="l")
plot(gamma1,type="l")