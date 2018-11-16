library(raster)
library(rgdal)
library(MASS)

# Load genetic data
gen.data <- data.frame(readxl::read_excel("~/Data for alnus integrative model/Genetic_Data_For_Integrative_Model.xlsx"))
coord.g <- subset(gen.data, select=c(lat, long))
x.g <- qnorm(pmax(gen.data$Cluster.1, gen.data$Cluster.2))
# x.g <- scale(subset(gen.data, select=max(c(Cluster.1, Cluster.2))))
n.g <- nrow(gen.data)
colnames(coord.g) <- c('lat', 'lon')
colnames(x.g) <- c('pr')

# Load sdm data
r <- raster("~/Data for alnus integrative model/FRU_lgmEM.asc.txt")
r.pts <- rasterToPoints(r, spatial = TRUE)
r.df <- data.frame(r.pts)
n.m <- nrow(r.df)
r.df['scaled'] <- (r.df$FRU_lgmEM.asc - min(r.df$FRU_lgmEM.asc))/(diff(range(r.df$FRU_lgmEM.asc)))
r.df$scaled <- qnorm(r.df$scaled)
# Choose k rows
k <- 1000
set.seed(1)
r.df <- r.df[sample(nrow(r.df), k), ]

coord.m <- subset(r.df, select=c(y,x))

x.m <- r.df$scaled
n.m <- nrow(r.df)
colnames(coord.m) <- c('lat', 'lon')
colnames(x.m) <- c('pr')
# r.pts
# # Verify if correctly loaded
# nrow(fru.sdm.table)
# ncol(fru.sdm.table)

plot(r)

# Load Pollen data
pol.data <- data.frame(readxl::read_excel('Pollen_Data_For_Inegrative_Model_v2.xlsx'))
coord.p <- subset(pol.data, select=c(Latitude, Longitude))
x.p <- qnorm(pol.data$Pollen.Data)
n.p <- nrow(pol.data)
colnames(coord.p) <- c('lat', 'lon')
colnames(x.p) <- c('pr')

# Union of all coordinates
n <- n.p + n.m + n.g
coords <- rbind(coord.m, coord.g, coord.p)
# Check any common coordinates
# dim(coords)
# dim(unique(coords))
xx <- rbind(x.m, x.g, x.p)
# x.m <- matrix(0, nrow = n, ncol = 1)

# x.m[1:n.m] =
# M, G, P matrices
M <- matrix(0, nrow = n.m, ncol = n)
M[, 1:n.m] <- diag(n.m)

G <- matrix(0, nrow = n.g, ncol = n)
G[ , (n.m+1):(n.m+n.g)] <- diag(n.g)

P <- matrix(0, nrow = n.p, ncol = n)
P[, (n.m+n.g+1):n] <- diag(n.p)


# Initial values
tt <- 2
tau <- 25
a <- 1
b <- 5
w <- 4
A <- 0.2
B <- 0.3
phi_sd <- 0.2

init0 <- rnorm(3, 0, sd=sqrt(tau))
init1 <- rnorm(3, 1, sd=sqrt(tau))

init.sigma2 <- 1/rgamma(4, shape = a, scale = b)

mu <- matrix(0, nrow = tt, ncol = 1)

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

phi <- matrix(0, nrow = tt, ncol = 1)

beta0[1] <- init0[1]
alpha0[1] <- init0[2]
gamma0[1] <- init0[3]

beta1[1] <- init1[1]
alpha1[1] <- init1[2]
gamma1[1] <- init1[3]


sigma2.m[1] <- init.sigma2[1]
sigma2.g[1] <- init.sigma2[2]
sigma2.p[1] <- init.sigma2[3]
sigma.2[1] <- init.sigma2[4]


mu[1] <- rnorm(1, 0, w)

x <- matrix(0, nrow = n, ncol = tt)

phi[1] <-  runif(1, A, B)
aa = 0.826
# Some helper functions
x.dist <- as.matrix(dist(coords))
minx = min(x.dist)
rangex = diff(range(x.dist))
x.dist <- apply(x.dist, MARGIN = 1, FUN = function(X) (X - minx)/rangex)
V <- function(phi) {return(exp(-x.dist/phi) + aa*diag(1065))}
Sigma.inv <- function(sigma.2, V) {return(sigma.2*solve(V))}

phi_target <- function(phi0, sigma2, x1, mu0){
  if (phi0 >= A && phi0 <= B){
    return((det(V(phi0))^-0.5) * exp(1/(2*sigma2) * (t(x1 - mu0) %*% solve(V(phi0)) %*% (x1 - mu0))))
  }
  return(0)
}

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
# Metropolis hastings algorithm
# MH <- function(old, sigma2, x1, mu0, sd) {
#   proposed <-  old + rnorm(1, mean=0, sd=sd)
#   ratio <- min(1, phi_target(proposed, sigma2, x1, mu0)/phi_target(old, sigma2, x1, mu0))
#   # print(ratio)
#   if (runif(1) <= ratio) { return(proposed) }
#   return(old)
# }

MH <- function(old, sigma2, x1, mu0) {
  proposed <-  old + rnorm(1, mean=0, sd=phi_sd)
  print(proposed)
  a = -1/(2*sigma2) * (t(x1 - mu0) %*% solve(V(proposed), (x1 - mu0)))
  b = -1/(2*sigma2) * (t(x1 - mu0) %*% solve(V(old), x1 - mu0))
  ratio <- min(1, sqrt(det(V(old)) / det(V(proposed))) * exp(a-b))
  if(runif(1) <= ratio) { return(proposed) }
  return(old)
}

for (i in 2:tt) {
# Try for one timestep
# i=2
  mm <- beta1[i-1]/sigma2.m[i-1] * t(M) %*% (x.m - beta0[i-1])
  gg <- alpha1[i-1]/sigma2.g[i-1] * t(G) %*% (x.g - alpha0[i-1])
  pp <- gamma1[i-1]/sigma2.p[i-1] * t(P) %*% (x.p - gamma0[i-1])
  sinv <- Sigma.inv(sigma.2[i-1], V(phi[i-1]))

  mu.tilda <- mu[i-1]*rowSums(sinv) + mm + gg + pp

  mm <- beta1[i-1]/sigma2.m[i-1] * (t(M) %*% M)
  gg <- alpha1[i-1]/sigma2.g[i-1] * (t(G) %*% G)
  pp <- gamma1[i-1]/sigma2.p[i-1] * (t(P) %*% P)

  Sigma.tilda <- makeSymm(solve(sinv + mm + gg + pp))

  # Sample x
  x[,i] <- mvrnorm(mu = Sigma.tilda %*% mu.tilda, Sigma = Sigma.tilda)

  # Sample mu
  k <- n/sigma.2[i-1] + 1/w
  mu[i,] <- rnorm(1, (rep(1,n) %*% sinv %*% x[,i]) / k, sqrt(1/k))

  # Sample sigma.2
  bb <- 0.5*t(x[,i] - mu[i,]) %*% solve(V(phi[i-1])) %*% (x[,i] - mu[i,]) + b
  sigma.2[i,] <- 1/rgamma(1, shape = n/2 + a, scale = bb)

  # Sample sigma2.m, sigma2.g, sigma2.p
  aa <- n.m/2 + a
  bb <- b + 0.5*sum((x.m - beta0[i-1] - beta1[i-1]*x[1:n.m,i])^2)
  sigma2.m[i] <- 1/rgamma(1,shape=aa,scale=bb)

  aa <- n.g/2 + a
  bb <- b + 0.5*sum((x.g - alpha0[i-1] - alpha1[i-1]*x[(n.m+1):(n.m+n.g),i])^2)
  sigma2.g[i] <- 1/rgamma(1,shape=aa,scale=bb)

  aa <- n.p/2 + a
  bb <- b + 0.5*sum((x.p - gamma0[i-1] - gamma1[i-1]*x[(n.m+n.g+1):n,i])^2)
  sigma2.p[i] <- 1/rgamma(1,shape=aa,scale=bb)

  # Sample beta1
  c.k <- 1/n.m * (t(x.m) %*% x[1:n.m,i])
  s.k <- 1/n.m * (t(x[,i]) %*% t(M) %*% M %*% x[,i])
  eta.k <- n.m * s.k / sigma2.m[i] + 1/tau

  mn <- n.m/sigma2.m[i] * (c.k - beta0[i-1]*mean(x[1:n.m,i])) + 1/tau

  beta1[i] <- rnorm(1, mean = mn / eta.k, sd = sqrt(1/eta.k))

  # Sample alpha1

  c.k <- 1/n.g * (t(x.g) %*% x[(n.m+1):(n.m + n.g),i])
  s.k <- 1/n.g * (t(x[,i]) %*% t(G) %*% G %*% x[,i])
  eta.k <- n.g * s.k / sigma2.g[i] + 1/tau

  mn <- n.g/sigma2.g[i] * (c.k - alpha0[i-1]*mean(x[(n.m+1):(n.m + n.g),i])) + 1/tau

  alpha1[i] <- rnorm(1, mean = mn / eta.k, sd = sqrt(1/eta.k))

  # Sample gamma1

  c.k <- 1/n.p * (t(x.p) %*% x[(n.m+n.g+1):n,i])
  s.k <- 1/n.p * (t(x[,i]) %*% t(P) %*% P %*% x[,i])
  eta.k <- n.p * s.k / sigma2.p[i] + 1/tau

  mn <- n.p/sigma2.p[i] * (c.k - gamma0[i-1]*mean(x[(n.m+n.g+1):n,i])) + 1/tau

  gamma1[i] <- rnorm(1, mean = mn / eta.k, sd = sqrt(1/eta.k))

  # Sample beta0
  rho.k <- n.m/sigma2.m[i] + 1/tau
  mn <- n.m/sigma2.m[i]*(mean(x.m) - beta1[i]*mean(x[1:n.m,i]))
  beta0[i] <- rnorm(1, mean = mn/rho.k, sd = 1/sqrt(rho.k))

  # Sample alpha0
  rho.k <- n.g/sigma2.g[i] + 1/tau
  mn <- n.g/sigma2.g[i]*(mean(x.g) - alpha1[i]*mean(x[(n.m+1):(n.m + n.g),i]))
  alpha0[i] <- rnorm(1, mean = mn/rho.k, sd = 1/sqrt(rho.k))

  # Sample gamma0
  rho.k <- n.p/sigma2.p[i] + 1/tau
  mn <- n.p/sigma2.p[i]*(mean(x.p) - gamma1[i]*mean(x[(n.m+n.g+1):n,i]))
  gamma0[i] <- rnorm(1, mean = mn/rho.k, sd = 1/sqrt(rho.k))

  # Sample phi by Metropolis Hastings
  phi[i] <- MH(phi[i-1], sigma.2[i], x[,i], mu[i], phi_sd)
}
