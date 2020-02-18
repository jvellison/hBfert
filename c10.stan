data {
matrix[30,3] X;
int N[30,40];
matrix[30,40] W;
matrix[30,30] cov;
}
transformed data {
matrix[30,40] logW;
logW = log(W);
}
parameters {
vector[3] beta;
vector<lower=0>[4] sigma;
vector<lower=0>[3] rho;
vector[30] eta;
matrix[30,39] epsilonrest;
}
transformed parameters {
matrix[30,40] logmu;
matrix[30,40] mu;
matrix[30,40] epsilon;
for (i in 1:30)
epsilon[i,(40-i)] = cov[i,]*eta;
for (i in 1:30) {
  epsilon[i,1:(39-i)] = epsilonrest[i,1:(39-i)];
  epsilon[i,(41-i):40] = epsilonrest[i,(40-i):39];
}
for (i in 1:40)
logmu[,i] = X*beta + logW[,i] + epsilon[,i];
mu = exp(logmu);
}
model {
for (i in 1:30)
N[i,1:(40-i)] ~ poisson(mu[i,1:(40-i)]);
eta ~ normal(0,sigma[1]);
epsilon[,1] ~ normal(0,sigma[2]);
for (j in 2:5)
  epsilon[,j] ~ normal(rho[1]*epsilon[,(j-1)],sigma[3]);
for (k in 6:40)
  epsilon[,k] ~ normal(rho[2]*epsilon[,(k-1)]+
  rho[3]*(epsilon[,(k-1)]+(0.033333)*(10*epsilon[,(k-1)]-epsilon[,(k-2)]-
  2*epsilon[,(k-3)]-3*epsilon[,(k-4)]-4*epsilon[,(k-5)])),sigma[4]);
beta ~ normal(0,30);
sigma ~ normal(0,0.25);
rho ~ normal(0,0.5);
}
