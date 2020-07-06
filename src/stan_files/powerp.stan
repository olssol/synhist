//
//  Power prior
//

data {
  //existing data
  int<lower  = 1>  S;
  real<lower = 0>  N[S];
  real<lower = 0>  YSUM[S];

  int<lower  = 1>  NCUR;
  int<lower  = 0>  YSUMCUR;

  real<lower = 0, upper = 1>  WEIGHTS[S];
}

parameters {
  real<lower = 0, upper = 1> theta;
}

model {
  // prior S=1: vague
  for (i in 1:S) {
    target += WEIGHTS[i] * beta_lpdf(theta | YSUM[i], N[i] - YSUM[i]);
  }

  //likelihood
  YSUMCUR ~ binomial(NCUR, theta);
}
