//
//  MIXTURE PRIOR
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
  real mix;
  // vague
  mix = 0;
  for (i in 1:S) {
    mix = mix + WEIGHTS[i] * exp(beta_lpdf(theta | YSUM[i], N[i] - YSUM[i]));
  }

  //prior
  target += log(mix);
  //likelihood
  YSUMCUR ~ binomial(NCUR, theta);
}
