functions {
  matrix interaction_rng(int S, vector income, vector age) {
    int N = rows(age);
    matrix[S, N] draws;
    for (s in 1:S) {
      real alpha = normal_rng(0, 2.5); // intercept
      real beta = normal_rng(0, 2);    // slope on income
      
      real gamma_0 = normal_rng(0, 1);
      real gamma_1 = normal_rng(0, 1);
      vector[N] gamma = gamma_0 + gamma_1 * income;
      
      vector[N] eta = alpha + beta * income + gamma .* age; // log-odds
      // same as: eta = alpha + beta * income + gamma_0 * age + gamma_1 * (income .* age);
      draws[s, ] = to_row_vector(bernoulli_logit_rng(eta));
    }
    return draws;
  }
}
