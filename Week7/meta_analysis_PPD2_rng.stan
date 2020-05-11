functions {
  matrix meta_analysis_PPD2_rng(int S, vector se) {
    int N = rows(se);
    matrix[S, N] draws;
    for (s in 1:S) {
      real mu = normal_rng(0, 1);
      real tau = fabs(cauchy_rng(0, 0.3));
      for (n in 1:N) {
        real epsilon = normal_rng(0, 1);
        real delta = mu + tau * epsilon;
        real y = normal_rng(delta, se[n]);
        draws[s, n] = y;
      }
    }
    return draws;
  }
}
