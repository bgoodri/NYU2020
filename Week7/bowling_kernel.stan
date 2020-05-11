functions { /* bowling_kernel.stan */
  real bowling_kernel(vector pi, vector a, 
                      int [ , ] x1_x2) {
    real log_like = 0; // categorical 
    real log_prior = dirichlet_lpdf(pi | a);
    for (n in 1:dims(x1_x2)[1]) {
      int x1 = x1_x2[n, 1];
      log_like += log(pi[x1 + 1]);
      if (x1 < 10) {  // not a strike
        int np1 = 10 - x1 + 1;
        vector[np1] pi_ = pi[1:np1] 
                        / sum(pi[1:np1]);
        int x2 = x1_x2[n, 2];
        log_like += log(pi_[x2 + 1]);
      }
    }
    return log_prior + log_like;
  }
}
