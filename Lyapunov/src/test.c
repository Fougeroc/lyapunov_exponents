#include "lyapunov_exponents.h"

typedef struct{
  int * entier;
}new_type;


int testtt(int ** tab)
{
  return(tab[1][1]);
}

void main() {
  size_t i;
  size_t test;
  test = 0;
  i = 1;
  ++i;
  if (test)
    printf("%i", i);
  size_t degree = 5;
  size_t perm[degree];
  permutation *sigma = new_permutation(degree);
  permutation *tau = new_permutation(degree);
  cyclic_permutation(1, sigma);
  inverse_permutation(sigma, tau);
  degree = 2;
  int gp_perm[6] = {0, 1, 1, 2, 2, 0};
  int gp_twin[6] = {5, 2, 1, 4, 3, 0};
  permutation *sigma_0 = new_permutation(degree);
  permutation *sigma_1 = new_permutation(degree);
  permutation *sigma_2 = new_permutation(degree);
  cyclic_permutation(0, sigma_0);
  cyclic_permutation(1, sigma_1);
  cyclic_permutation(1, sigma_2);
  generalized_permutation * gp = new_generalized_permutation(gp_perm, gp_twin, 3, 3);
  permutation * sig[3] = {sigma_0, sigma_1, sigma_2};
  quad_cyclic_cover *qcc = new_quad_cyclic_cover(gp, sig, degree, 3);
  print_quad_cyclic_cover(qcc);
  double theta[3];
  init_GS(3);
  lyapunov_exponents_H_plus(qcc, theta, 10000);
  print_quad_cyclic_cover(qcc);
}
