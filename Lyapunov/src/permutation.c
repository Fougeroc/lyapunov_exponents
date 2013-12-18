#include "lyapunov_exponents.h"

permutation *new_permutation(size_t degree)
{
	permutation *sigma = (permutation *) malloc(sizeof(permutation));
	sigma->perm = (size_t *) malloc(degree * sizeof(size_t));
	sigma->degree = degree;
	return sigma;
}

void perm_copy_table(permutation *perm, size_t* tab)
{
  memcpy(perm->perm, tab, perm->degree * sizeof(size_t));
}

void copy_permutation(permutation *perm, permutation *perm_buffer)
{
  /*Check the degrees*/
  if(perm->degree != perm_buffer->degree){
    fprintf(stderr, "The degree of the permutation and the buffer doesn't match !");
    return;
  }
  
  size_t i;
  
  for(i=0; i<perm->degree; ++i)
    perm_buffer->perm[i] = perm->perm[i];
}


void free_permutation(permutation ** sigma)
{
	if((*sigma) != NULL)
	{
		if((*sigma)->perm != NULL)
		{
			free((*sigma)->perm);
			(*sigma)->perm = NULL;
		}
		free(*sigma);
		*sigma = NULL;
	}
}

void print_permutation(permutation *sigma)
{
	size_t i;
	for(i=0; i < sigma->degree; ++i) printf(" %i", (sigma->perm)[i] + 1);
	printf("\n");
}



int check_permutation(permutation *sigma)
{
  size_t i, j;
  size_t seen;

  for(i=0; i < sigma->degree; ++i)
    {
      seen = 0;
      for(j=0; j < sigma->degree; ++j)
	{
	  if(sigma->perm[j] == i)
	    {
	      if(seen) 
		{
		  fprintf(stderr, "%i appears twice in the permutation", i);
		  return 1;
		}
	    }
	  else
	    seen = 1;
	  ++j;
	}
      if(!seen){
	fprintf(stderr, "%i doesn't appear in the permutation", i);
	return 1;
      }
    }
  for(i=0; i < sigma->degree; ++i)
    {
      if(sigma->perm[i] >= sigma->degree || sigma->perm[i] < 0)
	{
	  fprintf(stderr, "%i is not between 0 and %i", i, sigma->degree);
	  return 1;
	}
    }
  return 0;
}
	

void inverse_permutation(permutation *sigma, permutation *perm_buffer)
{ 
  /*Check the degrees*/
  if(sigma->degree != perm_buffer->degree){
    fprintf(stderr, "The degree of the permutation and the buffer doesn't match !");
    return;
  }

  size_t i;

  for(i=0; i < sigma->degree; ++i)
    perm_buffer->perm[sigma->perm[i]] = i;
}

void cyclic_permutation(int n, permutation *perm_buffer)
{
  size_t i;
  size_t degree = perm_buffer->degree;

  for(i=0; i < degree; ++i)
    perm_buffer->perm[i] = (i + n) % degree;
}

void perm_name(interval *inter, permutation *perm_buffer)
{
  size_t degree = inter->lab->sigma->degree;
  if(degree != perm_buffer->degree){
    fprintf(stderr, "The degree of the permutation and the buffer doesn't match !");
    return;
  }
  if(inter->orientation == 1)
    cyclic_permutation(0, perm_buffer);        //return identity
  else
    inverse_permutation(inter->lab->sigma, perm_buffer);
}

void perm_ident_rev(interval *inter, permutation *perm_buffer)
{
  if(inter->orientation == 1)
    copy_permutation(inter->lab->sigma, perm_buffer);
  else
    inverse_permutation(inter->lab->sigma, perm_buffer);
}

void perm_product(permutation *sigma, permutation *tau, permutation *perm_buffer)
{
  if(sigma->degree != tau->degree || sigma->degree != perm_buffer->degree)
    {
      fprintf(stderr, "Product of two permutations of different sets");
      return;
    }

  size_t i;
  for (i=0; i < sigma->degree; ++i)
    perm_buffer->perm[i] = tau->perm[sigma->perm[i]];
}