#include <stdlib.h>
#include "state.h"
#include "binhash.h"

sim_state_t* alloc_state(int n)
{
    sim_state_t* s = (sim_state_t*) calloc(1, sizeof(sim_state_t));
    s->n     = n;
    s->part  = (particle_t*) calloc(n, sizeof(particle_t));
    s->hash  = (particle_t**) calloc(HASH_SIZE, sizeof(particle_t*));
    return s;
}

void free_state(sim_state_t* s)
{
    free(s->hash);
    free(s->part);
    free(s);
}

// Function to compare particles for qsort
int compPart(const void *a, const void *b) {
  int ah = ((struct particle_t*) a)->hind;
  int bh = ((struct particle_t*) b)->hind;

  if (ah < bh) {
    return -1;
  }
  else if (ah == bh) {
    return 0;
  }
  else return 1;
}
