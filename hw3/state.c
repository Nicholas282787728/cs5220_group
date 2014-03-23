#include <stdlib.h>
#include "state.h"
#include "binhash.h"

sim_state_t* alloc_state(int n)
{
    sim_state_t* s = (sim_state_t*) calloc(1, sizeof(sim_state_t));
    s->n     = n;
    s->part  = (particle_t*) calloc(n, sizeof(particle_t));
    s->hash  = (particle_t**) calloc(HASH_SIZE, sizeof(particle_t*));
    s->proc_part = (particle_t**) calloc(MAX_PROC, sizeof(particle_t*));
    s->proc_hash = (particle_t***) calloc(MAX_PROC, sizeof(particle_t**));

    for (int i = 0; i < MAX_PROC; i++) {
      s->proc_part[i] = (particle_t*) calloc(n, sizeof(particle_t));
      s->proc_hash[i] = (particle_t**) calloc(HASH_SIZE, sizeof(particle_t*));
    }

    return s;
}

void free_state(sim_state_t* s)
{
    free(s->hash);
    free(s->part);
    for (int i = 0; i < MAX_PROC; i++) {
      free(s->proc_part[i]);
      free(s->proc_hash[i]);
    }
    free(s->proc_part);
    free(s->proc_hash);
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
