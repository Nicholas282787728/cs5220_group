#include <string.h>
#include <math.h>
#include <stdio.h>

#include "zmorton.h"
#include "binhash.h"

/*@q
 * ====================================================================
 */

/*@T
 * \subsection{Spatial hashing implementation}
 * 
 * In the current implementation, we assume [[HASH_DIM]] is $2^b$,
 * so that computing a bitwise of an integer with [[HASH_DIM]] extracts
 * the $b$ lowest-order bits.  We could make [[HASH_DIM]] be something
 * other than a power of two, but we would then need to compute an integer
 * modulus or something of that sort.
 * 
 *@c*/

#define HASH_MASK (HASH_DIM-1) // HASH_DIM is 0x10

unsigned particle_bucket(particle_t* p, float h)
{
    // We add 16 such that negative numbers are avoided. 
    unsigned ix = (p->x[0]/h + 16);
    unsigned iy = (p->x[1]/h + 16);
    unsigned iz = (p->x[2]/h + 16);
    return zm_encode(ix & HASH_MASK, iy & HASH_MASK, iz & HASH_MASK); // Hashes the last 4 digits
}

// Note: We check ALL buckets, even those that are weird... which hashing should take care of
void particle_neighborhood(unsigned* buckets, particle_t* p, float h)
{
//  unsigned ix = p->x[0]/h;
//  unsigned iy = p->x[1]/h;
//  unsigned iz = p->x[2]/h;
  unsigned ix = (p->x[0]/h + 16);
  unsigned iy = (p->x[1]/h + 16);
  unsigned iz = (p->x[2]/h + 16);
  unsigned x,y,z;

  int counter = 0;
  for (int i = -1; i < 2; i++) {
    for (int j = -1; j < 2; j++) {
      for (int k = -1; k < 2; k++) {
        x = ix + i;
        y = iy + j;
        z = iz + k;

        buckets[counter] = zm_encode(x & HASH_MASK,y & HASH_MASK,z & HASH_MASK);
        counter += 1;
      }
    }
  }
}

void hash_particles(sim_state_t* s, float h)
{

  // Unpack particles and hash
  particle_t* p = s->part;
  particle_t** hash = s->hash;
  int n = s->n;

  // First clear hashtable (TODO: Make this faster)
  for (int i = 0; i < HASH_SIZE; i++)
    hash[i] = NULL;

  // Loop through particles to hash
  for (int i = 0; i < n; i++) {
    // Hash using Z Morton
    int b = particle_bucket(&p[i], h);

    // Add particle to the start of the list of bin b
    p[i].next = hash[b];
    p[i].hind = b;
    hash[b] = &p[i];
  }

}

void hash_particles_proc(particle_t* p, particle_t** hash, int n, float h) {

  // First clear hashtable (TODO: Make this faster)
  for (int i = 0; i < HASH_SIZE; i++) {
    hash[i] = NULL;
  }

  // Loop through particles to hash
  for (int i = 0; i < n; i++) {
    // Hash using Z Morton
    int b = particle_bucket(&p[i], h);

    // Add particle to the start of the list of bin b
    p[i].next = hash[b];
    p[i].hind = b;
    hash[b] = &p[i];
  }

}
