#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

void hash_particles(sim_state_t* s,float h)
{

  // Unpack particles and hash
  particle_t* p = s->part;
  particle_t** hash = s->hash;
  int n = s->n;

  //printf("In hashing! Processor %d\n", omp_get_thread_num());

  // First clear hashtable (TODO: Make this faster)
//#pragma omp parallel for
  for (int i = 0; i < HASH_SIZE; i++)
    hash[i] = NULL;

  // Loop through particles to hash
//#pragma omp parallel for
//#pragma omp single // Debug this later
  
  for (int i = 0; i < n; i++) {
    // Hash using Z Morton
    int b = particle_bucket(&p[i], h);

    // Add particle to the start of the list of bin b
    // omp critical here
      p[i].next = hash[b];
      p[i].hind = b;
      hash[b] = &p[i];
  }
  

//printf("Reached end with processor %d\n", pInfo->proc );
}

// Bad idea
/*
void hash_particles_parallel(sim_state_t* s, proc_info* pInfo,float h)
{

  // Unpack particles and hash
  particle_t* p = s->part;
  particle_t** hash = s->hash;
  int n = s->n;

  //printf("In hashing! Processor %d\n", omp_get_thread_num());

  // First clear hashtable (TODO: Make this faster)
//#pragma omp parallel for
  for (int i = 0; i < HASH_SIZE; i++)
    hash[i] = NULL;

  // Loop through particles to hash
//#pragma omp parallel for
//#pragma omp single // Debug this later
 // #pragma omp parallel for 
  for (int i = pInfo->beg; i < pInfo->end; i++) {
    // Hash using Z Morton
    int b = particle_bucket(&p[i], h);

    // Add particle to the start of the list of bin b
     // #pragma omp critical
      p[i].next = hash[b];
      p[i].hind = b;
      hash[b] = &p[i];
  }
  

//printf("Reached end with processor %d\n", pInfo->proc );
}*/

void hash_particles_parallel(sim_state_t* s, proc_info* pInfo,float h, int hist[8][4096], int vec_count[8])
{

  // Unpack particles and hash
  particle_t* p = s->part;
  particle_t** hash = s->hash;
  particle_t*** hash_l = s->hash_l;
  //particle_t** hash_test = NULL;
  /*particle_t*** hash_old;
 
  
  hash_old = malloc(pInfo->nproc*sizeof(particle_t));
  for (int i = 0; i<pInfo->nproc;i++)
  {
      hash_old[i] = malloc(sizeof(hash));
  }
  particle_t** hash_test = NULL;    */
  //particle_t hash_old[pInfo->nproc][HASH_SIZE];

  /*particle_t*** hash_l;
  
  hash_l  = (particle_t***) calloc(pInfo->nproc, sizeof(particle_t**));
  for (int i = 0;i<pInfo->nproc;i++)
  {
  hash_l[i] = (particle_t**) calloc(HASH_SIZE, sizeof(particle_t*));
#pragma omp for
  for (int j = 0;j<HASH_SIZE;j++)
  {
  hash_l[i][j] = NULL;
  }
  }*/
  
//#pragma omp single
  //{
  for (int i = 0;i<pInfo->nproc;i++)
  {
  #pragma omp for    
  for (int j = 0;j<HASH_SIZE;j++)
  {
  hash_l[i][j] = NULL;
  }
  }
 //}
 /* #pragma omp for
  for (int i = 0; i < HASH_SIZE; i++)
  {
    hash[i] = NULL;
  }*/
  
  
  
  int n = s->n;
 /* int hist[pInfo->nproc][HASH_SIZE];


  memset(hist, -1, (pInfo->nproc)*HASH_SIZE*sizeof(hist[0][0]) );*/
  

  //printf("In hashing! Processor %d\n", omp_get_thread_num());

  // First clear hashtable (TODO: Make this faster)
/*#pragma omp for
  for (int i = 0; i < HASH_SIZE; i++)
  {
    hash[i] = NULL;
  }*/

  int b, i;
  int count =0;
  

//int* vec_count = (int*)calloc(pInfo->nproc,sizeof(int));

 /* int vec_count[pInfo->nproc];
  memset(vec_count, 0, (pInfo->nproc)*sizeof(int) );*/
  
  

//#pragma omp for private(i,b) shared()
  

#pragma omp for private(b)
  for (i = 0; i < n; i++) 
  {
    // Hash using Z Morton
    b = particle_bucket(&p[i], h);
    //printf("b is %d and processor %d \n", b, pInfo->proc);
    p[i].hind = b;
    p[i].next = hash_l[pInfo->proc][b];
    //p[i].next = hash[b];
   // printf("address is %d processor %d \n", hash_l[pInfo->proc][b],pInfo->proc);
   if(p[i].next == NULL)
    {
        hist[pInfo->proc][count] = i;
        count +=1;
       // printf("count is %d processor %d \n",count,pInfo->proc);
        vec_count[pInfo->proc] = count;
      //  printf("vec count is %d processor %d \n",vec_count[pInfo->proc],pInfo->proc);
        //hist[pInfo->proc][count] = i;
    }
    hash_l[pInfo->proc][b] = &p[i];
    
   //hash[b] = &p[i];
    }
  
/*#pragma omp single
  {
for (int i = 0; i<pInfo->nproc;i++)
{
    for (int j = 0;j<HASH_SIZE;j++)
    {
        if (hash_l[i][j] != NULL)
     printf("Hash is %d and i %d \n", hash_l[i][j], i);   
    }
}
  }*/
  
  #pragma omp for
  for (int i = 0; i < HASH_SIZE; i++)
  {
    hash[i] = hash_l[0][i];
  }
  
//#pragma omp barrier 
  

//#pragma omp single
  //{
for (int i = 1; i<pInfo->nproc;i++)
{
#pragma omp for
    for (int j = 0;j<vec_count[i];j++)
    {
        b = p[hist[i][j]].hind;
        //printf("next particle is %d and i %d \n", p[hist[i][j]].next, i);
        int local_i = i - 1;
        while (local_i >=0)
        {           
           // printf("local_i is %d \n", local_i);
            if (hash_l[local_i][b] != NULL)
            {
                //printf("local_i is %d and i %d \n", local_i, i);
                p[hist[i][j]].next = hash_l[local_i][b];
                local_i = -1;
            }
            hash[b] = hash_l[i][b];
            local_i -=1;
        }
        
    }
}     
 // }
 
 // exit(0);

  
 /* 
#pragma omp single 
 for (i = 0; i < n; i++) 
 {
   int b = particle_bucket(&p[i], h);
   hash[b] = &p[i];  
   //printf("i is %d address is %d \n", i, hash[b] );
 }*/
 // exit(0);
  
 /* #pragma omp single
for(int i =1;i< pInfo->nproc;i++)
  {
  for (int j = 0; j< HASH_SIZE; j++)
      {
          //hash_test = hash_old[i-1];
          //printf("j is: %lu \n", sizeof(hash_old)/sizeof(particle_t));
          //exit(0);
        if (hist[i][j] != -1)
          {
          b = p[hist[i][j]].hind;
          printf("processor %d hist[i][j] is %d \n", i, j);
          printf("j is %d \n", j);
          printf("b is %d \n", b);
          printf("address is %d \n", hash_l[i-1][b]);
          
             p[hist[i][j]].next = hash_l[i-1][b];
         }
      }
}
exit(0);
*/
  
 /*
#pragma omp barrier
  
#pragma omp single
for (int i=0; i<pInfo->nproc;i++)
{
    
    printf("vec count is %d \n",vec_count[i]);
    
}
   */

/*  
#pragma omp single
  {
  for(int i =0;i< pInfo->nproc;i++)
  {
      for (int j = 0; j< HASH_SIZE; j++)
      {
          hash[j] = hash_l[i][j];
      }
  }
  }
  */
  
    
 
/*
#pragma omp single
  {
  for(int i =1;i< pInfo->nproc;i++)
  {
      for (int j = 0; j< vec_count[i]; j++)
      {
          //hash_test = hash_old[i-1];
          //printf("j is: %lu \n", sizeof(hash_old)/sizeof(particle_t));
          //exit(0);
       // if (hist[i][j] != -1)
         // {
          int b = p[hist[i][j]].hind;
          
          //printf("j is %d \n", j);
         // printf("b is %d \n", b);
          //printf("address is %d \n", hash_l[i-1][b]);
          
         //printf("size p %d size hash %d \n", sizeof(p[hist[i][j]].next), sizeof(hash_l[i-1][b]));
      // printf("processor %d p %d index %d hash %d vec %d \n",i, p[hist[i][j]].next, hist[i][j], hash_l[i-1][b], vec_count[i]);
        // printf("hash %d \n",hash_l[i-1][b]);
          
          if( hash_l[i][b] != NULL )
          {
             p[hist[i][j]].next = hash_l[i-1][b];
          }
           //  hash_l[i-1][b] = &p[hist[i][j]];
             
            // hash_l[i-1][b] = &p[hist[i][j]];
        // }
      }
  
  }
  }
  */

/*    
#pragma omp single
  {
  for(int i =7;i>0;i--)
  {
      for (int j = 0; j< HASH_SIZE; j++)
      {
          if(hash[j] == NULL)
          {
          hash[j] = hash_l[i][j];
          }
      }
  }
  }
*/

    
    
 //exit(0);
  
/*#pragma omp barrier    
#pragma omp single  
  {
        int count_test = 0;
  for (int j = 0; j< vec_count[7]; j++)
  {
      //printf("p.next %d \n", p[hist[1][j]].next);
      if (p[hist[7][j]].next == NULL) 
      {
          count_test +=1;
          printf("vec_count %d count_test %d \n", vec_count[7], count_test);
      }
          
      }
  }
   
 
#pragma omp barrier    
exit(0);*/
  }
  
  
    
  

    
    
    
    
/*  
#pragma omp single
  for (int i = 0; i < n; i++)
  {

    // Add particle to the start of the list of bin b
     // #pragma omp critical
      int b = p[i].hind;
      p[i].next = hash[b];
      hash[b] = &p[i];
  }
*/  

  
  /*#pragma omp single 
  {
  for (int i = 0; i < n; i++) 
  {
    // Hash using Z Morton
    int b = particle_bucket(&p[i], h);
    p[i].hind = b;
    p[i].next = hash[b];
    hash[b] = &p[i];
  }
  
  }*/
  
//printf("Reached end with processor %d\n", pInfo->proc );


    //if(i == pInfo->end -1)
   // {
        //hash_old[pInfo->proc] = hash;
        //printf("size of hash is %d \n",HASH_SIZE*sizeof(particle_t) );
    
       // memcpy(&hash_old[pInfo->proc], hash, sizeof(hash));
   // }