#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <omp.h>

#include "vec3.h"
#include "io.h"
#include "params.h"
#include "state.h"
#include "binhash.h"
#include "interact.h"
#include "leapfrog.h"

/*@q
 * ====================================================================
 */

/*@T
 * \section{Initialization}
 *
 * We've hard coded the computational domain to a unit box, but we'd prefer
 * to do something more flexible for the initial distribution of fluid.
 * In particular, we define the initial geometry of the fluid in terms of an
 * {\em indicator function} that is one for points in the domain occupied
 * by fluid and zero elsewhere.  A [[domain_fun_t]] is a pointer to an
 * indicator for a domain, which is a function that takes two floats and
 * returns 0 or 1.  Two examples of indicator functions are a little box
 * of fluid in the corner of the domain and a circular drop.
 *@c*/
typedef int (*domain_fun_t)(float, float, float);

int box_indicator(float x, float y, float z)
{
  return (x < 0.5) && (y < 0.75) && (z < 0.5);
}

int circ_indicator(float x, float y, float z)
{
  float dx = (x-0.5);
  float dy = (y-0.5);
  float dz = (z-0.5);
  float r2 = dx*dx + dy*dy + dz*dz;
  return (r2 < 0.25*0.25*0.25);
}

int points_indicator(float x, float y, float z) {
  return (x >= 0.5 ||  x<=0.55) && (y == 0.5) && (z == 0.5);
}
/*@T
 *
 * The [[place_particles]] routine fills a region (indicated by the
 * [[indicatef]] argument) with fluid particles.  The fluid particles
 * are placed at points inside the domain that lie on a regular mesh
 * with cell sizes of $h/1.3$.  This is close enough to allow the
 * particles to overlap somewhat, but not too much.
 *@c*/
sim_state_t* place_particles(sim_param_t* param, 
    domain_fun_t indicatef)
{
  float h  = param->h;
  float hh = h/1.3;

  // Count mesh points that fall in indicated region.
  int count = 0;
  for (float x = 0; x < 1; x += hh)
    for (float y = 0; y < 1; y += hh)
      for (float z = 0; z < 1; z += hh)
        count += indicatef(x,y,z);

  // Populate the particle data structure
  sim_state_t* s = alloc_state(count);
  int p = 0;
  for (float x = 0; x < 1; x += hh) {
    for (float y = 0; y < 1; y += hh) {
      for (float z = 0; z < 1; z += hh) {
        if (indicatef(x,y,z)) {
          vec3_set(s->part[p].x, x, y, z);
          vec3_set(s->part[p].v, 0, 0, 0);
          s->part[p].ind = p; // Indexing the particles

          ++p;
        }
      }
    }
  }
  return s;
}

/*@T
 *
 * The [[place_particle]] routine determines the initial particle
 * placement, but not the desired mass.  We want the fluid in the
 * initial configuration to exist roughly at the reference density.
 * One way to do this is to take the volume in the indicated body of
 * fluid, multiply by the mass density, and divide by the number of
 * particles; but that requires that we be able to compute the volume
 * of the fluid region.  Alternately, we can simply compute the
 * average mass density assuming each particle has mass one, then use
 * that to compute the particle mass necessary in order to achieve the
 * desired reference density.  We do this with [[normalize_mass]].
 * 
 * @c*/
void normalize_mass(sim_state_t* s, proc_info* pInfo, sim_param_t* param)
{
  if (pInfo->proc == 0) {
    s->mass = 1; // Set mass with only one processor

    //hash_particles_parallel(s, pInfo, param->h);
    hash_particles(s, param->h); // Hashing with only one processor
  }

  float rho0 = param->rho0;
  float rho2s = 0;
  float rhos  = 0;
#pragma omp barrier // Barrier because need hashing information before computing density

  compute_density(s, pInfo, param);
#pragma omp barrier // Want all processor to finish updating their own densities

  //printf("Starting: %d\n", pInfo->proc);
#pragma omp parallel for reduction(+:rhos, rho2s)
  for (int i = 0; i < s->n; ++i) {
    rho2s += (s->part[i].rho)*(s->part[i].rho);
    rhos  += s->part[i].rho;
  }

#pragma omp single // Only one processor to update this
  {
    s->mass *= ( rho0*rhos / rho2s );
  }
}

sim_state_t* init_particles(sim_param_t* param)
{
  sim_state_t* s = place_particles(param, box_indicator);
  return s;
}

/*@T
 * \section{The [[main]] event}
 *
 * The [[main]] routine actually runs the time step loop, writing
 * out files for visualization every few steps.  For debugging
 * convenience, we use [[check_state]] before writing out frames,
 * just so that we don't spend a lot of time on a simulation that
 * has gone berserk.
 *@c*/

void check_state(sim_state_t* s, proc_info* pInfo)
{
  for (int i = pInfo->beg; i < pInfo->end; ++i) {
    float xi = s->part[i].x[0];
    float yi = s->part[i].x[1];
    float zi = s->part[i].x[2];

    //printf("%d: (%f %f %f) at with (%f, %f, %f) %d\n", i, xi, yi, zi, axi, ayi, azi, pInfo->proc);
    assert( xi >= 0 || xi <= 1 );
    assert( yi >= 0 || yi <= 1 );
    assert( zi >= 0 || zi <= 1 );
  }
}


/* Merge sorted lists A and B into list A.  A must have dim >= m+n */
void merge(particle_t* A, particle_t* B, int m, int n) 
{
  int i=0, j=0, k=0;
  int size = m+n;
  particle_t *C = (particle_t *)malloc(size*sizeof(particle_t));
  while (i < m && j < n) {
    if (A[i].hind <= B[j].hind) C[k] = A[i++];
    else C[k] = B[j++];
    k++;
  }
  if (i < m) for (int p = i; p < m; p++,k++) C[k] = A[p];
  else for (int p = j; p < n; p++,k++) C[k] = B[p];
  for( i=0; i<size; i++ ) A[i] = C[i];
  free(C);
}

void arraymerge(particle_t* a, int size, proc_info* pInfo)
{
  int N = pInfo->nproc; 
  int index[N+1];
  int i;

  // One could do this while in serial. The performance turned out to be
  // the same as parallel. 
  // #pragma omp single
  {
    while ( N>1 ) {

      for( i=0; i<N; i++ ) index[i]=i*size/N; index[N]=size;


      // This is an alternative way of writing the parallel for below
      /*if (pInfo->proc <N/2)
        {

        merge(a+index[pInfo->proc*2],a+index[pInfo->proc*2+1],index[pInfo->proc*2+1]
        -index[pInfo->proc*2],index[pInfo->proc*2+2]-index[pInfo->proc*2+1]); 
        }
#pragma omp barrier*/

#pragma omp for private(i)
      for( i=0; i<N; i+=2 ) 
      {
        merge(a+index[i],a+index[i+1],index[i+1]-index[i],index[i+2]-index[i+1]);
      }
      N /= 2;
    }
  }
}

void sort(sim_state_t* globalState, proc_info* pInfo)
{
    // The serial version
    //#pragma omp single // Implied barrier
         //qsort(globalState->part,globalState-> n, sizeof(particle_t), compPart);
            
           // Parallel version
           qsort(globalState->part+pInfo->beg, pInfo->end-pInfo->beg ,sizeof(particle_t),compPart);
           
           // Merging
           //#pragma omp barrier
           //if( pInfo->nproc >1 ) arraymerge(globalState->part, globalState->n, pInfo);
           
           
           /*#pragma omp single
           for (int i = 0; i<globalState->n;i++)
           {
               printf("particle %d bucket %d \n", i,globalState->part[i].hind);
           }*/
          
}

int main(int argc, char** argv)
{
  sim_param_t params;
  if (get_params(argc, argv, &params) != 0)
    exit(-1);

  // Create global
  sim_state_t* globalState = init_particles(&params);
  omp_set_num_threads(4);

#pragma omp parallel shared(globalState, params) 
  {
    int proc = omp_get_thread_num();
    int nproc = omp_get_num_threads();

    FILE* fp    = fopen(params.fname, "w");
    int nframes = params.nframes;
    int npframe = params.npframe;
    float dt    = params.dt;
    int n       = globalState->n;

    // Processor information and holder
    proc_info* pInfo = malloc(sizeof(proc_info)); 
    pInfo->proc = proc;
    pInfo->nproc = nproc;
    pInfo->beg = round((proc/(double)nproc)*n);
    pInfo->end = round(((proc+1)/(double)nproc)*n);
    pInfo->forceAccu = calloc(3*n, sizeof(float)); // Never used this...


    if (proc == 0) {
      printf("Running in parallel with %d processor\n", nproc);
    }

    normalize_mass(globalState, pInfo, &params);

    double t_start = omp_get_wtime();

    if (proc == 0) { // We only write for one processor
      write_header(fp, n, nframes, params.h);
      write_frame_data(fp, n, globalState, NULL);
    }

    if (proc == 0) {
      hash_particles(globalState, params.h);
    }
    //hash_particles_parallel(globalState, pInfo, params.h);

#pragma omp barrier // Need the hashing to be done

    compute_accel(globalState, pInfo, &params);

#pragma omp barrier
    leapfrog_start(globalState, pInfo, dt);
    check_state(globalState, pInfo);
    for (int frame = 1; frame < nframes; ++frame) {

      // We sort according to Z-Morton to ensure locality, need to implement paralle qsort
      if (frame % 5 == 0) {

        // Dividing into chunks of sorting each chunk
        // This alone turned out to better than sorting the entire array
        //qsort(globalState->part+pInfo->beg, pInfo->end-pInfo->beg ,sizeof(particle_t),compPart);
        // Sorting the array consisting of sorted chunks
        // This turned out to actually lower the performance. That's why
        // I commented it.
        // #pragma omp barrier
        //   if( pInfo->nproc >1 ) arraymerge(globalState->part, globalState->n, pInfo);
//#pragma omp barrier*/

        // Serial version
          
        /*#pragma omp single // Implied barrier
          qsort(globalState->part, n, sizeof(particle_t), compPart);*/
          sort(globalState, pInfo);
      }


#pragma omp barrier // Need sort to finish

    for (int i = 0; i < npframe; ++i) {
      if (proc == 0 && npframe % 4 == 0) { // Ammortize hashing cost
        hash_particles(globalState, params.h);        
      }

#pragma omp barrier
      compute_accel(globalState, pInfo, &params);
      leapfrog_step(globalState, pInfo, dt);
      check_state(globalState, pInfo);
#pragma omp barrier
    }

    if (proc == 0) {
      printf("Frame: %d of %d - %2.1f%%\n",frame, nframes, 
          100*(float)frame/nframes);
      write_frame_data(fp, n, globalState, NULL);
    }
  }

  double t_end = omp_get_wtime();

  if (proc == 0) {
    printf("Ran in %g seconds\n", t_end-t_start);
  }

  free(pInfo);
  fclose(fp);
}

free_state(globalState);
}

