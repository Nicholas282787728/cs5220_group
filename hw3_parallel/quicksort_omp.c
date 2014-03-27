
/* 	quicksort_omp: OpenMP parallel version of quick sort.
 *
 * 	usage: quicksort_omp Nelements [Nthreads]
 *	Nelements = # of integers to sort 
 *	Nthreads = # of OMP threads -- MUST BE A POWER OF 2! (Can be 1)
 *		   If unspecified then it uses all available cores 
 *
 *	Compile: gcc -O2 quicksort_omp.c -o quicksort_omp -lm -fopenmp -std=c99
 *
 *	Generates Nelements random integers and sorts them in ascending order.  
 *	Runs a check at the end to make sure the list is sorted correctly.
 *
 *	Algorithm: The C qsort() function is run on OMP threads to create
 *	a piecewise sorted list, and then those pieces are merged into
 *	a final sorted list.  Note that it requires a temporary array of
 *	(max) size Nelements for the merging.
 *	
 *	Benchmarks (YMMV): On an 8-core Xeon, sort time for 100 million 
 *	elements is 14.0s on 1 thread, 7.7s on 2, 4.6s on 4, and 3.6s on 8.
 *
 *	Romeel Dave' 5.April.2012
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#define VERBOSE 0

void srand48();
double drand48();

static int CmpInt(const void *a, const void *b)
{
	return ( *(int*)a - *(int*)b );
}

/* Merge sorted lists A and B into list A.  A must have dim >= m+n */
void merge(int A[], int B[], int m, int n) 
{
	int i=0, j=0, k=0;
	int size = m+n;
	int *C = (int *)malloc(size*sizeof(int));
	while (i < m && j < n) {
            if (A[i] <= B[j]) C[k] = A[i++];
            else C[k] = B[j++];
            k++;
	}
	if (i < m) for (int p = i; p < m; p++,k++) C[k] = A[p];
	else for (int p = j; p < n; p++,k++) C[k] = B[p];
	for( i=0; i<size; i++ ) A[i] = C[i];
	free(C);
}

/* Merges N sorted sub-sections of array a into final, fully sorted array a */
void arraymerge(int *a, int size, int *index, int N)
{
	int i;
	while ( N>1 ) {
	    for( i=0; i<N; i++ ) index[i]=i*size/N; index[N]=size;
#pragma omp parallel for private(i) 
	    for( i=0; i<N; i+=2 ) {
		if( VERBOSE ) fprintf(stderr,"merging %d and %d, index %d and %d (up to %d)\n",i,i+1,index[i],index[i+1],index[i+2]);
		merge(a+index[i],a+index[i+1],index[i+1]-index[i],index[i+2]-index[i+1]);
		if( VERBOSE ) for(int i=0; i<size; i++) fprintf(stderr,"after: %d %d\n",i,a[i]);
	    }
	    N /= 2;
	}
}

int main(int argc,char **argv)
{
	int i;
	if( argc != 2 && argc != 3 ) {
		fprintf(stderr,"usage: quicksort_omp Nelements [Nthreads]\n");
		return -1;
	}
// set up array to be sorted
	int size = atoi(argv[1]);
	int *a = (int *)malloc(size*sizeof(int));
	srand48(8675309);
	for(i=0; i<size; i++) a[i] = (int) (size*drand48());
// set up threads
	int threads = omp_get_max_threads();
	if( argc == 3 ) threads=atoi(argv[2]);
	omp_set_num_threads(threads);
	int *index = (int *)malloc((threads+1)*sizeof(int));
	for(i=0; i<threads; i++) index[i]=i*size/threads; index[threads]=size;

/* Main parallel sort loop */
	double start = omp_get_wtime();
#pragma omp parallel for private(i)
	for(i=0; i<threads; i++) qsort(a+index[i],index[i+1]-index[i],sizeof(int),CmpInt);
	double middle = omp_get_wtime();
/* Merge sorted array pieces */
	if( threads>1 ) arraymerge(a,size,index,threads);
	double end = omp_get_wtime();
	fprintf(stderr,"sort time = %g s, of which %g s is merge time\n",end-start,end-middle);

/* Check the sort -- output should never show "BAD: ..." */
	for(int i=1; i<size; i++) if( a[i-1]>a[i] ) fprintf(stderr,"BAD: %d %d %d\n",i,a[i-1],a[i]);
	return 0;
}

