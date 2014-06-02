/*

TODO
====

http://www.eecs.berkeley.edu/~carazvan/2013bootcamp/help/Theory.pdf

1 - Distribuera partiklar
2 - För alla egna partiklar: Fråga roten om alla partiklar som är runt partikeln
3 - Apply force
4 - move
5 - Gather
6 - Uppdatera celler ( i rot )
7 - repetera från 2

*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <set>
using std::vector;
using std::set;
typedef std::pair<int,int> pii;

vector<vector<set<unsigned int> > > M;

static int size; 
static double scale;


pii pos2grid(particle_t p) {
    return pii((int) (p.x*scale), (int) (p.y*scale));
}
pii up(pii p){
    return pii(p.first, p.second-1);
}
pii down(pii p){
    return pii(p.first, p.second+1);
}
pii left(pii p){
    return pii(p.first-1, p.second);
}
pii right(pii p){
    return pii(p.first+1, p.second);
}

pii upleft(pii p) {
    return up(left(p));
}

pii upright(pii p) {
    return up(right(p));
}

pii downleft(pii p) {
    return down(left(p));
}

pii downright(pii p) {
    return down(right(p));
}

void nf(particle_t *particles, pii pos, int i) {
#ifdef DEBUG        
        printf("set of particles around %i:\n",i);
#endif
    for (std::set<unsigned int>::iterator it = M[pos.first][pos.second].begin(); it != M[pos.first][pos.second].end(); ++it)
    {
#ifdef DEBUG        
        printf("%u\n",*it);
#endif

        if (i != (int) *it)
            apply_force(particles[i], particles[*it]);
    }
}


/**
 * main - where the fun begins!
 */
int main( int argc, char **argv )
{    
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    


    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    size = ceil(sqrt( 20 * n ));
    scale = double(size)/set_size(n);

    if (0 == rank) { // Setup this for root
        printf("integer size: %i\n", size);
        printf("Scale: %lf\n", scale);
        M.assign(size, vector<set<unsigned int> >(size));

        for (int i = 0; i < n; ++i ){
            pii p = pos2grid(particles[i]);
            M[p.first][p.second].insert(i);
        }
    }

    
    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
    
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    
    //
    //  allocate storage for local partition and neighbour buffer 
    //  (assumed that nlocal >> number of neighbours in our grid so the buffer doesn't overflow)
    //
    int nlocal = partition_sizes[rank];
    particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );
    particle_t *neighbours = (particle_t*) malloc( nlocal * sizeof(particle_t));
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );
    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {


        /*
        if (rank == 0) {
            for (int i = 0; i < size;++i) {
                particle_t tmp;
                MPI_Status status;
                MPI_Recv(&tmp, 1, PARTICLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            }
        } else {
            for(int i=0;i<nlocal;++i) {
                particle_t p = local[i];
                int n_neighbours;
                MPI_Send(&p, 1, PARTICLE, 0, 0, MPI_COMM_WORLD);
                MPI_Recv(&neighbours, )

            }
        }

        */

        
        MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
        
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        

        for( int i = 0; i < nlocal; i++ )
        {
            local[i].ax = local[i].ay = 0;
            for (int j = 0; j < nlocal; j++ )
                if (j != i)
                    apply_force( local[i], local[j] );
        }
        
        //
        //  move particles
        //
        for( int i = 0; i < nlocal; i++ )
            move( local[i] );

        
    }
    simulation_time = read_timer( ) - simulation_time;
    
    if( rank == 0 )
        printf( "n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time );
    
    //
    //  release resources
    //
    free( partition_offsets );
    free( partition_sizes );
    free( local );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
