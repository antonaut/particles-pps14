#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include "common.h"
#include <vector>
#include <set>

using std::vector;
using std::set;
typedef std::pair<int,int> pii;

//
//  global variables
//
int n, n_threads;
particle_t *particles;
FILE *fsave;
pthread_barrier_t barrier;

vector<vector<set<unsigned int> > > M;

static int size; 
static double scale;

pthread_mutex_t lock; // Used to lock access to M.


pii pos2grid(particle_t p){
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

void nf(pii pos, int i){
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



//
//  check that pthreads routine call was successful
//
#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}

//
//  This is where the action happens
//
void *thread_routine( void *pthread_id )
{
    int thread_id = *(int*)pthread_id;

    int particles_per_thread = (n + n_threads - 1) / n_threads;
    int first = min(  thread_id    * particles_per_thread, n );
    int last  = min( (thread_id+1) * particles_per_thread, n );
    
    //
    //  simulate a number of time steps
    //
    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //
        for( int i = first; i < last; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            particle_t p = particles[i];
            pii pos = pos2grid(p);
            
            nf(pos, i);

            if (pos.second != 0) {
                nf(up(pos),i);
                if (pos.first != 0)
                    nf(upleft(pos),i);
                if (pos.first != size-1)
                    nf(upright(pos), i);
            }

            if (pos.second < size-1){
                nf(down(pos),i);
                if (pos.first != 0)
                    nf(downleft(pos),i);
                if (pos.first != size-1)
                    nf(downright(pos), i);
            }

            if (pos.first != 0)
                nf(left(pos),i);
            if (pos.first != size-1)
                nf(right(pos), i);

        }
        pthread_barrier_wait( &barrier );
        
        //
        //  move particles & update grid
        //
        for( int i = first; i < last; i++ ) {
                    pii oldPos = pos2grid(particles[i]);
            
            move( particles[i] );
            
            pii newPos = pos2grid(particles[i]);

            // move patricle between cells
            if (oldPos != newPos){
/* CRITICAL SECTION ENTER */
                pthread_mutex_lock( &lock );
                M[oldPos.first][oldPos.second].erase(i);
                M[newPos.first][newPos.second].insert(i);
                pthread_mutex_unlock( &lock );
/* CRITICAL SECTION LEAVE */                

            }

        }
        pthread_barrier_wait( &barrier );
        
        //
        //  save if necessary
        //
        if( thread_id == 0 && fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    
    return NULL;
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //
    //  process command line
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-p <int> to set the number of threads\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    n = read_int( argc, argv, "-n", 1000 );
    n_threads = read_int( argc, argv, "-p", 2 );
    char *savename = read_string( argc, argv, "-o", NULL );
    
    //
    //  allocate resources
    //
    fsave = savename ? fopen( savename, "w" ) : NULL;

    particles = (particle_t*) malloc( n * sizeof(particle_t) );

    //set_size( n );
    size = ceil(sqrt( 20 * n ));
    scale = double(size)/set_size(n);
    printf("integer size: %i\n",size);
    printf("Scale: %lf\n",scale);
    init_particles( n, particles );

    M.assign(size, vector<set<unsigned int> >(size));

    for (int i = 0; i < n; ++i ){
        pii p = pos2grid(particles[i]);
        M[p.first][p.second].insert(i);
    }


    pthread_mutex_init( &lock, NULL );

    pthread_attr_t attr;
    P( pthread_attr_init( &attr ) );
    P( pthread_barrier_init( &barrier, NULL, n_threads ) );

    int *thread_ids = (int *) malloc( n_threads * sizeof( int ) );
    for( int i = 0; i < n_threads; i++ ) 
        thread_ids[i] = i;

    pthread_t *threads = (pthread_t *) malloc( n_threads * sizeof( pthread_t ) );
    
    //
    //  do the parallel work
    //
    double simulation_time = read_timer( );
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_create( &threads[i], &attr, thread_routine, &thread_ids[i] ) );
    
    thread_routine( &thread_ids[0] );
    
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_join( threads[i], NULL ) );
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, n_threads = %d, simulation time = %g seconds\n", n, n_threads, simulation_time );
    
    //
    //  release resources
    //
    P( pthread_barrier_destroy( &barrier ) );
    P( pthread_attr_destroy( &attr ) );
    P( pthread_mutex_destroy( &lock ) );
    free( thread_ids );
    free( threads );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}