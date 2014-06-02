//TODO: Antalet threads som parameter
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"
#include <set>
#include <omp.h>
//#define DEBUG
using std::vector;
using std::set;
typedef std::pair<int,int> pii;
vector<vector<set<unsigned int> > > M;
particle_t *particles;
static int size; 
static double scale;

//TODO: fixa skalning i pos2grid, återställ density i common.cpp, använd inte set_size för att ta reda på grid size.
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
//  benchmarking program
//
int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 100 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particles = (particle_t*) malloc( n * sizeof(particle_t) );
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

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //
#pragma omp parallel for default(none) shared(particles, size, n)
        for (int i = 0; i < n; ++i)
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


        /* OLD :D
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
                apply_force( particles[i], particles[j] );
        }
        */
        //
        //  move particles and update grid
        //
#pragma omp parallel for default(none) shared(particles, n, M)
        for( int i = 0; i < n; i++ ) {
            
            pii oldPos = pos2grid(particles[i]);
            
            move( particles[i] );
            
            pii newPos = pos2grid(particles[i]);

            // move patricle between cells
            if (oldPos != newPos){
#pragma omp critical
                {
                M[oldPos.first][oldPos.second].erase(i);
                M[newPos.first][newPos.second].insert(i);
                }
            }

        }


        
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
