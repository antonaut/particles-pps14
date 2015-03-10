/*

TODO
====

http://www.eecs.berkeley.edu/~carazvan/2013bootcamp/help/Theory.pdf

1 - Låt områded bestå av celler för att linearisera problemet
2 - Distribuera området i num_proc delar. Varje process har endast celler som ligger i sitt område.
3 - Apply force. Om grann-cellerna inte ägs av processen, be om dom från från de andra processerna via msg-passing
4 - move. Om en partikel åker in i en annan process område, skicka meddelande till den
5 - Tillbaks till 3

*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <set>
#include <algorithm>
#include <mpi/mpi.h>

using std::vector;
using std::set;
typedef std::pair<unsigned int, unsigned int> pii; // first: y second: x

// Grid of cells that contain the indices of particles for this process to process
// access: cells[y][x]
vector<vector<set<unsigned int> > > cells; 
// upper and lower bounds of the y-values of the area that the process is handling particles for
std::pair<unsigned int, unsigned int> cell_bounds; 
static double scale; // scale from grid position to particle position

pii pos2grid(const particle_t& p) {
    return pii(
			static_cast<unsigned int>(p.y*scale),
			static_cast<unsigned int>(p.x*scale)
			);
}
pii globalgrid2localgrid(const pii& p) {
	return pii(
			p.first - cell_bounds.first,
			p.second
			);
}
unsigned int pos2grid(double pos) {
	return static_cast<unsigned int>(pos*scale);
}

MPI_Datatype make_particle_datatype() {
	// six floats and a particle index
	int blockcounts[2] = {6,1};
	MPI_Datatype old_types[2] = {MPI_DOUBLE, MPI_UNSIGNED};
	MPI_Aint offsets[2], extent;
	MPI_Type_extent(MPI_DOUBLE, &extent);
	offsets[0] = 0;
	offsets[1] = 6 * extent;
	
	MPI_Datatype new_type;
	MPI_Type_struct(2, blockcounts, offsets, old_types, &new_type);
	MPI_Type_commit(&new_type);
	return new_type;
}

int main(int argc, char *argv[])
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
    int num_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &num_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL; // Savefile
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) ); // Array containing all particles
    
    // Define a PARTICLE datatype
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

	// Split the area into num_proc subareas. We go the easier route here and split on the vertical
    {
    	unsigned int cells_side = ceil(sqrt( 20 * n )); // square root of the total number of cells, a number that should be much larger than the number of particles
    	cells_side = cells_side - (cells_side%num_proc); // make divisible by num_proc
    	unsigned int per_proc = cells_side / num_proc;
    	cell_bounds.first = per_proc*rank;
		cell_bounds.second = cell_bounds.first + per_proc;
    	cells.assign(per_proc, std::vector<std::set<unsigned int> >(cells_side) );
    	scale = double(cells_side)/set_size(n); // Scale between cell index and particle position
    }

	//
	// Initialize particles in root
	//
    if (0 == rank) { 
    	// Let root print some info
        printf("Number of cells per process: %i * %i\n", cells.size(), cells.size());
        printf("Scale: %lf\n", scale);
        printf("Lower cell bound: %i Upper cell bound: %i\n", cell_bounds.first, cell_bounds.second);

        // Initialize particles
        init_particles( n, particles );
        
        // Distribute particles to their respective owning process
        // sort by y value
        std::sort(particles, particles + n, [](const particle_t& a, const particle_t& b) {
        	return a.y < b.y;
        });
    }
	
	printf("pid: %i Lower cell bound: %i Upper cell bound: %i\n", rank, cell_bounds.first, cell_bounds.second);
	
	// Broadcast the particles
	MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);
	
	// Place the relevant particles in the observed cells
	auto comp = [](const particle_t& p, double val) {
		return pos2grid(p.y) < val;
	};
	particle_t* begin = std::lower_bound(particles, particles+n, cell_bounds.first, comp);
	particle_t* end = std::lower_bound(particles, particles+n, cell_bounds.second, comp);
	fprintf(stderr, "pid: %i n: %u begin: %i end-1: %i begin.y: %lf end-1.y: %lf grid(begin): %u grid(end-1): %u\n", rank, n, begin-particles, (end-1)-particles, begin->y, (end-1)->y, pos2grid(begin->y), pos2grid((end-1)->y));
	for (particle_t* it = begin; it != end; ++it) {
		pii p = globalgrid2localgrid(pos2grid(*it));
		if (p.first >= cells.size()) fprintf(stderr, "Uh Oh!\n");
		cells[p.first][p.second].insert(it-particles);
	}
	
	// TODO: At the start of each iteration, get particles on at the top and bottom bounds from other processes
	// Then process own particles as usual. Then if a particle has left grid, send it to the corresponding process
	MPI_Finalize();
	return 0;
}