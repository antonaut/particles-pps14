/*

TODO
====

http://www.eecs.berkeley.edu/~carazvan/2013bootcamp/help/Theory.pdf

1 - Låt områded bestå av celler för att linearisera problemet
2 - Distribuera området i num_proc delar. Varje process har endast celler som ligger i sitt område.
3 - Apply force. Om grann-cellerna inte ägs av processen, be om dom från från de andra processerna via msg-passing
4 - move. Om en partikel åker in i en annan process område, skicka meddelande till den
5 - Tillbaks till 3

 *   Hela Området
 *    -----------------
 *    -  1            - Cell-Grid för process 1
 *    -               -
 *    -----------------
 *    -  2            - Cell-Grid för process 2
 *    -               -
 *    -----------------
 *    -  3            - Cell-Grid för process 3
 *    -               -
 *    -----------------
*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
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
static unsigned int num_cells_x, num_cells_y;
static double scale; // scale from grid position to particle position

enum MSG_TAG {
	LENDING_PARTICLES,
	GIVING_PARTICLES,
	COLLECT,
};

pii pos2grid(const particle_t& p) {
    return pii(
			static_cast<unsigned int>(p.y*scale),
			static_cast<unsigned int>(p.x*scale)
			);
}
pii globalgrid2localgrid(const pii& p) {
//	assert(p.first >= cell_bounds.first);
	return pii(
			p.first - cell_bounds.first,
			p.second
			);
}
unsigned int pos2grid(double pos) {
	return static_cast<unsigned int>(pos*scale);
}

pii up(pii p){
    return pii(p.first-1, p.second);
}
pii down(pii p){
    return pii(p.first+1, p.second);
}
pii left(pii p){
    return pii(p.first, p.second-1);
}
pii right(pii p){
    return pii(p.first, p.second+1);
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
template<typename Set>
void nf(Set& particles, particle_t& p1){
	for (auto& p2 : particles)
		apply_force(p1, p2);
}

void nf(pii pos, unsigned int p, particle_t* particles){
	for (const auto& p2 : cells[pos.first][pos.second]) {
		if (p != p2)
			apply_force(particles[p], particles[p2]);
	}
}

template<typename Type>
void RecvVarSize(std::vector<Type>& buffer, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status)
{
	MPI_Probe(source, tag, comm, status);
	int count;
	MPI_Get_count(status, MPI_BYTE, &count);
	buffer.resize(count / sizeof(Type));
	MPI_Recv(buffer.data(), buffer.size(), datatype, source, tag, comm, status);
}

#define sanity_check() for (const auto& row : cells) for (const auto& set : row) for (const auto& p : set) assert(pos2grid(particles[p]).first >= cell_bounds.first && pos2grid(particles[p]).first < cell_bounds.second);

struct Particle_w_index {
	unsigned int i;
	particle_t p;
};

MPI_Datatype make_particle_with_index_datatype() {
	// six floats and a particle index
//	int blockcounts[2] = {1,6};
//	
//	MPI_Datatype old_types[2] = {MPI_UNSIGNED, MPI_DOUBLE};
//	MPI_Aint offsets[2];
//
//	offsets[0] = offsetof(Particle_w_index, Particle_w_index::i);
//	offsets[1] = offsetof(Particle_w_index, Particle_w_index::p);
//	
//	fprintf(stderr, "offsets: %li %li\n", offsets[0], offsets[1], );
//	
//	MPI_Datatype new_type;
//	MPI_Type_create_struct(2, blockcounts, offsets, old_types, &new_type);
//	MPI_Type_commit(&new_type);
	MPI_Datatype new_type;
    MPI_Type_contiguous( sizeof(Particle_w_index), MPI_BYTE, &new_type );
    MPI_Type_commit( &new_type );
	int size;
	MPI_Type_size(new_type, &size);
	assert(size == sizeof(Particle_w_index));
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
	// Define a PARTICLE_WITH_INDEX datatype
	MPI_Datatype PARTICLE_WITH_INDEX = make_particle_with_index_datatype();

	// Split the area into num_proc subareas. We go the easier route here and split on the vertical
    {
    	unsigned int cells_side = ceil(sqrt( 20 * n )); // square root of the total number of cells, a number that should be much larger than the number of particles
    	cells_side = cells_side - (cells_side%num_proc); // make divisible by num_proc
    	unsigned int per_proc = cells_side / num_proc;
    	cell_bounds.first = per_proc*rank;
		cell_bounds.second = cell_bounds.first + per_proc;
    	cells.assign(per_proc, std::vector<std::set<unsigned int> >(cells_side) );
    	scale = double(cells_side)/set_size(n); // Scale between cell index and particle position
		num_cells_x = cells_side;
		num_cells_y = per_proc;
    }

	//
	// Initialize particles in root
	//
    if (0 == rank) { 
    	// Let root print some info
        printf("Number of cells per process: %lu * %lu\n", num_cells_x, num_cells_y);
        printf("Scale: %lf\n", scale);
        printf("Lower cell bound: %u Upper cell bound: %u\n", cell_bounds.first, cell_bounds.second);

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
	for (particle_t* it = begin; it != end; ++it) {
		pii p = globalgrid2localgrid(pos2grid(*it));
		cells[p.first][p.second].insert(it-particles);
	}
	sanity_check();
	//
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ ) 
//	for( int step = 0; step < 1; step++ ) 
	{
		
		sanity_check();
		//
		// Get data from neighbors
		//

		// asynchronously send particles at edge to neighbors, then receive the neighbors data
		std::vector<std::vector<particle_t>> top(num_cells_x), bottom(num_cells_x);
		{
		std::vector<particle_t> to_send_up, to_send_down;
		MPI_Request top_req, bot_req;
		if (rank > 0) {
			for (const auto& set : cells.front())
				for (const auto& p : set)
					to_send_up.push_back(particles[p]);
			MPI_Isend(to_send_up.data(), to_send_up.size(), PARTICLE, rank-1, LENDING_PARTICLES, MPI_COMM_WORLD, &top_req);
//			if (to_send_up.size()) fprintf(stderr, "pid %i sending up information about %lu particles\n", rank, to_send_up.size());
		}
		if (rank < num_proc-1) {
			to_send_down.clear();
			for (const auto& set : cells.back())
				for (const auto& p : set)
					to_send_down.push_back(particles[p]);
			MPI_Isend(to_send_down.data(), to_send_down.size(), PARTICLE, rank+1, LENDING_PARTICLES, MPI_COMM_WORLD, &bot_req);
//			if (to_send_down.size()) fprintf(stderr, "pid %i sending down information about %lu particles\n", rank, to_send_down.size());
		}
		
		// receive particles from neighbors
		// We get from at least one neighbor
		const unsigned int num_to_receive = ((rank==0) || (rank==num_proc-1)) ? 1 : 2;
		MPI_Status status;
		int count;
		for (unsigned int i = 0; i < num_to_receive; ++i) {
			MPI_Probe(MPI_ANY_SOURCE, LENDING_PARTICLES, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, PARTICLE, &count);
			std::vector<particle_t> loaned_particles(count);
			MPI_Recv(loaned_particles.data(), loaned_particles.size(), PARTICLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
//			if(loaned_particles.size()) fprintf(stderr, "pid %i receiving information about %lu particles\n", rank, loaned_particles.size());
			if (status.MPI_SOURCE == rank-1)
				for (const auto& p : loaned_particles)
					top[pos2grid(p.x)].push_back(p);
			else if (status.MPI_SOURCE == rank+1)
				for (const auto& p : loaned_particles)
					bottom[pos2grid(p.x)].push_back(p);
		}
		
		// Wait for neighbors to finish
		if (rank > 0)
			MPI_Wait(&top_req, &status);
		if (rank < num_proc-1)
			MPI_Wait(&bot_req, &status);
		}
		 
		sanity_check();
		//
        //  compute forces
        //
		// first for cells not on the edge
		{
			const unsigned int begin_y = rank == 0 ? 0 : 1;
			const unsigned int end_y = rank == num_proc-1 ? cells.size() : cells.size()-1;
//			for (unsigned int i = begin_y; i < end_y; ++i) {
			for (unsigned int i = 0; i < cells.size(); ++i) {
				for (unsigned int j = 0; j < cells[i].size(); ++j) {
					for (unsigned int pi : cells[i][j]) {

						particles[pi].ax = particles[pi].ay = 0;
						pii pos = {i,j};
						
						nf(pos, pi, particles);
						if (i > 0)
							nf(up(pos), pi, particles);
						if (i < end_y-1)
							nf(down(pos), pi, particles);
						if (j > 0) {
							nf(left(pos), pi, particles);
							if (i < end_y-1)
								nf(left(down(pos)), pi, particles);
							if (i > 0)
								nf(left(up(pos)), pi, particles);
						}
						if (j < cells[i].size()-1) {
							nf(right(pos), pi, particles);
							if (i < end_y-1)
								nf(right(down(pos)), pi, particles);
							if (i > 0)
								nf(right(up(pos)), pi, particles);
						}

					}
				}
			}
		}
		// then the top
		if (rank > 0) { // not the first process
			for (unsigned int x = 0; x < cells.front().size(); ++x) {
				for (int pi : cells.front()[x]) {
					particles[pi].ax = particles[pi].ay = 0;

					// FIXME: Do down as well!
					nf(top[x], particles[pi]);	// up
					if (x > 0)
						nf(top[x-1], particles[pi]); // up-left
					if (x < top.size()-1)
						nf(top[x+1], particles[pi]); // up-right
				}
			}
		}
		// then the bottom
		if (rank < num_proc-1) { // not the last process
			for (unsigned int x = 0; x < cells.back().size(); ++x) {
				for (int pi : cells.back()[x]) {
					particles[pi].ax = particles[pi].ay = 0;
					
					// FIXME: Do up as well!
					nf(bottom[x], particles[pi]);	// down
					if (x > 0)
						nf(bottom[x-1], particles[pi]); // down-left
					if (x < bottom.size()-1)
						nf(bottom[x+1], particles[pi]); // down-right
				}
			}
		}
		
//		sanity_check();
		for (const auto& row : cells) 
			for (const auto& set : row) 
				for (const auto& p : set) 
					assert(pos2grid(particles[p]).first >= cell_bounds.first && pos2grid(particles[p]).first < cell_bounds.second);
		//
        //  move particles and update grid
        //
		// keep track of the particles that change cells, as to not mutate the cells during iteration
		std::vector<std::pair<unsigned int, pii>> to_change;
		for (const auto& row : cells) {
			for (const auto& set : row) {
				for (const auto& p : set) {
					pii oldPos = pos2grid(particles[p]);
					move (particles[p]);
					pii newPos = pos2grid(particles[p]);
					if (oldPos != newPos)
						to_change.push_back({p, oldPos});
				}
			}
		}
		std::vector<Particle_w_index> to_send_down, to_send_up; //FIXME: Send to the correct process instead of just 'up' and 'down'
		for (auto p : to_change) {
			pii oldPos = globalgrid2localgrid(p.second);
			pii newPos = pos2grid(particles[p.first]);
			assert(oldPos.second >= 0 && oldPos.second < num_cells_x);
			assert(oldPos.first >= 0 && oldPos.first < num_cells_y);
			
			cells[oldPos.first][oldPos.second].erase(p.first);
			if (newPos.first >= cell_bounds.second) {
				to_send_down.push_back({p.first, particles[p.first]});
			}
			else if (newPos.first < cell_bounds.first) {
				to_send_up.push_back({p.first, particles[p.first]});
			}
			else {
				newPos = globalgrid2localgrid(newPos);
				assert(newPos.first < num_cells_y && newPos.first >= 0);
				cells[newPos.first][newPos.second].insert(p.first);
			}
		}
		
		sanity_check();
		
//		if (to_send_up.size()) fprintf(stderr, "pid %i had %lu particles going North\n", rank, to_send_up.size());
//		if (to_send_down.size()) fprintf(stderr, "pid %i had %lu particles going South\n", rank, to_send_down.size());
		
		// Do synchronous message passing in increasing rank order.
		
		// Send to the bottom
		if (rank < num_proc-1) {
			int ierr = MPI_Send(to_send_down.data(), to_send_down.size(), PARTICLE_WITH_INDEX, rank+1, GIVING_PARTICLES, MPI_COMM_WORLD);
			assert(ierr == MPI_SUCCESS);
//			if (to_send_down.size()) fprintf(stderr, "pid %i sent down %lu particles\n", rank, to_send_down.size());
		}
		// Receive from the top
		if (rank > 0) 
		{
			MPI_Status status;
			std::vector<Particle_w_index> recv_buf;
			RecvVarSize(recv_buf, PARTICLE_WITH_INDEX, rank-1, GIVING_PARTICLES, MPI_COMM_WORLD, &status);
//			if (recv_buf.size()) fprintf(stderr, "pid %i got %lu particles from above\n", rank, recv_buf.size());
			for (const auto& p : recv_buf) {
				pii pos = globalgrid2localgrid(pos2grid(p.p));
				assert(pos2grid(p.p).first < cell_bounds.second && pos2grid(p.p).first >= cell_bounds.first);
				assert(pos.first < num_cells_y);
				// update particle list
				particles[p.i] = p.p;
				// put in grid
				cells[pos.first][pos.second].insert(p.i);
			}
		}
		// Send to top
		if (rank > 0) {
			MPI_Send(to_send_up.data(), to_send_up.size(), PARTICLE_WITH_INDEX, rank-1, GIVING_PARTICLES, MPI_COMM_WORLD);
//			if (to_send_up.size()) fprintf(stderr, "pid %i sent down %lu particles\n", rank, to_send_up.size());
		}
		// Receive from bottom
		if (rank < num_proc-1) 
		{
			MPI_Status status;
			std::vector<Particle_w_index> recv_buf;
			RecvVarSize(recv_buf, PARTICLE_WITH_INDEX, rank+1, GIVING_PARTICLES, MPI_COMM_WORLD, &status);
//			if (recv_buf.size()) fprintf(stderr, "pid %i got %lu particles from below\n", rank, recv_buf.size());
			for (const auto& p : recv_buf) {
				pii pos = globalgrid2localgrid(pos2grid(p.p));
				assert(pos2grid(p.p).first < cell_bounds.second && pos2grid(p.p).first >= cell_bounds.first);
				assert(pos.first < num_cells_y);
				// update particle list
				particles[p.i] = p.p;
				// put in grid
				cells[pos.first][pos.second].insert(p.i);
			}
		}
		
		//
        //  save if necessary
        //
        if( savename && (step%SAVEFREQ) == 0 ) {
			// let 0 collect all particles and print it
			if (rank == 0) {
				MPI_Status status;
				std::vector<Particle_w_index> buf;
				for (unsigned int i = 0; i < num_proc-1; ++i) {
					RecvVarSize(buf, PARTICLE_WITH_INDEX, MPI_ANY_SOURCE, COLLECT, MPI_COMM_WORLD, &status);
					for (const auto& p : buf)
						particles[p.i] = p.p;
				}
				save( fsave, n, particles );
			// send particles to 0
			} else { 
				std::vector<Particle_w_index> to_send;
				for (const auto& row : cells)
				for (const auto& set : row)
				for (unsigned int i : set)
					to_send.push_back({i, particles[i]});
				MPI_Send(to_send.data(), to_send.size(), PARTICLE_WITH_INDEX, 0, COLLECT, MPI_COMM_WORLD);
			}
		}
		
	}
	
	simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
//    free( particles );
//    if( fsave )
//        fclose( fsave );
	
	MPI_Finalize();
	return 0;
}