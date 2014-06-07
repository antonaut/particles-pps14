// TODO
//
// [/] convert particle init and previous linear implementation
//
// [X] save: implement the save function
//
// [ ] algorithm:
// make a goroutine for every square on the grid
// make sure particles can be sent between squares
// compute forces by asking for your neighbour squares particles
package main

import (
	"flag"
	"fmt"
	"math"
	"math/rand"
	"os"
	"sync"
	"time"
)

var (
	nParticles int
	fileName   string
	size       int
	sizef      float64
	scale      float64
	particles  []Particle
	grid       [][]IntSet
	barrier    sync.WaitGroup
	firstSave  bool
	shouldSave bool
	file       *os.File
)

const (
	density  = 0.01
	mass     = 30
	cutoff   = 1
	min_r    = 0.01
	dt       = 0.00005
	NSTEPS   = 1000
	SAVEFREQ = 10
)

type IntSet struct {
	set map[int]struct{}
}

func (set *IntSet) Init() {
	set.set = make(map[int]struct{})
}

func (set *IntSet) Add(i int) bool {
	_, found := set.set[i]
	set.set[i] = *new(struct{})
	return !found //False if it existed already
}

type Particle struct {
	x, y, dx, dy, ax, ay float64
}

func (p *Particle) PositionOnGrid() (int, int) {
	return (int)(p.x * scale), (int)(p.y * scale)
}

func (p *Particle) ApplyForceOf(p2 *Particle) {
	dx := p2.x - p.x
	dy := p2.y - p.y
	r2 := dx*dx + dy*dy

	if r2 > cutoff*cutoff {
		return
	}

	r2 = math.Max(r2, min_r)
	r := math.Sqrt(r2)
	coef := (1 - cutoff/r) / r2 / mass
	p.ax += coef * dx
	p.ay += coef * dy
}

func (p *Particle) Move() {
	p.dx += p.ax * dt
	p.dy += p.ay * dt
	p.x += p.dx * dt
	p.y += p.dy * dt

	for p.x < 0 || p.x > sizef {
		if p.x < 0 {
			p.x = -p.x
		} else {
			p.x = 2*sizef - p.x
		}
		p.dx = -p.dx
	}

	for p.y < 0 || p.y > sizef {
		if p.y < 0 {
			p.y = -p.y
		} else {
			p.y = 2*sizef - p.y
		}
		p.dy = -p.dy
	}
}

func (p *Particle) SaveString() string {
	return fmt.Sprintf("%g %g\n", p.x, p.y)
}

func (p *Particle) String() (s string) {
	return fmt.Sprintf("{x:%g, y:%g, dx:%g, dy:%g, ax:%g, ay: %g", p.x, p.y, p.dx, p.dy, p.ax, p.ay)
}

func InitParticles() {

	sizef = math.Sqrt(density * float64(nParticles))

	rand.Seed(time.Now().UnixNano())

	sx := size
	sy := (nParticles + sx - 1) / sx

	shuffle := make([]int, nParticles)

	for i := 0; i < nParticles; i++ {
		shuffle[i] = i
	}

	for i := 0; i < nParticles; i++ {

		j := rand.Int() % (nParticles - i)
		k := shuffle[j]
		shuffle[j] = shuffle[nParticles-i-1]

		particles[i].x = sizef * float64(1+(k%sx)) / float64(1+sx)
		particles[i].y = sizef * float64(1+(k/sx)) / float64(1+sy)

		particles[i].dx = rand.Float64()*2 - 1
		particles[i].dy = rand.Float64()*2 - 1
	}
}

func InitGrid() {
	grid = make([][]IntSet, size)
	for i := 0; i < size; i++ {
		grid[i] = make([]IntSet, size)
		for j := 0; j < size; j++ {
			grid[i][j].Init()
		}
	}

	for i := 0; i < nParticles; i++ {
		x, y := particles[i].PositionOnGrid()
		grid[x][y].Add(i)
	}
}

func SimulateSquare(id int, particleCounter []chan int, particleStream []chan Particle) {
	fmt.Printf("Particle[%d] : %v\n", id, particles[id])
}

func save() {
	if !shouldSave {
		return
	}

	if firstSave {
		file.WriteString(fmt.Sprintf("%d %g\n", nParticles, sizef))
		firstSave = false
	}

	for i := 0; i < nParticles; i++ {
		_, err := file.WriteString(particles[i].SaveString())
		if err != nil {
			panic(err)
		}

	}

	err := file.Sync()
	if err != nil {
		panic(err)
	}
}

func main() {

	var err error

	// Setup flags
	flag.IntVar(&nParticles, "n", 1005, "<int> to set the number of particles")
	flag.StringVar(&fileName, "o", "", "<filename> to specify the output file name")
	flag.Parse()

	fmt.Printf("n: %d, o:%s\n", nParticles, fileName)

	// Calculate & setup global variables
	size = int(math.Ceil(math.Sqrt(float64(nParticles))))
	scale = float64(size) / math.Sqrt(density*float64(nParticles))
	particles = make([]Particle, nParticles)

	firstSave = true
	shouldSave = fileName != ""

	if shouldSave {
		file, err = os.Create(fileName)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Couldn't open file : %s\n", err)
			return
		}
	}

	// Setup particles and grid
	InitParticles()
	InitGrid()

	// Start timer
	start := time.Now()

	// Simulate particles
	for step := 0; step < NSTEPS; step++ {

		barrier.Add(nParticles)

		for i := 0; i < nParticles; i++ {
			go func(i int) {
				particles[i].ax = -0.1
				particles[i].ay = -0.1
				particles[i].Move()
				barrier.Done()
			}(i)
		}
		barrier.Wait()

		if step%SAVEFREQ == 0 {
			// Save if neccessary
			save()
		}
	}

	diffTime := time.Since(start)
	fmt.Printf("Time of simulation: %fs\n", diffTime.Seconds())
}
