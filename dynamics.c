

// pour malloc et pour les fichiers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



double** malloc_2d(size_t m, size_t n);
void free_2d(double** a);



/*
 * Allocation funcionts
 *
 * Atom coordinates will be provided in nanometers.
 * You will need to allocate two-dimensional arrays. 
 * The following functions should help you allocate and free 2-dimensional arrays:
 */
double** malloc_2d(size_t m, size_t n)
{
    double** a = malloc(m*sizeof(double*)); // memory allocation en 1 dimention
    if (a == NULL)
    {
        return NULL;
    }
    
    a[0] = malloc(n*m*sizeof(double));
    
    if (a[0] == NULL)
    {
        free(a);
        return NULL;
    }
    
    for (size_t i=1 ; i<m ; i++)
    {
        a[i] = a[i-1]+n;
    }

    return a;
}


void free_2d(double** a)
{
    free(a[0]);
    a[0] = NULL;
    free(a);
}



// 1 Describing the Atoms
// We will create functions to read the atomic data from an input file. 
// The input file follows this format:
// The first line contains the number of atoms (Natoms).
// Each subsequent line contains the x, y, and z coordinates followed by the mass of an atom.

size_t read_Natoms(FILE* input_file);
void read_molecule( FILE* input_file,
                    size_t Natoms,
                    double** coord,
                    double* mass);

void compute_distances( size_t Natoms,
                        double** coord,
                        double** distance);
                    
double V(double epsilon,
         double sigma,
         size_t Natoms,
         double** distance);

double T(size_t Natoms,
         double** velocity,
         double* mass);

void compute_acc(size_t Natoms,
                 double** coord,
                 double* mass,
                 double** distance,
                 double** acceleration);


 
// Create a function with the following prototype which reads the number of atoms from an opened input file:
// using C
size_t read_Natoms(FILE* input_file)
{
    size_t Natoms = 0;
    // The first line contains the number of atoms (Natoms).
 // size_t n = read(fd , &number_of_atoms , 1 * sizeof(int));
  //  size_t n = fread(&Natoms, sizeof(int), 1, input_file);
    fscanf(input_file, "%zu", &Natoms);
    return Natoms;
}


// Then, write a function with the following prototype:
// using C



void read_molecule( FILE* input_file, //!< @in read stream
                    size_t Natoms,    //!< @in number of atoms
                    double** coord,   //!< @var coordinates
                    double* mass)     //!< @var mass
{
    // The first line contains the number of atoms (Natoms).
    //  *  Each subsequent line contains the x, y, and z coordinates followed by the mass of an atom.
    for (size_t i = 0; i < Natoms; i++)
    {
        //double x = 0.0L;
        //double y = 0.0L;
        //double z = 0.0L;
        
        //fread(&x, sizeof(double), 1, input_file);
        //fread(&y, sizeof(double), 1, input_file);
        //fread(&z, sizeof(double), 1, input_file);
        
        //coord[i][0] = x;
        //coord[i][1] = y;
        //coord[i][2] = z;
        
        // fread(&(mass[i]), sizeof(double), 1, input_file);
    fscanf(input_file, "%lf %lf %lf %lf", &coord[i][0], &coord[i][1], &coord[i][2], &mass[i]);
    }
}


// using C
// Write another function which takes as input the array of coordinates
// and returns internuclear distances between each pair in a two-dimensional
// double array of size (Natoms × Natoms)
void compute_distances(size_t Natoms,
            double** coord,
            double** distance)
{

    for (size_t i = 0; i < Natoms; i++)
        for (size_t j = 0; j < Natoms; j++)
        {
            double dx = coord[i][0] - coord[j][0];
                    double dy = coord[i][1] - coord[j][1];
                    double dz = coord[i][2] - coord[j][2];
            distance[i][j] = sqrt(dx * dx + dy * dy + dz * dz); 

        }
}



// 2 The Lennard-Jones potential
// 
// Write a function which takes as input ϵ, σ, the number of atoms and the array of distances, 
// and computes the total potential energy
// where VLJ(r) is the Lennard-Jones potential VLJ(r)
// 
// We will now write a function to calculate the total potential energy of the system using the Lennard-Jones potential. 
// This function will take the following inputs: ϵ, σ, the number of atoms (Natoms), 
// and the array containing the distances between each pair of atoms.

// Throughout this tutorial, we will use Argon atoms (mass = 39.948 g/mol) with the following parameters:

// ϵ : 0.0661 j/mol
static const double EPSILON = 0.0661; // j/mol

// σ : 0.3345 nm
static const double SIGMA = 0.3345; // nm

double VLJ(double epsilon,
           double sigma,
           double r) 
{
    return 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));
}

// using C
double V(double epsilon,
         double sigma,
         size_t Natoms,
         double** distance)
{
    double LJ = 0;

    for (size_t i = 0; i < Natoms; i++)
        for (size_t j = i+1; j < Natoms; j++)
        {
            LJ += VLJ(epsilon, sigma, distance[i][j]);                
        }

    return LJ;
}


// 3 Computing the total energy
// 
// We will now write functions to compute the total kinetic energy and the total energy of the system.
// First, we will write a function to calculate the total kinetic energy of the system T. 
// This function will take as input the number of atoms, the array of velocities (3×Natoms) and masses (Natoms). 
// The velocities are initialized to zero. 

void init_to_zero_velocity(size_t Natoms,
                           double** velocity)
{
    for (size_t i = 0; i < Natoms; i++)
    {
        velocity[i][0] = 0.0;
        velocity[i][1] = 0.0;
        velocity[i][2] = 0.0;
    }
}

// The total kinetic energy T is given by T

// using C
double T(size_t Natoms,
         double** velocity,
         double* mass)
{
    double K = 0;

    for (size_t i = 0; i < Natoms; i++)
    {
        // velocity norm is sqrt(vix² + viy² + viz²)
    // velocity norm squared is vix² + viy² + viz²
        double velocity_squared_i = velocity[i][0] * velocity[i][0] + velocity[i][1] * velocity[i][1] + velocity[i][2] * velocity[i][2];
    K += mass[i]*velocity_squared_i;
    }

    return (1/2)*K;
}

// Next, we will write a function to compute the total energy of the system,
// which is the sum of the total kinetic energy T and the total potential energy V. 
// The total energy E is given by: E = T + V.

double E(double epsilon,
         double sigma,
         size_t Natoms,
         double** distance,
         double** velocity,
         double* mass)
{
    return T(Natoms, velocity, mass) +  V(epsilon, sigma, Natoms, distance);
}


// 4 Computing the acceleration
//
// The acceleration vector for each atom is given by:
// The analytical expression of the acceleration reads as:
//
// We need to write a function that computes the acceleration vector for each atom and stores it in a double precision array.

// using C
double U(double epsilon,
         double sigma,
         double r)
{
    return 24*(epsilon/r)*(pow(sigma/r,6)-2*pow(sigma/r,12));
}


void compute_acc(size_t Natoms,
                 double** coord,
                 double* mass,
                 double** distance,
                 double** acceleration)
{
        for (size_t i = 0; i < Natoms; i++)
            for (size_t k = 0; k < 3; k++)
	    {
	        double sum = 0;
        	for (size_t j = 1; j < Natoms; j++)
		{
		     sum += U(EPSILON, SIGMA, distance[i][j]) * (coord[i][k] - coord[j][k]) / distance[i][j];
		}
		acceleration[i][k] = -(1/mass[i]) * sum;
            }
}




int main()
{
    const char * fichier_molecules = "./inp.txt";
    
    FILE * input_file = fopen(fichier_molecules, "r");
    
    if (input_file == NULL)
    {
        return -1;
    }
    else
    {
        size_t Natoms = read_Natoms(input_file);

        fprintf(stdout, "Natoms = %zu\n", Natoms);

            
        double** coord = malloc_2d(Natoms, 3);
        double*  mass  = malloc(Natoms * sizeof(double));
        double** distance = malloc_2d(Natoms, Natoms);
        double** velocity = malloc_2d(Natoms, 3);
            
        // Read input file
        read_molecule(input_file, Natoms, coord, mass);

		// Display coor and mass
        for (size_t i = 0; i < Natoms; i++)
        {
            fprintf(stdout, "atom[%zu] = coord(%lf %lf %lf) mass=%lf\n", i, coord[i][0], coord[i][1], coord[i][2], mass[i]);
        }    
        fprintf(stdout, "\n");

		// Fill distance
        compute_distances(Natoms, coord, distance);
		
		// Display distance
        for (size_t i = 0; i < Natoms; i++)
	    for (size_t j = 0; j < Natoms; j++)
        {
            fprintf(stdout, "distance[%zu][%zu]=%lf\n", i, j, distance[i][j]);
        }    
        fprintf(stdout, "\n");

		// Init all velocities to zero
        init_to_zero_velocity(Natoms, velocity);

        double Vr = V(EPSILON, SIGMA, Natoms, distance);
        fprintf(stdout, "Vr=%lf\n", Vr);

		fprintf(stdout, "T=%lf\n", T(Natoms, velocity, mass));

        double Total_Energy = E(EPSILON, SIGMA, Natoms, distance, velocity, mass);
		fprintf(stdout, "E=%lf\n", Total_Energy);

        free(mass);
        free_2d(coord);
        fclose(input_file);
    }        
    
    
    return 0;
}





