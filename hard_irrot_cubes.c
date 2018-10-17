#include "mt19937.h" // Mersenne Twister (dsmft_genrand();)
#include "math_4d.h" // https://github.com/arkanis/single-header-file-c-libs/blob/master/math_3d.h
// #include <direct.h> // mkdir on Windows
#include <sys/types.h> // mkdir on Linux
#include <sys/stat.h> // mkdir on Linux
#include <unistd.h>
// #include <iostream> // C++
#include <stdio.h> // C
#include <math.h> // in Linux, make sure to gcc ... -lm, -lm stands for linking math library.
#include <string.h> // This is for C (strcpy, strcat, etc. ). For C++, use #include <string>
#include <assert.h>
#include <stdbool.h> // C requires this for (bool, true, false) to work

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#define NDIM 3
#define N 1000


/* Initialization variables */
// many have been made not constant, so that one can enter into the command line:
// modsim5.exe initial_packing_fraction Delta DeltaV BetaP initial_configuration output_folder
double packing_fraction = 0.5;
static double Delta = 0.1;  // delta and deltaV are dynamic, i.e. every output_steps steps,
static double DeltaV = 2.0; // they will be nudged a bit to keep the move and volume acceptance in between 0.4 and 0.6.
static double BetaP = 30;
char init_filename[] = "sc.txt";
char output_foldername[] = "datafolder_c_p=%4.1lf";

const int mc_steps = 20000;
const int output_steps = 100;
const double diameter = 1.0;

/* Simulation variables */
static double r[N][NDIM]; // position of center of cube // TODO: make vec3_t
static mat4_t m[N]; // rotation matrix of cube
static double box[NDIM]; // dimensions of box
// static double Normal[3][NDIM]; // the normal vector of an unrotated cube. Normal[0] is the normal in the x-dir, etc.
static vec3_t Normal[3];
static double Edge_Length = 1; // TODO: make write, and read. TODO: think about this, maybe just always = 1?
static int n_particles = 0;
static double ParticleVolume;
static double Phi = M_PI / 2.; // angle of slanted cube


/* Functions */
/// Returns a random number between low and high (including exactly low, but not exactly high).
double ran(double low, double high)
{
    return (high - low) * dsfmt_genrand() + low;
}

/// Scales the system with the scale factor
void scale(double scale_factor){
    for (int n = 0; n < n_particles; ++n) {
        for (int d = 0; d < NDIM; ++d)
            r[n][d] *= scale_factor;
    }
    for (int d = 0; d < NDIM; ++d)
        box[d] *= scale_factor;
}

/// Returns the distance between the nearest image of particle i and j.
/// We need to take into account the periodic boundary conditions,
/// we use the nearest image convention
double nearimgdist(int i, int j){
    double distsq = 0; // distance squared
    for (int dim = 0; dim < NDIM; dim++) {
        double tempdist; // we run over the dimensions and need a temporary variable to check nearest image

        tempdist = r[i][dim] - r[j][dim];
        if (tempdist > 0.5 * box[dim])
            tempdist -= box[dim];
        if (tempdist <= -0.5 * box[dim])
            tempdist += box[dim];

        distsq += pow(tempdist, 2.); // add the contribution of the nearest image to the squared distance
    }
    return sqrt(distsq);
}

/// checks if there is overlap between cubes along axis, between cubes i, j
/// Also the difference vector r2-r1 is given as it has already been calculated
/// and the (primary) direction of the normal is given by k = 0, 1, 2, 3, 4, 5 
/// to cut down on needed vertices. If axis is the nth normal (x=0, y=1, z=2)
/// of cube 1, k = n, if axis is the nth normal of cube 2, k = n + 3.
/// e.g. say k = 0, we only need vertex 0, 4 for cube 1.
bool is_collision_along_axis(vec3_t axis, int i, int j, vec3_t r2_r1, int k)
{
    
}

/// checks if cube i and cube j overlap using the separating axis theorem
bool is_overlap(int i, int j)
{
    vec3_t r1 = vec3(r[i][0], r[i][1], r[i][2]);
    vec3_t r2 = vec3(r[j][0], r[j][1], r[j][2]);
    vec3_t r2_r1 = v3_sub(r2, r1); // read as r2 - r1
    
    // If the cubes are more than their greatsphere apart, they couldn't possibly overlap.
    // Similarly, if they are less than their inscribed sphere apart, they couldn't possible NOT overlap.
    // TODO: make Phi-dependent.
    if(v3_length(r2_r1) > 1.73205080757 * Edge_Length) return false;
    if(v3_length(r2_r1) < Edge_Length) return true;

    // Now we use the separating axis theorem. Check for separation along all normals
    // and crossproducts between edges of the cubes. Only if along all these axes
    // we find no separation, we may conclude there is overlap.
    vec3_t axes[6 + 9]; // normals of r1, normals of r2, cross products between edges
    // TODO: ask/think about if there are really 9 axes due to cross products
    for (int k = 0; k < 3; k++) {
        axes[k] = m4_mul_pos(m[i], Normal[k]);
        axes[k + 3] = m4_mul_pos(m[j], Normal[k]);
    }

    bool is_collision = true;
    for (int k = 0; (k < 6) && is_collision; k++) { // TODO: check for cross products of edges between cubes
        is_collision = is_collision_along_axis(axes[k], i, j, r2_r1, k);
    }

    return is_collision;
}

/// Special distance function for move_particle().
/// It returns the distance between particle i and the particle at position trialr[3].
double nearimgdist_mov(double trialr[NDIM], int i)
{
    double distsq = 0; // distance squared
    for (int dim = 0; dim < NDIM; dim++) {
        double tempdist; // we run over the dimensions and need a temporary variable to check nearest image

        tempdist = trialr[dim] - r[i][dim];
        if (tempdist > 0.5 * box[dim])
            tempdist -= box[dim];
        if (tempdist <= -0.5 * box[dim])
            tempdist += box[dim];

        distsq += pow(tempdist, 2.); // add the contribution of the nearest image to the squared distance
    }
    return sqrt(distsq); //TODO: can check against distsq later, if sqrt is slow.
}

/// This function checks if the move- and volume-acceptances are too high or low,
/// and adjusts delta and deltaV accordingly.
void nudge_deltas(double mov, double vol)
{
    if(mov < 0.4) Delta  *= 0.9; // acceptance too low  --> decrease delta
    // if(mov > 0.6) Delta  *= 1.1; // acceptance too high --> increase delta
    if(vol < 0.4) DeltaV *= 0.9;
    if(vol > 0.6) DeltaV *= 1.1;
}

/// This function attempts to change the volume, returning 1 if succesful and 0 if not.
int change_volume(void)
{
    double dV = ran(-DeltaV, DeltaV);
    double oldvol = 1;
    for (int d = 0; d < NDIM; d++)
        oldvol *= box[d];

    double newvol = oldvol + dV;
    double scale_factor = pow(newvol / oldvol, 1. / NDIM); // the . of 1. is important, otherwise 1 / NDIM == 1 / 3 == 0

    // change the configuration
    scale(scale_factor);

    // TODO: now we need to check for overlaps, and reject if there are any
    bool is_collision = false;

    if (is_collision) { 
        scale(1./scale_factor); // move everything back
        return 0; // unsuccesful change
    }
    
    // Now that there are no collisions, we need to accept the change at a rate like the boltzmann factor.
    // Effectively the above code simulates a hard sphere potential (U = \infty if r < d, U = 0 else)
    double boltzmann = pow(M_E, -BetaP * dV) * pow(newvol / oldvol, n_particles);
    if(ran(0,1) < boltzmann) { // if the boltzmann factor > 1, the energy change is negative, and this move will always be accepted.
        return 1; // succesful change
    } else {
        scale(1./scale_factor); // move everything back
        return 0; // unsuccesful change
    }
}

/// This function reads the initial configuration of the system,
/// it trusts the user to supply data in the correct format.
void read_data(void)
{
    FILE* pFile = fopen(init_filename, "r");
    fscanf(pFile, "%i", &n_particles);

    // Gives a proper warning in the next line in main()
    if (n_particles > N) {
        n_particles = 0;
        return;
    }

    // Put the dimensions of the boundaries in box
    double templeft;
    double tempright;
    for (int i = 0; i < NDIM; i++) {
        fscanf(pFile, "%lf %lf", &templeft, &tempright);
        box[i] = tempright - templeft;
        if(box[i] < 0) {
            printf("Error: box dimension negative\n");
            n_particles = 0; // Gives a proper warning in the next line in main()
            return;
        }
    }

    // TODO: make write, and read edge length
    Edge_Length = 1;
    
    // TODO: make write, and read Phi.
    Phi = M_PI / 2.;

    // TODO: check if this is right
    Normal[0].x = cos(Phi); // normal on x-dir
    Normal[0].z = -1. * sin(Phi);
    Normal[1].y = 1.; // normal on y-dir
    Normal[2].z = 1.; // normal on z-dir

    // TODO: make proper particle_volume reading dependent on Phi (I think just this/cos(Phi))
    ParticleVolume = pow(Edge_Length, 3.);

    // Now load all particle positions into r
    double pos;
    for (int i = 0; i < n_particles; i++) {
        for (int dim = 0; dim < NDIM; dim++) {
            fscanf(pFile, "%lf", &pos); // Now pos contains what r[i][dim] should be
            r[i][dim] = pos;
            // m[i][dim][dim] = 1; // TODO: is this the right place and time and way to do this (initialize the rotation matrices)?
        }
    }

    fclose(pFile);
}

/// This function returns the jth vertex (j = 0, ... , 7) of cube number i:
/*        3----7
         /|   /|
z       1-+--5 |
| y     | 2--+-6
|/      |/   |/
 ----x  0----4 */
vec3_t get_vertex(int i, int j){
    vec3_t result = vec3(r[i][0], r[i][1], r[i][2]); // the center
    // TODO: make Phi != M_PI/2 compatible
    vec3_t offset = vec3(-Edge_Length / 2, -Edge_Length / 2, -Edge_Length / 2);
    if(j & 4) offset.x += Edge_Length; // x+ = 4, 5, 6, 7
    if(j & 2) offset.y += Edge_Length; // y+ = 2, 3, 6, 7
    if(j & 1) offset.z += Edge_Length; // z+ = 1, 3, 5, 7
    result = v3_add(result, offset);
    // TODO: add rotation
    return result;
}

/// This moves a random particle in a cube of volume (2 * delta)^3.
/// Note this gives particles a tendency to move to one of the 8 corners of the cube,
/// however it obeys detailed balance.^{[citation needed]}
int move_particle(void)
{
    // first choose the particle
    int index = (int) ran(0., n_particles); // goes from 0 to n_particles - 1

    // then choose the displacement from a cube with sides 2 * delta and store the attempted move in trialr
    double displacement[NDIM];
    double trialr[NDIM];
    for (int dim = 0; dim < NDIM; dim++) {
        displacement[dim] = ran(-Delta, Delta);

        // periodic boundary conditions happen here, fmod = floating point modulo. Since delta < box[dim],
        // the following expression will always return a positive number.
        trialr[dim] = fmod(r[index][dim] + displacement[dim] + box[dim], box[dim]);
    }

    // TODO: now check if this move is possible
    bool is_collision = false; //is_overlap(trialr, index);

    if (is_collision) {
        return 0; // unsuccesful move
    } else { // move, i.e. store trialmove in r
        for (int dim = 0; dim < NDIM; dim++) {
            r[index][dim] = trialr[dim];
        }
        return 1; // succesful move
    }
}

void write_data(int step) // TODO: how many decimal digits are needed? maybe 6 is too much.
{
    char buffer[128];
    strcpy(buffer, output_foldername);
    
    strcat(buffer, "/coords_step%07d.poly");
    char output_file[128];
    sprintf(output_file, buffer, step); // replace %07d with step and put in output_file.

    FILE* fp = fopen(output_file, "w");
    fprintf(fp, "%d\n", n_particles);
    fprintf(fp, "0.0\t0.0\t0.0\n"); // TODO: place the format in the readme file
    for (int d = 0; d < 9; ++d) { // dimensions of box
        if(d % 4 == 0){
            fprintf(fp, "%lf\t", box[d / 4]);
        } else {
            fprintf(fp, "0.000000\t");
        }
        if(d % 3 == 2) fprintf(fp, "\n");
    }
    for (int n = 0; n < n_particles; ++n) {
        for (int d = 0; d < NDIM; ++d)
            fprintf(fp, "%lf\t", r[n][d]); // the position of the center (TODO: check this: center of edge?)
        fprintf(fp, "%lf\t", Edge_Length);
        for (int d1 = 0; d1 < NDIM; d1++) {
            for (int d2 = 0; d2 < NDIM; d2++) {
                fprintf(fp, "%lf\t", m[n].m[d1][d2]); // TODO: maybe the order is the wrong way around?
            }
        }
        fprintf(fp, "10 %lf\n", Phi); // the visualizer wants color, apparently 10 is fine.
    }
    fclose(fp);
}

void set_packing_fraction(void)
{
    double volume = 1.0;
    for (int d = 0; d < NDIM; ++d)
        volume *= box[d];

    double target_volume = (n_particles * ParticleVolume) / packing_fraction;
    double scale_factor = pow(target_volume / volume, 1. / NDIM); // the . of 1. is important, otherwise 1 / NDIM == 1 / 3 == 0

    scale(scale_factor);
}




int main(int argc, char* argv[])
{
    // this first bit is for parsing the command line arguments
    /* if(argc == 7) {
        try {
            packing_fraction = std::stod(argv[1]);
            assert(packing_fraction < 1 && packing_fraction > 0);
            delta = std::stod(argv[2]);
            deltaV = std::stod(argv[3]);
            assert(delta > 0 && deltaV > 0);
            BetaP = std::stod(argv[4]);
            init_filename = argv[5];
            output_foldername = argv[6];
        } catch (int error) {
            printf("error %i: invalid parameters, \nusage: modsim5.exe PF delta deltaV BetaP initial_config_file output_folder\n", error);
            return 1;
        }
    } else {        
        // if the wrong # of arguments are given, display how to use the program
        printf("usage: modsim5.exe PF delta deltaV BetaP initial_config_file output_folder\n");
        printf("standard: %lf, %lf, %lf, %lf, %s, %s\n", packing_fraction, delta, deltaV, BetaP, init_filename.c_str(), output_foldername.c_str());
        return 1;
    } */

    /* if (NDIM == 3)
        particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if (NDIM == 2)
        particle_volume = M_PI * pow(diameter, 2.0) / 4.0;
    else {
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 2;
    } */

    read_data();

    if (n_particles == 0) {
        printf("Error: box dimensions, or n_particles = 0, or n_particles > %d\n", N);
        return 3;
    }

    sprintf(output_foldername, output_foldername, BetaP); // replace %4.1lf with BetaP
    // mkdir(output_foldername); // make the folder to store all the data in, if it already exists do nothing.
    mkdir(output_foldername, S_IRWXU | S_IRGRP | S_IROTH); // Linux needs me to set rights, this gives rwx to me and nobody else.
    set_packing_fraction();

    dsfmt_seed(1); // given seed for retestability
    
    printf("#Step \t Volume \t Move-acceptnce\t Volume-acceptance \n");

    double avg_vol = 0;
    int move_accepted = 0;
    int vol_accepted = 0;
    for (int step = 0; step < mc_steps; ++step) {
        for (int n = 0; n < n_particles; ++n) {
            if(ran(0, n_particles + 1) < n_particles){ // Have to do this for detailed balance
                move_accepted += move_particle();
            } else {
                vol_accepted += change_volume(); // DONE: ask if the order (move->volume) obeys detailed balance. IT DOESN'T!!!
            }
        }

        avg_vol += box[0] * box[1] * box[2];
        if (step % output_steps == 0) {
            double move_acceptance = (double)move_accepted / (n_particles * output_steps);
            double volume_acceptance = (double)vol_accepted / output_steps;

            printf("%d \t %lf \t %lf \t %lf \t delta = %lf, deltaV = %lf \n",
                step, box[0] * box[1] * box[2],
                move_acceptance,
                volume_acceptance,
                Delta, DeltaV);
            
            // Here is where delta and deltaV might get changed if necessary
            nudge_deltas(move_acceptance, volume_acceptance);
            // And reset for the next loop
            move_accepted = 0;
            vol_accepted = 0;
            write_data(step);
        }
    }

    // At the very end, print the average volume to the console.
    avg_vol /= mc_steps;
    printf("Average volume: %lf\n", avg_vol);

    return 0;
}