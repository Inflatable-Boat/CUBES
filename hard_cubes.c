#include "math_4d.h" // https://github.com/arkanis/single-header-file-c-libs/blob/master/math_3d.h
#include "mt19937.h" // Mersenne Twister (dsmft_genrand();)
// #include <direct.h> // mkdir on Windows
#include <sys/stat.h> // mkdir on Linux
#include <sys/types.h> // mkdir on Linux
#include <unistd.h>
// #include <iostream> // C++
#include <assert.h>
#include <stdbool.h> // C requires this for (bool, true, false) to work
#include <string.h> // This is for C (strcpy, strcat, etc. ). For C++, use #include <string>
// #include <math.h> // in "math_4d.h" // in Linux, make sure to gcc ... -lm, -lm stands for linking math library.
// #include <stdio.h> // in "math_4d.h" // C
#include <time.h> // time(NULL)

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
double packing_fraction = 0.4;
static double Delta = 0.05; // delta, deltaV, deltaR are dynamic, i.e. every output_steps steps,
static double DeltaR = 0.05; // they will be nudged a bit to keep
static double DeltaV = 2.0; // the move and volume acceptance in between 0.4 and 0.6.
static double BetaP = 60;
char init_filename[] = "sc7.txt";
char output_foldername[] = "datafolder_dir2=%4.1lf";

const int mc_steps = 20000;
const int output_steps = 100;
const double diameter = 1.0;

/* Simulation variables */
static vec3_t r[N]; // position of center of cube // TODO: make vec3_t
static mat4_t m[N]; // rotation matrix of cube
static double box[NDIM]; // dimensions of box
static vec3_t Normal[3]; // the normal vector of an unrotated cube. Normal[0] is the normal in the x-dir, etc.
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
void scale(double scale_factor)
{
    for (int n = 0; n < n_particles; ++n) {
        float* pgarbage = &(r[n].x);
        for (int d = 0; d < NDIM; ++d)
            *(pgarbage + d) *= scale_factor;
    }
    for (int d = 0; d < NDIM; ++d)
        box[d] *= scale_factor;
}

/// This function returns the offset of the jth vertex (j = 0, ... , 7)
/// of cube number i:     3----7
///                      /|   /|
///     z               1-+--5 |
///     | y             | 2--+-6
///     |/              |/   |/
///     0----x          0----4
vec3_t get_offset(int i, int j)
{
    // TODO: make Phi != M_PI/2 compatible
    vec3_t offset = vec3(-Edge_Length / 2, -Edge_Length / 2, -Edge_Length / 2);
    if (j & 4)
        offset.x += Edge_Length; // x+ = 4, 5, 6, 7
    if (j & 2)
        offset.y += Edge_Length; // y+ = 2, 3, 6, 7
    if (j & 1)
        offset.z += Edge_Length; // z+ = 1, 3, 5, 7

    offset = m4_mul_dir(m[i], offset);

    return offset;
}

/// Checks if there is overlap between cubes along axis, between cubes i, j.
/// Also the difference vector r2-r1 is given as it has already been calculated
bool is_collision_along_axis(vec3_t axis, int i, int j, vec3_t r2_r1)
{
    // TODO: I think the axis needs to be normalized:
    axis = v3_norm(axis);

    double min1, min2, max1, max2, temp;
    min1 = max1 = v3_dot(axis, get_offset(i, 0));
    min2 = max2 = v3_dot(axis, v3_add(r2_r1, get_offset(j, 0)));
    for (int n = 1; n < 7; n++) {
        temp = v3_dot(axis, get_offset(i, n));
        min1 = fmin(min1, temp);
        max1 = fmax(max1, temp);
        temp = v3_dot(axis, v3_add(r2_r1, get_offset(j, n)));
        min2 = fmin(min2, temp);
        max2 = fmax(max2, temp);
    }

    if (max1 < min2 || max2 < min1) {
        return false; // separation!
    } else {
        return true; // collision
    }
}

/// Checks if cube i and cube j overlap using the separating axis theorem
bool is_overlap(int i, int j)
{
    vec3_t r2_r1 = v3_sub(r[i], r[j]); // read as r2 - r1

    // We need to apply nearest image convention to r2_r1.
    // We use pointers to loop over the x, y, z members of the vec3_t type.
    float* pdist = &(r2_r1.x);
    for (int d = 0; d < NDIM; d++) {
        if (*(pdist + d) > 0.5 * box[d])
            *(pdist + d) -= box[d];
        if (*(pdist + d) < -0.5 * box[d])
            *(pdist + d) += box[d];
    }

    // If the cubes are more than their circumscribed sphere apart, they couldn't possibly overlap.
    // Similarly, if they are less than their inscribed sphere apart, they couldn't possible NOT overlap.
    // TODO: make Phi-dependent.
    if (v3_length(r2_r1) > 1.73205080757 * Edge_Length)
        return false;
    if (v3_length(r2_r1) < Edge_Length)
        return true;

    // Now we use the separating axis theorem. Check for separation along all normals
    // and crossproducts between edges of the cubes. Only if along all these axes
    // we find no separation, we may conclude there is overlap.
    vec3_t axes[6 + 9]; // 6 normals of r1 and r2, 9 cross products between edges
    for (int k = 0; k < 3; k++) {
        axes[k] = m4_mul_dir(m[i], Normal[k]);
        axes[k + 3] = m4_mul_dir(m[j], Normal[k]);
    }

    // Now load the cross products between edges
    vec3_t edges1[3], edges2[3]; // TODO: do this nicely in a for loop or sth.
    edges1[0] = v3_sub(get_offset(i, 0), get_offset(i, 1));
    edges1[1] = v3_sub(get_offset(i, 0), get_offset(i, 2));
    edges1[2] = v3_sub(get_offset(i, 0), get_offset(i, 4));
    edges2[0] = v3_sub(get_offset(j, 0), get_offset(j, 1));
    edges2[1] = v3_sub(get_offset(j, 0), get_offset(j, 2));
    edges2[2] = v3_sub(get_offset(j, 0), get_offset(j, 4));

    for (int k = 0; k < 9; k++) {
        axes[k + 6] = v3_cross(edges1[k / 3], edges2[k % 3]);
    }

    for (int k = 0; k < 15; k++)
        if (!is_collision_along_axis(axes[k], i, j, r2_r1))
            return false;
    // TODO: make smarter i.e. check only 2 points from first cube

    return true;
}

/// This function checks if the move- and volume-acceptances are too high or low,
/// and adjusts delta and deltaV accordingly.
void nudge_deltas(double mov, double vol, double rot)
{
    if (mov < 0.4)
        Delta *= 0.9; // acceptance too low  --> decrease delta
    if (mov > 0.6 && Delta < Edge_Length / 2.)
        Delta *= 1.1; // acceptance too high --> increase delta
    if (vol < 0.4)
        DeltaV *= 0.9;
    if (vol > 0.6)
        DeltaV *= 1.1;
    if (rot < 0.4)
        DeltaR *= 0.9;
    if (rot > 0.6 && DeltaR < M_PI / 4)
        DeltaR *= 1.1;
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

    // now we need to check for overlaps (only if scale_factor < 1), and reject if there are any
    bool is_collision = false;
    if (scale_factor < 1) {
        for (int i = 0; i < n_particles; i++) {
            for (int j = i + 1; j < n_particles; j++) {
                if (is_overlap(i, j)) {
                    is_collision = true;
                    i = j = n_particles; // i.e. break out of both loops
                }
            }
        }
    }

    if (is_collision) {
        scale(1. / scale_factor); // move everything back
        return 0; // unsuccesful change
    }

    // Now that there are no collisions, we need to accept the change at a rate like the boltzmann factor.
    // Effectively the above code simulates a hard sphere potential (U = \infty if r < d, U = 0 else)
    double boltzmann = pow(M_E, -BetaP * dV) * pow(newvol / oldvol, n_particles);
    if (ran(0, 1) < boltzmann) { // if the boltzmann factor > 1, the energy change is negative, and this move will always be accepted.
        return 1; // succesful change
    } else {
        scale(1. / scale_factor); // move everything back
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
        if (box[i] < 0) {
            printf("Error: box dimension negative\n");
            n_particles = 0; // Gives a proper warning in the next line in main()
            return;
        }
    }

    // TODO: make write, and read edge length / decide if it should be always 1
    Edge_Length = 1;

    // TODO: make write, and read Phi.
    Phi = M_PI / 2.;

    // now initialize the normals, put everything to zero first:
    for (int i = 0; i < 3; i++) {
        Normal[i].x = Normal[i].y = Normal[i].z = 0;
    }
    // TODO: check if this is right
    Normal[0].x = sin(Phi); // normal on x-dir
    Normal[0].z = -1. * cos(Phi);
    Normal[1].y = 1.; // normal on y-dir
    Normal[2].z = 1.; // normal on z-dir

    // TODO: make proper particle_volume reading dependent on Phi (I think just this/cos(Phi))
    ParticleVolume = pow(Edge_Length, 3.);

    // Now load all particle positions into r
    double pos;
    for (int i = 0; i < n_particles; i++) {
        for (int temp = 0; temp < 16; temp++)
            m[i].m[temp % 4][temp / 4] = 0; // TODO: is this right place/time/way?
        float* pgarbage = &(r[i].x);
        for (int d = 0; d < NDIM; d++) {
            fscanf(pFile, "%lf", &pos); // Now pos contains what r[i][dim] should be
            *(pgarbage + d) = pos;
            m[i].m[d][d] = 1; // TODO: is this the right place and time and way to do this (initialize the rotation matrices)?
        }
        m[i].m[3][3] = 1; // TODO: Is this necessary?
    }

    fclose(pFile);
}

/// This moves a random particle in a cube of volume (2 * delta)^3.
/// Note this gives particles a tendency to move to one of the 8 corners of the cube,
/// however it obeys detailed balance.^{[citation needed]}
/// returns 1 on succesful move, 0 on unsuccesful move.
int move_particle(void)
{
    // first choose the particle and remember its position
    int index = (int)ran(0., n_particles); // goes from 0 to n_particles - 1
    vec3_t memory = r[index];

    float* pgarbage = &(r[index].x);
    for (int d = 0; d < NDIM; d++) {
        *(pgarbage + d) += ran(-Delta, Delta);
        // periodic boundary conditions happen here, fmod = floating point modulo. Since delta < box[dim],
        // the following expression will always return a positive number.
        *(pgarbage + d) = fmodf(*(pgarbage + d) + box[d], box[d]);
    }

    // TODO: make cell structure in box so we don't need as many checks
    // if is_overlap(index, any other one) == true, it stops the loop,
    // as any collision results in an unsuccesful move.
    bool is_collision = false; // TODO: i < index + n_particles && (!is_collision) necessary?
    for (int i = index + 1; i < index + n_particles; i++) {
        if (is_overlap(index, i % n_particles)) {
            is_collision = true;
            break;
        }
    }

    if (is_collision) {
        r[index] = memory; // move back
        return 0; // unsuccesful move
    } else {
        return 1; // succesful move
    }
}

/// This rotates a random particle around a random axis
/// by a random angle \in [-DeltaR, DeltaR]
/// returns 1 on succesful rotation, 0 on unsuccesful rotation.
int rotate_particle(void)
{
    // first choose a random particle and remember its position
    int index = (int)ran(0., n_particles);

    // then choose a random axis (by picking points randomly in a ball)
    double x[3], dist;
    do {
        for (int i = 0; i < 3; i++)
            x[i] = ran(-1, 1);
        dist = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    } while (dist > 1); // TODO: dist = 0 problematic? can this happen at all?
    if(dist == 0) // TODO: remove this baby
        printf("\nholy shit the random number generator is exactly zero!\n");
    vec3_t rot_axis = vec3(x[0], x[1], x[2]); // the axis doesn't need to be normalized
    mat4_t rot_mx = m4_rotation(ran(-DeltaR, DeltaR), rot_axis);

    // rotate the particle
    m[index] = m4_mul(rot_mx, m[index]);
    // and check for collisions
    // TODO: make cell structure
    for (int i = index + 1; i < index + n_particles; i++) {
        if (is_overlap(index, i % n_particles)) {
            // rotate back
            m[index] = m4_mul(m4_invert_affine(rot_mx), m[index]);
            return 0; // unsuccesful rotation
        }
    }
    // if we arrive here, there is no overlap and the rotation is accepted!
    return 1;
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
    fprintf(fp, "0.0\t0.0\t0.0\n");
    for (int d = 0; d < 9; ++d) { // dimensions of box
        if (d % 4 == 0) {
            fprintf(fp, "%lf\t", box[d / 4]);
        } else {
            fprintf(fp, "0.000000\t");
        }
        if (d % 3 == 2)
            fprintf(fp, "\n");
    }
    for (int n = 0; n < n_particles; ++n) {
        float* pgarbage = &(r[n].x);
        for (int d = 0; d < NDIM; ++d)
            fprintf(fp, "%lf\t", *(pgarbage + d)); // the position of the center (TODO: check this: center of edge?)
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

    dsfmt_seed(1234);

    printf("#Step\tVolume\t acceptances\t\t\t deltas\n");

    double avg_vol = 0;
    int move_accepted = 0, vol_accepted = 0, rot_accepted = 0;
    int move_attempted = 0, vol_attempted = 0, rot_attempted = 0;
    for (int step = 0; step < mc_steps; ++step) {
        for (int n = 0; n < 2 * n_particles + 1; ++n) {
            // Have to randomize order of moves for detailed balance
            if (ran(0, 2 * n_particles + 1) < 2 * n_particles) {
                if (ran(0, 1) < 0.5) {
                    move_attempted++;
                    move_accepted += move_particle();
                } else {
                    rot_attempted++;
                    rot_accepted += rotate_particle();
                }
            } else {
                vol_attempted++;
                vol_accepted += change_volume(); // DONE: ask if the order (move->volume) obeys detailed balance. IT DOESN'T!!!
            }
        }

        avg_vol += box[0] * box[1] * box[2];
        if (step % output_steps == 0) {
            double move_acceptance = (double)move_accepted / move_attempted;
            double rotation_acceptance = (double)rot_accepted / rot_attempted;
            double volume_acceptance = (double)vol_accepted / vol_attempted;
            printf("%d\t%.3lf\t %.3lf\t%.3lf\t%.3lf\t %.3lf\t%.3lf\t%.3lf\n",
                step, box[0] * box[1] * box[2],
                move_acceptance,
                rotation_acceptance,
                volume_acceptance,
                Delta, DeltaR, DeltaV);

            // Here is where delta, deltaR, deltaV might get changed if necessary
            nudge_deltas(move_acceptance, volume_acceptance, rotation_acceptance);
            // And reset for the next loop
            move_attempted = rot_attempted = vol_attempted = 0;
            move_accepted = rot_accepted = vol_accepted = 0;
            write_data(step);
        }
    }

    // At the very end, print the average volume to the console.
    avg_vol /= mc_steps;
    printf("Average volume: %lf\n", avg_vol);

    return 0;
}