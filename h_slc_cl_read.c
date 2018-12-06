#include "math_5d.h" // https://github.com/arkanis/single-header-file-c-libs/blob/master/math_3d.h
#include "mt19937.h" // Mersenne Twister (dsmft_genrand();)
#ifdef _WIN32
#include <direct.h> // mkdir on Windows
#elif __linux__
#include <sys/stat.h> // mkdir on Linux
#include <sys/types.h> // mkdir on Linux
#endif
// #include <unistd.h> // sleep()
#include <stdbool.h> // C requires this for (bool, true, false) to work
#include <string.h> // This is for C (strcpy, strcat, etc. ). For C++, use #include <string>
// #include <math.h> // in "math_4d.h" // in Linux, make sure to gcc ... -lm, -lm stands for linking math library.
// #include <stdio.h> // in "math_4d.h" // C
// #include <iostream> // C++
#include <time.h> // time(NULL)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#define NDIM 3
#define N 8000 // 20^3 should be plenty for now
#define NC 16 // max number of cells per direction
#define MAXPC 64 // max number of particles per cell
// WARNING!, don't make NC much larger than 32, because
// WhichCubesInCell is an array of size NC^3, this is already a million at NC == 32.

/* Initialization variables */
// many have been made not constant, so that one can enter into the command line:
// a.exe mc_steps packing_fraction BetaP Phi
static int mc_steps;
static double packing_fraction;
static double BetaP;
static double Phi; // angle of slanted cube

const char labelstring[] = "v1_%02dpf%04.2lfp%04.1lfa%04.2lf";
// e.g. sl10pf0.50p08.0a1.25:
// 10 CubesPerDim, pack_frac 0.50, pressure 8.0, angle 1.25
const char usage_string[] = "usage: program.exe (r for read / c for create) \
(readfile / # cubes per dim) mc_steps packing_fraction BetaP Phi\n";

const int output_steps = 100;

/* Simulation variables */
// TODO: use malloc and pointers instead of global variables?
static double Delta = 0.05; // delta, deltaV, deltaR are dynamic, i.e. every output_steps steps,
static double DeltaR = 0.05; // they will be nudged a bit to keep
static double DeltaV = 2.0; // the move and volume acceptance in between 0.4 and 0.6.

static bool IsCreated; // Did we create a system or read a file

static vec3_t r[N]; // position of center of cube
static mat4_t m[N]; // rotation matrix of cube
static double CellLength; // The length of a cell
static int CellsPerDim; // number of cells per dimension
static int NumCubesInCell[NC][NC][NC]; // how many cubes in each cell
static int WhichCubesInCell[NC][NC][NC][MAXPC]; // which cubes in each cell
static int InWhichCellIsThisCube[N]; // (a,b,c) == NC*NC*a + NC*b + c
static double box[NDIM]; // dimensions of box
static vec3_t Normal[3]; // the normal vector of an unrotated cube. Normal[0] is the normal in the x-dir, etc.
static int n_particles = 0;
static int CubesPerDim;
static double ParticleVolume;
static double CosPhi; // cos and sin of Phi appear a lot, and are expensive to calculate.
static double SinPhi; // Since Phi doesn't change, it's faster to calculate only once.

/* Functions */
void print_celllists(void);

inline static double ran(double low, double high);
inline static int pos_mod_i(int a, int b);
inline static double pos_mod_f(double a, double b);
void scale(double scale_factor);

// initialization
int parse_commandline(int argc, char* argv[]);
void read_data2(char* init_file);
void create_system();
void set_packing_fraction(void);
void set_random_orientation(void);
void remove_overlap(void);
void initialize_cell_list(void);

// mc steps
// int move_particle(void);
int move_particle_cell_list(void);
int rotate_particle(void);
int change_volume(void);
void nudge_deltas(double mov, double vol, double rot);
void write_data(int step, FILE* fp_density, char datafolder_name[128]);

// collision detection
bool is_overlap_between(int i, int j);
bool is_collision_along_axis(vec3_t axis, int i, int j, vec3_t r2_r1);
vec3_t get_offset(int i, int j);
bool is_overlap_from(int index);
void update_cell_list(int index);
inline static void update_CellLength(void);
bool is_overlap(void);

/* Main */

int main(int argc, char* argv[])
{
    dsfmt_seed(time(NULL));

    if (parse_commandline(argc, argv)) { // if ... parsing fails
        printf(usage_string);
        return 1;
    };

    if (IsCreated) {
        set_packing_fraction();
        set_random_orientation(); // to stop pushing to rhombic phase
        remove_overlap(); // due to too high a packing fraction or rotation
    } else {
        if (is_overlap()) {
            printf("\n\n\tWARNING\n\n\nThe read file contains overlap.\n\
            Expanding system until no overlap left\n");
            remove_overlap();
        } else {
            printf("No overlap detected, continuing simulation.\n");
        }
    }
    initialize_cell_list();
    // DEBUG: print NumCubesInCell, and other cell list numbers.
    printf("celllength: %lf\nboxsize: %lf\ncellsperdim: %d\n\n", CellLength, box[0], CellsPerDim);
    for (int i = 0; i < CellsPerDim; i++) {
        for (int j = 0; j < CellsPerDim; j++) {
            for (int k = 0; k < CellsPerDim; k++) {
                printf("%2d ", NumCubesInCell[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
    char buffer[128] = "datafolder/";
    strcat(buffer, labelstring);
    char datafolder_name[128] = "";
    // replace all %d, %lf in buffer with values and put in datafolder_name
    sprintf(datafolder_name, buffer, CubesPerDim, packing_fraction, BetaP, Phi);

// make the folder to store all the data in, if it already exists do nothing.
#ifdef _WIN32
    mkdir("datafolder");
    mkdir(datafolder_name);
    mkdir("densities");
#elif __linux__
    // Linux needs me to set rights, this gives rwx to me and just r to all others.
    mkdir("datafolder", S_IRWXU | S_IRGRP | S_IROTH);
    mkdir(datafolder_name, S_IRWXU | S_IRGRP | S_IROTH);
    mkdir("densities", S_IRWXU | S_IRGRP | S_IROTH);
#elif
    printf("Please use Linux or Windows\n");
    exit(3);
#endif

    char buffer2[128] = "densities/";
    strcat(buffer2, labelstring);
    char density_filename[128] = "";
    // replace all %d, %lf in buffer2 with values and put in density_filename
    sprintf(density_filename, buffer2, CubesPerDim, packing_fraction, BetaP, Phi);

    FILE* fp_density = fopen(density_filename, "w");

    int mov_accepted = 0, vol_accepted = 0, rot_accepted = 0;
    int mov_attempted = 0, vol_attempted = 0, rot_attempted = 0;

    printf("#Step\tVolume\t acceptances\t\t\t deltas\n");
    for (int step = 0; step <= mc_steps; ++step) {
        for (int n = 0; n < 2 * n_particles + 1; ++n) {
            // Have to randomize order of moves to obey detailed balance
            int temp_ran = (int)ran(0, 2 * n_particles + 2);
            if (temp_ran < n_particles) {
                mov_attempted++;
                mov_accepted += move_particle_cell_list();
            } else if (temp_ran < 2 * n_particles) {
                rot_attempted++;
                rot_accepted += rotate_particle();
            } else {
                vol_attempted++;
                vol_accepted += change_volume();
            }
        }

        if (step % output_steps == 0) {
            double move_acceptance = (double)mov_accepted / mov_attempted;
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
            mov_attempted = rot_attempted = vol_attempted = 0;
            mov_accepted = rot_accepted = vol_accepted = 0;
            write_data(step, fp_density, datafolder_name);
            if (step % 10000 == 0) {
                if (is_overlap()) {
                    printf("Found overlap in this step!\n");
                }
            }
        }
    }

    fclose(fp_density); // densities/...

    return 0;
}

/* Functions implementation */

/// Returns a random number between low and high (including exactly low, but not exactly high).
inline static double ran(double low, double high)
{
    return (high - low) * dsfmt_genrand() + low;
}

/// returns a (mod b), nonnegative, given that a >= -b is always true
inline static int pos_mod_i(int a, int b)
{
    return (a + b) % b;
}

/// returns a (mod b), nonnegative, given that a >= -b is always true
inline static double pos_mod_f(double a, double b)
{
    return fmod(a + b, b);
}

/// Scales the system with the scale factor
void scale(double scale_factor)
{
    for (int n = 0; n < n_particles; ++n) {
        // We use pointers to loop over the x, y, z members of the vec3_t type.
        double* pgarbage = &(r[n].x);
        for (int d = 0; d < NDIM; ++d)
            *(pgarbage + d) *= scale_factor;
    }
    for (int d = 0; d < NDIM; ++d)
        box[d] *= scale_factor;
}

/// This function returns the offset of the jth vertex (j = 0, ... , 7)
/// from the center of cube number i:     3----7
/// (y points into the screen)           /|   /|
///     z                               1-+--5 |
///     | y                             | 2--+-6
///     |/                              |/   |/
///     0----x                          0----4
/// The angle Phi is the angle ∠104 in the picture above, i.e. the angle of
/// "the z-axis of the cube" with "the x-axis of the cube"
vec3_t get_offset(int i, int j)
{
    vec3_t offset = vec3(-0.5 * (1 + CosPhi), -0.5, -0.5 * SinPhi);
    if (j & 4) //   x+ = 4, 5, 6, 7
        offset.x += 1;
    if (j & 2) //   y+ = 2, 3, 6, 7
        offset.y += 1;
    if (j & 1) { // z+ = 1, 3, 5, 7
        offset.z += SinPhi;
        offset.x += CosPhi;
    }
    // offset = v3_muls(offset, Edge_Length); // Edge_Length == 1

    offset = m4_mul_dir(m[i], offset);

    return offset;
}

/// Checks if there is overlap between cubes along axis, between cubes i, j.
/// Also the difference vector r2-r1 is given as it has already been calculated
bool is_collision_along_axis(vec3_t axis, int i, int j, vec3_t r2_r1)
{
    // The axis doesn't need to be normalized, source:
    // https://gamedevelopment.tutsplus.com/tutorials/collision-detection-using-the-separating-axis-theorem--gamedev-169)
    // axis = v3_norm(axis);

    double min1, min2, max1, max2, temp;
    min1 = max1 = v3_dot(axis, v3_add(r2_r1, get_offset(i, 0)));
    min2 = max2 = v3_dot(axis, get_offset(j, 0));
    for (int n = 1; n < 8; n++) {
        temp = v3_dot(axis, v3_add(r2_r1, get_offset(i, n)));
        min1 = fmin(min1, temp);
        max1 = fmax(max1, temp);
        temp = v3_dot(axis, get_offset(j, n));
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
bool is_overlap_between(int i, int j)
{
    // the next line shouldn't be necessary!
    // if (i == j) return false; // don't check on overlap with yourself!
    vec3_t r2_r1 = v3_sub(r[i], r[j]); // read as r2 - r1

    // We need to apply nearest image convention to r2_r1.
    // We use pointers to loop over the x, y, z members of the vec3_t type.
    double* pdist = &(r2_r1.x);
    for (int d = 0; d < NDIM; d++) {
        if (*(pdist + d) > 0.5 * box[d])
            *(pdist + d) -= box[d];
        if (*(pdist + d) < -0.5 * box[d])
            *(pdist + d) += box[d];
    }

    // If the cubes are more than their circumscribed sphere apart, they couldn't possibly overlap.
    // Similarly, if they are less than their inscribed sphere apart, they couldn't possibly NOT overlap.
    double len2 = v3_dot(r2_r1, r2_r1); // sqrtf is slow so test length^2
    if (len2 > (3 + 2 * CosPhi)) // * Edge_Length * Edge_Length) // Edge_Length == 1
        return false;
    if (len2 < SinPhi * SinPhi)
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
    vec3_t edges1[3], edges2[3];
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
            return false; // found separation, no overlap!
    // TODO: make smarter e.g. check only 2 points from first cube

    // overlap on all axes ==> the cubes overlap
    return true;
}

/// This function checks if the move- and volume-acceptances are too high or low,
/// and adjusts delta and deltaV accordingly.
void nudge_deltas(double mov, double vol, double rot)
{
    if (mov < 0.3)
        Delta *= 0.9; // acceptance too low  --> decrease delta
    if (mov > 0.4 && Delta < 0.5) // Edge_Length == 1
        Delta *= 1.1; // acceptance too high --> increase delta
    if (vol < 0.1)
        DeltaV *= 0.9;
    if (vol > 0.2)
        DeltaV *= 1.1;
    if (rot < 0.3)
        DeltaR *= 0.9;
    if (rot > 0.4 && DeltaR < M_PI / 4)
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
    if (scale_factor < 1.0) {
        for (int i = 0; i < n_particles; i++) {
            if (is_overlap_from(i)) {
                is_collision = true;
                break;
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
        update_CellLength();
        return 1; // succesful change
    } else {
        scale(1. / scale_factor); // move everything back
        return 0; // unsuccesful change
    }
}

/// This function reads the initial configuration of the system,
/// it reads data in the same format as it outputs data.
/// It also initializes SinPhi, CosPhi, Normal, ParticleVolume, r, and m (rot. mx.).
void read_data2(char* init_file)
{
    FILE* pFile = fopen(init_file, "r");
    if (NULL == pFile) {
        printf("file not found: %s\n", init_file);
        exit(1);
    }
    // read n_particles
    if (!fscanf(pFile, "%d", &n_particles)) {
        printf("failed to read num of particles\n");
        exit(2);
    }
    if (n_particles < 1 || n_particles > N) {
        printf("num particles %d, go into code and change N (max nparticles)\n", n_particles);
        exit(2);
    }
    double garbagef; // garbage (double) float

    // three zeroes which do nothing
    if (!fscanf(pFile, "%lf %lf %lf", &garbagef, &garbagef, &garbagef)) {
        printf("failed to read three zeroes\n");
        exit(2);
    }

    for (int d = 0; d < 9; ++d) { // dimensions of box
        if (d % 4 == 0) {
            if (!fscanf(pFile, "%lf\t", &(box[d / 4]))) {
                printf("failed to read dimensions of box\n");
                exit(2);
            }
        } else {
            if (!fscanf(pFile, "%lf\t", &garbagef)) {
                printf("failed to read dimensions of box\n");
                exit(2);
            }
        }
    }

    // now read the particles
    bool rf = true; // read flag. if false, something went wrong
    for (int n = 0; n < n_particles; ++n) {
        // We use pointers to loop over the x, y, z members of the vec3_t type.
        double* pgarbage = &(r[n].x);
        for (int d = 0; d < NDIM; ++d) // the position of the center of cube
            rf = rf && fscanf(pFile, "%lf\t", (pgarbage + d));
        rf = rf && fscanf(pFile, "%lf\t", &garbagef); // Edge_Length == 1
        for (int d1 = 0; d1 < NDIM; d1++) {
            for (int d2 = 0; d2 < NDIM; d2++) {
                rf = rf && fscanf(pFile, "%lf\t", &(m[n].m[d1][d2]));
            }
        }
        rf = rf && fscanf(pFile, "%lf", &garbagef);
        if (n == 0) {
            rf = rf && fscanf(pFile, "%lf\n", &Phi); // only read Phi once
        } else {
            if (EOF == fscanf(pFile, "%lf\n", &garbagef)) {
                printf("reached end of read file unexpectedly\n");
                exit(2);
            }
        }
    }
    if (!rf) { // if read flag false, something went wrong.
        printf("failed to read (one of the) particles\n");
        exit(2);
    }
    fclose(pFile);

    SinPhi = sin(Phi);
    CosPhi = cos(Phi);
    ParticleVolume = SinPhi; // Edge_Length == 1

    // now initialize the normals, put everything to zero first:
    for (int i = 0; i < 3; i++) {
        Normal[i].x = Normal[i].y = Normal[i].z = 0;
    }
    Normal[0].x = SinPhi; // normal on x-dir
    Normal[0].z = -CosPhi;
    Normal[1].y = 1.; // normal on y-dir
    Normal[2].z = 1.; // normal on z-dir
}

/// This function creates the initial configuration of the system.
/// It also initializes SinPhi, CosPhi, Normal, ParticleVolume, r, and m (rot. mx.).
void create_system()
{
    // initialize n_particles
    n_particles = CubesPerDim * CubesPerDim * CubesPerDim;
    if (n_particles > N) {
        printf("num particles too large, go into code and change N (max nparticles)\n");
        exit(2);
    }

    // initialize box
    for (int d = 0; d < NDIM; d++) {
        box[d] = CubesPerDim; // this will be changed in set_packingfraction
    }

    // initialize the particle positions (r) on a simple cubic lattice
    {
        int index = 0;
        for (int i = 0; i < CubesPerDim; i++) {
            for (int j = 0; j < CubesPerDim; j++) {
                for (int k = 0; k < CubesPerDim; k++) {
                    r[index++] = vec3(i + 0.5, j + 0.5, k + 0.5);
                }
            }
        }
    }

    SinPhi = sin(Phi);
    CosPhi = cos(Phi);
    ParticleVolume = SinPhi; // Edge_Length == 1

    // now initialize the normals, put everything to zero first:
    for (int i = 0; i < 3; i++) {
        Normal[i].x = Normal[i].y = Normal[i].z = 0;
    }
    Normal[0].x = SinPhi; // normal on x-dir
    Normal[0].z = -CosPhi;
    Normal[1].y = 1.; // normal on y-dir
    Normal[2].z = 1.; // normal on z-dir

    // now initialize the rotation matrices
    for (int i = 0; i < n_particles; i++) {
        for (int temp = 0; temp < 16; temp++)
            m[i].m[temp % 4][temp / 4] = 0; // everything zero first
        for (int d = 0; d < NDIM; d++) {
            m[i].m[d][d] = 1; // 1 on the diagonal
        }
    }
}

/// Initializes CellsPerDim, CellLength,
/// NumCubesInCell[][][], WhichCubesInCell[][][][] and InWhichCellIsThisCube[]
void initialize_cell_list(void)
{
    // first initialize the lists to zero
    for (int i = 0; i < NC; i++)
        for (int j = 0; j < NC; j++)
            for (int k = 0; k < NC; k++) {
                NumCubesInCell[i][j][k] = 0;
                for (int l = 0; l < MAXPC; l++) {
                    WhichCubesInCell[i][j][k][l] = -1;
                }
            }

    // the minimum cell length is the diameter of circumscribed sphere
    double min_cell_length = sqrt(3 + 2 * CosPhi);
    // the minimum box size is (the volume of a maximally packed box)^(1/3)
    double min_box_size = pow(n_particles * ParticleVolume, 1. / 3.);

    // CellsPerDim must not be too large, so that
    // box[0] / CellsPerDim = CellLength >= min_cell_length
    // therefore, rounding CellsPerDim down is okay.
    CellsPerDim = (int)(min_box_size / min_cell_length);
    // NC is the maximum number of cells
    CellsPerDim = (CellsPerDim < NC) ? CellsPerDim : NC;

    update_CellLength();
    // Update CellLength every volume change, so that
    // all particles stay in the same cell when the system is scaled

    // now assign each cube to the cell they are in,
    // and count how many cubes are in each cell.
    for (int i = 0; i < n_particles; i++) {
        int x = r[i].x / CellLength;
        int y = r[i].y / CellLength;
        int z = r[i].z / CellLength;
        // add particle i to WhichCubesInCell at the end of the list
        // and add one to the counter of cubes of this cell (hence the ++)
        WhichCubesInCell[x][y][z][NumCubesInCell[x][y][z]++] = i;
        // and keep track of in which cell this cube is
        InWhichCellIsThisCube[i] = NC * NC * x + NC * y + z;
    }
}

/// Must be called every succesful volume change, and during initialization
inline static void update_CellLength(void)
{
    CellLength = box[0] / CellsPerDim;
}

/// This function returns if cube number index overlaps, using cell lists
bool is_overlap_from(int index)
{
    bool is_collision = false;
    int cell = InWhichCellIsThisCube[index];
    // convert cell number to x, y, z coordinates
    int x = cell / (NC * NC);
    int y = (cell / NC) % NC;
    int z = cell % NC;
    // loop over all neighbouring cells
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            for (int k = -1; k < 2; k++) {
                // now loop over all cubes in this cell, remember periodic boundary conditions
                int loop_x = pos_mod_i(x + i, CellsPerDim);
                int loop_y = pos_mod_i(y + j, CellsPerDim);
                int loop_z = pos_mod_i(z + k, CellsPerDim);
                int num_cubes = NumCubesInCell[loop_x][loop_y][loop_z];
                for (int cube = 0; cube < num_cubes; cube++) {
                    int index2 = WhichCubesInCell[loop_x][loop_y][loop_z][cube];
                    // if checking your own cell, do not check overlap with yourself
                    if (index == index2) {
                        continue; // TODO: is this efficient enough?
                    }

                    if (is_overlap_between(index, index2)) {
                        is_collision = true;
                        // and break out of all loops
                        cube = N;
                        i = j = k = 2;
                    }
                }
            }
        }
    }
    return is_collision;
}

/// This moves a random particle in a cube of volume (2 * delta)^3.
/// Note this gives particles a tendency to move to one of the 8 corners of the cube,
/// however it obeys detailed balance.^{[citation needed]}
/// It checks for overlap using the cell lists.
/// returns 1 on succesful move, 0 on unsuccesful move.
int move_particle_cell_list(void)
{
    // first choose the cube and remember its position
    int index = (int)ran(0., n_particles); // goes from 0 to n_particles - 1
    vec3_t r_old = r[index];

    // move the cube
    // We use pointers to loop over the x, y, z members of the vec3_t type.
    double* pgarbage = &(r[index].x);
    for (int d = 0; d < NDIM; d++) {
        *(pgarbage + d) += ran(-Delta, Delta);
        // periodic boundary conditions happen here, in pos_mod_f. Since delta < box[dim],
        // the following expression will always return a positive number.
        // *(pgarbage + d) = fmodf(*(pgarbage + d) + box[d], box[d]);
        *(pgarbage + d) = pos_mod_f(*(pgarbage + d), box[d]);
        // TODO: maybe make faster by only checking on boundary cells
    }

    update_cell_list(index);

    // and check for overlaps
    if (is_overlap_from(index)) {
        r[index] = r_old; // move back
        update_cell_list(index); // and re-update the cell list. // TODO: make more efficient
        return 0; // unsuccesful move
    } else {
        // remember to update (Num/Which)CubesInCell and InWhichCellIsThisCube
        return 1; // succesful move
    }
}

/// This function updates the cell list. At this point, r[index] contains the new
/// (accepted) position, and we need to check if is in the same cell as before.
/// If not, update (Num/Which)CubesInCell and InWhichCellIsThisCube.
void update_cell_list(int index)
{
    int cell_old = InWhichCellIsThisCube[index];
    int x_old = cell_old / (NC * NC);
    int y_old = (cell_old / NC) % NC;
    int z_old = cell_old % NC;
    vec3_t r_new = r[index];
    int x_new = r_new.x / CellLength;
    int y_new = r_new.y / CellLength;
    int z_new = r_new.z / CellLength;
    if (x_new == CellsPerDim || y_new == CellsPerDim || z_new == CellsPerDim) {
        // DEBUG
        printf("new coordinate is exactly(ish) the box size\n");
        printf("cube %d moved to (%lf, %lf, %lf), but boxsize is %lf.\n", index, r_new.x, r_new.y, r_new.z, box[0]);
        printf("CellLength is %lf, so xyz is (%d, %d, %d).\n", CellLength, x_new, y_new, z_new);
        printf("celllength: %lf\nboxsize: %lf\ncellsperdim: %d\n\n", CellLength, box[0], CellsPerDim);
        for (int i = 0; i < CellsPerDim; i++) {
            for (int j = 0; j < CellsPerDim; j++) {
                for (int k = 0; k < CellsPerDim; k++) {
                    printf("%2d ", NumCubesInCell[i][j][k]);
                }
                printf("\n");
            }
            printf("\n");
        }

        // exit(6); // yeah so this actually happens for floats.
        x_new = x_new % CellsPerDim;
        y_new = y_new % CellsPerDim;
        z_new = z_new % CellsPerDim;
        // ENDEBUG
    }
    if (x_old == x_new && y_old == y_new && z_old == z_new) {
        return; // still in same box, don't have to change anything
    } else {
        // update in which cell this cube is
        InWhichCellIsThisCube[index] = NC * NC * x_new + NC * y_new + z_new;
        // update WhichCubesInCell, first check at what index the moved cube was
        int cube = 0;
        while (index != WhichCubesInCell[x_old][y_old][z_old][cube]) {
            cube++;
            if (cube >= MAXPC) {
                printf("infinite loop in update_cell_list\n");
                // DEBUG
                // printf("in cell %d %d %d, trying to find %d:\n", x_old, y_old, z_old, index);
                // for (int i = 0; i < NC; i++) {
                //     printf("%d ", WhichCubesInCell[x_old][y_old][z_old][i]);
                // }
                // printf("num of particles: %d\n", NumCubesInCell[x_old][y_old][z_old]);
                exit(5); // not so infinite anymore eh
                // ENDEBUG
            }
        }
        // now cube contains the index in the cell list

        // add the cube to the new cell and add one to the counter (hence ++)
        WhichCubesInCell[x_new][y_new][z_new][NumCubesInCell[x_new][y_new][z_new]++] = index;
        // and remove the cube in the old cell by replacing it with the last
        // in the list and remove one from the counter (hence --)
        int last_in_list = WhichCubesInCell[x_old][y_old][z_old][--NumCubesInCell[x_old][y_old][z_old]];
        WhichCubesInCell[x_old][y_old][z_old][cube] = last_in_list;
        return;
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
    } while (dist > 1 || dist == 0); // dist == 0 doesn't seem to happen
    vec3_t rot_axis = vec3(x[0], x[1], x[2]); // the axis doesn't need to be normalized
    mat4_t rot_mx = m4_rotation(ran(-DeltaR, DeltaR), rot_axis);

    // rotate the particle
    m[index] = m4_mul(rot_mx, m[index]);

    // and check for overlaps
    if (is_overlap_from(index)) {
        // overlap, rotate back
        m[index] = m4_mul(m4_invert_affine(rot_mx), m[index]);
        return 0; // unsuccesful rotation
    } else {
        // no overlap, rotation is accepted!
        // since the cube didn't move, there is no need to update the cell lists
        return 1; // succesful rotation
    }
}

void write_data(int step, FILE* fp_density, char datafolder_name[128])
{
    fprintf(fp_density, "%lf\n", n_particles * ParticleVolume / (box[0] * box[1] * box[2]));
    fflush(fp_density); // write the densities everytime we have one, otherwise it waits for ~400 lines

    char buffer[128];
    strcpy(buffer, datafolder_name);
    strcat(buffer, "/coords_step%07d.poly");

    char datafile[128];
    sprintf(datafile, buffer, step); // replace %07d with step and put in output_file.
    // char datafile[128];
    // strcpy(datafile, datafolder_name);
    // strcat(datafile, "/coords_step%07d.poly");
    // sprintf(datafile, datafile, step); // replace %07d with step and put in output_file.
    // THIS GOES WRONG
    
    FILE* fp = fopen(datafile, "w");
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
        // We use pointers to loop over the x, y, z members of the vec3_t type.
        double* pgarbage = &(r[n].x);
        for (int d = 0; d < NDIM; ++d)
            fprintf(fp, "%lf\t", *(pgarbage + d)); // the position of the center of cube
        fprintf(fp, "1\t"); // Edge_Length == 1
        for (int d1 = 0; d1 < NDIM; d1++) {
            for (int d2 = 0; d2 < NDIM; d2++) {
                fprintf(fp, "%lf\t", m[n].m[d1][d2]);
            }
        }
        fprintf(fp, "10 %lf\n", Phi); // 10 is for slanted cubes.
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

/// After putting the (slanted) cubes at the right position,
/// rotate each cube 90 degrees around a random axis-aligned axis a few times
void set_random_orientation(void)
{
    // first rotate every particle a number of times along one of the box' axes
    for (int i = 0; i < n_particles; i++) {
        for (int j = 0; j < 12; j++) { // TODO: is this random enough? legit?
            vec3_t axis = vec3(1, 0, 0);
            int random = (int)ran(0, 6);
            if (random < 2) {
                axis = vec3(0, 0, 1);
            } else if (random < 4) {
                axis = vec3(0, 1, 0);
            } else {
                axis = vec3(1, 0, 0);
            }
            if (random & 1) {
                axis = v3_muls(axis, -1);
            }
            m[i] = m4_mul(m4_rotation(M_PI / 2, axis), m[i]);
        }
    }
}

/// returns true if there is overlap in the system, else false.
bool is_overlap(void)
{
    for (int i = 0; i < n_particles; i++) {
        for (int j = i + 1; j < n_particles; j++) {
            if (is_overlap_between(i, j)) {
                return true;
            }
        }
    }
    return false;
}

/// reduces initial packing fraction if there is overlap
void remove_overlap(void)
{
    while (is_overlap()) {
        scale(1.01); // make the system bigger and check for collisions again
        printf("initial packing fraction lowered to %lf\n",
            n_particles * ParticleVolume / (box[0] * box[1] * box[2]));
    }
}

/// Put parsing the commandline in a function.
/// If something goes wrong, return != 0
int parse_commandline(int argc, char* argv[])
{
    if (argc != 7) {
        printf("need 6 arguments:\n");
        return 3;
    }
    if (EOF == sscanf(argv[3], "%d", &mc_steps)) {
        printf("reading mc_steps has failed\n");
        return 1;
    };
    if (EOF == sscanf(argv[4], "%lf", &packing_fraction)) {
        printf("reading packing_fraction has failed\n");
        return 1;
    };
    if (EOF == sscanf(argv[5], "%lf", &BetaP)) {
        printf("reading BetaP has failed\n");
        return 1;
    };
    if (EOF == sscanf(argv[6], "%lf", &Phi)) {
        printf("reading Phi has failed\n");
        return 1;
    };
    if (mc_steps < 100) {
        printf("mc_steps > 99\n");
        return 2;
    }
    if (packing_fraction <= 0 || packing_fraction > 1) {
        printf("0 < packing_fraction <= 1\n");
        return 2;
    }
    if (BetaP <= 0 || BetaP >= 100) {
        printf("0 < BetaP < 100\n");
        return 2;
    }
    if (Phi <= 0 || Phi > M_PI / 2) {
        printf("0 < Phi < 1.57079632679\n");
        return 2;
    }

    // read data from a file or create a system with CubesPerDim cubes per dimension
    if (strcmp(argv[1], "read") == 0 || strcmp(argv[1], "r") == 0) {
        printf("reading file %s...\n", argv[2]);
        read_data2(argv[2]);
        CubesPerDim = (int)pow(n_particles, 1. / 3.);
        IsCreated = false;
    } else if (strcmp(argv[1], "create") == 0 || strcmp(argv[1], "c") == 0) {
        CubesPerDim = atoi(argv[2]);
        if (CubesPerDim < 4 || CubesPerDim > 41) {
            printf("3 < CubesPerDim < 42, integer!\n%s", usage_string);
            return 1;
        }
        create_system(CubesPerDim);
        IsCreated = true;
        printf("creating system with %d cubes...\n", n_particles);
    } else {
        printf("error reading first argument: %s\n", argv[1]);
        return 1;
    }

    return 0; // no exceptions, run the program
}

/// DEBUG method
/// CellsPerDim, CellLength,
/// NumCubesInCell[][][], WhichCubesInCell[][][][] and InWhichCellIsThisCube[]
void print_celllists(void)
{
    printf("boxsize: %lf", box[0]);
    printf("\tCellLength: %lf", CellLength);
    printf("\tCellsPerDim: %d\n", CellsPerDim);

    printf("InWhichCellIsThisCube:\n");
    for (int i = 0; i < n_particles; i++) {
        int cell = InWhichCellIsThisCube[i];
        printf("%d: %d, %d, %d\n", i, cell / (NC * NC), (cell / NC) % NC, cell % NC);
    }

    printf("\nNumCubesInCell:\n");
    for (int i = 0; i < CellsPerDim; i++) {
        for (int j = 0; j < CellsPerDim; j++) {
            for (int k = 0; k < CellsPerDim; k++) {
                printf("%d, %d, %d: %d\n", i, j, k, NumCubesInCell[i][j][k]);
            }
        }
    }

    printf("\nWhichCubesInCell:\n");
    for (int i = 0; i < CellsPerDim; i++) {
        for (int j = 0; j < CellsPerDim; j++) {
            for (int k = 0; k < CellsPerDim; k++) {
                printf("%d, %d, %d: ", i, j, k);
                for (int l = 0; l < NumCubesInCell[i][j][k]; l++) {
                    printf("%d ", WhichCubesInCell[i][j][k][l]);
                }
                printf("\n");
            }
        }
    }
}