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
#ifndef INVFPI
#define INVFPI 0.07957747154594766788444188168625718101
#endif

#define NDIM 3
#define N 8000 // 20^3 should be plenty for now
#define NC 16 // max number of cells per direction
#define MAXPC 64 // max number of particles per cell
// WARNING!, don't make NC much larger than 32, because
// WhichCubesInCell is an array of size NC^3, this is already a million at NC == 32.

typedef struct {
    double re;
    double im;
} compl_t; // complex number struct
typedef struct {
    double nz;
    double si;
    double co;
    int n;
} bndT; // TODO: ask Frank wtf this is
typedef struct {
    int n;
    bndT* bnd;
} blistT; // list of bndTs
blistT* blist;

/* Tunable variables */
const double BONDLENGTHSQ = 1.55*1.55;

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
const char usage_string[] = "usage: program.exe CubesPerDim \
mc_steps packing_fraction BetaP Phi\n";

const int output_steps = 100;

/* Simulation variables */
// TODO: use malloc and pointers instead of global variables?
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
// functions from Frank's order parameter code (crystal.c)
double dotprod(compl_t* vec1, compl_t* vec2, int l);
float plgndr(int l, int m, double x);
double gammln(double xx);
double facs(int l, int m);
void order(int l, bndT* bnd, compl_t* res1, compl_t* res2);
compl_t* calc_order(void);

inline static int pos_mod_i(int a, int b);

// initialization
int parse_commandline(int argc, char* argv[]);
void read_data2(char* init_file, int is_not_first_time);
void initialize_cell_list(void);

// collision detection
bool is_overlap_between(int i, int j);
bool is_collision_along_axis(vec3_t axis, int i, int j, vec3_t r2_r1);
vec3_t get_offset(int i, int j);
bool is_overlap_from(int index);
inline static void update_CellLength(void);
bool is_overlap(void);

void print_celllists(void);

/* Main */

int main(int argc, char* argv[])
{
    dsfmt_seed(time(NULL));

    if (parse_commandline(argc, argv)) { // if ... parsing fails
        printf(usage_string);
        return 1;
    };

    char buffer[128] = "datafolder/";
    strcat(buffer, labelstring);
    char datafolder_name[128] = "";
    // replace all %d, %lf in buffer with values and put in datafolder_name
    sprintf(datafolder_name, buffer, CubesPerDim, packing_fraction, BetaP, Phi);
    printf("datafolder_name: %s\n", datafolder_name);

    char outputfile_name[128] = "";
    strcpy(outputfile_name, datafolder_name);
    strcat(outputfile_name, "/order.txt");
    FILE* fp_order = fopen(outputfile_name, "w");

    strcat(buffer, "/coords_step%07d.poly");

    for (int step = 0; step <= mc_steps; step += 100) {
        char datafile_name[128] = "";
        // replace all %d, %lf in buffer with values and put in datafile_name
        sprintf(datafile_name, buffer, CubesPerDim, packing_fraction, BetaP, Phi, 0);
        printf("datafile_name: %s\n", datafile_name);

        read_data2(datafolder_name, 0);
        if (is_overlap()) {
            printf("\n\n\tWARNING\n\n\nThe read file contains overlap.\nExiting.\n");
            exit(4);
        }

        compl_t*order = calc_order();
        fprintf(fp_order, "%lf ", 1);

        free(order);
    }

    fclose(fp_order);

    return 0;
}

/* Functions implementation */

/// returns a (mod b), nonnegative, given that a >= -b is always true
inline static int pos_mod_i(int a, int b)
{
    return (a + b) % b;
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

/// This function reads the initial configuration of the system,
/// it reads data in the same format as it outputs data.
/// It also initializes SinPhi, CosPhi, Normal, ParticleVolume, r, and m (rot. mx.).
void read_data2(char* init_file, int is_not_first_time)
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

    if (!is_not_first_time) { // only initialize these once
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

    initialize_cell_list();
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
                        continue;
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

/// Put parsing the commandline in a function.
/// If something goes wrong, return != 0
int parse_commandline(int argc, char* argv[])
{
    if (argc != 6) {
        printf("need 5 arguments:\n");
        return 3;
    }
    if (EOF == sscanf(argv[1], "%d", &CubesPerDim)) {
        printf("reading CubesPerDim has failed\n");
        return 1;
    };
    if (EOF == sscanf(argv[2], "%d", &mc_steps)) {
        printf("reading mc_steps has failed\n");
        return 1;
    };
    if (EOF == sscanf(argv[3], "%lf", &packing_fraction)) {
        printf("reading packing_fraction has failed\n");
        return 1;
    };
    if (EOF == sscanf(argv[4], "%lf", &BetaP)) {
        printf("reading BetaP has failed\n");
        return 1;
    };
    if (EOF == sscanf(argv[5], "%lf", &Phi)) {
        printf("reading Phi has failed\n");
        return 1;
    };
    if (CubesPerDim < 1 || CubesPerDim > 20) {
        printf("0 < CubesPerDim < 21\n");
        return 2;
    }
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

    return 0; // no exceptions, run the program
}

/************************************************
 *             DOTPROD
 * Dot product
 ***********************************************/
double dotprod(compl_t* vec1, compl_t* vec2, int l)
{
    double res = 0;
    int m;
    for (m = -l; m <= l; m++)
        res += (*(vec1 + m + l)).re * (*(vec2 + m + l)).re + (*(vec1 + m + l)).im * (*(vec2 + m + l)).im;
    return res;
}

/************************************************
 *             PLGNDR
 * Legendre polynomial
 ***********************************************/
float plgndr(int l, int m, double x)
{
    double fact, pll = 0.0, pmm, pmmp1, somx2;
    int i, ll;
    if (m < 0 || m > l || fabs(x) > 1.0)
        printf("Bad arguments in routine plgndr %i %i %f\n", l, m, fabs(x));
    pmm = 1.0;
    if (m > 0) {
        somx2 = sqrt((1.0 - x) * (1.0 + x));
        fact = 1.0;
        for (i = 1; i <= m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m)
        return pmm;
    else {
        pmmp1 = x * (2 * m + 1) * pmm;
        if (l == (m + 1))
            return pmmp1;
        else {
            for (ll = m + 2; ll <= l; ll++) {
                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return pll;
        }
    }
}

/************************************************
 *             GAMMLN
 * Log of the gamma function
 ***********************************************/
double gammln(double xx)
{
    double x, y, tmp, ser;
    static double cof[6] = { 76.18009172947146, -86.50532032941677,
        24.01409824083091, -1.231739572450155,
        0.1208650973866179e-2, -0.5395239384953e-5 };
    int j;
    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++)
        ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

/************************************************
 *             FACS
 * Calculate factorials
 ***********************************************/
double facs(int l, int m)
{
    static double* fac_table[14] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }; //max l=10
    int a, b;
    if (fac_table[l] == NULL) {
        fac_table[l] = malloc(sizeof(double) * (2 * l + 1));
        for (a = 0; a < 2 * l + 1; a++) {
            b = a - l;
            fac_table[l][a] = exp(gammln(l - b + 1) - gammln(l + b + 1));
        }
    }
    return fac_table[l][m + l];
}

/************************************************
 *             ORDER
 * Calculate q for a pair of particles
 ***********************************************/
void order(int l, bndT* bnd, compl_t* res1, compl_t* res2)
{
    double fc, p, f, s, r, sp, spp, c, cp, cpp;
    double z;
    int m = 0;
    z = bnd->nz;

    //for m=0
    p = plgndr(l, 0, z);
    fc = facs(l, 0);
    f = sqrt((2 * l + 1) * INVFPI * fc);
    r = p * f;
    (res1 + 0)->re += r;
    (res1 + 0)->im += 0;
    (res2 + 0)->re += r * minpow(l);
    (res2 + 0)->im += 0;
    s = 0;
    sp = 0;
    c = 1;
    cp = 0;

    for (m = 1; m <= l; m++) {
        //positive m
        p = plgndr(l, m, z);
        fc = facs(l, m);
        f = sqrt((2 * l + 1) * INVFPI * fc);
        r = p * f;
        cpp = cp;
        cp = c;
        if (m == 1)
            c = bnd->co;
        else
            c = 2.0 * bnd->co * cp - cpp; //some cosine tricks

        spp = sp;
        sp = s;
        if (m == 1)
            s = bnd->si;
        else
            s = 2.0 * bnd->co * sp - spp; //some sine tricks

        (res1 + m)->re += r * c;
        (res1 + m)->im += r * s;
        (res2 + m)->re += r * c;
        (res2 + m)->im += r * s;

        //negative m
        r *= minpow(m);
        (res1 - m)->re += r * c;
        (res1 - m)->im += -r * s;
        (res2 - m)->re += r * c;
        (res2 - m)->im += -r * s;

        //printf ("Test: %d, %lf, %lf (cumu: %lf, %lf) \n", m, r*c,-r*s, (res1-m)->re, (res1-m)->im);
    }
}

/************************************************
 *             CALC_ORDER
 * Calculate q_4 for all particles
 ***********************************************/
compl_t* calc_order(void)
{
    int i, j, m;
    compl_t *q1, *q2;
    const int l = 4;
    double temp;
    compl_t* orderp = (compl_t*)malloc(sizeof(compl_t) * n_particles * (l * 2 + 1));
    memset(orderp, (int)0.0, sizeof(compl_t) * n_particles * (l * 2 + 1));
    for (i = 0; i < n_particles; i++) { // TODO: ask Frank if particlestocount == n_part
        q1 = (orderp + i * (2 * l + 1) + l);
        // for (j = 0; j < blist[i].n; j++) { // TODO I assume this is looping over neighbours
        //     if (blist[i].bnd[j].n > i) {
        //         q2 = (orderp + blist[i].bnd[j].n * (2 * l + 1) + l);
        //         order(l, &(blist[i].bnd[j]), q1, q2);
        //     }
        // }
        for (j = i + 1; j < n_particles; j++) {
            double dist2 = v3_dot(r[i], r[j]);
            if (dist2 < BONDLENGTHSQ) {
                q2 = (orderp + j * (2 * l + 1) + l);
                order(l, argh, q1, q2);
            }
        }
    }
    //normalize vector
    for (i = 0; i < n_particles; i++) {
        temp = 1.0 / sqrt(dotprod(orderp + i * (2 * l + 1), orderp + i * (2 * l + 1), l));
        for (m = -l; m <= l; m++) {
            (*(orderp + i * (2 * l + 1) + m + l)).re *= temp;
            (*(orderp + i * (2 * l + 1) + m + l)).im *= temp;
        }
        //end
    }
    return orderp;
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