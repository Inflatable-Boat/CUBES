#include "math_3d.h" // https://github.com/arkanis/single-header-file-c-libs/blob/master/math_3d.h
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
// to be read from the command line
static int StepNumber; // the stepnumber the init file should be read from

static double Phi; // angle of slanted cube

const char usage_string[] = "usage: program.exe datafolder/<input> StepNumber datafolder/<output>\n";

/* Simulation variables */
static vec3_t r[N]; // position of center of cube
static mat4_t m[N]; // rotation matrix of cube
static double box[NDIM]; // dimensions of box
static vec3_t Normal[3]; // the normal vector of an unrotated cube. Normal[0] is the normal in the x-dir, etc.
static int n_particles = 0;
// static int CubesPerDim;
static double ParticleVolume;
static double CosPhi; // cos and sin of Phi appear a lot, and are expensive to calculate.
static double SinPhi; // Since Phi doesn't change, it's faster to calculate only once.

/* Functions */
inline static double ran(double low, double high);

int parse_commandline(int argc, char* argv[]);
void read_data(char* init_file);
void write_data(char buffer[], int expected_order);

/* Main */

int main(int argc, char* argv[])
{
    dsfmt_seed(time(NULL));

    if (parse_commandline(argc, argv)) { // if ... parsing fails
        printf(usage_string);
        return 1;
    };

    // orient the cubes all the same way and write this file to ...100.poly
    for (int i = 0; i < n_particles; i++) {
        m[i] = m4_identity();
    }
    write_data(argv[3], 1);

    // now rotate each particle around random box edge axes a few times,
    // this should still give perfect symmetry. write this to ...200.poly
    for (int i = 0; i < n_particles; i++) {
        for (int j = 0; j < 12; j++) {
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
    write_data(argv[3], 2);

    // now orient each cubes randomly and write this file to ...000.poly
    for (int i = 0; i < n_particles; i++) {
        for (int rot = 0; rot < 24; rot++) {
            // choose a random axis (by picking points randomly in a ball)
            double x[3], dist;
            do {
                for (int i = 0; i < 3; i++)
                    x[i] = ran(-1, 1);
                dist = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
            } while (dist > 1 || dist == 0); // dist == 0 doesn't seem to happen
            vec3_t rot_axis = vec3(x[0], x[1], x[2]); // the axis doesn't need to be normalized
            mat4_t rot_mx = m4_rotation(ran(-M_PI, M_PI), rot_axis);

            // rotate the particle
            m[i] = m4_mul(rot_mx, m[i]);
        }
    }
    write_data(argv[3], 0);

    return 0;
}

/* Functions implementation */

/// Returns a random number between low and high (including exactly low, but not exactly high).
inline static double ran(double low, double high)
{
    return (high - low) * dsfmt_genrand() + low;
}

/// This function reads the initial configuration of the system,
/// it reads data in the same format as it outputs data.
/// It also initializes SinPhi, CosPhi, Normal, ParticleVolume, r, and m (rot. mx.).
void read_data(char* init_folder)
{
    char init_file_with_path[128] = "datafolder/";
    strcat(init_file_with_path, init_folder);
    strcat(init_file_with_path, "/coords_step%07d.poly");
    // printf("init_folder: %s\ninit_file_with_path: %s\n", init_folder, init_file_with_path);

    char final_init_folder[255] = "";
    sprintf(final_init_folder, init_file_with_path, StepNumber);
    // printf("final: %s\n", final_init_folder);

    FILE* pFile = fopen(final_init_folder, "r");
    if (NULL == pFile) {
        printf("file not found: %s\n", final_init_folder);
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

void write_data(char buffer[], int expected_order)
{
    char init_file_with_path[128] = "datafolder/";
    strcat(init_file_with_path, buffer);
    strcat(init_file_with_path, "/coords_step%07d.poly");

    char buffer2[128] = "";
    sprintf(buffer2, init_file_with_path, 100 * expected_order);
    
    FILE* fp = fopen(buffer2, "w");
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

/// Put parsing the commandline in a function.
/// If something goes wrong, return != 0
int parse_commandline(int argc, char* argv[])
{
    if (argc != 4) {
        printf("need 3 arguments:\n");
        return 3;
    }
    if (EOF == sscanf(argv[2], "%d", &StepNumber)) {
        printf("reading StepNumber has failed\n");
        return 1;
    };

    printf("reading file datafolder/%s/coords_step%07d.poly...\n", argv[1], StepNumber);
    read_data(argv[1]);

    // argv[3] is output_datafolder
    char datafolder_name[128] = "datafolder/";
    strcat(datafolder_name, argv[3]);
#ifdef _WIN32
    mkdir("v12");
    mkdir(datafolder_name);
#elif __linux__
    // Linux needs me to set rights, this gives rwx to me and just r to all others.
    mkdir("v12", S_IRWXU | S_IRGRP | S_IROTH);
    mkdir(datafolder_name, S_IRWXU | S_IRGRP | S_IROTH);
#elif
    printf("Please use Linux or Windows\n");
    exit(3);
#endif


    return 0; // no exceptions, run the program
}