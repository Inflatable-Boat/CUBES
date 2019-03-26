#ifndef V_21HEADER
#define V_21HEADER

#include "math_3d.h" // https://github.com/arkanis/single-header-file-c-libs/blob/master/math_3d.h

#define NDIM 3
#define N 8000 // 20^3 should be plenty for now
#define NC 16 // max number of cells per direction
#define MAXPC 64 // max number of particles per cell
// WARNING!, don't make NC much larger than 32, because
// WhichCubesInCell is an array of size MAXPC * NC^3, this is already 2 million at NC == 32.

typedef struct {
    vec3_t r[N]; // position of center of cube
    mat4_t m[N]; // rotation matrix of cube
    double CellLength; // The length of a cell
    int NumCubesInCell[NC][NC][NC]; // how many cubes in each cell
    int WhichCubesInCell[NC][NC][NC][MAXPC]; // which cubes in each cell
    int InWhichCellIsThisCube[N]; // (a,b,c) == NC*NC*a + NC*b + c
    double box[NDIM]; // dimensions of box
    int clust_size; // size of largest cluster
    double energy; // bias potential energy
    vec3_t offsets[8]; // the offsets of the cube, they depend only on phi
} system_t;

void copy_system(system_t* dest, system_t* src)
{
    for (int i = 0; i < NDIM; i++) {
        dest->box[i] = src->box[i];
    }
    for (int i = 0; i < N; i++) {
        dest->r[i] = src->r[i];
        dest->m[i] = src->m[i];
        dest->InWhichCellIsThisCube[i] = src->InWhichCellIsThisCube[i];
    }
    for (int i = 0; i < NC; i++) {
        for (int j = 0; j < NC; j++) {
            for (int k = 0; k < NC; k++) {
                dest->NumCubesInCell[i][j][k] = src->NumCubesInCell[i][j][k];
                for (int l = 0; l < MAXPC; l++) {
                    dest->WhichCubesInCell[i][j][k][l] = src->WhichCubesInCell[i][j][k][l];
                }
            }
        }
    }
}

#endif // V_21HEADER