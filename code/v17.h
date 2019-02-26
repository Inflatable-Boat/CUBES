#ifndef V_17HEADER
#define V_17HEADER

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
} system_t;

void copy_system(system_t dest, system_t src)
{
    printf("\ncopying:\nr[0] source = %lf, %lf, %lf\nr[0] dest = %lf, %lf, %lf\n",
        src.r[0].x, src.r[0].y, src.r[0].z, dest.r[0].x, dest.r[0].y, dest.r[0].z);
    for (int i = 0; i < NDIM; i++) {
        dest.box[i] = src.box[i];
    }
    for (int i = 0; i < N; i++) {
        dest.r[i] = src.r[i];
        dest.m[i] = src.m[i];
        dest.InWhichCellIsThisCube[i] = src.InWhichCellIsThisCube[i];
    }
    printf("done copying:\nr[0] source = %lf, %lf, %lf\nr[0] dest = %lf, %lf, %lf\n",
        src.r[0].x, src.r[0].y, src.r[0].z, dest.r[0].x, dest.r[0].y, dest.r[0].z);
    for (int i = 0; i < NC; i++) {
        for (int j = 0; j < NC; j++) {
            for (int k = 0; k < NC; k++) {
                dest.NumCubesInCell[i][j][k] = src.NumCubesInCell[i][j][k];
                for (int l = 0; l < MAXPC; l++) {
                    dest.WhichCubesInCell[i][j][k][l] = src.WhichCubesInCell[i][j][k][l];
                }
            }
        }
    }
}

#endif // V_17HEADER