#include "math_3d.h"
#include "21.h"
#include <malloc.h>
// #include <math.h>
// #include <stdbool.h>
// #include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <time.h>
// #include <unistd.h>
#ifdef _WIN32
#include <direct.h> // mkdir on Windows
#elif __linux__
#include <sys/stat.h> // mkdir on Linux
#include <sys/types.h> // mkdir on Linux
#endif

#define MAX_PART_CELL 40
#define MAX_NEIGHBORS 50
#define INVFPI 0.07957747154594766788444188168625718101
#define M_PI 3.14159265358979323846
// #define N 8000
// #define NDIM 3

//tunable parameters:
double bndLength = 1.55; //distance cutoff for bonds
double bnd_cuttoff = 0.6; //Order to be called a correlated bond
int nbnd_cuttoff = 4; //Number of correlated bonds for a crystalline particle
double obnd_cuttoff = 0.0; //Order to be in the same cluster (0.0 for all touching clusters 0.9 to see defects)

//////////////////////////////////////////////////////////////// my additions
// const char usage_string[] = "usage: program.exe datafolder read_per output_per
// (t/transl or o/orient or b/both) (sl_norm or unsl_norm or edge)\n
// output_per = 0 means don't save snapshots\n";
// static double Phi;
// static int output_per;
// static int read_per;
// static int n_particles;
// static double ruud_box[3];
// static vec3_t ruud_r[N]; // position of center of cube
// static mat4_t ruud_m[N]; // rotation matrix of cube
// static vec3_t Normal[3];
// static double SinPhi;
// static double CosPhi;
// static double ParticleVolume;
// typedef enum { transl = 1,
//     orient = 2,
//     axes_slanted_normals = 4,
//     axes_unslanted_normals = 8,
//     axes_edges = 16 } order_mode_enum_t;
// order_mode_enum_t order_mode;
// static char datafolder_name[256];
extern int n_particles;
extern vec3_t Normal[3];
extern double CosPhi;
extern double SinPhi;
extern system_t* sim;
int order_mode; // transl = 1, sl = 2, unsl = 3, edge = 4

//////////////////////////////////////////////////////////////// end of my additions

typedef struct {
    double re;
    double im;
} compl_t;
typedef struct {
    double nz;
    double si;
    double co;
    int n;
} bnd_t;
typedef struct {
    int n;
    bnd_t* bnd;
} blist_t;
blist_t* blist;
int* numconn;
double xsize, ysize, zsize;

double bndLengthSq;

int n_cells;
double avr_d, max_d, min_cell;
int nx, ny, nz;

double percentage;
int maxsize;
int numclus;

typedef struct {
    double x;
    double y;
    double z;
    double xhalf;
    double yhalf;
    double zhalf;
    double oneoverx;
    double oneovery;
    double oneoverz;
    double min;
} tbox_t;
tbox_t box;

typedef struct {
    // double x;
    // double y;
    // double z;
    double d;
    double r;
    int cell;
    char c;
    char str[256];
} tpart_t;
tpart_t* part;

typedef struct scell {
    double x1;
    double y1;
    double z1;
    double x2;
    double y2;
    double z2;
    int buren[27];
    int n;
    int particles[MAX_PART_CELL];
} tcells_t;
tcells_t* cells;

double sqr(double x)
{
    return x * x;
}

/************************************************
 *             PLGNDR
 * Legendre polynomial
 ***********************************************/
float plgndr(int l, int m, double x)
{
    double fact, pll = 0.0, pmm, pmmp1, somx2;
    int i, ll;
    if (m < 0 || m > l || fabs(x) > 1.0) {
        printf("Bad arguments in routine plgndr %i %i %f\n", l, m, fabs(x));
        x *= 0.99999;
    }
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
 * Calculate factorials (l - m)! / (l + m)!
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
 *             MINPOW
 * Powers of -1
 ***********************************************/
double minpow(int m)
{
    if ((m & 1) == 1)
        return -1.0;
    else
        return 1.0;
}

/************************************************
 *             ORIENT_ORDER
 * Calculate i for a particle along one axis
 ***********************************************/
void orient_order(int l, int i, compl_t* res, int axis)
{
    vec3_t dir; // = v3_norm(m4_mul_dir(ruud_m[i], Normal[axis]));
    // order_mode: transl = 1, sl = 2, unsl = 3, edge = 4
    if (order_mode == 2) {
        // take the normal of the slanted cube and rotate it
        dir = m4_mul_dir(sim->m[i], Normal[axis]);
    } else if (order_mode == 3) {
        // we need to construct the new vector. we do it using pointers:
        dir = vec3(0, 0, 0);
        // now set the axis-th member of dir to 1 in the ugliest way possible
        *(&(dir.x) + axis) = 1.;
        // now we have constructed the vector, rotate it
        dir = m4_mul_dir(sim->m[i], dir);
    } else if (order_mode == 4) {
        // we need to construct an edge of the cube
        if (axis == 0) {
            dir = vec3(1, 0, 0);
        } else if (axis == 1) {
            dir = vec3(0, 1, 0);
        } else if (axis == 2) {
            dir = vec3(CosPhi, 0, SinPhi);
        } else {
            printf("this message shouldn't be visible 2\n");
            exit(6);
        }
        // now we have constructed the vector, rotate it
        dir = m4_mul_dir(sim->m[i], dir);
    } else {
        printf("this message shouldn't be visible 1\n");
        exit(5);
    }
    double fc, p, f, s, r, sp, spp, c, cp, cpp;
    double z;
    int m = 0;
    // z = bnd->nz;
    z = dir.z;
    double si, co; // sin and cosine of azimuthal angle
    if (dir.x == 0 && dir.y == 0) {
        si = 0;
        co = 1;
    } else {
        double dxy = 1.0 / sqrt(dir.x * dir.x + dir.y * dir.y);
        si = dir.y * dxy;
        co = dir.x * dxy;
    }

    //for m=0
    p = plgndr(l, 0, z);
    fc = facs(l, 0);
    f = sqrt((2 * l + 1) * INVFPI * fc);
    r = p * f;
    (res + 0)->re += r;
    (res + 0)->im += 0;
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
            c = co;
        else
            c = 2.0 * co * cp - cpp; //some cosine tricks

        spp = sp;
        sp = s;
        if (m == 1)
            s = si;
        else
            s = 2.0 * co * sp - spp; //some sine tricks

        (res + m)->re += r * c;
        (res + m)->im += r * s;

        //negative m
        r *= minpow(m);
        (res - m)->re += r * c;
        (res - m)->im += -r * s;

        //printf ("Test: %d, %lf, %lf (cumu: %lf, %lf) \n", m, r*c,-r*s, (res1-m)->re, (res1-m)->im);
    }
}

/************************************************
 *             TRANSL_ORDER
 * Calculate q for a pair of particles
 ***********************************************/
void transl_order(int l, bnd_t* bnd, compl_t* res1, compl_t* res2)
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
 *             DOTPROD
 * Dot product
 ***********************************************/
// expects pointer to the start of vectors
double dotprod(compl_t* vec1, compl_t* vec2, int l)
{
    double res = 0;
    int m;
    for (m = -l; m <= l; m++)
        res += (*(vec1 + m + l)).re * (*(vec2 + m + l)).re + (*(vec1 + m + l)).im * (*(vec2 + m + l)).im;
    return res;
}

/************************************************
 *             COORDS2CELL
 * Find cell associated with particle position
 ***********************************************/
int coords2cell(double x, double y, double z)
{
    int i = ((int)(nx * (x + box.xhalf) * box.oneoverx) % nx) * ny * nz + ((int)(ny * (y + box.yhalf) * box.oneovery) % ny) * nz + ((int)(nz * (z + box.zhalf) * box.oneoverz) % nz);
    return i;
}

/************************************************
 *             UPDATE_NBLISTP
 * Find neighbors of a particle
 ***********************************************/
void update_nblistp(int p)
{
    int i, j, c, cellp, id;
    bnd_t* bnd;
    double d, dxy, dx, dy, dz;
    cellp = part[p].cell;
    blist[p].n = 0;
    for (i = 26; i >= 0; i--) {
        c = cells[cellp].buren[i];
        for (j = cells[c].n - 1; j >= 0; j--) {
            id = cells[c].particles[j];
            if (id != p) {
                dx = sim->r[p].x - sim->r[id].x;
                dy = sim->r[p].y - sim->r[id].y;
                dz = sim->r[p].z - sim->r[id].z;

                if (dx > box.xhalf)
                    dx -= box.x;
                else if (dx < -box.xhalf)
                    dx += box.x;
                if (dy > box.yhalf)
                    dy -= box.y;
                else if (dy < -box.yhalf)
                    dy += box.y;
                if (dz > box.zhalf)
                    dz -= box.z;
                else if (dz < -box.zhalf)
                    dz += box.z;

                d = dx * dx + dy * dy + dz * dz;

                if (d < bndLengthSq) {
                    bnd = &(blist[p].bnd[blist[p].n]);
                    blist[p].n++;
                    //          printf ("Test: %d, %d, %d, %d, %d\n", p, id, blist[p].n, i, j);
                    bnd->n = id;
                    if (id > p) {
                        d = 1.0 / sqrt(d);
                        //            bnd->nx = dx*d;
                        //            bnd->ny = dy*d;
                        bnd->nz = dz * d;
                        if (dx == 0 && dy == 0) {
                            bnd->si = 0;
                            bnd->co = 1;
                        } else {
                            dxy = 1.0 / sqrt(dx * dx + dy * dy);
                            bnd->si = dy * dxy;
                            bnd->co = dx * dxy;
                        }
                    }
                }
            }
        }
    }
}

/************************************************
 *             UPDATE_NBLIST
 * Make neighbor list for all particles
 ***********************************************/
void update_nblist(void)
{
    int p;
    for (p = n_particles - 1; p >= 0; p--) {
        update_nblistp(p);
    }
}

/************************************************
 *             INIT_NBLIST
 * Initialize the neighbor list
 ***********************************************/
// mallocs blist, numconn, blist[i].bnd
void init_nblist(void)
{
    int p;
    blist = (blist_t*)malloc(sizeof(blist_t) * n_particles);
    numconn = (int*)malloc(sizeof(int) * n_particles);
    for (p = n_particles - 1; p >= 0; p--) {
        blist[p].bnd = (bnd_t*)malloc(sizeof(bnd_t) * MAX_NEIGHBORS);
        numconn[p] = 0;
        update_nblistp(p);
    }
}

/************************************************
 *             INIT_CELLS
 * Initialize the cell list
 ***********************************************/
// mallocs cells
void init_cells(void)
{
    n_cells = 0;
    int x, y, z, i, j, k;
    int dx, dy, dz, a, b, c;
    // double cellsize = max_d;
    double cellsize = bndLength;
    // if (cellsize < bndLength)
    //     cellsize = bndLength;
    double ruud_min_cell_size = sqrt(3 + 2 * CosPhi);
    if (cellsize < ruud_min_cell_size)
        cellsize = ruud_min_cell_size;
    nx = floor(box.x / cellsize);
    ny = floor(box.y / cellsize);
    nz = floor(box.z / cellsize);
    min_cell = box.x / nx;
    if (box.y / ny < min_cell)
        min_cell = box.y / ny;
    if (box.z / nz < min_cell)
        min_cell = box.z / nz; //set min_cell
    // printf("cells: %d, %d, %d (%lf)\n", nx, ny, nz, cellsize);
    n_cells = nx * ny * nz;
    // printf ("Cells: %d (%lf)\n", n_cells, cellsize);
    cells = (tcells_t*)malloc(n_cells * sizeof(tcells_t));
    for (x = 0; x < nx; x++)
        for (y = 0; y < ny; y++)
            for (z = 0; z < nz; z++) {
                i = x * ny * nz + y * nz + z;
                cells[i].x1 = x * box.x / nx - box.xhalf;
                cells[i].y1 = y * box.y / ny - box.yhalf;
                cells[i].z1 = z * box.z / nz - box.zhalf;
                cells[i].x2 = (x + 1) * box.x / nx - box.xhalf;
                cells[i].y2 = (y + 1) * box.y / ny - box.yhalf;
                cells[i].z2 = (z + 1) * box.z / nz - box.zhalf;
                cells[i].n = 0;
                k = 0;
                for (a = -1; a < 2; a++)
                    for (b = -1; b < 2; b++)
                        for (c = -1; c < 2; c++) { //adding one self as a neighbor now (handy or not??)
                            dx = a;
                            dy = b;
                            dz = c;
                            if (x + dx < 0)
                                dx = nx - 1;
                            if (x + dx > nx - 1)
                                dx = -nx + 1;
                            if (y + dy < 0)
                                dy = ny - 1;
                            if (y + dy > ny - 1)
                                dy = -ny + 1;
                            if (z + dz < 0)
                                dz = nz - 1;
                            if (z + dz > nz - 1)
                                dz = -nz + 1;
                            j = (x + dx) * ny * nz + (y + dy) * nz + (z + dz);
                            //      printf("%i   %i %i %i  %i %i %i\n",j,dx,dy,dz, x+dx,y+dy,z+dz);
                            cells[i].buren[k] = j;
                            k++;
                        }
            }
    for (i = 0; i < n_particles; i++) {
        j = coords2cell(sim->r[i].x, sim->r[i].y, sim->r[i].z);
        if (cells[j].n + 1 >= MAX_PART_CELL) {
            printf("ERROR: Too many particles in a cell!\n");
            exit(7);
        }
        cells[j].particles[cells[j].n] = i;
        part[i].cell = j;
        cells[j].n++;
    }
}

/************************************************
 *             GETLINE_NUMBER
 * Read a line from a file
 * Ignore lines starting with #
 ***********************************************/
void getline_number(char* str, int n, FILE* f)
{
    int comment = 1;
    while (comment) {
        char* res = fgets(str, n, f);
        if (!res)
            return;
        if (str[0] != '#')
            comment = 0;
    }
}

/************************************************
 *             INIT_BOX
 * Initialize box parameters
 ***********************************************/
void init_box(double x_box, double y_box, double z_box)
{
    box.x = x_box;
    box.y = y_box;
    box.z = z_box;
    box.min = box.x;
    if (box.min > box.y)
        box.min = box.y;
    if (box.min > box.z)
        box.min = box.z;
    box.xhalf = box.x * 0.5;
    box.yhalf = box.y * 0.5;
    box.zhalf = box.z * 0.5;
    box.oneoverx = 1.0 / box.x;
    box.oneovery = 1.0 / box.y;
    box.oneoverz = 1.0 / box.z;
}

/************************************************
 *             CONVERT_DATA
 * converts data from my type to this file's weird type
 ***********************************************/
// mallocs part,
// indirect malloc cells, blist, numconn, blist[i].bnd
int convert_data()
{
    // n_part = particlestocount = n_particles;
    part = malloc(n_particles * sizeof(tpart_t));
    for (int n = 0; n < n_particles; n++) {
        // part[n].x = sim->r[n].x; // just use sim->r
        // part[n].y = sim->r[n].y;
        // part[n].z = sim->r[n].z;
        part[n].c = 'a';
        part[n].d = 1.0;
    }
    init_box(sim->box[0], sim->box[1], sim->box[2]);

    bndLengthSq = bndLength * bndLength;
    init_cells();
    init_nblist();
    return 0;
}

/************************************************
 *             CALC_ORIENT_ORDER
 * Calculate i_4 for all particles
 ***********************************************/
// mallocs (return value) orderp
compl_t* calc_orient_order(void)
{
    compl_t *q1, *q2, *q3;
    const int l = 4;
    // allocate two spaces too many to make place for q2, q3
    compl_t* orderp = (compl_t*)malloc(sizeof(compl_t) * (n_particles + 2) * (l * 2 + 1));
    memset(orderp, (int)0.0, sizeof(compl_t) * (n_particles + 2) * (l * 2 + 1));
    for (int i = 0; i < n_particles; i++) {
        q1 = (orderp + i * (2 * l + 1) + l);
        q2 = (orderp + (i + 1) * (2 * l + 1) + l); // temporary place for q2, q3
        q3 = (orderp + (i + 2) * (2 * l + 1) + l); // will be overwritten next loop

        // the folling gives order vectors i_{4,m} corresponding to the three particle axes
        // and stores them in q1, q2, q3. Next we have to add q2, q3 to q1 and normalize
        orient_order(l, i, q1, 0);
        orient_order(l, i, q2, 1);
        orient_order(l, i, q3, 2);
        for (int m = -l; m <= l; m++) {
            (q1 + m)->re += (q2 + m)->re + (q3 + m)->re;
            (q1 + m)->im += (q2 + m)->im + (q3 + m)->im;
        }
        //normalize vector
        double temp = 1.0 / sqrt(dotprod(q1 - l, q1 - l, l));
        for (int m = -l; m <= l; m++) {
            (q1 + m)->re *= temp;
            (q1 + m)->im *= temp;
        }
    }
    return orderp;
}

/************************************************
 *             CALC_TRANSL_ORDER
 * Calculate q_4 for all particles
 ***********************************************/
// mallocs (return value) orderp
compl_t* calc_transl_order(void)
{
    int i, j, m;
    compl_t *q1, *q2;
    const int l = 4;
    double temp;
    compl_t* orderp = (compl_t*)malloc(sizeof(compl_t) * n_particles * (l * 2 + 1));
    memset(orderp, (int)0.0, sizeof(compl_t) * n_particles * (l * 2 + 1));
    for (i = 0; i < n_particles; i++) {
        q1 = (orderp + i * (2 * l + 1) + l);
        for (j = 0; j < blist[i].n; j++) {
            if (blist[i].bnd[j].n > i) {
                q2 = (orderp + blist[i].bnd[j].n * (2 * l + 1) + l);
                transl_order(l, &(blist[i].bnd[j]), q1, q2);
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

/************************************************
 *             CALC_CONN
 * Find crystalline bonds between particles
 ***********************************************/
// mallocs (return value) conn
int* calc_conn(compl_t* orderp) //calculates "connected" particles
{
    int i, j, z, np = 0;
    const int l = 4;
    int* conn = malloc(sizeof(int) * n_particles);
    for (i = 0; i < n_particles; i++) {
        z = 0;
        for (j = 0; j < blist[i].n; j++) {
            if (dotprod(orderp + i * (2 * l + 1), orderp + blist[i].bnd[j].n * (2 * l + 1), l) > bnd_cuttoff) {
                z++;
                //        printf ("Test %d, %d\n", i, particlestocount);
            }
        }
        if (z >= nbnd_cuttoff) {
            np++;
            conn[i] = 1;
        } else {
            conn[i] = 0;
        }
        numconn[i] = z;
    }
    return conn;
}

/************************************************
 *             SAVE_CLUSS_DATA
 * Output order data: every 100 steps output 6 lines:
 * frac X, # clusters, max cluster size
 * # of neighbours histogram, 0-50
 * # of crystalline neighbours histogram, 0-50
 * q_4.q_4 of all bonds histogram
 * q_4.q_4 of all bonds, averaged for each cube histogram
 * all cluster sizes
 ***********************************************/
void save_cluss_data(int step, int* cluss, int* size, int big, int nn, int mode, compl_t* orderp)
{
    int i;
    const int l = 4;

    char buffer[256] = ""; // the coordinates filename
    char buffern[256] = ""; // the other data filename
    strcat(buffer, datafolder_name);
    strcat(buffern, datafolder_name);
    if (mode & transl) {
        strcat(buffer, "/clust_transl");
        strcat(buffern, "/v22_transl");
    } else if (mode & orient) {
        strcat(buffer, "/clust_orient");
        strcat(buffern, "/v22_orient");
        if (order_mode & axes_slanted_normals) {
            strcat(buffer, "_sl_normals");
            strcat(buffern, "_sl_normals");
        } else if (order_mode & axes_unslanted_normals) {
            strcat(buffer, "_unsl_normals");
            strcat(buffern, "_unsl_normals");
        } else if (order_mode & axes_edges) {
            strcat(buffer, "_edges");
            strcat(buffern, "_edges");
        } else {
            printf("invalid axes order_mode in save_cluss: %d", order_mode);
            printf(" mode: %d\n", mode);
            exit(4);
        }
    } else {
        printf("invalid order_mode in save_cluss: %d", order_mode);
        printf(" mode: %d\n", mode);
        exit(3);
    }
    strcat(buffer, "6");
    strcat(buffern, "6");
    strcat(buffer, "_coords%07d.poly");

    char fn[256] = "";
    sprintf(fn, buffer, step);

    if (!step) { // if first step
        remove(buffern); // remove the data-for-each-step-file
    }

    FILE* data_file;
    if ((data_file = fopen(buffern, "a")) == NULL) {
        printf("couldn't open clust_file buffern = %s\n", buffern);
        exit(-2);
    }
    fprintf(data_file, "%lf  %d  %d\n", percentage, numclus, maxsize);

    // now make and print the histogram of number of neighbours to file
    // note: prints fraction of cubes with this amount of neighbours
    {
        // make local scope to avoid clutter of random variables
        // make a histogram of number of neighbours:
        double nb_hist[MAX_NEIGHBORS + 1];
        memset(nb_hist, 0, (MAX_NEIGHBORS + 1) * sizeof(double));
        for (int i = 0; i < n_particles; i++) {
            // blist[i].n is the # of neighbours of particle i
            nb_hist[blist[i].n]++;
        }
        double oneovern_particles = 1.0 / n_particles;
        for (int i = 0; i <= MAX_NEIGHBORS; i++) {
            nb_hist[i] *= oneovern_particles;
        }
        for (int i = 0; i <= MAX_NEIGHBORS; i++) {
            fprintf(data_file, "%lf ", nb_hist[i]);
        }
        fprintf(data_file, "\n");
    }

    // now make and print the histogram of number of crystalline neighbours to file
    // note: prints fraction of cubes with this amount of neighbours
    {
        // make local scope to avoid clutter of random variables
        // make a histogram of number of neighbours:
        double xtal_nb_hist[MAX_NEIGHBORS + 1];
        memset(xtal_nb_hist, 0, (MAX_NEIGHBORS + 1) * sizeof(double));
        for (int i = 0; i < n_particles; i++) {
            // numconn[i] is the # of crystalline neighbours of particle i
            xtal_nb_hist[numconn[i]]++;
        }
        double oneovern_particles = 1.0 / n_particles;
        for (int i = 0; i <= MAX_NEIGHBORS; i++) {
            xtal_nb_hist[i] *= oneovern_particles;
        }
        for (int i = 0; i <= MAX_NEIGHBORS; i++) {
            fprintf(data_file, "%lf ", xtal_nb_hist[i]);
        }
        fprintf(data_file, "\n");
    }

    // now we make and print the |q_4|^2 histogram.
    // note: prints fraction of cubes with this |q_4|^2
    {
        const int NBINS = 100;
        double q4_hist[NBINS];
        memset(q4_hist, 0, NBINS * sizeof(double));
        int count = 0;
        for (int i = 0; i < n_particles; i++) {
            // printf("particle i: %d\n", i);
            for (int j = 0; j < blist[i].n; j++) {
                // printf("j = %d ", j);
                if (blist[i].bnd[j].n > i) {
                    // printf("yes blist.bnd.n > %d\n", i);
                    double order = dotprod(orderp + i * (2 * l + 1),
                        orderp + blist[i].bnd[j].n * (2 * l + 1), l);
                    // printf("order = %lf\n", order);
                        
                    if (order > 1.00001 || order < -1.00001) {
                        printf("unexpected order dotprod between %d and %d: %lf\nExiting.\n", i, blist[i].bnd[j].n, order);
                        exit(7);
                    }
                    if (order >= 1)
                        order = 0.999999; // to avoid rounding errors in the next step
                    if (order <= -1)
                        order = -0.999999;
                    q4_hist[(int) (NBINS * (order + 1) * 0.5)]++;
                    count++;
                }
            }
        }
        double oneovercount = 1;
        if (count)
            oneovercount = 1.0 / count;
        for (int i = 0; i < NBINS; i++) {
            q4_hist[i] *= oneovercount;
        }
        fprintf(data_file, "%d ", count);
        for (int i = 0; i < NBINS; i++) {
            fprintf(data_file, "%lf ", q4_hist[i]);
        }
        fprintf(data_file, "\n");
    }

    // now we make and print the |q_4|^2 histogram, averaged per particle.
    // note: prints fraction of cubes with this |q_4|^2
    {
        const int NBINS = 100;
        double avg_q4_hist[NBINS];
        memset(avg_q4_hist, 0, NBINS * sizeof(double));
        for (int i = 0; i < n_particles; i++) {
            // printf("particle i: %d\n", i);
            double order = 0;
            for (int j = 0; j < blist[i].n; j++) {
                // printf("j = %d ", j);
                order += dotprod(orderp + i * (2 * l + 1),
                    orderp + blist[i].bnd[j].n * (2 * l + 1), l);
                // printf("order = %lf\n", order);
            }
            if (blist[i].n) {
                order /= blist[i].n; // take the average
            }
            if (order > 1.00001 || order < -1.00001) {
                printf("unexpected order of particle %d: %lf\nExiting.\n", i, order);
                exit(7);
            }
            if (order >= 1)
                order = 0.999999; // to avoid rounding errors in the next step
            if (order <= -1)
                order = -0.999999;
            avg_q4_hist[(int) (NBINS * (order + 1) * 0.5)]++;
        }
        double oneovern_particles = 1.0 / n_particles;
        for (int i = 0; i < NBINS; i++) {
            avg_q4_hist[i] *= oneovern_particles;
        }
        for (int i = 0; i < NBINS; i++) {
            fprintf(data_file, "%lf ", avg_q4_hist[i]);
        }
        fprintf(data_file, "\n");
    }

    // now make and print all clustersizes
    for (int i = 1; i < nn; i++) {
        fprintf(data_file, "%d ", size[i]);
    }
    fprintf(data_file, "\n");

    fclose(data_file);

    if (output_per == 0) { // 0 means no save
        return;
    }
    FILE* clust_file;
    if ((step % (output_per)) == 0) {
        if ((clust_file = fopen(fn, "w")) == NULL) {
            printf("couldn't open clust_file fn = %s\n", fn);
            exit(-2);
        }
    } else {
        return;
    }
    // printf("step %d, saving to file %s\n", step, fn);
    printf("step %d, ", step);
    fprintf(clust_file, "%d\n0 0 0\n", n_particles);
    fprintf(clust_file, "%lf 0 0\n0 %lf 0\n0 0 %lf\n", box.x, box.y, box.z);

    // printf("defining sorta ");
    int* sorta = malloc(sizeof(int) * nn);
    // printf("and rank\n");
    int* rank = malloc(sizeof(int) * nn);
    // printf("looping over sorta\n");
    for (i = 0; i < nn; i++)
        sorta[i] = i;

    int cmpr(const void* a, const void* b)
    {
        return -size[*((int*)a)] + size[*((int*)b)];
    }

    // printf("sorting sorta[] if nn = %d >= 2\n", nn);
    // if (nn >= 2)
    //     qsort(sorta, nn, sizeof(int), &cmpr);
    // printf("sorting sorta[]\n");
    qsort(sorta, nn, sizeof(int), &cmpr);

    for (i = 0; i < nn; i++) {
        int bha;
        for (bha = 0; sorta[bha] != i; bha++)
            ;
        rank[i] = bha;
    }

    // printf("saving to files\n");
    for (i = 0; i < n_particles; i++) {
        int rnk = 0;
        if (cluss[i] && size[cluss[i]] > 2) {
            rnk = rank[cluss[i]];
            if (rnk > 0)
                rnk = ((rnk - 1) % 25) + 1;
            part[i].d = 1.0;
        } else {
            part[i].d = 0.1;
        }
        fprintf(clust_file, "%lf %lf %lf %lf ", sim->r[i].x, sim->r[i].y, sim->r[i].z, part[i].d);
        for (int d1 = 0; d1 < NDIM; d1++) {
            for (int d2 = 0; d2 < NDIM; d2++) {
                fprintf(clust_file, "%lf ", ruud_m[i].m[d1][d2]);
            }
        }
        // 10 = slanted cube, phi = angle, rnk = color
        fprintf(clust_file, "10 %lf %d\n", Phi, rnk);
    }
    fclose(clust_file);
    free(sorta);
    free(rank);
}

/************************************************
 *             CALC_CLUSTERS
 * Find clusters of bonded particles
 ***********************************************/
void calc_clusters(int* conn, compl_t* orderp, int mode)
{
    int* cluss = malloc(sizeof(int) * n_particles);
    int* size = malloc(sizeof(int) * n_particles);
    int i, cs, cn = 1, big = 0;
    // int bc = -1;
    const int l = 4;
    // printf("definition setcluss\n");
    void setcluss(int pn)
    {
        cluss[pn] = cn;
        int jj;
        for (jj = 0; jj < blist[pn].n; jj++) {
            int tmp = blist[pn].bnd[jj].n;
            if (conn[tmp] != 0 && cluss[tmp] == 0 && dotprod(orderp + pn * (2 * l + 1), orderp + tmp * (2 * l + 1), l) > obnd_cuttoff) {
                //0.9 gives nice results 0.6 gives all touching nuclei as one big nuclei
                cs++;
                setcluss(tmp);
            }
        }
    }

    // printf("initializing cluss[] and size[] to zero\n");
    for (i = 0; i < n_particles; i++) {
        cluss[i] = 0;
        size[i] = 0;
    }

    // printf("loop to setcluss(i) and initialize cn\n");
    for (i = 0; i < n_particles; i++) {
        cs = 0;
        if (conn[i] == 1 && cluss[i] == 0) {
            cs++; //has at least size 1;
            setcluss(i);
            size[cn] = cs;
            if (cs > big) {
                big = cs;
                // bc = cn;
            }
            cn++;
        }
    }

    // printf("calculating average cluster size\n");
    //calculate average cluster size
    /* int tcs = 0;
    for (i = 0; i < cn; i++) {
        tcs += size[i];
    }

    percentage = tcs / (double)n_particles; */
    maxsize = big;
    numclus = cn - 1;

    // final point where data starts getting freed
    free(cluss);
    free(size);
}

/************************************************
 *             MAIN
 ***********************************************/
int biggest_cluster_size_and_order(int what_order)
{
    order_mode = what_order; // transl = 1, sl = 2, unsl = 3, edge = 4
    compl_t* order = NULL;
    int* connections = NULL;

    convert_data(); // mallocs part, indirectly cells, blist, numconn, blist[i].bnd

    if (order_mode == 1) { // translational order
        order = calc_transl_order(); // mallocs order
        connections = calc_conn(order); // mallocs connections
    } else if (order_mode >= 2) { // orientational order
        order = calc_orient_order(); // mallocs order
        connections = calc_conn(order); // mallocs connections
    }
    calc_clusters(connections, order, order_mode); // and save if output_per > 0
    free(order);
    free(connections);
    free(part);
    free(cells);
    free(numconn);

    for (int i = 0; i < n_particles; i++)
        free(blist[i].bnd);
    free(blist);

    return maxsize;
}