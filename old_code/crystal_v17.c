#include "math_3d.h"
#include <malloc.h>
// #include <math.h>
// #include <stdbool.h>
// #include <stdio.h>
#include <stdlib.h>
// #include <string.h>
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
#define N 8000
#define NDIM 3

//tunable parameters:
double bndLength = 1.55; //distance cutoff for bonds
double bnd_cuttoff = 0.6; //Order to be called a correlated bond
int nbnd_cuttoff = 4; //Number of correlated bonds for a crystalline particle
double obnd_cuttoff = 0.0; //Order to be in the same cluster (0.0 for all touching clusters 0.9 to see defects)

//////////////////////////////////////////////////////////////// my additions
// extern const char labelstring[];
// extern int CubesPerDim;
// extern int mc_steps;
// extern double Phi;
extern int n_particles;
// extern double box[NDIM];
// extern vec3_t rr[N]; // position of center of cube
// extern mat4_t mm[N]; // rotation matrix of cube
// extern double SinPhi;
extern double CosPhi;
// extern double ParticleVolume;
extern system_t sim;
//////////////////////////////////////////////////////////////// end of my additions

typedef struct {
    double re;
    double im;
} compl_t;
//typedef struct {double nx; double ny; double nz; double si; double co ;int n;}bndT_t;
typedef struct {
    double nz;
    double si;
    double co;
    int n;
} bndT_t;
typedef struct {
    int n;
    bndT_t* bnd;
} blistT_t;
// typedef struct {
//     int n;
//     int* nb;
// } nlistT;
blistT_t* blist;
int* numconn;
double xsize, ysize, zsize;

// int dontsave = 0;
double bndLengthSq;
// int n_part;
char* filename;
// char*

int n_cells;
double avr_d, max_d, min_cell;
int nx, ny, nz;

double percentage;
int maxsize;
int numclus;
// int particlestocount;

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
tbox_t tbox;

typedef struct {
    double x;
    double y;
    double z;
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
 *             DIST_PART
 * Distance between two particles
 ***********************************************/
double dist_part(int i, int j)
{
    double dx = part[i].x - part[j].x;
    double dy = part[i].y - part[j].y;
    double dz = part[i].z - part[j].z;
    if (dx > tbox.xhalf)
        dx = dx - tbox.x;
    else if (dx < -tbox.xhalf)
        dx = dx + tbox.x;
    if (dy > tbox.yhalf)
        dy = dy - tbox.y;
    else if (dy < -tbox.yhalf)
        dy = dy + tbox.y;
    if (dz > tbox.zhalf)
        dz = dz - tbox.z;
    else if (dz < -tbox.zhalf)
        dz = dz + tbox.z;
    return dx * dx + dy * dy + dz * dz;
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
 *             ORDER
 * Calculate q for a pair of particles
 ***********************************************/
void order(int l, bndT_t* bnd, compl_t* res1, compl_t* res2)
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
    int i = ((int)(nx * (x + tbox.xhalf) * tbox.oneoverx) % nx) * ny * nz + ((int)(ny * (y + tbox.yhalf) * tbox.oneovery) % ny) * nz + ((int)(nz * (z + tbox.zhalf) * tbox.oneoverz) % nz);
    return i;
}

/************************************************
 *             UPDATE_NBLISTP
 * Find neighbors of a particle
 ***********************************************/
void update_nblistp(int p)
{
    int i, j, c, cellp, id;
    bndT_t* bnd;
    double d, dxy, dx, dy, dz;
    cellp = part[p].cell;
    blist[p].n = 0;
    for (i = 26; i >= 0; i--) {
        c = cells[cellp].buren[i];
        for (j = cells[c].n - 1; j >= 0; j--) {
            id = cells[c].particles[j];
            if (id != p) {
                dx = part[p].x - part[id].x;
                dy = part[p].y - part[id].y;
                dz = part[p].z - part[id].z;

                if (dx > tbox.xhalf)
                    dx -= tbox.x;
                else if (dx < -tbox.xhalf)
                    dx += tbox.x;
                if (dy > tbox.yhalf)
                    dy -= tbox.y;
                else if (dy < -tbox.yhalf)
                    dy += tbox.y;
                if (dz > tbox.zhalf)
                    dz -= tbox.z;
                else if (dz < -tbox.zhalf)
                    dz += tbox.z;

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
void init_nblist(void)
{
    int p;
    blist = (blistT_t*)malloc(sizeof(blistT_t) * n_particles);
    numconn = (int*)malloc(sizeof(int) * n_particles);
    for (p = n_particles - 1; p >= 0; p--) {
        blist[p].bnd = (bndT_t*)malloc(sizeof(bndT_t) * MAX_NEIGHBORS);
        numconn[p] = 0;
        update_nblistp(p);
    }
}

/************************************************
 *             INIT_CELLS
 * Initialize the cell list
 ***********************************************/
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
    nx = floor(tbox.x / cellsize);
    ny = floor(tbox.y / cellsize);
    nz = floor(tbox.z / cellsize);
    min_cell = tbox.x / nx;
    if (tbox.y / ny < min_cell)
        min_cell = tbox.y / ny;
    if (tbox.z / nz < min_cell)
        min_cell = tbox.z / nz; //set min_cell
    // printf("cells: %d, %d, %d (%lf)\n", nx, ny, nz, cellsize);
    n_cells = nx * ny * nz;
    // printf ("Cells: %d (%lf)\n", n_cells, cellsize);
    cells = (tcells_t*)malloc(n_cells * sizeof(tcells_t));
    for (x = 0; x < nx; x++)
        for (y = 0; y < ny; y++)
            for (z = 0; z < nz; z++) {
                i = x * ny * nz + y * nz + z;
                cells[i].x1 = x * tbox.x / nx - tbox.xhalf;
                cells[i].y1 = y * tbox.y / ny - tbox.yhalf;
                cells[i].z1 = z * tbox.z / nz - tbox.zhalf;
                cells[i].x2 = (x + 1) * tbox.x / nx - tbox.xhalf;
                cells[i].y2 = (y + 1) * tbox.y / ny - tbox.yhalf;
                cells[i].z2 = (z + 1) * tbox.z / nz - tbox.zhalf;
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
        j = coords2cell(part[i].x, part[i].y, part[i].z);
        if (cells[j].n + 1 >= MAX_PART_CELL) {
            printf("ERROR: Too many particles in a cell!\n");
            exit(666);
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
 * Initialize tbox parameters
 ***********************************************/
void init_box(double x_box, double y_box, double z_box)
{
    tbox.x = x_box;
    tbox.y = y_box;
    tbox.z = z_box;
    tbox.min = tbox.x;
    if (tbox.min > tbox.y)
        tbox.min = tbox.y;
    if (tbox.min > tbox.z)
        tbox.min = tbox.z;
    tbox.xhalf = tbox.x * 0.5;
    tbox.yhalf = tbox.y * 0.5;
    tbox.zhalf = tbox.z * 0.5;
    tbox.oneoverx = 1.0 / tbox.x;
    tbox.oneovery = 1.0 / tbox.y;
    tbox.oneoverz = 1.0 / tbox.z;
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
    for (i = 0; i < n_particles; i++) {
        q1 = (orderp + i * (2 * l + 1) + l);
        for (j = 0; j < blist[i].n; j++) {
            if (blist[i].bnd[j].n > i) {
                q2 = (orderp + blist[i].bnd[j].n * (2 * l + 1) + l);
                order(l, &(blist[i].bnd[j]), q1, q2);
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
                //        printf ("Test %d, %d\n", i, n_particles);
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
 *             CALC_CLUSTERS
 * Find clusters of bonded particles
 ***********************************************/
void calc_clusters(int* conn, compl_t* orderp)
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
    // int tcs = 0;
    // for (i = 0; i < cn; i++) {
    //     tcs += size[i];
    // }

    // percentage = tcs / (double)n_particles;
    maxsize = big;
    // numclus = cn - 1;
    // printf("%i clusters, Max size: %i, %% crystalline particles: %f\n", numclus, big, percentage);

    free(cluss);
    free(size);
}

/************************************************
 *             MAIN
 ***********************************************/
/// returns the size of the larges cluster
int crystal(void)
{
    // int f = 0, i = 0;
    compl_t* order = NULL;
    int* connections = NULL;

    init_box(sim.box[0], sim.box[1], sim.box[2]);
    part = (tpart_t*)malloc(n_particles * sizeof(tpart_t));

    for (int n = 0; n < n_particles; n++) {
        part[n].x = sim.r[n].x;
        part[n].y = sim.r[n].y;
        part[n].z = sim.r[n].z;
        part[n].c = 'a';
        part[n].d = 1.0;
    }

    bndLengthSq = bndLength * bndLength;
    init_cells();
    init_nblist();
    // end load particles

    order = calc_order();
    connections = calc_conn(order);
    calc_clusters(connections, order);
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