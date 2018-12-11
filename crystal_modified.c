#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>

#define MAX_PART_CELL 40
#define MAX_NEIGHBORS 50
#define INVFPI 0.07957747154594766788444188168625718101
#define M_PI 3.14159265358979323846
#define N 8000

//tunable parameters:
double bndLength = 1.55; //distance cutoff for bonds
double bnd_cuttoff = 0.6; //Order to be called a correlated bond
int nbnd_cuttoff = 4; //Number of correlated bonds for a crystalline particle
double obnd_cuttoff = 0.0; //Order to be in the same cluster (0.0 for all touching clusters 0.9 to see defects)

//////////////////////////////////////////////////////////////// my additions
const char labelstring[] = "v1_%02dpf%04.2lfp%04.1lfa%04.2lf";
// e.g. sl10pf0.50p08.0a1.25:
// 10 CubesPerDim, pack_frac 0.50, pressure 8.0, angle 1.25
const char usage_string[] = "usage: program.exe CubesPerDim \
mc_steps packing_fraction BetaP Phi\n";
static int CubesPerDim;
static int mc_steps;
static double packing_fraction;
static double BetaP;
static double Phi;
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

int dontsave = 0;
double bndLengthSq;
int n_part;
char* filename;

int n_cells;
double avr_d, max_d, min_cell;
int nx, ny, nz;

double percentage;
int maxsize;
int numclus;
int particlestocount;

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
    if (dx > box.xhalf)
        dx = dx - box.x;
    else if (dx < -box.xhalf)
        dx = dx + box.x;
    if (dy > box.yhalf)
        dy = dy - box.y;
    else if (dy < -box.yhalf)
        dy = dy + box.y;
    if (dz > box.zhalf)
        dz = dz - box.z;
    else if (dz < -box.zhalf)
        dz = dz + box.z;
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
void order(int l, bndT_t* bnd, compl_t*res1, compl_t*res2)
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
double dotprod(compl_t*vec1, compl_t*vec2, int l)
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
    for (p = particlestocount - 1; p >= 0; p--) {
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
    blist = (blistT_t*)malloc(sizeof(blistT_t) * n_part);
    numconn = (int*)malloc(sizeof(int) * n_part);
    for (p = n_part - 1; p >= 0; p--) {
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
    double cellsize = max_d;
    if (cellsize < bndLength)
        cellsize = bndLength;
    nx = floor(box.x / cellsize);
    ny = floor(box.y / cellsize);
    nz = floor(box.z / cellsize);
    min_cell = box.x / nx;
    if (box.y / ny < min_cell)
        min_cell = box.y / ny;
    if (box.z / nz < min_cell)
        min_cell = box.z / nz; //set min_cell
    printf("cells: %d, %d, %d (%lf)\n", nx, ny, nz, cellsize);
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
    for (i = 0; i < particlestocount; i++) {
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
 *             LOAD
 * Read in coordinate file
 ***********************************************/
int load(tpart_t** pp, int* np, int close)
{
    static FILE* fp = NULL;
    static double rstart = 0;
    particlestocount = 0;
    //  double xa,xb,ya,yb,za,zb;
    if (close && fp) {
        fclose(fp);
        return -1;
    }
    int n;
    char str[256];
    tpart_t* p;
    if (fp == NULL)
        fp = fopen(filename, "r");
    if (fp) {
        getline_number(str, 256, fp);
        if (sscanf(str, "%i\n", np) != 1) {
            if (sscanf(str, "&%i\n", np) != 1)
                return -1; //// new snapshot
        }

        getline_number(str, 256, fp);
        /*    sscanf(str,"%le %le\n",&xa,&xb);
    getline_number(str,256,fp);
    sscanf(str,"%le %le\n",&ya,&yb);
    getline_number(str,256,fp);
    sscanf(str,"%le %le\n",&za,&zb);*/

        int lengths = sscanf(str, "%lf %lf %lf\n", &xsize, &ysize, &zsize);
        if (lengths != 3) {
            printf("Problem! Boxsize unknown.\n");
            fclose(fp);
            exit(666);
        }
        init_box(xsize, ysize, zsize);

        p = (tpart_t*)malloc(*np * sizeof(tpart_t));

        double t1, t2, t3, t4;
        // int ti;
        avr_d = 0.0;
        max_d = 0.0;
        char type;
        for (n = 0; n < *np; n++) {
            getline_number(str, 256, fp);
            int r;
            r = sscanf(str, "%c %le %le %le %le %[^\n]s\n", &type, &t1, &t2, &t3, &t4, p[n].str);

            if (r >= 4) {
                p[n].x = t1;
                p[n].y = t2;
                p[n].z = t3;
                p[n].c = type;
            } else {
                printf("Error reading coords.in\n");
                exit(666);
            }
            if (r >= 5) {
                if (rstart == 0) {
                    particlestocount++;
                    rstart = t4;
                } else if (t4 == rstart)
                    particlestocount++;
                else
                    rstart = -1;
                p[n].d = t4;
                if (rstart > 0)
                    avr_d += p[n].d;
            } else {
                p[n].d = 1.0;
                particlestocount = *np;
                avr_d++;
            }
            if (max_d < p[n].d)
                max_d = p[n].d;
            //      p[n].x-=(xa+xb)*0.5;
            //      p[n].y-=(ya+yb)*0.5;
            //      p[n].z-=(za+zb)*0.5;
            //      if(p[n].x > xb-xa || p[n].x < xa-xb) printf("ERROR: Particle outside box\n");
            //      if(p[n].y > yb-ya || p[n].y < ya-yb) printf("ERROR: Particle outside box\n");
            //      if(p[n].z > zb-za || p[n].z < za-zb) printf("ERROR: Particle outside box\n");
            p[n].r = p[n].d * 0.5;
            /*    for(n=0;n<n_part;n++){
      p[n].x -= (xb+xa)/2.0;
      p[n].y -= (yb+ya)/2.0;
      p[n].z -= (zb+za)/2.0;
    }*/
        }

        avr_d /= (double)particlestocount;
        *pp = p;
    } else {
        printf("Error coord file not found!\n");
        exit(666);
    }
    //  printf ("Particles: %d (%lf, %lf)\n", particlestocount, avr_d, max_d);
    //  particlestocount = *np;
    return 1; ////
}

/************************************************
 *             LOAD_PARTICLES
 * Read file and initialize stuff
 ***********************************************/
int load_particles()
{
    if (load(&part, &n_part, 0) == -1) //// if new snapshot
        return -1; //// new snapshot 
    //  bndLength=1.4*avr_d;
    bndLengthSq = bndLength * bndLength;
    init_cells();
    init_nblist();
    return 0; //// not end of file, continue
}

/************************************************
 *             CALC_ORDER
 * Calculate q_4 for all particles
 ***********************************************/
compl_t*calc_order(void)
{
    int i, j, m;
    compl_t*q1, *q2;
    const int l = 4;
    double temp;
    compl_t*orderp = (compl_t*)malloc(sizeof(compl_t) * n_part * (l * 2 + 1));
    memset(orderp, (int)0.0, sizeof(compl_t) * n_part * (l * 2 + 1));
    for (i = 0; i < particlestocount; i++) {
        q1 = (orderp + i * (2 * l + 1) + l);
        for (j = 0; j < blist[i].n; j++) {
            if (blist[i].bnd[j].n > i) {
                q2 = (orderp + blist[i].bnd[j].n * (2 * l + 1) + l);
                order(l, &(blist[i].bnd[j]), q1, q2);
            }
        }
    }
    //normalize vector
    for (i = 0; i < particlestocount; i++) {
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
int* calc_conn(compl_t*orderp) //calculates "connected" particles
{
    int i, j, z, np = 0;
    const int l = 4;
    int* conn = malloc(sizeof(int) * n_part);
    for (i = 0; i < particlestocount; i++) {
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
 *             SAVE_CLUSS
 * Output snapshot in which clusters are colored
 * Fluid particles printed small
 ***********************************************/
void save_cluss(int* cluss, int* size, int big, int nn)
{
    int i;
    static int first = 1;

    char fn[100] = "clus";
    strcat(fn, filename);
    char dfn[100] = "nclus";
    strcat(dfn, filename);
    FILE* file;
    FILE* datafile;
    if (first) {
        file = fopen(fn, "w");
        datafile = fopen(dfn, "w");
        //printf ("Writing to %s\n", fn);
        first = 0;
    } else {
        file = fopen(fn, "a");
        datafile = fopen(dfn, "a");
    }
    fprintf(file, "&%i \n", n_part);
    //   printf("%f %f\n",-box.xhalf,box.xhalf);
    //   printf("%f %f\n",-box.yhalf,box.yhalf);
    //   printf("%f %f\n",-box.zhalf,box.zhalf);
    fprintf(file, "%.12lf %.12lf %.12lf\n", xsize, ysize, zsize);

    int* sorta = malloc(sizeof(int) * nn);
    int* rank = malloc(sizeof(int) * nn);
    for (i = 0; i < nn; i++)
        sorta[i] = i;

    int cmpr(const void* a, const void* b)
    {
        return -size[*((int*)a)] + size[*((int*)b)];
    }

    qsort(sorta, nn, sizeof(int), &cmpr);

    for (i = 0; i < nn; i++) {
        int bha;
        for (bha = 0; sorta[bha] != i; bha++)
            ;
        rank[i] = bha;
    }

    if (dontsave == 0) {
        for (i = 0; i < n_part; i++) {
            if (cluss[i] && size[cluss[i]] > 2) {
                int rnk = rank[cluss[i]];
                if (rnk > 0)
                    rnk = ((rnk - 1) % 25) + 1;
                char ch = 'a' + rnk;
                //          char ch = 'a' + (numconn[i]);
                fprintf(file, "%c %lf %lf %lf %lf ", ch, part[i].x, part[i].y, part[i].z, part[i].d);
                fprintf(file, "%s\n", part[i].str);
            } else {
                char ch = 'a';
                //          char ch = 'a' + (numconn[i]);
                fprintf(file, "%c %lf %lf %lf %lf ", ch, part[i].x, part[i].y, part[i].z, part[i].d / 10);
                fprintf(file, "%s\n", part[i].str);
            }
        }
    }
    fprintf(datafile, "%lf  %d  %d\n", percentage, numclus, maxsize);
    fclose(file);
    fclose(datafile);
    free(sorta);
    free(rank);
}

/************************************************
 *             CALC_CLUSTERS
 * Find clusters of bonded particles
 ***********************************************/
void calc_clusters(int* conn, compl_t*orderp)
{
    int* cluss = malloc(sizeof(int) * n_part);
    int* size = malloc(sizeof(int) * n_part);
    int i, cs, cn = 1, big = 0;
    // int bc = -1;
    const int l = 4;
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

    for (i = 0; i < n_part; i++) {
        cluss[i] = 0;
        size[i] = 0;
    }

    for (i = 0; i < particlestocount; i++) {
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

    //calculate average cluster size
    int tcs = 0;
    for (i = 0; i < cn; i++) {
        tcs += size[i];
    }

    save_cluss(cluss, size, big, cn);
    percentage = tcs / (double)n_part;
    maxsize = big;
    numclus = cn - 1;
    printf("% i clusters, Max size: %i, Percentage of crystalline particles %f\n", numclus, big, percentage);

    free(cluss);
    free(size);
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
 *             MAIN
 ***********************************************/
int main(int argc, char** argv)
{
    int f = 0, i = 0;
    compl_t*order = NULL;
    int* connections = NULL;

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
    strcat(outputfile_name, "/order.txt"); // TODO: put in separate directory in main?
    FILE* fp_order = fopen(outputfile_name, "w");

    strcat(buffer, "/coords_step%07d.poly");

    for (int step = 0; step <= mc_steps; step += 100) {
        char datafile_name[128] = "";
        // replace all %d, %lf in buffer with values and put in datafile_name
        sprintf(datafile_name, buffer, CubesPerDim, packing_fraction, BetaP, Phi, step);
        printf("datafile_name: %s\n", datafile_name);

        read_data2(datafile_name, 0);
        if (is_overlap()) {
            printf("\n\n\tWARNING\n\n\nThe read file contains overlap.\nExiting.\n");
            exit(4);
        }

        // compl_t*order = calc_order();
        fprintf(fp_order, "%lf ", 1.0);

        // free(order);
    }

    fclose(fp_order);

    // while (load_particles() != -1) {
    //     printf("Frame %d: ", f);
    //     f++;
    //     order = calc_order();

    //     connections = calc_conn(order);
    //     calc_clusters(connections, order);
    //     free(order);
    //     free(connections);
    //     free(part);
    //     free(cells);
    //     free(numconn);
    //     for (i = 0; i < n_part; i++)
    //         free(blist[i].bnd);
    //     free(blist);
    // }
    // if (i == 0)
    //     exit(999);
    // load(NULL, NULL, 1);
    return 0;
}