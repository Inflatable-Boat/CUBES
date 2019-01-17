#include "math_3d.h"
#include "mt19937.h"
#include <stdio.h>
#include <stdlib.h> // malloc(), free
#include <time.h> // time()
#include <stdbool.h>
#include "test2.c"

inline static double ran(double low, double high)
{
    return (high - low) * dsfmt_genrand() + low;
}

/* void write(char datafolder_name[128])
{
    char buffer[128];
    strcpy(buffer, datafolder_name);

    strcat(buffer, "/coords_step%07d.poly");
    char output_file[128];
    sprintf(output_file, buffer, 100); // replace %07d with step and put in output_file.

    printf("output_file = %s\n", output_file);

    char datafile[128] = "";
    strcpy(datafile, datafolder_name);
    strcat(datafile, "/coords_step%07d.poly");
    sprintf(datafile, datafile, 100); // replace %07d with step and put in output_file.

    printf("output_file = %s\n", output_file);
}

const char labelstring[] = "v1_%02dpf%04.2lfp%04.1lfa%04.2lf"; */


int main(int argc, char* argv[])
{
    dsfmt_seed(time(NULL));
    /* vec3_t aap = vec3(2, 3, 4);
    float *noot = &aap.x;
    *(noot + 1) = 5;
    for (int i = 0; i < 3; i++) {
        printf("%3.3lf\n", *(noot + i));
    } */

    /*     mat4_t aap = mat4(1,2,3,0,-1,2,-3,0,5,4,3,0,0,0,0,0);
    vec3_t noot = vec3(1, -1, 2);
    mat4_t mies = m4_rotation(0.1, noot);

    m4_print(aap); printf("\n");

    aap = m4_mul(mies, aap);

    m4_print(aap);printf("\n");

    aap = m4_mul(m4_invert_affine(mies), aap);

    m4_print(aap); */

    /*     vec3_t* aap = malloc(10 * sizeof(vec3_t));
    if(!aap){
        printf("mem alloc failed");
        return 1;
    }
    printf("mem alloc succeeded\n");
    for (int i = 0; i < 10; i++) {
        *(aap + i) = vec3(i, 10 - i, (i - 3) * (i - 3));
    }
    printf("init succeeded\n");

    float* pv3;
    for (int i = 0; i < 10; i++) {
        pv3 = &(aap + i)->x; // == &((aap + i)->x);
        for (int d = 0; d < 3; d++) {
            printf("%f ", *(pv3 + d));
        }
        printf("\n");
    }

    printf("size = %ld\n", sizeof(vec3_t));
    printf("aap     = %ld\naap + 1 = %ld\n", aap, aap + 1);

    free(aap); */
/*     int aap = 0;
    double noot = 0.1 / aap;
    if(noot < 0.4)
        printf("inf < 0.4\n");
    if(noot > 0.6)
        printf("inf > 0.6\n"); */
/*     for(int random = 0; random < 6; random++)
        if(random & 1)
            printf("%d\n", random); //1 3 5*/

    /* int aap = 2;
    int noot = 2;
    int mies = 3;

    if(aap == noot == mies)
        printf("ja");
    if(aap == noot == mies - 1)
        printf("nee"); // doesn't print because ((2 == 2) == true) != 2?? */
    /* for (int i = -3; i < 3; i++) {
        printf("%d\n", ((i + 3) % 3));
    } */

    /* const int N = 38;
    int aap[N][N][N][N]; // N = 38 is the max for which it works, 38^4 < 2^21
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                for (int l = 0; l < N; l++) {
                    aap[i][j][k][l] = i * N * N * N + j * N * N + k * N + l;
                }
            }
        }
    }
    for (int i = 0; i < N * N; i++) {
        printf("%d ", aap[0][0][i/N][i%N]);
        if((i + 1 % N) != 0)
            printf("\n");
    } */

    /* while(true)
    { // geen fucking warning ech nie????
        printf("ech nie\n");
    } */

    /* int aap = 1;
    printf("%d %d   ", aap++, aap++); // 2 1
    printf("%d\n", aap); // 3 */

    /* int aap[3];
    printf("%d", aap[3000]); // segmentation fault (core dumped) */

    /* FILE* fp = fopen("read_file_test.txt", "r");
    int a=0,b=0;
    double e=0,f=0;
    int garbage;
    fscanf(fp, "%d %d %d %d", &a, &b, &garbage, &garbage);
    fscanf(fp, "%lf %lf %lf %lf", &e, &f, &garbage, &garbage);
    printf("%d %d %lf %lf\n", a, b, e, f);
    fclose(fp); */

    // printf("%d\n", (int)(14.015651999999999/14.015652));
    /* const char string[] = "this is a \
        multiline string, does\
                an enter or tab  or two appear?\
            hello\n"; // no enter, yes tab
    printf(string); */

    /* if (strcmp(argv[1], "3") == 0) {
        printf("aap\n");
    } */
    
    /* char string[64] = "";
    strcat(string, "dit en zo: %02d %4.2lf\n"); // dit gaat helemaal naar de tering
    sprintf(string, string, 8, 612321423234315675675676523.3324328);
    printf("%s", string); */

    // char buffer[128] = "datafolder/";
    // strcat(buffer, labelstring);
    // // char datafolder_name[128] = "";
    // // replace %4.1lf with packing_fraction and BetaP and Phi
    // printf("%s\n", buffer);
    // sprintf(buffer, buffer, // undefined behaviour, sometimes works, sometimes not
    // /* CubesPerDim */ 8,
    // /* packing_fraction */ 0.5,
    // /* BetaP */ 12.,
    // /* Phi */ 1.25);

    // printf("%s\n", buffer);
    // write(buffer);
    
    /* int aap = 1, noot = 0, mies = -1;
    printf("%d, %d, %d becomes\n%d, %d, %d\n", aap, noot, mies, ~aap, ~noot, ~mies);
    for (int i = 31; i >= 0; i--) {
        if (mies & (1 << i)) printf("1");
        else printf("0");
    }
    printf("\n"); */

    /* char aap[32] = "aap";
    char noot[32] = "noot baby";
    strcpy(noot, aap);
    printf("noot now contains %s\n", noot); // prints aap
    printf("looping over noot, though:\n");
    for (int i = 0; i < 32; i++) {
        printf("%c", noot[i]); // prints aap baby
    }
    printf("\n"); */

    /* char aap[32];
    double noot[32];
    printf("%ld, %ld\n", sizeof aap, sizeof noot); // prints 32, 256 */

    /* double noot = 0, mies = 0;
    int aap = sscanf("a 3.0 4.0 5.0\n", "%lf %lf", &noot, &mies);
    printf("read %d arguments, namely:\n1st: %lf,\n2nd: %lf\n", aap, noot, mies); */
    // prints 0 arguments, namely 0 and 0.

    /* double aap = 0.1;
    printf("aap is %lf\n", aap); // 0.1
    printf("aap is now %lf\n", (aap+=0.1) - 0.1); // 0.1
    printf("aap is %lf\n", aap); // 0.2 */
    
    /* aap = 3;
    add_one_to_aap();
    printf("%d\n", aap); */


    return 0;
}



