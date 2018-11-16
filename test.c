#include "math_4d.h"
#include "mt19937.h"
#include <stdio.h>
#include <stdlib.h> // malloc(), free
#include <time.h> // time()
#include <stdbool.h>

inline static double ran(double low, double high)
{
    return (high - low) * dsfmt_genrand() + low;
}

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


    return 0;
}



