#include "math_4d.h"
#include <stdio.h>
#include <time.h> // time()
#include <stdlib.h> // malloc(), free


int main(int argc, char* argv[])
{
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
    double aap;
    sscanf(argv[1], "%lf", &aap);
    printf("%lf\n", aap);

    return 0;
}