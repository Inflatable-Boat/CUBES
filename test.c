#include "math_4d.h"
#include <stdio.h>
#include <time.h>    // time()


int main()
{
    /* vec3_t aap = vec3(2, 3, 4);
    float *noot = &aap.x;
    *(noot + 1) = 5;
    for (int i = 0; i < 3; i++) {
        printf("%3.3lf\n", *(noot + i));
    } */

    mat4_t aap = mat4(1,2,3,0,-1,2,-3,0,5,4,3,0,0,0,0,0);
    vec3_t noot = vec3(1, -1, 2);
    mat4_t mies = m4_rotation(0.1, noot);

    m4_print(aap); printf("\n");

    aap = m4_mul(mies, aap);

    m4_print(aap);printf("\n");

    aap = m4_mul(m4_invert_affine(mies), aap);

    m4_print(aap);

    return 0;
}