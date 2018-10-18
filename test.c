#include "math_4d.h"
#include <stdio.h>


int main()
{
    mat4_t aap;
    mat4_t noot = mat4(
        1, 2, 3, 4,
        0, 5, 0, 0,
        0, 0, 6, 0,
        7, 8, 9, 1
    );

    m4_print(aap);
    printf("\n");
    for (int i = 0; i < 16; i++) {
        printf("%lf ", aap.m[i % 4][i / 4]);
        if(i % 4 == 3) printf("\n");
    }
    m4_print(noot);

    return 0;
}