#include "math_4d.h"
#include <stdio.h>
#include <time.h>    // time()


int main()
{
    vec3_t aap = vec3(2, 3, 4);
    float *noot = &aap.x;
    *(noot + 1) = 5;
    for (int i = 0; i < 3; i++) {
        printf("%3.3lf\n", *(noot + i));
    }

    return 0;
}