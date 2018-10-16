#include "math_3d.h"
#include <stdio.h>

static vec3_t aap[3];

vec3_t dikke_swag_functie(vec3_t noot)
{
    vec3_t result = v3_add(noot, vec3(0, 1, 2));
    return result;
}

int main()
{
    aap[0] = vec3(3, 6, -1);
    aap[1] = vec3(3, -3, -1);
    printf("%lf\n", v3_dot(aap[0], aap[1]));
    aap[2] = dikke_swag_functie(vec3(3,2,1));
    printf("%lf %lf %lf\n", aap[2].x,aap[2].y,aap[2].z);

    printf("%d %d %d %d\n", 1 & 2, 2 & 6, 1 | 7, 1 | 4);

    return 0;
}