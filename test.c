#include <stdio.h>
#include "math_3d.h"

static vec3_t aap[3];

int main(){
    aap[2].x = 5;
    for(int i = 0; i<3; i++){
        printf("%lf ", aap[i].x);
    }

    return 0;
}