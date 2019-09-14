#ifndef RNDM_H_
#define RNDM_H_
#include "vec3.h"

vec3 random_in_unit_sphere() {
    vec3 p;
    do
    {
        p = 2.0 * vec3((float)rand() / (float)RAND_MAX, (float)rand() / (float)RAND_MAX, (float)rand() / (float)RAND_MAX) 
        - vec3(1, 1, 1);
    } while (dot(p,p) >= 1.0);
    return p;
}

inline vec3 random_cosine_direction() {
    float r1 = drand();
    float r2 = drand();
    float z = sqrt(1 - r2);
    float phi = 2 * 3.1415926535 * r1;
    float x = cos(phi) * 2 * sqrt(r2);
    float y = sin(phi) * 2 * sqrt(r2);
    return vec3(x, y, z);
}

#endif