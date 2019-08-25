#ifndef RAY_H_
#define RAY_H_
#include "vec3.h"

class ray
{
public:
    ray() {}
    ray(const vec3 &org, const vec3 &dir, float time = 0.0) :
        _org(org),
        _dir(unit_vector(dir)),
        _time(time) {}

    vec3  origin()    const { return _org; }
    vec3  direction() const { return _dir; }
    float time()      const { return _time;}

    vec3 point_at_parameter(const float t) const { return _org + t * _dir; }

private:
    vec3  _org;
    vec3  _dir;
    float _time;
};

#endif