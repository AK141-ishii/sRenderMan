#ifndef PDF_H_
#define PDF_H_

#include "vec3.h"
#include "onb.h"
#include "rndm.h"

class pdf {
public:
    virtual float value (const vec3 & direction ) const = 0;
    virtual vec3 generate() const = 0 ;
};

class cosine_pdf : public pdf {
public:
    cosine_pdf(const vec3 &w) { uvw.build_from_w(w); }
    virtual float value(const vec3 &direction) const {
        float cosine = dot(unit_vector(direction), uvw.w());
        if (cosine > 0)
            return cosine / 3.1415926535;
        else
            return 0.0;
    }
    virtual vec3 generate() const {
        return uvw.local(random_cosine_direction());
    }
    onb uvw;
};

class hitable_pdf : public pdf {
public:
    hitable_pdf(hitable *p_, const vec3& origin_) : hit_ptr(p_), origin(origin_) {};
    virtual float value (const vec3 & direction ) const {
        return hit_ptr->pdf_value(origin, direction);
    }
    virtual vec3 generate() const {
        return hit_ptr->random(origin);
    }

    vec3 origin;
    hitable *hit_ptr;
};

class mixture_pdf : public pdf {
public:
    mixture_pdf(pdf *p0, pdf *p1) { p[0] = p0; p[1] = p1; }
    virtual float value(const vec3 &direction) const {
        return 0.5 * p[0]->value(direction) + 0.5 * p[1]->value(direction);
    }
    virtual vec3 generate() const {
        if (drand() < 0.5) {
            return p[0]->generate();
        }
        else
            return p[1]->generate();
    }
    pdf *p[2];
};

#endif