#ifndef PRIMITIVE_H_
#define PRIMITIVE_H_

#include "vec3.h"
#include "material.h"
#include "hitable.h"
#include "hitable_list.h"
#include "aabb.h"
#include <float.h>
#include <math.h>

/* SPHERE */

void get_sphere_uv(const vec3 &p, float &u, float &v) {
    float phi = atan2(p.z(), p.x());
    float theta = asin(p.y());
    const float pi = 3.1415926535;
    u = 1 - (phi + pi) / (2 * pi);
    v = (theta + pi / 2) / pi;
}

class sphere : public hitable {
public:
    sphere(){}
    sphere(vec3 cen, float r, material *mat) : center(cen), radius(r), mat_ptr(mat) {};
    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec)const;
    virtual bool bounding_box(float t0, float t1, aabb& box) const ;
    vec3 center;
    float radius;
    material *mat_ptr;
};

bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
    // D>=0 @ (r.ori + t * r.dir - center)^2 = R^2
    vec3 oc = r.origin() - center;
    float a = dot(r.direction(), r.direction());
    float b = dot(oc, r.direction());
    float c = dot(oc,oc) - radius*radius;
    float discriminant = b*b - a*c;
    if(discriminant > 0){
        // smaller t : near to origin
        float temp = (-b - sqrt(b*b-a*c))/a;
        if(temp < t_max && temp > t_min){
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            get_sphere_uv((rec.p - center) / radius, rec.u, rec.v);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;
            return true;
        }
        // larger t : far from origin
        temp = (-b + sqrt(b*b - a*c))/a;
        if(temp < t_max && temp > t_min){
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            get_sphere_uv((rec.p - center) / radius, rec.u, rec.v);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
}

bool sphere::bounding_box(float t0, float t1, aabb& box) const {
    box = aabb(center-vec3(radius, radius, radius),
            center + vec3(radius, radius, radius));
    return true;
}


class moving_sphere: public hitable {
    public:
    moving_sphere() {}
    moving_sphere(vec3 cen0, vec3 cen1, float t0, float t1, float r, material *m)
    : center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), mat_ptr(m) {};
    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec)const;
    virtual bool bounding_box(float t0, float t1, aabb& box) const ;
    vec3 center(float time) const;
    vec3 center0, center1;
    float time0, time1;
    float radius;
    material *mat_ptr;
};

vec3 moving_sphere::center(float time) const{
    return center0 + ((time - time0) / (time1 - time0))*(center1 - center0);
}

bool moving_sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec)const{
    // D>=0 @ (r.ori + t * r.dir - center)^2 = R^2
    vec3 oc = r.origin() - center(r.time());
    float a = dot(r.direction(), r.direction());
    float b = dot(oc, r.direction());
    float c = dot(oc,oc) - radius*radius;
    float discriminant = b*b - a*c;
    if(discriminant > 0){
        // smaller t : near to origin
        float temp = (-b - sqrt(b*b-a*c))/a;
        if(temp < t_max && temp > t_min){
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center(r.time())) / radius;
            rec.mat_ptr = mat_ptr;
            return true;
        }
        // larget t : far from origin
        temp = (-b + sqrt(b*b - a*c))/a;
        if(temp < t_max && temp > t_min){
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center(r.time())) / radius;
            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
}

bool moving_sphere::bounding_box(float t0, float t1, aabb& box)const {
    aabb box0(center(t0)- vec3(radius, radius, radius), 
            center(t0) + vec3(radius, radius, radius));
    aabb box1(center(t1) - vec3(radius, radius, radius),
            center(t1)+vec3(radius, radius, radius));
    box = surrounding_box(box0, box1);
    return true;
}

/* RECT */

class xy_rect: public hitable {
    public:
        xy_rect() { };
        xy_rect(float _x0, float _x1, float _y0, float _y1, float _k, material *mat)
           : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat) { };
        virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
        virtual bool bounding_box(float t0, float t1, aabb& box) const {
            box = aabb(vec3(x0, y0, k-0.0001), vec3(x1, y1, k+0.0001));
            return true;
        }
        material *mp;
        float x0, x1, y0, y1, k;
};

bool xy_rect::hit(const ray& r, float t0, float t1, hit_record& rec) const {
    float t = (k - r.origin().z()) / r.direction().z();
    if (t < t0 || t > t1) return false;
    float x = r.origin().x() + t * r.direction().x();
    float y = r.origin().y() + t * r.direction().y();
    if (x < x0 || x > x1 || y < y0 || y > y1) return false;
    rec.u = (x-x0) / (x1-x0);
    rec.v = (y-y0) / (y1-y0);
    rec.t = t;
    rec.mat_ptr = mp;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(0,0,1);
    return true;
}

class yz_rect: public hitable {
    public:
        yz_rect() { };
        yz_rect(float _y0, float _y1, float _z0, float _z1, float _k, material *mat)
            : y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat){};
        virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
        virtual bool bounding_box(float t0, float t1, aabb& box) const {
            box = aabb(vec3(k-0.0001, y0, z0), vec3(k+0.0001, y1, z1));
            return true;
        }
        material *mp;
        float y0, y1, z0, z1, k;
};

bool yz_rect::hit(const ray &r, float t0, float t1, hit_record &rec) const {
    float t = (k-r.origin().x()) / r.direction().x();
    if (t < t0 || t > t1) return false;
    float y = r.origin().y() + t * r.direction().y();
    float z = r.origin().z() + t * r.direction().z();
    if (y < y0 || y > y1 || z < z0 || z > z1) return false;
    rec.u = (y-y0) / (y1-y0);
    rec.v = (z-z0) / (z1-z0);
    rec.t = t;
    rec.mat_ptr = mp;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(1,0,0);
    return true;
}

class zx_rect: public hitable {
    public:
        zx_rect() { };
        zx_rect(float _z0, float _z1, float _x0, float _x1, float _k, material *mat)
            : z0(_z0), z1(_z1), x0(_x0), x1(_x1), k(_k), mp(mat) { };
        virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
        virtual bool bounding_box(float t0, float t1, aabb& box) const {
            box = aabb(vec3(x0, k - 0.0001, z0), vec3(x1, k + 0.0001, z1));
            return true;
        }
        virtual float pdf_value(const vec3 &origin, const vec3 &v) const {
            hit_record rec;
            if(this->hit(ray(origin, v), 0.0001, FLT_MAX, rec)) {
                float area = (x1 - x0) * (z1 - z0);
                float distance_squared = rec.t * rec.t * v.squared_length();
                float cosine = fabs(dot(v, rec.normal) / v.length());
                return distance_squared / (cosine * area);
            }
            else return 0.0;
        }
        virtual vec3 random(const vec3 &origin) {
            vec3 random_point = vec3(x0 + drand() * (x1 - x0), k, z0 + drand() * (z1 - z0));
            return random_point - origin;
        }
        material *mp;
        float z0, z1, x0, x1, k;
};

bool zx_rect::hit(const ray&r, float t0, float t1, hit_record& rec) const {
    float t = (k-r.origin().y()) / r.direction().y();
    if (t < t0 || t > t1) return false;
    float z = r.origin().z() + t * r.direction().z();
    float x = r.origin().x() + t * r.direction().x();
    if (x < x0 || x > x1 || z < z0 || z > z1) return false;
    rec.u = (z-z0) / (z1-z0);
    rec.v = (x-x0) / (x1-x0);
    rec.t = t;
    rec.mat_ptr = mp;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(0,1,0);
    return true;
}


/* BOX */
class box: public hitable {
public:
    box() {}
    box(const vec3 &p0, const vec3 &p1, material *ptr);
    virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
    virtual bool bounding_box(float t0, float t1, aabb& box) const {
        box = aabb(pmin, pmax);
        return true;
    }
    vec3 pmin, pmax;
    hitable *list_ptr;
};

box::box(const vec3& p0, const vec3& p1, material *ptr) {
    pmin = p0;
    pmax = p1;
    hitable **list = new hitable*[6];
    list[0] = new xy_rect(pmin.x(), pmax.x(), pmin.y(), pmax.y(), pmax.z(), ptr);
    list[1] = new flip_normals(new xy_rect(pmin.x(), pmax.x(), pmin.y(), pmax.y(), pmin.z(), ptr));
    list[2] = new yz_rect(pmin.y(), pmax.y(), pmin.z(), pmax.z(), pmax.x(), ptr);
    list[3] = new flip_normals(new yz_rect(pmin.y(), pmax.y(), pmin.z(), pmax.z(), pmin.x(), ptr));
    list[4] = new zx_rect(pmin.z(),pmax.z(), pmin.x(), pmax.x(), pmax.y(), ptr);
    list[5] = new flip_normals(new zx_rect(pmin.z(),pmax.z(), pmin.x(), pmax.x(), pmin.y(), ptr));
    list_ptr = new hitable_list(list, 6);
}

bool box::hit(const ray& r, float t0, float t1, hit_record& rec) const {
    return list_ptr->hit(r, t0, t1, rec);
}


/* VOLUME */

class constant_medium : public hitable {
    public:
        constant_medium(hitable *b, float d, texture *a):boundary(b), density(d) {
            phase_function = new isotropic(a);
        }
        virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
        virtual bool bounding_box(float t0, float t1, aabb& box) const {
            return boundary->bounding_box(t0, t1, box);
        }
        hitable *boundary;
        float density;
        material *phase_function;
};

bool constant_medium::hit(const ray &r, float t_min, float t_max, hit_record &rec) const {
    bool db = (drand() < 0.00001);
    db = false;
    hit_record rec1, rec2;
    if (boundary->hit(r, -FLT_MAX, FLT_MAX, rec1)) {
        if (boundary->hit(r, rec1.t + 0.0001, FLT_MAX, rec2)) {
            if (db) fprintf(stderr, "\nt0 t1 %f %f\n", rec1.t, rec2.t);
            if (rec1.t < t_min) rec1.t = t_min;
            if (rec2.t > t_max) rec2.t = t_max;
            if (rec1.t >= rec2.t) return false;
            if (rec1.t < 0) rec1.t = 0;
            float distance_inside_boundary = (rec2.t - rec1.t) * r.direction().length();
            float hit_distance = -(1.0 / density) * log(drand());
            if (hit_distance < distance_inside_boundary) {
                rec.t = rec1.t + hit_distance / r.direction().length();
                rec.p = r.point_at_parameter(rec.t);
                rec.normal = vec3(1, 0, 0); // arbitrary
                rec.mat_ptr = phase_function;
                if (db) fprintf(stderr, "rec.t = %f\n", rec.t);
                if (db) fprintf(stderr, "hit_distance = %f\n", hit_distance);
                if (db) std::cerr << "rec.p = " << rec.p << "\n";
                return true;
            }
        }
    }
    return false;
}

/* TUBE */
class tube : public hitable {
public:
    tube(){}
    tube(vec3 b, vec3 e, float r, material *mat) : begin(b), end(e), radius(r), mat_ptr(mat){};
    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec)const;

    vec3 begin, end;
    float radius;
    material *mat_ptr;
};

bool tube::hit(const ray& r, float t_min, float t_max, hit_record& rec)const{

    return false;
}



#endif
