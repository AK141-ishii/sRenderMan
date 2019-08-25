#include <iostream>
#include <random>
#include "texture.h"
#include "primitive.h"
#include "hitable_list.h"
#include "camera.h"
#include "material.h"
#include "scene.h"
#include <float.h>

vec3 color(const ray &r, hitable *world, int depth) {
    hit_record rec;
    if (world->hit(r, 0.001, FLT_MAX, rec)) {
        ray scattered;
        vec3 attenuation;
        vec3 emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return emitted + attenuation * color(scattered, world, depth + 1);
        else
            return emitted;
        
    }
    else 
        return vec3(0, 0, 0);
}


int main() {
    int nx = 1200;
    int ny = 800;
    int ns = 100;
    std::cout << "P3\n" << nx << " " << ny << "\n255\n";

    hitable *world = cornell_box();

    vec3 lookfrom(278, 278, -800);
    vec3 lookat(278, 278, 0);
    float dist_to_focus = 10.0;
    float aperture = 0.0;
    float vfov = 40.0;

    camera cam(lookfrom, lookat, vec3(0,1,0), vfov, float(nx)/float(ny), aperture, dist_to_focus, 0.0, 1.0);

    std::cerr << "RENDERING START" << std::endl;

    for (int j = ny - 1; j >= 0; j--) {

        if((j+1)%100==0){// report
            std::cerr << "Round " << ny - j - 1 << "/" << ny <<"\t" ;
        } if(j%10==0){ std::cerr << "■"; }

        for (int i = 0; i < nx; i++) {
            vec3 col(0, 0, 0);
            for (int s = 0; s < ns; s++) {
                float u = float(i + (float)rand()/(float)RAND_MAX) / float(nx);
                float v = float(j + (float)rand()/(float)RAND_MAX) / float(ny);
                ray r = cam.get_ray(u, v);
                col += color(r, world, 0);
            }
            col /= float(ns);
            col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
            int ir = int(255.99 * col[0]);
            int ig = int(255.99 * col[1]);
            int ib = int(255.99 * col[2]);

            std::cout << ir << " " << ig << " " << ib << "\n";
        }
        if(j%100==0){// report
            std::cerr << " done;" << std::endl;
        }

    }

    return 0;
}