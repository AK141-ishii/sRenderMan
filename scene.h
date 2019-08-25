#ifndef SCENE_H_
#define SCENE_H_
#include <random>
#include "hitable_list.h"
#include "primitive.h"
#include "material.h"
#include "texture.h"
#include "bvh.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


hitable *random_scene() {
    int n=5000;
    hitable **list = new hitable*[n+1];
    texture *checker = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)), new constant_texture(vec3(0.9, 0.9, 0.9)));
    list[0] = new sphere(vec3(0,-1000,0),1000,new lambertian(checker));
    int i=1;
    for(int a=-11;a<11;a++){
        for(int b=-11;b<11;b++){
            float choose_mat = drand();
            vec3 center(a+0.9*drand(),0.2,b+0.9*drand());
            if((center-vec3(4,0.2,0)).length() > 0.9){
                if(choose_mat < 0.8) { //diffuse
                    list[i++] = new moving_sphere(center, center + vec3(0.5 * drand(), 0, 0.5 * drand()), 0.0, 1.0, 0.2, new lambertian(new constant_texture(vec3(drand(), drand(), drand()))));
                }
                else if(choose_mat < 0.95){ //metal
                    list[i++] = new sphere(center, 0.2, new metal(vec3(
                        0.5*(1+drand()), 0.5*(1+drand()), 0.5*(1+drand())),
                        0.5*drand()));
                }
                else {
                    list[i++] = new sphere(center, 0.2, new dielectric(1.5)); 
                }
            }
        }
    }
    list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
    list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(new constant_texture(vec3(0.4, 0.2, 0.1))));
    list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
    return new bvh_node(list,i, 0.0, 1.0);
}

hitable *two_sphere() {
    texture *checker = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)), new constant_texture(vec3(0.9, 0.9, 0.9)));
    int n = 50;
    hitable **list = new hitable*[n+1];
    list[0] = new sphere(vec3(0, -10, 0), 10, new lambertian(checker));
    list[1] = new sphere(vec3(0, 10, 0), 10, new lambertian(checker));

    return new hitable_list(list, 2);
}

hitable *two_perlin_spheres() {
    texture *pertext = new noise_texture(5);
    hitable **list = new hitable*[2];
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
    list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian(pertext));
    return new hitable_list(list, 2);
}

hitable *earth_with_metal() {
    int nx, ny, nn;
    unsigned char * tex_data = stbi_load("4096px-Nasa_land_ocean_ice_8192.jpg", &nx, &ny, &nn, 0);
    material *mat = new lambertian(new image_texture(tex_data, nx, ny));
    hitable **list = new hitable*[7];
    list[0] = new sphere(vec3(-8, 0, 0), 2, new metal(vec3(0.5, 0.5, 0.5), 0));
    list[1] = new sphere(vec3(8, 0, 0), 2, new metal(vec3(0.5, 0.5, 0.5), 0));
    list[2] = new sphere(vec3(0, 10, 0), 2, new metal(vec3(0.5, 0.5, 0.5), 0));
    list[3] = new sphere(vec3(0, 0, -8), 2, new metal(vec3(0.5, 0.5, 0.5), 0));
    list[4] = new sphere(vec3(0, 0, 8), 2, new metal(vec3(0.5, 0.5, 0.5), 0));

    list[5] = new sphere(vec3(0, 3, 0), 3, mat);
    list[6] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(new constant_texture(vec3(0.8, 0.8, 0.8))));
    return new hitable_list(list, 7);
}

hitable * simple_light() {
    texture *pertext = new noise_texture(4);
    hitable **list = new hitable*[4];
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
    list[1] = new sphere(vec3(0,2,0), 2, new lambertian(pertext));
    list[2] = new sphere(vec3(0,7,0), 2, new diffuse_light(new constant_texture(vec3(4,4,4))));
    list[3] = new xy_rect(3, 5, 1, 3, -2, new diffuse_light(new constant_texture(vec3(4, 4, 4))));
    return new hitable_list(list, 4);
}

hitable * cornell_box() {
    hitable **list = new hitable *[8];
    int i =0;
    material *red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
    material *white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
    material *green = new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
    material *light = new diffuse_light(new constant_texture(vec3(15, 15, 15)));
    list[i++] = new flip_normals(new yz_rect(0, 555, 0, 555, 555, green));
    list[i++] = new yz_rect(0, 555, 0, 555, 0, red);
    list[i++] = new zx_rect(227, 332, 213, 343, 554, light);
    list[i++] = new flip_normals(new zx_rect(0, 555, 0, 555, 555, white));
    list[i++] = new zx_rect(0, 555, 0, 555, 0, white);
    list[i++] = new flip_normals(new xy_rect(0, 555, 0, 555, 555, white));
    list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 165, 165), white), -18), vec3(130, 0, 65));
    list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 330, 165), white), 15), vec3(265, 0, 295));
    return new hitable_list(list, i);
}


#endif