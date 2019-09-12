#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <random>
#include <float.h>
#include <time.h>
#include "texture.h"
#include "primitive.h"
#include "hitable_list.h"
#include "camera.h"
#include "material.h"
#include "scene.h"

vec3 color(const ray &r, hitable *world, int depth) {
    hit_record rec;
    if (world->hit(r, 0.001, FLT_MAX, rec)) {

        if(rec.t==0)return vec3(0, 0, 0);

        ray scattered;
        vec3 emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
        float pdf;
        vec3 albedo;
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, albedo, scattered, pdf)) {
            vec3 on_light = vec3(213 + drand() * (343 - 213), 554, 227 + drand() * (332 - 227));
            vec3 to_light = on_light - rec.p;
            float distance_squared = to_light.squared_length();
            to_light.make_unit_vector();
            if (dot(to_light, rec.normal) < 0)
                return emitted;
            float light_area = (343-213)*(332-227);
            float light_cosine = fabs(to_light.y());
            if(light_cosine < 0.000001)
                return emitted;
            pdf = distance_squared / (light_cosine * light_area);
            scattered = ray(rec.p, to_light, r.time());

            return emitted + albedo * rec.mat_ptr->scattering_pdf(r, rec, scattered) * color(scattered, world, depth + 1) / pdf;
        }
        else
            return emitted;
    }
    else 
        return vec3(0, 0, 0);
}

class framebuffer {
public:
    static constexpr size_t kBytesPerPixel = 3;

    framebuffer(size_t width, size_t height):
        width_(width),
        height_(height),
        data_((uint8_t*)malloc(width * height * kBytesPerPixel)) {}

    ~framebuffer() { free(data_); }

    void set_pixel(size_t row,
                    size_t col,
                    size_t r,
                    size_t g,
                    size_t b ) {
        const size_t idx = kBytesPerPixel * (row * width_ + col);
        data_[idx + 0] = b;
        data_[idx + 1] = g;
        data_[idx + 2] = r;
    }

    void save(const char* file_path) const {
        FILE* fptr = fopen(file_path, "wb");
        assert(fptr);
        putc(0,fptr);
        putc(0, fptr);
        putc(2, fptr); /* uncompressed RGB */
        putc(0, fptr);
        putc(0, fptr);
        putc(0, fptr);
        putc(0, fptr);
        putc(0, fptr);
        putc(0, fptr);
        putc(0, fptr); /* X origin */
        putc(0, fptr);
        putc(0, fptr); /* y origin */
        putc((width_ & 0x00FF), fptr);
        putc((width_ & 0xFF00) / 256, fptr);
        putc((height_ & 0x00FF), fptr);
        putc((height_ & 0xFF00) / 256, fptr);
        putc(24, fptr); /* 24 bit bitmap */
        putc(0, fptr);
        fwrite(data_, kBytesPerPixel, width_ * height_, fptr);
        fclose(fptr);
    }

    size_t width(){return width_;}
    size_t height(){return height_;}

private:
    uint8_t *data_;
    size_t  width_;
    size_t  height_;
};

void cornell_box(hitable **scene, camera **cam, float aspect)
{
    int i = 0;
    hitable ** list = new hitable*[8];
    material *red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
    material *white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
    material *green = new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
    material *light = new diffuse_light(new constant_texture(vec3(15, 15, 15)));

    list[i++] = new flip_normals(new yz_rect(0, 555, 0, 555, 555, green));
    list[i++] = new yz_rect(0, 555, 0, 555, 0, red);
    list[i++] = new flip_normals(new zx_rect(227, 332, 213, 343, 554, light));
    list[i++] = new flip_normals(new zx_rect(0, 555, 0, 555, 555, white));
    list[i++] = new zx_rect(0, 555, 0, 555, 0, white);
    list[i++] = new flip_normals(new xy_rect(0, 555, 0, 555, 555, white));
    list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 165, 165), white), -18), vec3(130, 0, 65));
    list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 330, 165), white), 15), vec3(265, 0, 295));

    *scene = new hitable_list(list, i);

    vec3 lookfrom(278, 278, -800);
    vec3 lookat (278,278,0);
    float dist_to_focus = 10.0;
    float aperture = 0.0;
    float vfov = 40.0;
    *cam = new camera(lookfrom, lookat, vec3(0,1,0), vfov, aspect, aperture, dist_to_focus, 0.0, 1.0);
}

int main(int argc, char *argv[]) {
    int nx = 1200;
    int ny = 800;
    int ns = 10;

    std::fstream fs;
    fs.open("tmp.ppm", std::ios::out);
    fs << "P3\n" << nx << " " << ny << "\n255\n";

    hitable *world;
    camera *cam;
    cornell_box(&world, &cam, (float)nx/(float)ny);

//    vec3 lookfrom(278, 278, -500);
//    vec3 lookat(278, 278, 0);
//    float dist_to_focus = 10.0;
//    float aperture = 0.0;
//    float vfov = 40.0;
//
//    camera cam(lookfrom, lookat, vec3(0,1,0), vfov, float(nx)/float(ny), aperture, dist_to_focus, 0.0, 1.0);

    std::cerr << "RENDERING START" << std::endl;
    time_t begin_time = time(NULL);
    time_t mid_time_before = begin_time;
    time_t mid_time_current;



    for (int j = ny - 1; j >= 0; j--) {

        if((j+1)%100==0){// report
            std::cerr << "Round " << ny - j - 1 << "/" << ny <<"\t" ;
        } if(j%10==0){ std::cerr << "■"; }

        for (int i = 0; i < nx; i++) {
            vec3 col(0, 0, 0);
            for (int s = 0; s < ns; s++) {
                float u = float(i + (float)rand()/(float)RAND_MAX) / float(nx);
                float v = float(j + (float)rand()/(float)RAND_MAX) / float(ny);
                ray r = cam->get_ray(u, v);
                col += color(r, world, 0);
            }
            col /= float(ns);
            col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
            col[0] = (col[0] > 1) ? 1.0f : col[0];
            col[1] = (col[1] > 1) ? 1.0f : col[1];
            col[2] = (col[2] > 1) ? 1.0f : col[2];
            int ir = int(255.99 * col[0]);
            int ig = int(255.99 * col[1]);
            int ib = int(255.99 * col[2]);

            fs << ir << " " << ig << " " << ib << "\n";
        }
        if(j%100==0){// report
            mid_time_current = time(NULL);
            std::cerr << " done | " << difftime(mid_time_current, mid_time_before) << " sec" << std::endl;
            mid_time_before = mid_time_current;
        }

    }
    int spending_time = (int)difftime(mid_time_current, begin_time);
    std::cerr << "TIME : " << (spending_time / 60) << " min " << (spending_time % 60) << " sec" << std::endl;

    fs.close();

    return 0;
}
