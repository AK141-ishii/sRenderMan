#include <stdio.h>
#include <assert.h>
#include <stdint.h>
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


int main() {
    int nx = 1200;
    int ny = 800;
    int ns = 100;
    std::cout << "P3\n" << nx << " " << ny << "\n255\n";

    hitable *world = final();

    vec3 lookfrom(278, 278, -500);
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
            col[0] = (col[0] > 1) ? 1.0f : col[0];
            col[1] = (col[1] > 1) ? 1.0f : col[1];
            col[2] = (col[2] > 1) ? 1.0f : col[2];
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
