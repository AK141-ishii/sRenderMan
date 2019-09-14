#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <memory>
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
#include "pdf.h"

inline vec3 de_nan(const vec3& v) {
    vec3 temp = v;
    if(!(temp[0] == temp[0])) temp[0] = 0;
    if(!(temp[1] == temp[1])) temp[1] = 0;
    if(!(temp[2] == temp[2])) temp[2] = 0;
    return temp;
}

vec3 color(const ray &r, hitable *world, hitable *light_shape, int depth) {
    hit_record hrec;
    if (world->hit(r, 0.001, FLT_MAX, hrec)) {

        if (hrec.t < 0.001 || !(hrec.t == hrec.t) || !(hrec.v == hrec.v) || !(hrec.u == hrec.u))
            return vec3(0, 0, 0);

        scatter_record srec;
        vec3 emitted = hrec.mat_ptr->emitted(r, hrec, hrec.u, hrec.v, hrec.p);
        if (depth < 50 && hrec.mat_ptr->scatter(r, hrec, srec)) {
            if(srec.is_specular) {
                return srec.attenuation * color(srec.specular_ray, world, light_shape, depth + 1);
            }
            else {
                hitable_pdf plight(light_shape, hrec.p);
                mixture_pdf p(&plight, srec.pdf_ptr);
                ray scattered = ray(hrec.p, p.generate(), r.time());
                float pdf_val = p.value(scattered.direction());
                delete(srec.pdf_ptr);
                return emitted + srec.attenuation * hrec.mat_ptr->scattering_pdf(r, hrec, scattered) * color(scattered, world, light_shape, depth + 1) / pdf_val;
            }
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

int main(int argc, char *argv[]) {
    int nx = 1200;
    int ny = 800;
    int ns = 10;

    std::fstream fs;
    fs.open("tmp.ppm", std::ios::out);
    fs << "P3\n" << nx << " " << ny << "\n255\n";

    hitable *world;
    camera *cam;
    float aspect = float(nx) / float(ny);
    hitable_list *hlist;

    light_through_glass(&world, &cam, aspect, &hlist);

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
                float u = float(i + drand()) / float(nx);
                float v = float(j + drand()) / float(ny);
                ray r = cam->get_ray(u, v);
                col += de_nan(color(r, world, hlist, 0));
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
