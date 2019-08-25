#Makefile

SRC :=
SRC += main.cpp
SRC += vec3.h
SRC += ray.h
SRC += primitive.h
SRC += hitable.h
SRC += hitable_list.h
SRC += camera.h
SRC += material.h
SRC += scene.h
SRC += aabb.h
SRC += bvh.h
SRC += texture.h
SRC += perlin.h

a.out: ${SRC}
	${CXX} -g3 -fopenmp -o $@ main.cpp

.PHONY: clean
clean:
	${RM} a.out

.PHONY: out
out: a.out
	./a.out > tmp.ppm

# ::: AK141 note :::
#
# $@ : Target Name (such as DBT.o...)
#
# $^ : File Name list
# $< : First File Name (such as DBT.cpp...)
#
# ${CXX} = C++ compiler (like g++...)
# ${RM}  = rm -f
#
# .PHONY: <- No dependency & No output
# .SUFFIXED: <- I don't know, but it may be useful.
