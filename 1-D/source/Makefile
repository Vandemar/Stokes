HOME = /home/cig/workspace/stokes/1-D
CPP = g++
CPPFLAGS= -g

all: advect

advect: 
	${CPP} ${CPPFLAGS} -I ${HOME}/include -L /usr/local/lib -o advect main.cpp /usr/local/lib/libboost_program_options.so.1.55.0 geometry.cpp  

clean:
	rm advect
