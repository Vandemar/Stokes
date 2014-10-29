CPP = g++
CPPFLAGS= -g

all: advect

advect: 
	${CPP} ${CPPFLAGS} -I /home/clear-skies/workspace/stokes/1-D/include/ -L /usr/local/lib -o advect main.cpp /usr/local/lib/libboost_program_options.so geometry.cpp  

clean:
	rm advect
