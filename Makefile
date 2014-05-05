CC = gcc
CFLAGS = -g -Wall -O2 #-pedantic -ansi -W -Wmissing-prototypes -Wstrict-prototypes  -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -fshort-enums -fno-common 
LDLIBS = -lgsl -lgslcblas -lm

MAIN = main

libs_for_intelc = -L/media/SPACE/intel/lib/intel64 -mkl -lgsl
headers_for_intelc = -I/media/SPACE/intel/mkl/include -I/media/SPACE/intel/include/intel64

ifeq ($(CC),icc)
	LDLIBS = $(libs_for_intelc)
	CPPFLAGS = $(headers_for_intelc)
	MAIN = main-intel
else ifeq ($(CC),clang)
	MAIN = main-clang
endif

main.o:	main.c funcs.h

Reg_cc.o: Reg_cc.c funcs.h

Reg_ss.o: Reg_ss.c funcs.h

Img_cs.o: Img_cs.c funcs.h

Img_sc.o: Img_sc.c funcs.h

Reg_c0.o: Reg_c0.c funcs.h

Img_c0.o: Img_c0.c Expi.o funcs.h

Img_ss.o: Img_ss.c Expi.o funcs.h

Img_cc.o: Img_cc.c Expi.o funcs.h

Expi.o:	Expi.c funcs.h

Reg_s0.o: Reg_s0.c funcs.h

Img_s0.o: Img_s0.c funcs.h

Reg_sc.o: Reg_sc.c funcs.h

integs.o: integs.c funcs.h

red_gen.o: red_gen.c funcs.h

cp_gen.o: cp_gen.c funcs.h

mat_file.o: mat_file.c funcs.h

evol.o:	evol.c funcs.h

hamil.o: hamil.c 

station.o: station.c 

entropy.o: entropy.c

objects = main.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Img_ss.o Img_cc.o Reg_s0.o Img_s0.o Reg_sc.o red_gen.o integs.o Expi.o evol.o cp_gen.o mat_file.o hamil.o station.o entropy.o

main: $(objects)
	$(CC) -o $(MAIN) $(objects) $(LDLIBS) $(LDFLAGS)

.PHONY: clean

clean: 
	rm -f $(objects) REDFIELD_MATRIX CP_MATRIX EVOLUTION.dat

