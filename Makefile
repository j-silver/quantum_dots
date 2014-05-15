CC = gcc
CFLAGS = -O2 -march=core2 -Wall -ansi -Wpointer-arith -Wcast-qual -Wcast-align -Wshadow -Wconversion -Wmissing-prototypes -Wstrict-prototypes -fno-common #-g -W -pedantic -Wwrite-strings -Wnested-externs -fshort-enums 
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

hamil.o: hamil.c funcs.h

station.o: station.c funcs.h

entropy.o: entropy.c funcs.h

current.o: current.c funcs.h

write.o: write.c funcs.h

objects = main.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Img_ss.o Img_cc.o Reg_s0.o Img_s0.o Reg_sc.o red_gen.o integs.o Expi.o evol.o cp_gen.o mat_file.o hamil.o station.o entropy.o write.o 

data_files = REDFIELD_MATRIX CP_MATRIX RED-EVOLUTION.dat RED-CURRENT.dat RED-ENTROPY.dat CP-EVOLUTION.dat CP-CURRENT.dat CP-ENTROPY.dat INTEGRALS.dat RED-STAT-CURR-T.dat CP-STAT-CURR-T.dat RED-STAT-CURR-O.dat CP-STAT-CURR-O.dat

data: main
	./main

main: $(objects)
	$(CC) -o $(MAIN) $(objects) $(LDLIBS) $(LDFLAGS)



#
# Stationary currents as functions of T/D and Omega/D
#
current_tdel.o: current_tdel.c funcs.h

current_tdel:	current_tdel.o red_gen.o cp_gen.o station.o integs.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Img_ss.o Img_cc.o Reg_s0.o Img_s0.o Reg_sc.o Expi.o
	$(CC) -o current_tdel current_tdel.o red_gen.o cp_gen.o station.o integs.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Img_ss.o Img_cc.o Reg_s0.o Img_s0.o Reg_sc.o Expi.o $(LDLIBS) $(LDFLAGS)

current_omegad.o: current_omegad.c funcs.h

current_omegad:	current_omegad.o red_gen.o cp_gen.o station.o integs.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Img_ss.o Img_cc.o Reg_s0.o Img_s0.o Reg_sc.o Expi.o
	$(CC) -o current_omegad current_omegad.o red_gen.o cp_gen.o station.o integs.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Img_ss.o Img_cc.o Reg_s0.o Img_s0.o Reg_sc.o Expi.o $(LDLIBS) $(LDFLAGS)

#
# Evolution of the Redfield and CP dynamics
#
red_evol.o: red_evol.c funcs.h initial.h

red_evol: red_evol.o evol.o mat_file.o station.o entropy.o
	$(CC) -o red_evol red_evol.o evol.o mat_file.o station.o entropy.o $(LDLIBS) $(LDFLAGS) 

cp_evol.o: cp_evol.c funcs.h initial.h

cp_evol: cp_evol.o evol.o mat_file.o station.o entropy.o
	$(CC) -o cp_evol cp_evol.o evol.o mat_file.o station.o entropy.o $(LDLIBS) $(LDFLAGS)


.PHONY: clean

clean: 
	rm -f $(objects) $(data_files) main current_tdel current_tdel.o current_omegad current_omegad.o red_evol.o cp_evol.o red_evol cp_evol

