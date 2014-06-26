CC = gcc-4.9
CFLAGS = -O2 -march=core2 -Wall -ansi -Wpointer-arith -Wcast-qual -Wcast-align -Wshadow -Wconversion -Wmissing-prototypes -Wstrict-prototypes -fno-common -Wnested-externs -Wfloat-equal -fstack-protector -Wstack-protector # -g -W -pedantic -Wwrite-strings  -fshort-enums 
LDLIBS = -lgsl -lgslcblas -lm

libs_for_intelc = -L/media/SPACE/intel/lib/intel64 -mkl -lgsl
headers_for_intelc = -I/media/SPACE/intel/mkl/include -I/media/SPACE/intel/include/intel64

ifeq ($(CC),icc)
	CFLAGS = -O2 -march=core2 -Wall -ansi -Wpointer-arith -Wcast-qual -Wshadow -Wconversion -Wmissing-prototypes -Wstrict-prototypes -fno-common -Wfloat-equal 
	LDLIBS = $(libs_for_intelc)
	CPPFLAGS = $(headers_for_intelc)
	MAIN = main-intel
else ifeq ($(CC),clang)
	MAIN = main-clang
else ifeq ($(CC),gcc-4.9)
	CFLAGS += -fdiagnostics-color
endif


#
# main programs
#
# Evolution of the Redfield and CP dynamics (state, current and entropy)
#
red_evol.o: red_evol.c funcs.h initial.h

red_evol: red_evol.o evol.o mat_file.o entropy.o total_current.o
	$(CC) -o red_evol red_evol.o evol.o mat_file.o entropy.o total_current.o $(LDLIBS) $(LDFLAGS) 

cp_evol.o: cp_evol.c funcs.h initial.h

cp_evol: cp_evol.o evol.o mat_file.o entropy.o total_current.o
	$(CC) -o cp_evol cp_evol.o evol.o mat_file.o entropy.o total_current.o $(LDLIBS) $(LDFLAGS)

	
#
# Stationary currents as functions of T/D and Omega/D
#
current_tdel.o: current_tdel.c funcs.h

current_tdel: tdel_objects
	$(CC) -o current_tdel $(TDEL_OBJECTS) $(LDLIBS) $(LDFLAGS)


current_omegad.o: current_omegad.c funcs.h

current_omegad:	omegad_objects
	$(CC) -o current_omegad $(OMEGAD_OBJECTS) $(LDLIBS) $(LDFLAGS)

# 
# Total current
#
total_current.o: total_current.c funcs.h initial.h


#
# Asymptotic current and generators matrices
#
asymptotic.o: asymptotic.c funcs.h initial.h

asymptotic: asymptotic_objects
	$(CC) -o asymptotic $(ASYMPTOTIC_OBJECTS) $(LDLIBS) $(LDFLAGS)


#
# Positivity check 
#
sample.o: sample.c funcs.h

r0dot.o: r0dot.c funcs.h

polar.o: polar.c funcs.h

sample: sample.o r0dot.o polar.o mat_file.o
	$(CC) -o sample sample.o r0dot.o polar.o mat_file.o $(LDLIBS) $(LDFLAGS)

# Average
#
average.o: average.c funcs.h initial.h

average: average_objects
	$(CC) -o average $(AVERAGE_OBJECTS) $(LDLIBS) $(LDFLAGS)

#
# Object files
#
stationary.o: stationary.c funcs.h 

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

entropy.o: entropy.c funcs.h

current.o: current.c funcs.h

write.o: write.c funcs.h


ASYMPTOTIC_OBJECTS = asymptotic.o stationary.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Img_ss.o Img_cc.o Reg_s0.o Img_s0.o Reg_sc.o red_gen.o integs.o Expi.o evol.o cp_gen.o mat_file.o hamil.o entropy.o write.o 

OMEGAD_OBJECTS = stationary.o current_omegad.o red_gen.o cp_gen.o integs.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Img_ss.o Img_cc.o Reg_s0.o Img_s0.o Reg_sc.o Expi.o

TDEL_OBJECTS = stationary.o current_tdel.o red_gen.o cp_gen.o integs.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Img_ss.o Img_cc.o Reg_s0.o Img_s0.o Reg_sc.o Expi.o

SAMPLE_OBJECTS = sample.o r0dot.o polar.o mat_file.o

EVOL_OBJECTS = red_evol.o cp_evol.o total_current.o

AVERAGE_OBJECTS = average.o evol.o entropy.o mat_file.o

DATA_FILES = REDFIELD_MATRIX CP_MATRIX RED-EVOLUTION.dat RED-CURRENT.dat RED-ENTROPY.dat CP-EVOLUTION.dat CP-CURRENT.dat CP-ENTROPY.dat INTEGRALS.dat RED-STAT-CURR-T.dat CP-STAT-CURR-T.dat RED-STAT-CURR-O.dat CP-STAT-CURR-O.dat CP_STATIONARY.dat RED_STATIONARY.dat POS_VIOLATIONS CP-ENTROPY-PROD.dat RED-ENTROPY-PROD.dat

asymptotic_objects: $(ASYMPTOTIC_OBJECTS)

omegad_objects: $(OMEGAD_OBJECTS)
	
tdel_objects: $(TDEL_OBJECTS) 
	
sample_objects: $(SAMPLE_OBJECTS)

average_objects: $(AVERAGE_OBJECTS)


.PHONY: clean

clean_backups:
	rm -f *~

clean: 
	rm -f $(ASYMPTOTIC_OBJECTS) $(OMEGAD_OBJECTS) $(TDEL_OBJECTS) $(SAMPLE_OBJECTS) $(DATA_FILES) $(EVOL_OBJECTS) $(AVERAGE_OBJECTS) asymptotic current_omegad current_tdel red_evol cp_evol sample average

