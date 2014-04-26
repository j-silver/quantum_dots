CC = gcc
CFLAGS = -g -Wall

main.o:	main.c
	$(CC) $(CFLAGS) -c main.c

Reg_cc.o: Reg_cc.c
	$(CC) $(CFLAGS) -c Reg_cc.c

Reg_ss.o: Reg_ss.c
	$(CC) $(CFLAGS) -c Reg_ss.c

Img_cs.o: Img_cs.c
	$(CC) $(CFLAGS) -c Img_cs.c

Img_sc.o: Img_sc.c
	$(CC) $(CFLAGS) -c Img_sc.c

Reg_c0.o: Reg_c0.c
	$(CC) $(CFLAGS) -c Reg_c0.c

Img_c0.o: Img_c0.c Expi.o
	$(CC) $(CFLAGS) -c Img_c0.c

Img_ss.o: Img_ss.c Expi.o
	$(CC) $(CFLAGS) -c Img_ss.c

Img_cc.o: Img_cc.c Expi.o
	$(CC) $(CFLAGS) -c Img_cc.c

Expi.o:	Expi.c
	$(CC) $(CFLAGS) -c Expi.c

Reg_s0.o: Reg_s0.c
	$(CC) $(CFLAGS) -c Reg_s0.c

Img_s0.o: Img_s0.c
	$(CC) $(CFLAGS) -c Img_s0.c

Reg_sc.o: Reg_sc.c
	$(CC) $(CFLAGS) -c Reg_sc.c

integs.o: integs.c
	$(CC) $(CFLAGS) -c integs.c

red_gen.o: red_gen.c
	$(CC) $(CFLAGS) -c red_gen.c

main:	main.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Img_ss.o Img_cc.o  Reg_s0.o Img_s0.o Reg_sc.o red_gen.o integs.o 
	$(CC) -lgsl -lgslcblas -lm -o main main.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Img_c0.o Expi.o Img_ss.o Img_cc.o  Reg_s0.o Img_s0.o Reg_sc.o red_gen.o integs.o

clean: 
	rm -f main main.o Reg_cc.o Reg_ss.o Img_cs.o Img_sc.o Reg_c0.o Reg_c0.o Img_c0.o Expi.o Img_ss.o Img_cc.o  Img_s0.o Reg_s0.o Reg_sc.o red_gen.o integs.o

