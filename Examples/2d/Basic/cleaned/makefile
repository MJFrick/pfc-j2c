CC = h5pcc
IDIR = ../include
ODIR = ../obj

CFLAGS = -I./include -L/usr/local/lib #-lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas

LIBS = -lfftw3_mpi -lfftw3 -lm -lgsl -lgslcblas

main: obj/state.o obj/main.o obj/io.o obj/dynamics.o obj/wisdom.o obj/error.o
	$(CC) obj/* $(CFLAGS) $(LIBS) -o main

obj/main.o: src/main.c
	$(CC) -c src/main.c $(CFLAGS) -o obj/main.o

obj/state.o: src/state.c include/pfc_header.h
	$(CC) -c src/state.c $(CFLAGS) -o obj/state.o

obj/io.o: src/io.c include/pfc_header.h
	$(CC) -c src/io.c $(CFLAGS) -o obj/io.o

obj/dynamics.o: src/dynamics.c include/pfc_header.h
	$(CC) -c src/dynamics.c $(CFLAGS) -o obj/dynamics.o

obj/wisdom.o: src/wisdom.c include/pfc_header.h
	$(CC) -c src/wisdom.c $(CFLAGS) -o obj/wisdom.o

obj/error.o: src/error.c include/pfc_header.h
	$(CC) -c src/error.c $(CFLAGS) -o obj/error.o

.PHONY: clean

clean:
	rm -f data/Data.h5
	rm -f obj/*
	rm -f main
