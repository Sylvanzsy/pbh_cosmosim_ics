f90invoke = gfortran -g -fbacktrace

CC = gcc

OBJS = random_object.o kdtree.o pbh_ic.o

OPTF = -DTEST

# General rules
#.SUFFIXES:
.SUFFIXES: .exe .o .F90 .f90 .mod .f .c

.F90.o:
	$(f90invoke) -c $(OPTF) $< -fbounds-check -o $*.o $(F90) 

.f90.o:
	$(f90invoke) -c $(OPTF) $< -fbounds-check -o $*.o $(F90) 

.f.o:
	$(f90invoke) -c $(OPTF) $< -fbounds-check -o $*.o $(F90)

.c.o:
	$(CC) -c $< -fbounds-check -o $*.o

all:	pbh_ic.exe
#all:	test_cfe.exe

clean:
	\rm -f ./*.o *.o *.mod *.exe

pbh_ic.o: random_object.o
pbh_ic.o: kdtree.o

pbh_ic.exe: $(OBJS)
	$(f90invoke)  $(OBJS) -o pbh_ic.exe

#test_cfe.exe: $(OBJS)
#	$(f90invoke) $(OBJS) -o test_cfe.exe -lm

#	$(CC) $(OBJS) -o test_cfe.exe -lgfortran -lm
