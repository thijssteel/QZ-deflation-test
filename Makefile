include make.inc


LIBOBJS = ./obj/dhgeqzv1.o ./obj/shgeqzv1.o ./obj/zhgeqzv1.o

all: lib tests

lib: $(LIBOBJS)

tests: ./bin/test_accuracy.out ./bin/test_inf.out ./bin/test_special_matrix.out

./obj/dhgeqzv1.o : ./src/dhgeqzv1.f make.inc
	$(FC) $(FFLAGS) $(INCFLAG)./include -o $@ -I./include -c $<

./obj/shgeqzv1.o : ./src/shgeqzv1.f make.inc
	$(FC) $(FFLAGS) $(INCFLAG)./include -o $@ -I./include -c $<

./obj/util.o : ./tests/util.f90 make.inc
	$(FC) $(FFLAGS) $(INCFLAG)./include -o $@ -I./include -c $<

./bin/%.out: ./tests/%.f90 ./obj/util.o $(LIBOBJS) make.inc
	$(FC) $(FFLAGS) $(INCFLAG)./include -o $@ $< ./obj/util.o $(LIBOBJS) $(LIBS)