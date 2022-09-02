CFLAGS = -O
CONTOBJ = connect.o contree.o
EXELIST = kin3Dcont kin2Dcont

.c.o:
	cc $(CFLAGS) -c $*.c

all: $(EXELIST)

kin2Dcont: kin2Dcont.o $(CONTOBJ) nrmat.o readpoints.o
	cc -o $@ $(CFLAGS) kin2Dcont.o $(CONTOBJ) nrmat.o readpoints.o -lm

kin3Dcont: kin3Dcont.o $(CONTOBJ) nrmat.o readpoints.o
	cc -o $@ $(CFLAGS) kin3Dcont.o $(CONTOBJ) nrmat.o readpoints.o -lm

clean:
	-rm -f *.ckp *.o

spotless:
	-rm -f *.ckp *.o $(EXELIST)

install: probe
	echo "no install"

# DO NOT DELETE THIS LINE -- make depend uses it
connect.o: connect.h contree.h connect.c
kin3Dcont.o: connect.h contree.h nrmat.h readpoints.h kin3Dcont.c
kin2Dcont.o: connect.h contree.h nrmat.h readpoints.h kin2Dcont.c
contree.o: contree.h ctdata.h contree.c
nrmat.o: nrmat.h nrmat.c
readpoints.o: connect.h contree.h nrmat.h readpoints.h readpoints.c
# DO NOT DELETE THIS 2nd LINE -- make depend uses it
