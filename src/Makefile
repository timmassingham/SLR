CC = gcc
CFLAGS = -O4 -funroll-loops -DNDEBUG -fomit-frame-pointer -std=gnu99
#CFLAGS = -O4 -funroll-loops -DNDEBUG -DBLAS -DLAPACK -fomit-frame-pointer -std=gnu99
#CFLAGS = -static -O3 -funroll-loops -DNDEBUG -std=gnu99 -Wall
#CFLAGS =  -DNDEBUG -pg  -Wall -DBLAS -DLAPACK -std=gnu99
#CFLAGS = -Wall -g -DWARNINGS -std=gnu99 -pedantic
#CFLAGS = -O -funroll-loops  -std=gnu99 -DNDEBUG
LD = ld
BLASLIB = atlas
LDFLAGS =  -lm
INCFLAGS = 
LIBDICTDIR = ../libdict-0.2.1
#LDFLAGS = -lm

objects = like.o tree.o data.o rng.o model.o  bases.o codonmodel.o gencode.o utility.o matrix.o optimize_new.o spinner.o options.o tree_data.o linemin.o gamma.o statistics.o mystring.o nucmodel.o root.o vec.o
lineobj = brent/fminbr.o


Slr: slr.o $(objects) $(lineobj) Meschach/meschach.a  libdict
	gcc  $(INCFLAGS) $(CFLAGS) -o ../bin/$@ slr.o $(objects) $(lineobj) $(LDFLAGS) Meschach/meschach.a $(LIBDICTDIR)/libdict.a 

libdict:
	cd $(LIBDICTDIR) && make
.c.o:
	$(CC) $(CFLAGS) $(INCFLAGS) -o $@ -c $<
.f.o:
	g77 -O3 -o $@ -c $<
clean:
	rm -f *~ *.o *.obj tn/*.o tn/*.obj *.obj brent/*.o brent/*.obj Meschach/*.o Meschach/*.obj *.exe *.lib Meschach/*.a && cd $(LIBDICTDIR) && make clean

Meschach/meschach.a:
	cd Meschach && make basic
