CC = gcc
CFLAGS = -O4 -funroll-loops -DNDEBUG -fomit-frame-pointer -std=gnu99
#CFLAGS = -static -O3 -funroll-loops -DNDEBUG -std=gnu99 -Wall
#CFLAGS = -Wall -g -DWARNINGS -std=gnu99 -pedantic
#CFLAGS = -O -funroll-loops  -std=gnu99 -DNDEBUG
LD = ld
LDFLAGS =  -lm 
INCFLAGS = -I/System/Library/Frameworks/vecLib.framework/Headers/
LIBDICTDIR = ../libdict-0.2.1
#LDFLAGS = -lm

objects = like.o tree.o data.o rng.o model.o  bases.o codonmodel.o gencode.o utility.o matrix.o optimize_new.o spinner.o options.o tree_data.o linemin.o gamma.o statistics.o mystring.o nucmodel.o root.o vec.o brent.o


Slr: slr.o $(objects) libdict
	gcc  -o ../bin/$@ slr.o $(objects) $(LDFLAGS) $(LIBDICTDIR)/libdict.a -lblas -llapack 

libdict:
	cd $(LIBDICTDIR) && make
.c.o:
	$(CC) $(CFLAGS) $(INCFLAGS) -o $@ -c $<

clean:
	rm -f *~ *.o *.obj *.obj *.exe *.lib && cd $(LIBDICTDIR) && make clean

