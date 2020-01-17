# specify the include directory (containing all .h files)
IDIR =../include
# specify compiler ... 
CC=g++
# and compiler flags (-I$(IDIR): include the directory; -O3: pretty fast optimization).
CFLAGS=-I$(IDIR) -I "C:/Program Files/eigen-3.3.7"

# specify object directory to hide them from plain sight
ODIR=obj

# If the project depends on libraries specify them here (used in the linker)
LIBS=-lm 

# All user defined headers 
_DEPS = Particle.h PBC.h RateEngine.h Site.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# All objects that need to be created (for every .c-file we include a .o-file here)
_OBJ = main.o Particle.o PBC.o RateEngine.o Site.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# Check for dependecies
$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# Compile the program
kmc: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

# If you need a fresh build call clean
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 