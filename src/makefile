# specify the include directory (containing all .h files)
IDIR =../include
# specify compiler ... 
CC=g++
# and compiler flags (-I$(IDIR): include the directory; -O3: pretty fast optimization).
CFLAGS=-I$(IDIR) -I"C:/Program Files/eigen-3.3.7" -I"C:/Program Files/boost_1_72_0" -O3

# specify object directory to hide them from plain sight
ODIR=obj

# If the project depends on libraries specify them here (used in the linker)
LIBS=-lm 

# All user defined headers 
_DEPS = Particle.h PBC.h RateEngine.h Site.h KmcRun.h NextEventList.h OutputManager.h EnumNames.h RandomEngine.h Site.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# All objects that need to be created (for every .c-file we include a .o-file here)
_OBJ = main.o RateEngine.o Site.o NextEventList.o OutputManager.o RandomEngine.o KmcRun.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# Check for dependecies
$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# Compile the program
kmc: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS) -DNDEBUG

# If you need a fresh build call clean
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 