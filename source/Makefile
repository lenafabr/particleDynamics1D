OBJS = keys.o readkey.o genutil.o mt19937.o inputparams.o particleutils.o partdynamics.o fusionutils.o
TESTOBJS = keys.o readkey.o genutil.o mt19937.o inputparams.o particleutils.o 

FC = gfortran
FLAGS = -O3 
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid -fcheck=all -Wno-tabs
LIBS =  -lblas -llapack
OSGFLAGS = -static
OSGLIBS =  -lblas -llapack -lgfortran -lm

# define directories
PDIR = ./profile
ODIR = .

# executable name
EXE = ../partdynamics1D.exe
OSGEXE=../partdynamics1D_osg.exe
TESTEXE = ../testrun.exe
PROFEXE = ../partdynamics1D_prof.exe

.SUFFIXES: .o .f90 

# general suffix rules

$(PDIR)/%.o : tools/%.f90
	mkdir -p $(@D)
	$(FC) $(FLAGS) $(DEBUGFLAGS) -pg -c -o $@ $<

$(ODIR)/%.o : tools/%.f90
	mkdir -p $(@D)
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c -o $@ $<

$(PDIR)/%.o : %.f90
	mkdir -p $(@D)
	$(FC) $(FLAGS) $(DEBUGFLAGS) -pg -c -o $@ $<

$(ODIR)/%.o : %.f90
	mkdir -p $(@D)
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c -o $@ $<

default: $(addprefix $(ODIR)/,$(OBJS) main.o) # recipe for ODIR/OBJS
	$(FC) $(FLAGS) $(DEBUGFLAGS) $^ -o $(EXE) $(LIBS)

profile: $(addprefix $(PDIR)/,$(OBJS) main.o) # recipe for PDIR/OBJS
	$(FC) $(FLAGS) $(DEBUGFLAGS) -pg $^ -o $(PROFEXE) $(LIBS)

osg: $(OBJS) main.o
	$(FC) $(FLAGS) $(OSGFLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(OSGEXE) $(OSGLIBS)

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

clean: 
	rm *.o *.mod $(ODIR)/*.o $(ODIR)/*.mod $(PDIR)/*.o $(PDIR)/*.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
$(ODIR)/readkey.o: $(addprefix $(ODIR)/,inputparams.o keys.o genutil.o)
$(ODIR)/genutil.o: $(addprefix $(ODIR)/,mt19937.o)
$(ODIR)/partdynamics.o: $(addprefix $(ODIR)/,keys.o mt19937.o particleutils.o fusionutils.o)
$(ODIR)/fusionutils.o: $(addprefix $(ODIR)/,keys.o mt19937.o particleutils.o)
$(ODIR)/particleutils.o: $(addprefix $(ODIR)/,keys.o genutil.o mt19937.o)

$(PDIR)/readkey.o: $(addprefix $(ODIR)/,inputparams.o keys.o genutil.o)
$(PDIR)/genutil.o: $(addprefix $(ODIR)/,mt19937.o)
$(PDIR)/partdynamics.o: $(addprefix $(ODIR)/,keys.o mt19937.o particleutils.o fusionutils.o)
$(PDIR)/fusionutils.o: $(addprefix $(ODIR)/,keys.o mt19937.o particleutils.o)
$(PDIR)/particleutils.o: $(addprefix $(ODIR)/,keys.o genutil.o mt19937.o)