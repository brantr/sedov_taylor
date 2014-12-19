EXEC   = sedov_taylor

OPTIMIZE =  -O2

OPT   = $(OPTIMIZE)

OBJS   = main.o routines.o interpolation.o sedov_taylor.o rk_int.o bisection_root_finding.o

CC     = g++

INCL   = routines.h interpolation.h sedov_taylor.h rk_int.h bisection_root_finding.h

LIBS   = -lm -lgsl -lgslcblas

CFLAGS = $(OPT)

$(EXEC): $(OBJS) 
	 $(CC) $(OBJS) $(LIBS) -o $(EXEC)  

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

