EXEC   = sedov_taylor

OPTIMIZE =  -O2

OPT   = $(OPTIMIZE)

OBJS   = main.o routines.o sedov_taylor.o rk_int.o bisection_root_finding.o

CC     = g++

INCL   = sedov_taylor.h routines.h rk_int.h bisection_root_finding.h

LIBS   = -lm -lgsl -lgslcblas

CFLAGS = $(OPT)

$(EXEC): $(OBJS) 
	 $(CC) $(OBJS) $(LIBS) -o $(EXEC)  

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

