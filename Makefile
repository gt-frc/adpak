FCOMP = gfortran
FCFLAGS = 
PROGRAM =  adpak
SRCS = dutilits.f adpak.f main.f
OBJECTS = $(SRCS:.f=.o)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^

%.o: %.f08
	$(FCOMP) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM)
