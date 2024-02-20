MAKEFLAGS += --no-builtin-rules --no-builtin-variables

FC := pgf90 -acc
FCFLAGS := -acc
LD := $(FC)

RM := rm -f

SRCS := utilities.f90 \
		porosity.f90 \
		surfacer.f90 \
		main.f90
PROG := namac
OBJS := $(addsuffix .o, $(SRCS))

.PHONY: all clean
all: $(PROG)

$(PROG): $(OBJS)
	$(LD) -o $@ $^

$(OBJS): %.o: %
	$(FC) -c -o $@ $<

$(OBJS): $(MAKEFILE_LIST)

clean:
	$(RM) $(filter %.o, $(OBJS)) $(wildcard *.mod) $(PROG$(PROG): $(OBJS)
