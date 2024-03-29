HOST = $(shell hostname)
BANG   =  $(shell expr match `hostname` ccom-bang)
BANG-COMPUTE   =  $(shell expr match `hostname` compute)


ifneq ($(BANG), 0)
PUB = /share/class/public/cse160-wi13
else
ifneq ($(BANG-COMPUTE), 0)
PUB = /share/class/public/cse160-wi13
else
PUB = /Path/To/Your/Own/machine/here
endif
endif
include $(PUB)/Arch/arch.gnu-4.7_c++11

# XTRAFLAGS += -fsource-asm -S

# OpenMP options only
# If you want to add dyamic scheduling, set dyn=1 on make command line
# Also enables compilation with openmp
ifeq ($(dyn),1)
C++FLAGS    +=  -fopenmp -DDYN
CFLAGS      +=  -fopenmp -DDYN
endif

# Sets the chunk size for dynamic scheduling
# Also enables compilation with openmp, and dynamic scheduling
ifdef chunk
C++FLAGS    +=  -fopenmp -DDYN -DCHUNK=$(chunk)
CFLAGS      +=  -fopenmp -DDYN -DCHUNK=$(chunk)
endif
# End openMP options


# If you want to add symbol table information for gdb/cachegrind
# specify gdb=1 on the "make" command line
ifeq ($(gdb), 1)
        CFLAGS += -g
        LDFLAGS += -g
        C++FLAGS += -g
endif   


# If you want to add symbol table information for gprof
# specify gprof=1 on the "make" command line
ifeq ($(gprof), 1)
        CFLAGS += -g -pg
        C++FLAGS += -g -pg
        LDFLAGS += -g -pg
endif

# If you want to compile for single precision,
# specify single=1 on the "make" command line
ifeq ($(single), 1)
else
    C++FLAGS += -D_DOUBLE
    CFLAGS += -D_DOUBLE
endif



# If you want to use restrict pointers, make restrict=1
# This applies to the hand code version
ifeq ($(restrict), 1)
    C++FLAGS += -D__RESTRICT
    CFLAGS += -D__RESTRICT
ifneq ($(CARVER), 0)
    C++FLAGS += -restrict
    CFLAGS += -restrict
endif
endif

WARNINGS        += -Wall -pedantic
C++FLAGS 	+= $(WARNINGS)

C++LINK		= $(C++)

C++FLAGS += $(OPTIMIZATION)

TARGETS = nbody

all:	$(TARGETS)

SOBJ = nbody.o common.o cmdLine.o particlesHelper.o Plotting.o Box.o Bins.o

nbody: $(SOBJ)
	$(C++LINK) $(C++FLAGS) -o $@ $(LDLIBS) $(SOBJ)

# ===========================
#
.PHONY: clean
clean:
	rm -f *.o $(TARGETS)

