HOST = $(shell hostname)
BANG   =  $(shell expr match `hostname` ccom-bang)
BANG-COMPUTE   =  $(shell expr match `hostname` compute)



ifneq ($(BANG), 0)
PUB     = /share/class/public/cse160-wi13
include $(PUB)/Arch/arch.gnu-4.7_c++11.generic
else
ifneq ($(BANG-COMPUTE), 0)
PUB     = /share/class/public/cse160-wi13
include $(PUB)/Arch/arch.gnu-4.7_c++11.generic
endif
endif

#
# Also add debugging output controlled via
# conditional compilation with the DEBUG macro
# Add symbol table information for gdb/cachegrind
ifeq ($(debug), 1)
        CFLAGS += -g -DDEBUG
        LDFLAGS += -g
        C++FLAGS += -g -DDEBUG
endif   


# Add symbol table information for gprof
ifeq ($(gprof), 1)
        CFLAGS += -g -pg
        C++FLAGS += -g -pg
        LDFLAGS += -g -pg
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




app:		apf

OBJECTS = apf.o solve.o Plotting.o cmdLine.o util.o Report.o
ifneq ($(mpi),1)
OBJECTS += Timer.o
endif

apf:	        $(OBJECTS) 
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJECTS)  $(LDLIBS)

.PHONY: clean
clean:	
	$(RM) *.o apf;
	$(RM) core;
