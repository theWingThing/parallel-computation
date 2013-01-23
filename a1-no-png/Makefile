PUB = .
# include $(PUB)/Arch/arch.gnu-4.7_c++11
include $(PUB)/Arch/arch.gnu.pthread

# XTRAFLAGS += -fsource-asm -S

#
# You may override the optimization level
# set in the "arch" file by uncommenting the
# next line and modifying the optimization level as needed
#
# OPTIMIZATION    =  -O0

# Add symbol table information for gdb/cachegrind
ifeq ($(dbg), 1)
        CFLAGS += -g
        LDFLAGS += -g
 C++FLAGS += -g
endif


# Add symbol table information for gprof
ifeq ($(gprof), 1)
        CFLAGS += -g -pg
        C++FLAGS += -g -pg
        LDFLAGS += -g -pg
endif



# INCLUDES += -I$(PUB)/include
# INCLUDES += -DNO_FREETYPE

# LDLIBS += -L$(PUB)/lib 

UTIL    = Timer.o util.o splot.o

OBJ     = main.o pmdb.o smdb.o Report.o cmdLine.o $(UTIL)

mdb:		$(OBJ)
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

.PHONY: clean
clean:
	$(RM) *.o mdb;
	$(RM) core;
	$(RM) cscope.*
	$(RM) tags

cleanall:
	$(RM) *.o mdb *.png *.rgb log;
	$(RM) core;
	$(RM) cscope.*
	$(RM) tags
