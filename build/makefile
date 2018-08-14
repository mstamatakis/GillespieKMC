COMPILER = gfortran

OUTPUT := ./GillespieKMC.x

FLAGS := -O2 -ffree-line-length-1024
LIBS := 

SOURCE = constants_module parser_module process_details_module execution_queue_class execution_queue_binary_heap_class

EXECUTABLE = on-lattice

OBJS = $(SOURCE:%=%.o)
ALL_OBJS = $(SOURCE:%=%.o) $(EXECUTABLE:%=%.o)

.PHONY : all
all: $(ALL_OBJS) $(EXECUTABLE)

.PHONY : clean
clean:
	@echo "Cleaning directory"
	@\rm -rf *.o *.mod *.exe $(OUTPUT)

%.o: ../%.f90
	@echo "Compiling: " $*
	$(COMPILER) -c  ${FLAGS} $< -o $@

%:  %.o $(ALL_OBJS)
	@echo "    Linking program:   " $@
	$(COMPILER) $(FLAGS) $(ALL_OBJS) $(LIBS) -o $(OUTPUT)

