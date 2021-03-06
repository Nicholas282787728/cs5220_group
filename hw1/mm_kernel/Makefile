# ---
# Platform-dependent configuration
#
# If you have multiple platform-dependent configuration options that you want
# to play with, you can put them in an appropriately-named Makefile.in.
# For example, the default setup has a Makefile.in.icc and Makefile.in.gcc.

PLATFORM=icc
include Makefile.in.$(PLATFORM)

.PHONY:	all
all:	ktimer ktimer-simple ktimer-doublerainbow ktimer-doublerainbow_4x4 ktimer-bindel

.PHONY: run
run:
	make clean
	make
	csub ./ktimer

# ---
# Rules to build the drivers

ktimer: $(OBJS) kdgemm_doublerainbow.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)

ktimer-%: $(OBJS) kdgemm_%.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS)
    
# --
# Rules to build object files

ktimer.o: ktimer.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $<

%.o: %.c
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(CPPFLAGS) $<

%.s: %.c
	$(CC) -S $(CFLAGS) $(OPTFLAGS) $(CPPFLAGS) $<

%.o: %.f
	$(FC) -c $(FFLAGS) $<

.PHONY: clean
clean:
	rm -f ktimer ktimer-* *.o
