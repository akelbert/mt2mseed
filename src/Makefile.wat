#
# Wmake File - for Watcom's wmake
# Use 'wmake -f Makefile.wat'

.BEFORE
	@set INCLUDE=.;$(%watcom)\H;$(%watcom)\H\NT
	@set LIB=.;$(%watcom)\LIB386

cc     = wcc386
cflags = -zq
lflags = OPT quiet OPT map LIBRARY ..\libmseed\libmseed.lib
cvars  = $+$(cvars)$- -DWIN32

BIN = ..\nims2mseed.exe

INCS = -I..\libmseed

all: $(BIN)

$(BIN):	nims2mseed.obj readNIMSbin.obj
	wlink $(lflags) name $(BIN) file {nims2mseed.obj readNIMSbin.obj}

# Source dependencies:
nims2mseed.obj:		readNIMSbin.h nims2mseed.c
readNIMSbin.obj:	readNIMSbin.h readNIMSbin.c

# How to compile sources:
.c.obj:
	$(cc) $(cflags) $(cvars) $(INCS) $[@ -fo=$@

# Clean-up directives:
clean:	.SYMBOLIC
	del *.obj *.map $(BIN)
