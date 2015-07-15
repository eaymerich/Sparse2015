#----------------------------------------------------------
# This section detects if we are compiling on Windows.
# If not, we assume a *nix environment.
# Code adapted from:
#	http://stackoverflow.com/questions/4058840/makefile-that-distincts-between-windows-and-unix-like-systems
#
ifeq ($(OS),Windows_NT)
	RM=del /Q
	EXT=.exe
else
#	ifeq ($(shell,uname),Linux)
	RM=rm -f
	EXT=
#	endif
endif

TARGET=diag$(EXT) csr$(EXT) tjds$(EXT) oski$(EXT) csc$(EXT) sss$(EXT)
FLAGS= -O3
COMMON_DEPS=mmio.c mmio.h clock.h coo.h

all: $(TARGET)

csc$(EXT): csc.c csc.h $(COMMON_DEPS)
	gcc $(FLAGS) -o csc$(EXT) csc.c mmio.c

oski$(EXT): oski.c oski.h $(COMMON_DEPS)
	gcc $(FLAGS) -o oski$(EXT) oski.c mmio.c

csr$(EXT): csr.c csr.h $(COMMON_DEPS)
	gcc $(FLAGS) -o csr$(EXT) csr.c mmio.c
	
tjds$(EXT): tjds.c tjds.h $(COMMON_DEPS)
	gcc $(FLAGS) -o tjds$(EXT) tjds.c mmio.c

sss$(EXT): sss.c sss.h $(COMMON_DEPS)
	gcc $(FLAGS) -o sss$(EXT) sss.c mmio.c
	
diag$(EXT): diag.c mmio.c mmio.h coo.h
	gcc $(FLAGS) -o diag$(EXT) diag.c mmio.c
	
clean:
	$(RM) $(TARGET)
