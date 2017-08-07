NTHREADS = 2

# Compiler and its basic options:
CC = g++
BASIOPTS = -march=native -m64 -cpp -fmax-errors=1

# Warning-related options:
BASIOPTS += -Wall -Wno-maybe-uninitialized -Wno-attributes -Wno-unused-result

# Option FAST (default):
FAST_OPTS = -Ofast

# Option DEBUG:
DEBUG_OPTS = --coverage -g -Og -fstack-check -fsanitize=null -fbounds-check -Ddebug

# Checks chosen option:
ifeq ($(DEBUG), 1)
  OPTS = $(BASIOPTS) $(DEBUG_OPTS)
else
  OPTS = $(BASIOPTS) $(FAST_OPTS)
endif

SRCDIR = ./src
OBJDIR = $(SRCDIR)/obj
BINDIR = ./bin
TSTDIR = ./test

LIBS = -lemdee -lm

.PHONY: all clean test

.DEFAULT_GOAL := all

# Phony targets:

all: $(BINDIR)/simtop

clean:
	rm -rf $(OBJDIR) $(BINDIR) *.gcda *.gcno

test: $(BINDIR)/simtop
	$(BINDIR)/simtop $(NTHREADS) $(TSTDIR)/tip3p_sample

$(BINDIR)/simtop: $(OBJDIR)/simtop.o $(OBJDIR)/tops.o
	mkdir -p $(BINDIR)
	$(CC) $(OPTS) -o $@ $^ $(LIBS)

$(OBJDIR)/simtop.o: $(SRCDIR)/simtop.cc $(OBJDIR)/tops.o
	$(CC) $(OPTS) -c -o $@ $<

$(OBJDIR)/tops.o: $(SRCDIR)/tops.cc $(SRCDIR)/tops.h $(SRCDIR)/vecmat3.h
	mkdir -p $(OBJDIR)
	$(CC) $(OPTS) -c -o $@ $<
