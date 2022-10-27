SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:

MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

WGSIZE = 64
EXTRA_FLAGS =

MACHINE = $(shell uname -m)

ifeq ($(MACHINE), x86_64)
# Intel platforms benefit more from unrolling, specially Skylake and later
WGSIZE = 256
endif

VI = 4
ifeq ($(shell expr $(WGSIZE) \< 32), 1)
   VI = 1
endif

CFLAGS_CHAPEL = --warnings --fast --no-ieee-float --mllvm "-force-vector-interleave=$(VI)"
CFLAGS = $(CFLAGS_CHAPEL) -s WGSIZE=$(WGSIZE)
CHPL = chpl

# -------
EXE = bude


.PHONY: all $(EXE) clean

all: $(EXE)

$(EXE): Bude.chpl Atom.chpl Configuration.chpl Context.chpl FFParams.chpl Helper.chpl VecPoseInner.chpl
	$(CHPL) $(CFLAGS) Bude.chpl -o $@ $(EXTRA_FLAGS)

clean:
	rm -f $(EXE)
