SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:

MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

WGSIZE = 64

MACHINE = $(shell uname -m)

ifeq ($(MACHINE), x86_64)
# Intel platforms benefit more from unrolling, specially Skylake and later
WGSIZE = 256
endif

VI = 4
ifeq ($(shell expr $(WGSIZE) \< 32), 1)
   VI = 1
endif

CHPLFLAGS = --warnings --fast --no-ieee-float --mllvm "-force-vector-interleave=$(VI)" -s WGSIZE=$(WGSIZE)
CHPL = chpl

# -------
EXE = bude


.PHONY: all $(EXE) clean

all: $(EXE)

$(EXE): Bude.chpl
	$(CHPL) $(CHPLFLAGS) Bude.chpl -o $@

clean:
	rm -f $(EXE)
