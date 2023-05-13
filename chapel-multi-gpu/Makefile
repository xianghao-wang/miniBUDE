SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# -------

OPTIONS = 
CHPL_FLAGS = --fast --no-ieee-float $(OPTIONS)
CHPL = chpl

# -------
EXE = bude
.PHONY: all $(EXE) clean

all: $(EXE)

$(EXE): Bude.chpl
	$(CHPL) $(CHPL_FLAGS) Bude.chpl -o $@

clean:
	rm -f $(EXE)