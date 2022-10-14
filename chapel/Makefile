# bude_options = --no-ieee-float --fast -s WGSIZE=256 --mllvm="-force-vector-interleave=4"
# vec_options = -g --mllvm "--pass-remarks-analysis=loop-vectorize --pass-remarks=loop-vectorize -force-vector-interleave=4" 
# extra_ops = --copy-propagation --dead-code-elimination --loop-invariant-code-motion --inline-iterators --privatization --optimize-range-iteration

# bude:
# 	chpl $(bude_options) Bude.chpl -o bude

# vec:
# 	chpl $(bude_options) $(vec_options) Bude.chpl -o bude 2>vec_msg

# .PHONY: clean
# clean:
# 	rm -f bude

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:

MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# -------
MACHINE = $(shell uname -m)

ifeq ($(MACHINE), x86_64)
ARCHFLAG = march

# On Skylake and Cascade Lake, 256-bit vectors are used by default, but 512-bit is more beneficial
AVX512   =  $(if $(filter skylake% cascadelake% %avx512, $(ARCH)),yes,)

# Intel platforms benefit more from unrolling, specially Skylake and later
WGSIZE = 256

else
# The target CPU is specificed differently on x86 and on aarch64
# https://community.arm.com/developer/tools-software/tools/b/tools-software-ides-blog/posts/compiler-flags-across-architectures-march-mtune-and-mcpu
ARCHFLAG = mcpu
endif

WGSIZE = 64
ARCH   = native
AVX512 =
EXTRA_FLAGS =

CC = chpl

CFLAGS_CHAPEL = --warnings --fast --no-ieee-float --mllvm "-force-vector-interleave=4"
CFLAGS = $(CFLAGS_CHAPEL) -s WGSIZE=$(WGSIZE)


# -------
EXE = bude

all: $(EXE)

$(EXE): Bude.chpl Atom.chpl Configuration.chpl Context.chpl FFParams.chpl Helper.chpl VecPoseInner.chpl
	$(CC) $(CFLAGS) Bude.chpl -o $@ $(EXTRA_FLAGS)

clean:
	rm -f $(EXE)
