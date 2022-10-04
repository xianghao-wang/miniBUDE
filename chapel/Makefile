bude_options = --no-ieee-float --fast
vec_options = -g --mllvm "--pass-remarks-analysis=loop-vectorize --pass-remarks=loop-vectorize" 

bude:
	chpl $(bude_options) Bude.chpl -o bude

vec:
	chpl $(bude_options) $(vec_options) Bude.chpl -o bude 2>vec_msg

.PHONY: clean
clean:
	rm -f bude

