bude_options = --no-ieee-float --fast
vec_options = -g --mllvm "--pass-remarks-analysis=loop-vectorize --pass-remarks=loop-vectorize" 
extra_ops = --copy-propagation --dead-code-elimination --loop-invariant-code-motion --inline-iterators --privatization --optimize-range-iteration

bude:
	chpl $(bude_options) Bude.chpl -o bude

vec:
	chpl $(bude_options) $(vec_options) Bude.chpl -o bude 2>vec_msg

.PHONY: clean
clean:
	rm -f bude

