# miniBUDE Chapel with GPU Support

This is an implementation of miniBUDE using Chapel with GPU support.

## Building

### Prerequisites
 * Chapel >= 1.30

### Compile
The compiling is simply using Makefile.
```
make 
```

## Running
```
./bude -w=64
```
This implementation sets the number of blocks on GPU by passing `--wgsize`. The `-n` and `-i` parameters are available. Run bude `-h` for a help message.