# miniBUDE Chapel

This is an implementation of miniBUDE using Chapel, targetting either host or GPUs.

## Building

### Prerequisites
 * Chapel >= 1.30

### Compile
The compiling is simply using Makefile.
```
make 
```

Specifying work group size per thread
```
make WGSIZE=64
```

## Running
```
./bude -w=64 --ngpu=2
```
This implementation sets the number of blocks on GPU by passing `--wgsize` or `-w` and the number of devices by `--ngpu`. The `-n` and `-i` parameters are available. Run bude `-h` for a help message.
* NOTE: `WGSIZE` for building and `--wgsize` for running are different 2 concepts.

```
./bude --host
```
This make the program target host.