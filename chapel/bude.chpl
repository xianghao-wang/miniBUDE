use Config;

record ParameterConfig {
  const natlig: int;
  const natpro: int;
  const ntypes: int;
  const nposes: int;
  
  // TODO: Pointer in Chapel

  const deckDir: string;
  const iterations: int;

  proc loadParameters(args: [] string) {

  }
}

// Read command line arguments
// config const iterations = DEFAULT_ITERS
// config const numposes = DEFAULT_NPOSES
// config const 

proc main(args: [] string) {
  writeln(DATA_DIR);
}
