module Bude {
    use Context;

    var params: context;

    proc main(args: [] string) {
      // Load context
      params = new context(args);
      params.load();

      // Show meta-information
      writeln("");
      writeln("Poses     :", params.nposes);
      writeln("Iterations:", params.iterations);
      writeln("Ligands   :", params.natlig);
      writeln("Proteins  :", params.natpro);
      writeln("Deck      :", params.deckDir);
    }
}