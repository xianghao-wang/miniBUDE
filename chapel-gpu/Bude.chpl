module Bude {
  use Time;
  use IO;
  use GPU;
  use CTypes;
  use ArgumentParser;

  param DEFAULT_ITERS = 8;
  param DEFAULT_NPOSES = 65536;
  param REF_NPOSES = 65536;
  param DATA_DIR = "../data/bm1";
  param FILE_LIGAND = "/ligand.in";
  param FILE_PROTEIN = "/protein.in";
  param FILE_FORCEFIELD = "/forcefield.in";
  param FILE_POSES = "/poses.in";
  param FILE_REF_ENERGIES = "/ref_energies.out";

  record atom {
    var x, y, z: real(32);
    var aType: int(32);
  }

  record ffParams {
    var hbtype: int(32);
    var radius: real(32);
    var hphb: real(32);
    var elsc: real(32);
  }

  class params {
    var deckDir: string;
    var iterations: int;
    var natlig, natpro, ntypes, nposes: int;

    var proteinDomain: domain(1);
    var ligandDomain: domain(1);
    var forcefieldDomain: domain(1);
    var posesDomain: domain(2);

    var protein: [proteinDomain] atom; 
    var ligand: [ligandDomain] atom;
    var forcefield: [forcefieldDomain] ffParams;
    var poses: [posesDomain] real(32);

    proc init() {  }

    proc load(args: [] string) throws {
      // Parsing command-line parameters
      var parser = new argumentParser();
      var iterationsArg = parser.addOption(name="iterations"
        , opts=["-i", "--iteartions"]
        , defaultValue=DEFAULT_ITERS: string
        , valueName="I"
        , help="Repeat kernel I times (default: " + DEFAULT_ITERS: string + ")");
      var deckArg = parser.addOption(name="deck"
        , opts=["--deck"]
        , defaultValue=DATA_DIR
        , valueName="DECK"
        , help="Use the DECK directory as input deck (default: " + DATA_DIR + ")");
      var numposesArg = parser.addOption(name="numposes"
        , opts=["-n", "--numposes"]
        , defaultValue=DEFAULT_NPOSES: string
        , valueName="N"
        , help="Compute energies for N poses (default: " + DEFAULT_NPOSES: string + ")");
      parser.parseArgs(args);

      // Store these parameters
      try {
        this.iterations = iterationsArg.value(): int;
        if (this.iterations < 0) then throw new Error();
      } catch {
        writeln("Invalid number of iterations");
      }

      try {
        this.nposes = numposesArg.value(): int;
        if (this.nposes < 0) then throw new Error();
      } catch {
        writeln("Invalid number of poses");
      }
      
      this.deckDir = deckArg.value(); 
    }
  }

  proc main(args: [] string) {
    var context = new unmanaged params(); // Unmanaged preventing deallocating class during computation
    try! context.load(args);

    delete context;
  }
}