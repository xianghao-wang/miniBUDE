module Bude {
  use Time;
  use IO;
  use GPU;
  use CTypes;
  use ArgumentParser;
  use Path;
  use CKernel;
  use GKernel;
  use ChplConfig;

  // Program context parameters
  param DEFAULT_ITERS = 8;
  param DEFAULT_NPOSES = 65536;
  param DEFAULT_WGSIZE = 64;
  param REF_NPOSES = 65536;
  param DATA_DIR = "../data/bm1";
  param FILE_LIGAND = "/ligand.in";
  param FILE_PROTEIN = "/protein.in";
  param FILE_FORCEFIELD = "/forcefield.in";
  param FILE_POSES = "/poses.in";
  param FILE_REF_ENERGIES = "/ref_energies.out";

  // Energy evaluation parameters
  const CNSTNT: real(32) = 45.0;
  const HBTYPE_F: real(32) = 70.0;
  const HBTYPE_E: real(32) = 69.0;
  const HARDNESS: real(32) = 38.0;
  const NPNPDIST: real(32) = 5.5;
  const NPPDIST: real(32) = 1.0;

  // Configurations
  config param NUM_TD_PER_THREAD: int(32) = 4; // Work per core

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
    var iterations: int(32);
    var natlig, natpro, ntypes, nposes: int(32);
    var wgsize: int;

    var proteinDomain: domain(1,int(32));
    var ligandDomain: domain(1,int(32));
    var forcefieldDomain: domain(1,int(32));
    var posesDomain: domain(2,int(32));

    var protein: [proteinDomain] atom; 
    var ligand: [ligandDomain] atom;
    var forcefield: [forcefieldDomain] ffParams;
    var poses: [posesDomain] real(32);

    proc init() { }

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
      var wgsizeArg = parser.addOption(name="wgsize"
        , opts=["-w", "--wgsize"]
        , defaultValue=DEFAULT_WGSIZE: string
        , valueName="W"
        , help="Set the number of blocks to W (default: " + DEFAULT_WGSIZE: string + ")");
      parser.parseArgs(args);

      // Store these parameters
      try {
        this.iterations = iterationsArg.value(): int(32);
        if (this.iterations < 0) then throw new Error();
      } catch {
        writeln("Invalid number of iterations");
        exit(1);
      }

      try {
        this.nposes = numposesArg.value(): int(32);
        if (this.nposes < 0) then throw new Error();
      } catch {
        writeln("Invalid number of poses");
        exit(1);
      }

      try {
        this.wgsize = wgsizeArg.value(): int;
        if (this.wgsize < 0) then throw new Error();
      } catch {
        writeln("Invalid number of blocks");
        exit(1);
      }
      
      this.deckDir = deckArg.value(); 
      
      // Load data
      var length: int;
      var aFile: file;
      var reader: fileReader();
      
      // Read ligand
      aFile = openFile(this.deckDir + FILE_LIGAND, length);
      this.natlig = (length / c_sizeof(atom)): int(32);
      this.ligandDomain = {0..<this.natlig};
      reader = aFile.reader(kind=iokind.native, region=0..);
      reader.read(this.ligand);
      reader.close(); aFile.close();

      // Read protein
      aFile = openFile(this.deckDir + FILE_PROTEIN, length);
      this.natpro = (length / c_sizeof(atom)): int(32);
      this.proteinDomain = {0..<this.natpro};
      reader = aFile.reader(kind=iokind.native, region=0..);
      reader.read(this.protein);
      reader.close(); aFile.close();

      // Read force field
      aFile = openFile(this.deckDir + FILE_FORCEFIELD, length);
      this.ntypes = (length / c_sizeof(ffParams)): int(32);
      this.forcefieldDomain = {0..<this.ntypes};
      reader = aFile.reader(kind=iokind.native, region=0..);
      reader.read(this.forcefield);
      reader.close(); aFile.close();

      // Read poses
      this.posesDomain = {0..<6:int(32), 0..<this.nposes};
      aFile = openFile(this.deckDir + FILE_POSES, length);
      const reader_tmp =  aFile.reader(kind=iokind.native, region=0.., locking=false);
      var available = (length / (6 * c_sizeof(real(32)): int)): int(32);
      var current = 0:int(32);
      while (current < this.nposes) {
        var fetch = this.nposes - current;
        if (fetch > available) then fetch = available;

        for i in posesDomain.dim(0) {
          const base = i*available*c_sizeof(real(32)):int;
          const amount = fetch*c_sizeof(real(32)):int;
          reader_tmp.seek(region=base..base+amount);
          reader_tmp.read(this.poses(i, current..));
        }
        current += fetch;
      }
      this.nposes = current;
      reader.close(); aFile.close();
    }
  }

  var context: params = new params();

  proc main(args: [] string) {
    try! context.load(args);

    // Show meta-information
    writeln("");
    writeln("Poses     : ", context.nposes);
    writeln("Iterations: ", context.iterations);
    writeln("Ligands   : ", context.natlig);
    writeln("Proteins  : ", context.natpro);
    writeln("Deck      : ", context.deckDir);

    // Compute
    var energies: [0..<context.nposes] real(32);
    compute(energies);

    // Validate
    var length: int;
    const ref_energies = openFile(context.deckDir+FILE_REF_ENERGIES, length);
    var e: real(32);
    var diff: real(32);
    var maxdiff: real(32) = -100.0;
    var n_ref_poses = context.nposes;
    if (context.nposes > REF_NPOSES) {
      writeln("Only validating the first ", REF_NPOSES, " poses");
      n_ref_poses = REF_NPOSES;
    }
    var reader = try! ref_energies.reader();
    for i in 0..<n_ref_poses {
      try! reader.read(e);
      if (abs(e) < 1.0 && abs(energies(i)) < 1.0) {
        continue;
      }
      diff = abs(e - energies(i)) / e;
      if (diff > maxdiff) {
        maxdiff = diff;
      }
    }
    writef("\nLargest difference was %{.###}%%.\n\n", 100 * maxdiff);
  }

  proc compute(results: [] real(32)) {

    if (CHPL_GPU == "nvidia") {
      writeln("\nRunning Chapel on ", here.gpus.size, (if here.gpus.size > 1 then " GPUs" else " GPU"));
      gkernel(context, results);
    } else if (CHPL_GPU == "none") {
      writeln("\nRunning Chapel");
      ckernel(context, results);
    } else {
      writeln("\n" + CHPL_GPU + " is not supported.");
      exit(1);
    }
    
  } // main

  

  proc openFile(fileName: string, ref length: int): file {
    try {
      const aFile = open(fileName, ioMode.r);
      length = aFile.size;
      return aFile;
    } catch {
      try! stderr.writeln("Failed to open '", fileName, "'");
      exit(0);
    }
  }

  proc timestampMS() {
    return timeSinceEpoch().totalSeconds() * 1000;
  }

  proc printTimings(timeMS: real(64)) {
    const ms = timeMS / context.iterations;
    const runtime = ms * 1e-3;

    const ops_per_wg = NUM_TD_PER_THREAD * 27 + context.natlig * (2 + NUM_TD_PER_THREAD * 18 + context.natpro * (10 + NUM_TD_PER_THREAD * 30)) + NUM_TD_PER_THREAD;
    const total_ops = ops_per_wg * (context.nposes / NUM_TD_PER_THREAD);
    const flops = total_ops / runtime;
    const gflops = flops / 1e9;

    const interactions = 1.0 * context.nposes * context.natlig * context.natpro;
    const interactions_per_sec = interactions / runtime;

    // Print stats
    writef("- Total time:     %7.3dr ms\n", timeMS);
    writef("- Average time:   %7.3dr ms\n", ms);
    writef("- Interactions/s: %7.3dr billion\n", (interactions_per_sec / 1e9));
    writef("- GFLOP/s:        %7.3dr\n", gflops);
  }
}