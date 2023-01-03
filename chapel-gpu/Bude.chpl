module Bude {
  use IO;
  use Time;
  use AutoMath;

  config param WGSIZE = 4,
    DEFAULT_ITERS = 8,
    DEFAULT_NPOSES = 65536,
    REF_NPOSES = 65536,
    DATA_DIR = "../data/bm1",
    FILE_LIGAND = "/ligand.in",
    FILE_PROTEIN = "/protein.in",
    FILE_FORCEFIELD = "/forcefield.in",
    FILE_POSES = "/poses.in",
    FILE_REF_ENERGIES = "/ref_energies.out",
    ATOM_SIZE = 16,
    FFPARAMS_SIZE = 16;

  // Energy evaluation parameters
  const CNSTNT: real(32) = 45.0;
  const HBTYPE_F: real(32) = 70.0;
  const HBTYPE_E: real(32) = 69.0;
  const HARDNESS: real(32) = 38.0;
  const NPNPDIST: real(32) = 5.5;
  const NPPDIST: real(32) = 1.0;

  const WORK_GROUP = 0..<WGSIZE;

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

  record context {
    var iterations: int = DEFAULT_ITERS;
    var nposes: int = DEFAULT_NPOSES;
    var deckDir: string = DATA_DIR;

    // Domains for arrays
    var natlig: int;
    var ligandDom: domain(1);
    var natpro: int;
    var proteinDom: domain(1);
    var ntypes: int;
    var forcefieldDom: domain(1);
    var posesDom: domain(2);
    
    var ligand: [ligandDom] atom;
    var protein: [proteinDom] atom;
    var forcefield: [forcefieldDom] ffParams;
    var poses: [posesDom] real(32);

    proc init() { }

    /* Initialize context with arguments */
    proc init(args: [] string) { }

    /* Load input data */
    proc load() { }

    /* Helper functions */
    
    /* Load data from file to array */
    proc loadData(aFile: file, ref A: [] ?t, size: int) { }

    /* Load data with length of the given bytes */
    proc loadDataPiece(aFile: file, ref A: ?t, base: int, offset: int) {  }
  
    /* Convert the give string to integer */
    proc parseInt(ref x: int, s: string): int { }

    /* Open file with given mode */
    proc openFile(parent: string, child: string, mode: iomode, ref length: int): file { }
  }

  var params: context;

  proc main(args: [] string) {
    params = new context(args);
    params.load();

    // Show meta-information
    writeln("");
    writeln("Poses     : ", params.nposes);
    writeln("Iterations: ", params.iterations);
    writeln("Ligands   : ", params.natlig);
    writeln("Proteins  : ", params.natpro);
    writeln("Deck      : ", params.deckDir);

    // Compute
    var energiesChapel: [0..<params.nposes] real(32);
    compute(energiesChapel);

    // TODO: validation
  }

  /* Compute results with loaded data */
  proc compute(results: [] real(32)) {

  }

  /* Core computing function */
  proc fasten_main(
    natlig: int,
    natpro: int,
    protein: [] atom,
    ligand: [] atom,
    transforms: [] real(32),
    results: [] real(32),
    forcefield: [] ffParams,
    group: int) { }

  /* Get current time */
  proc timestamp: real(64) {
    return getCurrentTime(unit=TimeUnits.milliseconds);
  }

  /* Print timing */
  proc printTimings(start: real(64), end: real(64)) { }
}