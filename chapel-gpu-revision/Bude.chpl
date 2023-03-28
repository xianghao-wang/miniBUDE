module Bude {
  use Context;

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

  // Debug only
  config const DEBUG = false;

  // Energy evaluation parameters
  const CNSTNT: real(32) = 45.0;
  const HBTYPE_F: real(32) = 70.0;
  const HBTYPE_E: real(32) = 69.0;
  const HARDNESS: real(32) = 38.0;
  const NPNPDIST: real(32) = 5.5;
  const NPPDIST: real(32) = 1.0;
  
  config param NUM_TD_PER_THREAD = 2;

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

    // Debug information
    if (DEBUG) {
      const posesSize = params.nposes * 16 / 1000.0;
      const ligSize = params.natlig * 16 / 1000.0;
      const proSize = params.natpro * 16 / 1000.0;

      writeln("");
      writeln("===== DEBUG START =====");
      writeln("Memory Usage");
      writeln("Poses     : ", posesSize, " KB");
      writeln("Ligands   : ", ligSize, " KB");
      writeln("Proteins  : ", proSize, " KB");
      writeln("Total     : ", posesSize + ligSize + proSize, " KB");
      writeln("Locale");
      writeln("Poses     : ", params.forcefield.locale);
      writeln("Ligands   : ", params.ligand.locale);
      writeln("Proteins  : ", params.protein.locale);
      writeln("=====  DEBUG END  =====");
    }

    // Compute
    var energiesChapel: [0..<params.nposes] real(32);
    compute(energiesChapel);
  }

  proc compute(results: [] real(32)) {
    on here.gpus[0] {
      
    }
  }
}