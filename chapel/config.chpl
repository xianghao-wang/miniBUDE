module Config {
  const WGSIZE = 4,
        DEFAULT_ITERS = 8,
        DEFAULT_NPOSES = 65536,
        REF_NPOSES = 65536,
        DATA_DIR = "../data/bm1",
        FILE_LIGAND = "/ligand.in",
        FILE_PROTEIN = "/protein.in",
        FILE_FORCEFIELD = "/forcefield.in",
        FILE_POSES = "/poses.in",
        FILE_REF_ENERGIES = "/ref_energies.out";

  record FFParams {
    var hbtype: int;
    var radius: real;
    var hphb: real;
    var elsc: real;
  }

  record Atom {
    var x: real;
    var y: real;
    var z: real;
    var aType: int;
  }
}