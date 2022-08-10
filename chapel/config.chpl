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
    var hbtype: int(32);
    var radius: real(32);
    var hphb: real(32);
    var elsc: real(32);
  }

  record Atom {
    var x, y, z: real(32);
    var aType: int(32);
  }
}