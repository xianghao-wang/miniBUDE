module BudeGPU {
  use Bude;

  /* Compute results with loaded data */
  proc compute(device: locale, results: [] real(32)) {
    writeln("\nRunning Chapel GPU\n");

    on device {
      // Copy data to GPU
      var protein = params.protein;
      var ligand = params.ligand;
      var forcefield = params.forcefield;
      var buffer: [0..<params.nposes] real(32);

      
    }
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
}