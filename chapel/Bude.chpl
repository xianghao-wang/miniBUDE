module Bude {
  use IO;
  use Time;
  use Context;
  use Helper;
  use Atom;
  use FFParams;
  use VecPoseInner;
  use Configuration;
  use AutoMath;
  use IO.FormattedIO;

  var params: context;
  var resultDom: domain(1);
  var results: [resultDom] real(32);

  proc main(args: [] string) {
    // Load context
    params = new context(args);
    params.load();

    // Show meta-information
    writeln("");
    writeln("Poses     : ", params.nposes);
    writeln("Iterations: ", params.iterations);
    writeln("Ligands   : ", params.natlig);
    writeln("Proteins  : ", params.natpro);
    writeln("Deck      : ", params.deckDir);

    
    var energiesChapel: [0..<params.nposes] real(32);
    // Compute
    compute(energiesChapel);

    // Validate energies
    var length: int;
    const ref_energies = openFile(params.deckDir, FILE_REF_ENERGIES, iomode.r, length);
    var e: real(32);
    var diff: real(32);
    var maxdiff: real(32) = -100.0;
    var n_ref_poses = params.nposes;
    if (params.nposes > REF_NPOSES) {
      writeln("Only validating the first ", REF_NPOSES, " poses");
      n_ref_poses = REF_NPOSES;
    }

    var reader = try! ref_energies.reader();
    for i in dom0(n_ref_poses) {
      try! reader.read(e);
      if (abs(e) < 1.0 && abs(energiesChapel(i)) < 1.0) {
        continue;
      }

      diff = abs(e - energiesChapel(i)) / e;
      if (diff > maxdiff) {
        maxdiff = diff;
      }
    }

    writef("\nLargest difference was %{.###}%%.\n\n", 100 * maxdiff);
  }


  proc compute(ref results: [] real(32)) {
    writeln("\nRunning Chapel");

    var buffer: [dom0(params.nposes)] real(32);    
    // Copy data
    var poses = params.poses;
    var protein = params.protein;
    var ligand = params.ligand;
    var forcefield = params.forcefield;

    // Warm-up
    forall group in 0..<params.nposes/WGSIZE {
      fasten_main(params.natlig, params.natpro, protein, ligand,
                  poses, buffer, forcefield, group);
    }

    // Core part of computing
    const start = timestamp();
    for itr in 0..<params.iterations {
      forall group in dom0(params.nposes / WGSIZE) {
        fasten_main(params.natlig, params.natpro, protein, ligand, poses, buffer, forcefield, group);
      }
    }
    const end = timestamp();

    // Copy to result
    results = buffer;

    printTimings(start, end);
  }

  proc timestamp(): real(64) {
    return getCurrentTime(unit=TimeUnits.milliseconds);
  }

  proc printTimings(start: real(64), end: real(64)) {
    const ms = (end - start) / params.iterations;
    const runtime = ms * 1e-3;

    const ops_per_wg = WGSIZE * 27 + params.natlig * (2 + WGSIZE * 18 + params.natpro * (10 + WGSIZE * 30)) + WGSIZE;
    const total_ops = ops_per_wg * (params.nposes / WGSIZE);
    const flops = total_ops / runtime;
    const gflops = flops / 1e9;

    const total_finsts = 25.0 * params.natpro * params.natlig * params.nposes;
    const finsts = total_finsts / runtime;
    const gfinsts = finsts / 1e9;

    const interactions = 1.0 * params.nposes * params.natlig * params.natpro;
    const interactions_per_sec = interactions / runtime;

    // Print stats
    writef("- Total time:     %7.3dr ms\n", end - start);
    writef("- Average time:   %7.3dr ms\n", ms);
    writef("- Interactions/s: %7.3dr billion\n", (interactions_per_sec / 1e9));
    writef("- GFLOP/s:        %7.3dr\n", gflops);
    writef("- GFInst/s:       %7.3dr\n", gfinsts);
  }
}