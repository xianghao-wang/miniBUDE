module Bude {
  use Time;
  use Context;
  use Helper;
  use Atom;
  use FFParams;
  use VecPoseInner;
  use Configuration;

  var params: context;
  var resultDom: domain(1);
  var results: [resultDom] real(32);

  proc main(args: [] string) {
    // Load context
    params = new context(args);
    params.load();
    resultDom = {0..(params.nposes-1)};

    // Show meta-information
    writeln("");
    writeln("Poses     :", params.nposes);
    writeln("Iterations:", params.iterations);
    writeln("Ligands   :", params.natlig);
    writeln("Proteins  :", params.natpro);
    writeln("Deck      :", params.deckDir);

    
    var energiesChapel: [dom0(params.nposes)] real(32);
    // Compute
    compute(energiesChapel);
  }


  proc compute(ref results: [] real(32)) {
    writeln("\nRunning Chapel");

    var poses: [params.poses.domain] real(32);
    var protein: [params.protein.domain] atom;
    var ligand: [params.ligand.domain] atom;
    var forcefield: [params.forcefield.domain] ffParams;
    // Note: array is initialised when created
    var buffer: [dom0(params.nposes)] real(32);
    
    forall i in poses.domain {
      poses[i] = params.poses[i];
    }

    forall i in protein.domain {
      protein[i] = params.protein[i];
    }

    forall i in ligand.domain {
      ligand[i] = params.ligand[i];
    }

    forall i in forcefield.domain {
      forcefield[i] = params.forcefield[i];
    }

    // Warm-up
    forall group in dom0(params.nposes/WGSIZE) {
      fasten_main(params.natlig, params.natpro, protein, ligand,
                  poses, buffer, forcefield, group);
    }

    // Core part of computing
    // timestamp: start
    const start = timestamp();

    // NOTE: SIMD
    for itr in 1..params.iterations {
      forall group in dom0(params.nposes / WGSIZE) {
        fasten_main(params.natlig, params.natpro, protein, ligand, poses, buffer, forcefield, group);
      }
    }

    // timestamp: end
    const end = timestamp();

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
    writeln("- Total time:     ", (end-start), " ms");
    writeln("- Average time:   ", ms, " ms");
    writeln("- Interactions/s: ", (interactions_per_sec / 1e9), " billion");
    writeln("- GFLOP/s:        ", gflops);
    writeln("- GFInst/s:       ", gfinsts);
  }
}