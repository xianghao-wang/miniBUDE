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
    // TODO: restirct keyword

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

    writeln("warmup");
    // TODO: infinite loop noted
    // Warm-up
    coforall group in dom0(params.nposes/WGSIZE) {
      fasten_main(params.natlig, params.natpro, protein, ligand,
                  poses, buffer, forcefield, group);
    }

    // Core part of computing
    // timestamp: start
    const start = timestamp();
    writeln("start");

    // NOTE: SIMD
    // coforall itr in 1..params.iterations {
      forall group in dom0(params.nposes / WGSIZE) {
        fasten_main(params.natlig, params.natpro, protein, ligand, poses, buffer, forcefield, group);
      }
    // }

    // timestamp: end
    const end = timestamp();
    writeln("end");

    results = buffer;
    writeln(end - start);
  }

  proc timestamp(): real(64) {
    return getCurrentTime(unit=TimeUnits.milliseconds);
  }
}