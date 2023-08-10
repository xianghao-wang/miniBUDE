module Bude {
  use Time;
  use IO;
  use GPU;
  use Math;
  use CTypes;
  use ArgumentParser;
  use Path;
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
      reader = aFile.reader(kind=iokind.native, region=0.., locking=false);
      var available = (length / (6 * c_sizeof(real(32)): int)): int(32);
      var current = 0:int(32);
      while (current < this.nposes) {
        var fetch = this.nposes - current;
        if (fetch > available) then fetch = available;

        for i in posesDomain.dim(0) {
          const base = i*available*c_sizeof(real(32)):int;
          const amount = fetch*c_sizeof(real(32)):int;
          reader.seek(region=base..base+amount);
          reader.read(this.poses(i, current..));
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
  } // main

  proc compute(results: [] real(32)) {
    if (CHPL_GPU == "nvidia") {
      writeln("\nRunning Chapel on ", here.gpus.size, (if here.gpus.size > 1 then " GPUs" else " GPU"));
      gpukernel(context, results);
    } else if (CHPL_GPU == "none") {
      writeln("\nRunning Chapel");
      cpukernel(context, results);
    } else {
      writeln("\n" + CHPL_GPU + " is not supported.");
      exit(1);
    }
  } 

  proc gpukernel(context: params, results: [] real(32)) {
    const ngpu = here.gpus.size;
    var times: [0..<ngpu] real;
    coforall (gpu, gpuID) in zip(here.gpus, here.gpus.domain) do on gpu {
      const iterations = context.iterations: int(32);
      const nposes = (context.nposes / ngpu) : int(32);
      const natlig = context.natlig: int(32);
      const natpro = context.natpro: int(32);
      const wgsize = context.wgsize: int(32);

      const protein = context.protein;
      const ligand = context.ligand;
      const forcefield = context.forcefield;
      const poses: [0..<6:int(32), 0..<nposes] real(32) = context.poses[.., gpuID*nposes..<(gpuID+1)*nposes];
      var buffer: [0..<nposes] real(32);

      times[gpuID] = timestampMS();
      for i in 0..<iterations {
        foreach ii in 0..<nposes/NUM_TD_PER_THREAD {
          __primitive("gpu set blockSize", wgsize);
          const ind = ii * NUM_TD_PER_THREAD;
          var etot: NUM_TD_PER_THREAD * real(32);
          var transform: NUM_TD_PER_THREAD * (3 * (4 * real(32)));

          for param jj in 0..<NUM_TD_PER_THREAD {
            const ix = ind + jj;
            // Compute transformation matrix
            const sx = sin(poses(0, ix));
            const cx = cos(poses(0, ix));
            const sy = sin(poses(1, ix));
            const cy = cos(poses(1, ix));
            const sz = sin(poses(2, ix));
            const cz = cos(poses(2, ix));
            transform(jj)(0)(0) = cy*cz;
            transform(jj)(0)(1) = sx*sy*cz - cx*sz;
            transform(jj)(0)(2) = cx*sy*cz + sx*sz;
            transform(jj)(0)(3) = poses(3, ix);
            transform(jj)(1)(0) = cy*sz;
            transform(jj)(1)(1) = sx*sy*sz + cx*cz;      
            transform(jj)(1)(2) = cx*sy*sz - sx*cz;
            transform(jj)(1)(3) = poses(4, ix);
            transform(jj)(2)(0) = -sy;
            transform(jj)(2)(1) = sx*cy;
            transform(jj)(2)(2) = cx*cy;
            transform(jj)(2)(3) = poses(5, ix);
            etot[jj] = 0.0;
          } // for jj in 0..<NUM_TD_PER_THREAD

          for il in 0..<natlig {
            const l_atom = ligand[il];
            const l_params = forcefield[l_atom.aType];
            const lhphb_ltz = l_params.hphb < 0.0;
            const lhphb_gtz = l_params.hphb > 0.0;

            // Transform ligand atom
            var lpos_x: NUM_TD_PER_THREAD * real(32);
            var lpos_y: NUM_TD_PER_THREAD * real(32);
            var lpos_z: NUM_TD_PER_THREAD * real(32);

            for param jj in 0..<NUM_TD_PER_THREAD {
              lpos_x[jj] = transform(jj)(0)(3)
                + l_atom.x * transform(jj)(0)(0)
                + l_atom.y * transform(jj)(0)(1)
                + l_atom.z * transform(jj)(0)(2);

              lpos_y[jj] = transform(jj)(1)(3)
                + l_atom.x * transform(jj)(1)(0)
                + l_atom.y * transform(jj)(1)(1)
                + l_atom.z * transform(jj)(1)(2);

              lpos_z[jj] = transform(jj)(2)(3)
                + l_atom.x * transform(jj)(2)(0)
                + l_atom.y * transform(jj)(2)(1)
                + l_atom.z * transform(jj)(2)(2);
            } // foreach jj in 0..<NUM_TD_PER_THREAD

            for ip in 0..< natpro {
              const p_atom = protein[ip];
              const p_params = forcefield[p_atom.aType];

              const radij = p_params.radius + l_params.radius;
              const r_radij = 1.0 / radij;

              const elcdst = if
                p_params.hbtype == HBTYPE_F && l_params.hbtype == HBTYPE_F
                then 4.0: real(32)
                else 2.0: real(32);

              const elcdst1 = if
                p_params.hbtype == HBTYPE_F && l_params.hbtype == HBTYPE_F
                then 0.25: real(32)
                else 0.5: real(32);

              const type_E = p_params.hbtype == HBTYPE_E || l_params.hbtype == HBTYPE_E;
              const phphb_ltz = p_params.hphb <  0;
              const phphb_gtz = p_params.hphb >  0;
              const phphb_nz  = p_params.hphb != 0;

              const p_hphb = p_params.hphb 
                * if phphb_ltz && lhphb_gtz then -1.0: real(32) else 1.0: real(32);

              const l_hphb = l_params.hphb 
                * if phphb_gtz && lhphb_ltz then -1.0: real(32) else 1.0: real(32);

              const distdslv =
                if phphb_ltz
                then (
                  if lhphb_ltz
                  then NPNPDIST
                  else NPPDIST
                ) else (
                  if lhphb_ltz
                  then NPPDIST
                  else -max(real(32))
                );

              const r_distdslv = 1.0 / distdslv;
              const chrg_init = l_params.elsc * p_params.elsc;
              const dslv_init = p_hphb + l_hphb; 
              
              for param jj in 0..<NUM_TD_PER_THREAD {
                const x = lpos_x[jj] - p_atom.x;
                const y = lpos_y[jj] - p_atom.y;
                const z = lpos_z[jj] - p_atom.z;
                const distij = sqrt(x*x + y*y + z*z);

                const distbb = distij - radij;
                const zone1 = distbb < 0.0;

                etot[jj] += (1.0 - distij * r_radij) * (if zone1 then 2.0*HARDNESS else 0.0);

                // Calculate formal and dipole charge interactions
                var chrg_e =
                  chrg_init * (
                    if zone1 
                    then 1.0: real(32)
                      else 1.0: real(32) - distbb * elcdst1
                  ) * (
                    if distbb < elcdst 
                    then 1.0: real(32)
                    else 0.0: real(32)
                  );
                
                var neg_chrg_e = -abs(chrg_e);
                chrg_e = if type_E then neg_chrg_e else chrg_e;
                etot[jj] += chrg_e * CNSTNT;

                const coeff = 1.0 - distbb * r_distdslv;
                var dslv_e = dslv_init 
                  * if distbb < distdslv && phphb_nz then 1.0: real(32) else 0.0: real(32);

                dslv_e *= if zone1 then 1.0: real(32) else coeff;
                etot[jj] += dslv_e;                  
              } // foreach jj in 0..<NUM_TD_PER_THREAD
            } // foreach ip in 0..< context.natpro
          } // foreach il in 0..<context.natlig
          for param jj in 0..<NUM_TD_PER_THREAD {
            buffer[ind+jj] = etot[jj] * 0.5;
          }
        } // foreach ii in 0..<context.nposes/NUM_TD_PER_THREAD
        results[gpuID*nposes..<(gpuID+1)*nposes] = buffer;
      } // for iter in 0..<iterations
      times[gpuID] = timestampMS() - times[gpuID];
    }

    printTimings(max reduce times);
  }

  proc cpukernel(context: params, results: [] real(32)) {
    var buffer: [0..<context.nposes] real(32); 
    var poses = context.poses;
    var protein = context.protein;
    var ligand = context.ligand;
    var forcefield = context.forcefield;

    const natlig = context.natlig: int(32);
    const natpro = context.natpro: int(32);
    const nposes = context.nposes: int(32);

    // Warm-up
    forall group in 0..<nposes/NUM_TD_PER_THREAD {
      fasten_main(natlig, natpro, protein, ligand,
                  poses, buffer, forcefield, group: int(32));
    }

    // Core part of computing
    const start: real = timestampMS();
    for itr in 0..<context.iterations {
      forall group in 0..<nposes / NUM_TD_PER_THREAD {
        fasten_main(natlig, natpro, protein, ligand,
                  poses, buffer, forcefield, group: int(32));
      }
    }
    const end: real = timestampMS();

    // Copy to result
    results = buffer;

    printTimings(end - start);
  }

  private proc fasten_main(
    natlig: int(32),
    natpro: int(32),
    protein: [] atom,
    ligand: [] atom,
    transforms: [] real(32),
    results: [] real(32),
    forcefield: [] ffParams,
    group: int(32)) {

    var transform: [0..<3:int(32), 0..<4:int(32), 0..<NUM_TD_PER_THREAD] real(32) = noinit;
    var etot: [0..<NUM_TD_PER_THREAD] real(32) = noinit;

    // Compute transformation matrix
    foreach i in 0..<NUM_TD_PER_THREAD {
      const ix = group*NUM_TD_PER_THREAD + i;
      const sx = sin(transforms(0, ix));
      const cx = cos(transforms(0, ix));
      const sy = sin(transforms(1, ix));
      const cy = cos(transforms(1, ix));
      const sz = sin(transforms(2, ix));
      const cz = cos(transforms(2, ix));
      transform(0, 0, i) = cy*cz;
      transform(0, 1, i) = sx*sy*cz - cx*sz;
      transform(0, 2, i) = cx*sy*cz + sx*sz;
      transform(0, 3, i) = transforms(3, ix);
      transform(1, 0, i) = cy*sz;
      transform(1, 1, i) = sx*sy*sz + cx*cz;      
      transform(1, 2, i) = cx*sy*sz - sx*cz;
      transform(1, 3, i) = transforms(4, ix);
      transform(2, 0, i) = -sy;
      transform(2, 1, i) = sx*cy;
      transform(2, 2, i) = cx*cy;
      transform(2, 3, i) = transforms(5, ix);

      etot[i] = 0.0;
    }
    
    foreach il in 0..<natlig {
      const l_atom = ligand[il];
      const l_params = forcefield[l_atom.aType];
      const lhphb_ltz = l_params.hphb < 0.0;
      const lhphb_gtz = l_params.hphb > 0.0;

      // Transform ligand atom
      var lpos_x: [0..<NUM_TD_PER_THREAD] real(32) = noinit;
      var lpos_y: [0..<NUM_TD_PER_THREAD] real(32) = noinit;
      var lpos_z: [0..<NUM_TD_PER_THREAD] real(32) = noinit;

      foreach l in 0..<NUM_TD_PER_THREAD {
        lpos_x[l] = transform(0, 3, l)
          + l_atom.x * transform(0, 0, l)
          + l_atom.y * transform(0, 1, l)
          + l_atom.z * transform(0, 2, l);

        lpos_y[l] = transform(1, 3, l)
          + l_atom.x * transform(1, 0, l)
          + l_atom.y * transform(1, 1, l)
          + l_atom.z * transform(1, 2, l);

        lpos_z[l] = transform(2, 3, l)
          + l_atom.x * transform(2, 0, l)
          + l_atom.y * transform(2, 1, l)
          + l_atom.z * transform(2, 2, l);
      }

      foreach ip in 0..<natpro {
        const p_atom = protein(ip);
        const p_params = forcefield(p_atom.aType);

        const radij = p_params.radius + l_params.radius;
        const r_radij = 1.0 / radij;

        const elcdst = if
          p_params.hbtype == HBTYPE_F && l_params.hbtype == HBTYPE_F
          then 4.0: real(32)
          else 2.0: real(32);

        const elcdst1 = if
          p_params.hbtype == HBTYPE_F && l_params.hbtype == HBTYPE_F
          then 0.25: real(32)
          else 0.5: real(32);

        const type_E = p_params.hbtype == HBTYPE_E || l_params.hbtype == HBTYPE_E;
        const phphb_ltz = p_params.hphb <  0;
        const phphb_gtz = p_params.hphb >  0;
        const phphb_nz  = p_params.hphb != 0;

        const p_hphb = p_params.hphb 
          * if phphb_ltz && lhphb_gtz then -1.0: real(32) else 1.0: real(32);

        const l_hphb = l_params.hphb 
          * if phphb_gtz && lhphb_ltz then -1.0: real(32) else 1.0: real(32);

        const distdslv =
          if phphb_ltz
          then (
            if lhphb_ltz
            then NPNPDIST
            else NPPDIST
          ) else (
            if lhphb_ltz
            then NPPDIST
            else -max(real(32))
          );

        const r_distdslv = 1.0 / distdslv;
        const chrg_init = l_params.elsc * p_params.elsc;
        const dslv_init = p_hphb + l_hphb; 

        foreach l in 0..<NUM_TD_PER_THREAD {
          // Calculate distance between atoms
          const x = lpos_x(l) - p_atom.x;
          const y = lpos_y(l) - p_atom.y;
          const z = lpos_z(l) - p_atom.z;
          const distij = sqrt(x * x + y * y + z* z); 

          // Calculate the sum of the sphere radii
          const distbb = distij - radij;
          const zone1 = distbb < 0.0: real(32);

          // Calculate steric energy
          etot[l] += (1.0 - distij * r_radij)
            * if zone1 then 2.0: real(32) * HARDNESS else 0.0: real(32);

          // Calculate formal and dipole charge interactions
          var chrg_e =
            chrg_init * (
              if zone1 
              then 1.0: real(32)
                else 1.0: real(32) - distbb * elcdst1
            ) * (
              if distbb < elcdst 
              then 1.0: real(32)
              else 0.0: real(32)
            );
          
          var neg_chrg_e = -abs(chrg_e);
          chrg_e = if type_E then neg_chrg_e else chrg_e;
          etot[l] += chrg_e * CNSTNT;

          const coeff = 1.0 - distbb * r_distdslv;
          var dslv_e = dslv_init 
            * if distbb < distdslv && phphb_nz then 1.0: real(32) else 0.0: real(32);

          dslv_e *= if zone1 then 1.0: real(32) else coeff;
          etot[l] += dslv_e;
        }
      }
    }

    results[group * NUM_TD_PER_THREAD..<(group + 1) * NUM_TD_PER_THREAD] = 0.5 : real(32) * etot;
  }

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
