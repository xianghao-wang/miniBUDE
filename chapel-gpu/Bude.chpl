module Bude {
  use Context;
  use Time;
  use IO;
  use Commons;
  use GpuDiagnostics;
  use GPU;

  param NUM_TD_PER_THREAD = 1;
  param DEFAULT_ITERS = 8;
  param DEFAULT_NPOSES = 65536;
  param REF_NPOSES = 65536;
  param DATA_DIR = "../data/bm1";
  param FILE_LIGAND = "/ligand.in";
  param FILE_PROTEIN = "/protein.in";
  param FILE_FORCEFIELD = "/forcefield.in";
  param FILE_POSES = "/poses.in";
  param FILE_REF_ENERGIES = "/ref_energies.out";
  param ATOM_SIZE = 16;
  param FFPARAMS_SIZE = 16;
  param WGSIZE = 16;

  // Energy evaluation parameters
  const CNSTNT: real(32) = 45.0;
  const HBTYPE_F: real(32) = 70.0;
  const HBTYPE_E: real(32) = 69.0;
  const HARDNESS: real(32) = 38.0;
  const NPNPDIST: real(32) = 5.5;
  const NPPDIST: real(32) = 1.0;
  

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

    // Compute
    var energiesChapel: [0..<params.nposes] real(32);
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
    for i in 0..<n_ref_poses {
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

  proc compute(results: [] real(32)) {
    on here.gpus[0] {
      // Copy data to device
      var protein = params.protein;
      var ligand = params.ligand;
      var forcefield = params.forcefield;
      var poses = params.poses;
      var buffer: [0..<params.nposes] real(32);

      const iterations = params.iterations;
      const nposes = params.nposes;
      const natlig = params.natlig;
      const natpro = params.natpro;

      const start = timestamp();
      for itr in 0..<iterations {
        // fasten_main
        foreach ii in 0..<nposes/NUM_TD_PER_THREAD {
          const ind = ii * NUM_TD_PER_THREAD;
          var etot: NUM_TD_PER_THREAD * real(32);
          var transform: NUM_TD_PER_THREAD * (3 * (4 * real(32)));

          for jj in 0..<NUM_TD_PER_THREAD {
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

              for jj in 0..<NUM_TD_PER_THREAD {
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
                
                for jj in 0..<NUM_TD_PER_THREAD {
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
              } // foreach ip in 0..< params.natpro
          } // foreach il in 0..<params.natlig
          for jj in 0..<NUM_TD_PER_THREAD {
            // assertOnGpu();
            buffer[ind+jj] = etot[jj] * 0.5;
          }
        } // foreach ii in 0..<params.nposes/NUM_TD_PER_THREAD
        // Copy results from device to host
      } // for itr in 0..<params.iterations

      const end = timestamp();
      printTimings(start, end);

      results = buffer;
    } // on  
  } // main

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