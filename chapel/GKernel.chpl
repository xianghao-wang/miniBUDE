module GKernel {
  use Bude;
  use Math;

  proc gkernel(context: params, results: [] real(32)) {
    const ngpu = here.gpus.size;
    var times: [0..<ngpu] real;
    coforall (gpu, gpuID) in zip(here.gpus, here.gpus.domain) with (ref times) do on gpu {
      const iterations = context.iterations: int(32);
      const nposes = (context.nposes / ngpu) : int(32);
      const natlig = context.natlig: int(32);
      const natpro = context.natpro: int(32);
      const wgsize = context.wgsize: int(32);

      const protein = context.protein;
      const ligand = context.ligand;
      const forcefield = context.forcefield;
      const poses: [{0:int(32)..<6:int(32), 0..<nposes}] real(32) = context.poses[{0..<6, gpuID*nposes..<(gpuID+1)*nposes}];
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
}