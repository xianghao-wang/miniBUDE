module CKernel {
  use Bude;
  use Math;

  proc ckernel(context: params, results: [] real(32)) {
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

  proc fasten_main(
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
}