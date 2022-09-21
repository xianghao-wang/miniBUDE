module VecPoseInner {
  use Configuration;
  use Helper;
  use AutoMath;
  use Atom;
  use FFParams;

  // TODO: macro in chapel
  // Configurations
  param CNSTNT: real(32) = 45.0;
  param HBTYPE_F: real(32) = 70.0;
  param HBTYPE_E: real(32) = 69.0;
  param HARDNESS: real(32) = 38.0;
  param NPNPDIST: real(32) = 5.5;
  param NPPDIST: real(32) = 1.0;

  proc fasten_main(
    natlig: int,
    natpro: int,
    ref protein: [] atom,
    ref ligand: [] atom,
    ref transforms: [] real(32),
    ref results: [] real(32),
    ref forcefield: [] ffParams,
    group: int) {

    const wg_ran = 0..WGSIZE-1;

    var transform: [0..2, 0..3, wg_ran] real(32);
    var etot: [wg_ran] real(32);

    for l in wg_ran {
      const ix = group * WGSIZE + l;

      // Compute transformation matrix
      const sx = sin(transforms(0, ix));
      const cx = cos(transforms(0, ix));
      const sy = sin(transforms(1, ix));
      const cy = cos(transforms(1, ix));
      const sz = sin(transforms(2, ix));
      const cz = cos(transforms(2, ix));

      transform(0, 0, l) = cy*cz;
      transform(0, 1, l) = sx*sy*cz - cx*sz;
      transform(0, 2, l) = cx*sy*cz + sx*sz;
      transform(0, 3, l) = transforms(3, ix);
      transform(1, 0, l) = cy*sz;
      transform(1, 1, l) = sx*sy*sz + cx*cz;
      transform(1, 2, l) = cx*sy*sz - sx*cz;
      transform(1, 3, l) = transforms(4, ix);
      transform(2, 0, l) = -sy;
      transform(2, 1, l) = sx*cy;
      transform(2, 2, l) = cx*cy;
      transform(2, 3, l) = transforms(5, ix);

      etot[l] = 0.0;
    }

    {
      // Loop over ligand atoms
      var il = 0;
      do {
        const l_atom = ligand[il];
        const l_params = forcefield[l_atom.aType];
        const lhphb_ltz = l_params.hphb < 0;
        const lhphb_gtz = l_params.hphb > 0;

        var lpos_x: [wg_ran] real(32);
        var lpos_y: [wg_ran] real(32);
        var lpos_z: [wg_ran] real(32);

        
        for l in wg_ran {
          lpos_x(l) = transform(0, 3, l)
            + l_atom.x * transform(0, 0, l)
            + l_atom.y * transform(0, 1, l)
            + l_atom.z * transform(0, 2, l);
          lpos_y(l) = transform(1, 3, l)
            + l_atom.x * transform(1, 0, l)
            + l_atom.y * transform(1, 1, l)
            + l_atom.z * transform(1, 2, l);
          lpos_z(l) = transform(2, 3, l)
            + l_atom.x * transform(2, 0, l)
            + l_atom.y * transform(2, 1, l)
            + l_atom.z * transform(2, 2, l);
        }

        // Loop over protein atoms
        var ip = 0;
        do {
          const p_atom = protein(ip);
          const p_params = forcefield(p_atom.aType);

          const radij: real(32) = p_params.radius + l_params.radius;
          const r_radij: real(32) = 1.0 / radij;

          const elcdst: real(32) = if
            p_params.hbtype == HBTYPE_F && l_params.hbtype == HBTYPE_F
            then 4.0: real(32)
            else 2.0: real(32);

          const elcdst1: real(32) = if
            p_params.hbtype == HBTYPE_F && l_params.hbtype == HBTYPE_F
            then 0.25: real(32)
            else 0.5: real(32);

          const type_E = 
            p_params.hbtype == HBTYPE_E ||
            l_params.hbtype == HBTYPE_E;

          const phphb_ltz = p_params.hphb <  0: real(32);
          const phphb_gtz = p_params.hphb >  0: real(32);
          const phphb_nz  = p_params.hphb != 0: real(32);

          const p_hphb = p_params.hphb * (
            if phphb_ltz && lhphb_gtz
              then -1.0: real(32)
              else 1.0: real(32)
          );
          const l_hphb = l_params.hphb * (
            if phphb_gtz && lhphb_ltz
              then -1.0: real(32)
              else 1.0: real(32)
          );
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
          const r_distdslv: real(32) = 
            1.0 / distdslv;

          const chrg_init: real(32) =
            l_params.elsc * p_params.elsc;

          const dslv_init: real(32) =
            p_hphb + l_hphb; 

          // NOTE: SIMD v.s. data parallelism
          for l in wg_ran {
            // Calculate distance between atoms
            const x: real(32) = lpos_x(l) - p_atom.x;
            const y: real(32) = lpos_y(l) - p_atom.y;
            const z: real(32) = lpos_z(l) - p_atom.z;
            const distij: real(32) = sqrt(x * x + y * y + z* z); 

            // Calculate the sum of the sphere radii
            const distbb: real(32) = distij - radij;

            const zone1: bool = distbb < 0.0;

            // Calculate steric energy
            etot[l] = etot[l] + (1.0 - distij * r_radij) * (
              if zone1 then 2.0 * HARDNESS else 0.0
            );

            // Calculate formal and dipole charge interactions
            var chrg_e: real(32) =
              chrg_init * (
                if zone1 then 1.0: real(32) else 1.0: real(32) - distbb * elcdst1
              ) * (
                if distbb < elcdst then 1.0: real(32) else 0.0: real(32)
              );
            var neg_chrg_e: real(32) = -abs(chrg_e);
            chrg_e = if type_E then neg_chrg_e else chrg_e;
            etot[l] = etot[l] + chrg_e * CNSTNT;

            const coeff:real (32) = 1.0: real(32) - distbb * r_distdslv;
            var dslv_e: real(32) = dslv_init * (
              if distbb < distdslv && phphb_nz
              then 1.0: real(32)
              else 0.0: real(32)
            );

            dslv_e *= if zone1 then 1.0: real(32) else coeff;
            etot[l] = etot[l] + dslv_e;
          }
          ip += 1;
        } while (ip < natpro);
        il += 1;
      } while (il < natlig);

      for l in wg_ran {
        results[group * WGSIZE + l] = etot[l] * 0.5;
      }
    }
  }

}