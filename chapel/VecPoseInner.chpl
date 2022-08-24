module VecPoseInner {
  use Configuration;
  use Helper;
  use AutoMath;
  use Atom;
  use FFParams;

  // TODO: macro in chapel

  const CNSTNT = 45.0,
        HBTYPE_F = 70,
        HBTYPE_E = 69,
        HARDNESS = 38.0,
        NPNPDIST = 5.5,
        NPPDIST = 1.0;

  proc fasten_main(
    natlig: int,
    natpro: int,
    ref protein: [] atom,
    ref ligand: [] atom,
    ref transforms: [] real(32),
    ref results: [] real(32),
    ref forcefield: [] ffParams,
    group: int) {
                    
    var transform: [{0..2, 0..3, 0..(WGSIZE-1)}] real(32);
    var etot: [{0..(WGSIZE-1)}] real(32);

    // TODO: SIMD vectorisation
    // forall l in dom0(WGSIZE) {
    //   const ix = group * WGSIZE + l;

    //   // Compute transformation matrix
    //   const sx = 
    // }
  }

}