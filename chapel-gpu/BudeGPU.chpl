module BudeGPU {
    use Bude;

    proc compute(device: locale, results: [] real(32)) {
        // Copy data from host to device
        on loc {
            var protein = params.protein;
            var ligand = params.ligand;
            var forcefield = params.forcefield;
            var buffer: [0..<params.nposes] real(32);
            var poses = params.poses;

            forall ii in params.iterations {
                foreach jj in 0..<WGSIZE {
                    fasten_main(
                        params.natlig,
                        params.natpro,
                        params.protein,
                        params.ligand,
                        poses,
                        buffer,
                        forcefield,
                        jj
                    );
                }
            }
        }
    }
}