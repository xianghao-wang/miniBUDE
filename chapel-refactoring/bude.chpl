module Bude {
    var natlig: int;
    var natpro: int;
    var ntypes: int;
    var nposes: int = DEFAULT_NPOSES;
    
    var proteinDom: domain(1);
    var protein: [proteinDom] atom;

    var ligandDom: domain(1);
    var ligand: [ligandDom] atom;

    var forcefieldDom: domain(1);
    var forcefield: [forcefieldDom] ffParams;

    var posesDom: domain(1);
    var poses: [posesDom] ([6] real(32));

    var deckDir: string = DATA_DIR;
    var iterations: int = DEFAULT_ITERS;
}