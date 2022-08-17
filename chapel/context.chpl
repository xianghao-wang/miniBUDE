module Context {
  use Configuration;
  use Helper;
  use Atom;
  use FFParams;
  use IO;
  use FileSystem;

  record context {
    var iterations: int = DEFAULT_ITERS;
    var nposes: int = DEFAULT_NPOSES;
    var deckDir: string = DATA_DIR;

    // Domains for arrays
    var natlig: int;
    var ligandDom: domain(1);
    var natpro: int;
    var proteinDom: domain(1);
    var ntypes: int;
    var forcefieldDom: domain(1);
    var posesDom: domain(2);
    
    var ligand: [ligandDom] atom;
    var protein: [proteinDom] atom;
    var forcefield: [forcefieldDom] ffParams;
    var poses: [posesDom] real(32);

    proc init() {
      
    }

    proc init(args: [] string) {
      /* Load command-line arguements */
      const argc = args.size;
      var arg: string;

      var t_deckDir = DATA_DIR;

      var i = 1;
      while i < argc {
        arg = args[i];
        if arg == "--iterations" || arg == "-i" {
          if i + 1 >= argc || parseInt(this.iterations, args[i+1]) < 0 {
            writeln("Invalid number of iterations");
            exit(1);
          }
          i += 1;
        } else if arg == "--numposes" || arg == "-n" {
          if i + 1 >= argc || parseInt(this.nposes, args[i+1]) < 0 {
            writeln("Invalid number of poses");
            exit(1);
          }
          i += 1;
        } else if arg == "--help" || arg == "-h" {
          writeln("");
          writeln("Usage: ./bude [OPTIONS]");
          writeln("Options:");
          writeln("  -h  --help               Print this message");
          writeln("  -i  --iterations I       Repeat kernel I times (default: ", DEFAULT_ITERS, ")");
          writeln("  -n  --numposes   N       Compute energies for N poses (default: ", DEFAULT_NPOSES, ")");
          writeln("      --deck       DECK    Use the DECK directory as input deck (default: ", DATA_DIR, ")");
          writeln("");
        } else if arg == "--deck" {
          if (i + 1 >= argc) {
            writeln("Invalid deck");
            exit(1);
          }
          t_deckDir = args[i + 1];
          i += 1;
        } else {
          writeln("Unrecognized argument '", arg, "' (try '--help')\n");
          exit(1);
        }
        i += 1;
      }

      this.deckDir = t_deckDir;
      var length: int;
      var aFile: file;

      /* init ligand array */
      aFile = openFile(this.deckDir, FILE_LIGAND, iomode.r, length);
      this.natlig = length / ATOM_SIZE;
      this.ligandDom = {0..(this.natlig-1)};

      /* init protein array */
      aFile = openFile(this.deckDir, FILE_PROTEIN, iomode.r, length);
      this.natpro = length / ATOM_SIZE;
      this.proteinDom = {0..(this.natpro-1)};

      /* init forcefield array */
      aFile = openFile(this.deckDir, FILE_FORCEFIELD, iomode.r, length); 
      this.ntypes = length / FFPARAMS_SIZE;
      this.forcefieldDom = {0..(this.ntypes-1)};

      /* init poses array */
      this.posesDom = {
        0..5,
        (0..this.nposes-1)
      };
    }

    proc load() {
      var length: int;
      var aFile: file;
      
      /* load ligand */
      aFile = openFile(this.deckDir, FILE_LIGAND, iomode.r, length);
      loadData(aFile, this.ligand, ATOM_SIZE);

      /* load protein */
      aFile = openFile(this.deckDir, FILE_PROTEIN, iomode.r, length);
      loadData(aFile, this.protein, ATOM_SIZE);

      /* load forcefields */
      aFile = openFile(this.deckDir, FILE_FORCEFIELD, iomode.r, length);
      loadData(aFile, this.forcefield, FFPARAMS_SIZE);

      /* load poses */
      aFile = openFile(this.deckDir, FILE_POSES, iomode.r, length);
      var available = length / (6 * 4);
      var cur_poses = 0, fetch, address = 0;
      while (cur_poses < this.nposes) {
        fetch = this.nposes - cur_poses;
        if (fetch > available) {
          fetch = available;
        }
        address = 0;
        for i in 0..5 {
          address += i * available * 4;
          for j in 0..(fetch-1) {
            loadDataPiecie(aFile, this.poses(i, cur_poses+j), address+4*j, 4);
          }
        }
        address = 0; // rewind
        cur_poses += fetch;
      }

      this.nposes = cur_poses;
    }
  }
}
