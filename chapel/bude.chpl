module Bude {
    use Config;
    use Helper;
    use Atom;
    use FFParams;

    record ParameterConfig {
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

    proc init(args: [] string) {

      /* load command-line parameters */
      const argc = args.size;
      var arg: string;
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
          this.deckDir = arg;
          i += 1;
        } else {
          writeln("Unrecognized argument '", arg, "' (try '--help')\n");
          exit(1);
        }
        i += 1;
      }

      /* Load file content */
      var aFile: file;
      var length: int;
      aFile = openFile()
    }
}