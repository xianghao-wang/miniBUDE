module Helper {
  use IO;
  use FileSystem;
  use Atom;

  /* Open a file */
  proc openFile(parent: string, child: string, mode: iomode, ref length: int): file {
    const name = parent + child;
    var aFile: file;

    try {
      aFile = open(name, mode);
      length = aFile.size;
    } catch {
      try {
        stderr.writeln("Failed to open '", name, "'");
        exit(0);
      } catch {
        exit(0);
      }
    }

    return aFile;
  }

  /* Load data from file to record array */
  proc loadData(aFile: file, A: [] ?t, size: int) {
    const n = A.size;

    for i in 0..(n - 1) {
      try {
        var r = aFile.reader(kind=iokind.native, start=i * size, end=(i + 1) * size);
        r.read(A[i]);
        r.close();
      } catch {
        // TODO: fix error
        // stderr.writeln("Failed to read '", aFile.path, "'");
        exit(0);
      }
    }
  }

  /* Load data piece */
  proc loadDataPiecie(aFile: file, ref A: ?t, base: int, offset: int) {
    try {
      var r = aFile.reader(kind=iokind.native, start=base, end=base+offset);
      r.read(A);
      r.close();
    } catch {
      // TODO: fix error
      try! stderr.writeln("Failed to read '", aFile.path, "'");
      exit(0);
    }
  }

  proc parseInt(ref x: int, s: string): int {
    try {
      x = s: int;
    } catch {
      return -1;
    }
    return x;
  }

  proc dom0(exclusiveUpper: int): domain(1) {
    return {0..(exclusiveUpper-1)};
  }
}