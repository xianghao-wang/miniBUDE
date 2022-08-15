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
        exit(0);
      }
    }
  }
}