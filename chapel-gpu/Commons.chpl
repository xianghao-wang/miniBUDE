module Commons {
  use IO;

  /* Some helper functions */

  /* Load data from file to array */
  proc loadData(aFile: file, ref A: [] ?t, size: int) {
    const n = A.size;
    var readChannel = try! aFile.reader(kind=iokind.native, region=0..n*size);
    try! readChannel.read(A);
    try! readChannel.close();
  }

  /* Load data with length of the given bytes */
  proc loadDataPiece(aFile: file, ref A: ?t, base: int, offset: int) {
    var r = try! aFile.reader(kind=iokind.native, region=base..base+offset);
    try! r.read(A);
    try! r.close();
  }

  /* Convert the give string to integer */
  proc parseInt(ref x: int, s: string): int {
    try {
      x = s: int;
    } catch {
      return -1;
    }
    return x;
  }

  /* Open file with given mode */
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
}