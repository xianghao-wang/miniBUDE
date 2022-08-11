module Helper {
    use FileSystem;
    use IO;
    
    proc parseInt(ref x: int, s: string): int {
      try {
        x = s: int;
      } catch {
        return -1;
      }
      return x;
    }

    proc openFile(parent: string, child: string, mode: iomode, ref length: int): file {
      const name = parent + child;
      const aFile: file;

      try {           
        aFile = open(name, mode);
      } catch {
        fprintf(stderr, "Failed to open '", name, "'\n");
        exit(1);
      }

      length = aFile.getFileSize();
      return aFile;
    }


}