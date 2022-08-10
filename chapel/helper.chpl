module Helper {
    proc parseInt(ref x: int, s: string): int {
        try {
            x = s: int;
        } catch {
            return -1;
        }

        return x;
    }
}