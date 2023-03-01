use AutoMath;

const N: int = 32;


proc test1() {
    var C: [0..<N] real(32);
    on here.gpus[0] {
        var A: [0..<N] real(32);
        var B: [0..<N] real(32) = noinit;

        foreach i in 0..<N {
            B[i] = sqrt(A[i]);
        }
    }

}

proc main() {
    test1();
}