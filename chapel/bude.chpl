module Bude {
    use Context;

    var params: context;

    proc main(args: [] string) {
        params = new context(args);
    }
}