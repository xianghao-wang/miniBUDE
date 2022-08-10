record Atom {
  var x: real;
  var y: real;
  var z: real;
  var aType: int;
}

record FFParams {
  var hbtype: int;
  var radius: real;
  var hphb: real;
  var elsc: real;
}

record ParameterConfig {
  const natlig: int;
  const natpro: int;
  const ntypes: int;
  const nposes: int;

  proc loadParameters(args: [] string) {

  }
}

proc main(args: [] string) {

}
