# Usage
# ./bench_single folder_name

n=65536
i=16
WG=(64 128 256 512)

folder=$(pwd)
cd $1
make clean
make

echo benchmarking $1

out=${folder}/$(basename $(pwd)).txt
echo "" > $out
for wgsize in ${WG[@]}
do
  t=$n
  for s in {1..7}
  do
    echo "" >> $out
    echo ./bude -i $i -n $n --wgsize $wgsize >> $out
    ./bude -i $i -n $n --wgsize $wgsize >> $out
    echo "####################################################" >> $out
    n=`expr $n \* 2`
  done
  n=$t
done