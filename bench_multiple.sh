# Usage
# ./bench_multiple folder_name max_device_numbers

n=4194304
i=16
WG=(64 128 256 512)

folder=$(pwd)
cd $1
make clean
make

echo benchmarking $1 on multiple GPUs

out=${folder}/$(basename $(pwd))_multi.txt
echo "" > $out
for wgsize in ${WG[@]}
do
  for (( ngpu=1; ngpu<=$2; ngpu++ ))
  do
    echo "" >> $out
    echo ./bude -i $i -n $n --wgsize $wgsize --ngpu $ngpu >> $out
    ./bude -i $i -n $n --wgsize $wgsize --ngpu $ngpu >> $out
    echo "####################################################" >> $out
  done
done