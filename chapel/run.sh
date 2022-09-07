echo "Compiling Bude..."
chpl Bude.chpl --fast
echo "Compiling complete!"
echo "Running Bude..."
./Bude
echo "Bude cleaned!"
rm ./Bude