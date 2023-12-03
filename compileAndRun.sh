rm -R build
cmake -S . -B build
cd build/
make

cd ..
cd app/
./flowBetweenFlatPlates
