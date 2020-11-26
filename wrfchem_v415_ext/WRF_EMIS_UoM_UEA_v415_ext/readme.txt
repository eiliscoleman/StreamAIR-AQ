This is the emission script from OIVIND. I have modified and added other modules for processing
naei and TNO emissions for CRIMECH and CBMZ UK speciations. 

To run emission script on archer:

module swap PrgEnv-cray PrgEnv-gnu
module load cray-netcdf
make clean
make -f makefile_new
./main

