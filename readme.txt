binddir = "where is executable file" (like../nfs/amino-home/zhng/zhanglab_programs/EMCMR/)
datadir = "where is the fold 'CaDensity'"(like../nfs/amino-home/zhng/zhanglab_programs/EMCMR/)


how to run this program?
./UseDensityx4_align input_PDB_file(like..T0914_dom2.pdb) input_MRC_file(like..T0914_dom2.mrc) density_map_resolation(like..6.0) sampling_rate(like..0.0)  date_static_file(binddir) PDB_name contact_map_file_fold

how to complie this program?
/opt/rh/devtoolset-7/root/usr/bin/g++ MC_EMSD33_00303_32_12_partE_test_50.cpp CaDensityx1_2.cpp SplineInterp.cpp GeometryTools.cpp readpdb.cpp GetVdwRadius.cpp xray_scattering.cpp Rotbuilder.cpp randomx.cpp -o UseDensityx4_align -I/(the location of data file)../CaDensity -L/(the location of data file)../CaDensity/fourier -lfourierx -O3

Our data set are available at http://182.92.205.0/zhangbiao/EMCMR-bm.zip

