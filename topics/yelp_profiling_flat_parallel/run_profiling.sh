CURRENT_DIR=$(pwd)
cd ../
OUTPUTFILENAME=profiling_yelp_8_cores_1 ./ompsub_patterned_name.sh -n 8 julia -p 8 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_8_cores_2 ./ompsub_patterned_name.sh -n 8 julia -p 8 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_8_cores_3 ./ompsub_patterned_name.sh -n 8 julia -p 8 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_7_cores_1 ./ompsub_patterned_name.sh -n 7 julia -p 7 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_7_cores_2 ./ompsub_patterned_name.sh -n 7 julia -p 7 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_7_cores_3 ./ompsub_patterned_name.sh -n 7 julia -p 7 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_6_cores_1 ./ompsub_patterned_name.sh -n 6 julia -p 6 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_6_cores_2 ./ompsub_patterned_name.sh -n 6 julia -p 6 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_6_cores_3 ./ompsub_patterned_name.sh -n 6 julia -p 6 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_5_cores_1 ./ompsub_patterned_name.sh -n 5 julia -p 5 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_5_cores_2 ./ompsub_patterned_name.sh -n 5 julia -p 5 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_5_cores_3 ./ompsub_patterned_name.sh -n 5 julia -p 5 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_4_cores_1 ./ompsub_patterned_name.sh -n 4 julia -p 4 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_4_cores_2 ./ompsub_patterned_name.sh -n 4 julia -p 4 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_4_cores_3 ./ompsub_patterned_name.sh -n 4 julia -p 4 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_3_cores_1 ./ompsub_patterned_name.sh -n 3 julia -p 3 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_3_cores_2 ./ompsub_patterned_name.sh -n 3 julia -p 3 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_3_cores_3 ./ompsub_patterned_name.sh -n 3 julia -p 3 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_2_cores_1 ./ompsub_patterned_name.sh -n 2 julia -p 2 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_2_cores_2 ./ompsub_patterned_name.sh -n 2 julia -p 2 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_2_cores_3 ./ompsub_patterned_name.sh -n 2 julia -p 2 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_1_cores_1 ./ompsub_patterned_name.sh -n 1 julia -p 1 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_1_cores_2 ./ompsub_patterned_name.sh -n 1 julia -p 1 driver_profile_enron.jl
OUTPUTFILENAME=profiling_yelp_1_cores_3 ./ompsub_patterned_name.sh -n 1 julia -p 1 driver_profile_enron.jl

condor_wait profiling_yelp_8_cores_1.l
condor_wait profiling_yelp_8_cores_2.l
condor_wait profiling_yelp_8_cores_3.l
condor_wait profiling_yelp_7_cores_1.l
condor_wait profiling_yelp_7_cores_2.l
condor_wait profiling_yelp_7_cores_3.l
condor_wait profiling_yelp_6_cores_1.l
condor_wait profiling_yelp_6_cores_2.l
condor_wait profiling_yelp_6_cores_3.l
condor_wait profiling_yelp_5_cores_1.l
condor_wait profiling_yelp_5_cores_2.l
condor_wait profiling_yelp_5_cores_3.l
condor_wait profiling_yelp_4_cores_1.l
condor_wait profiling_yelp_4_cores_2.l
condor_wait profiling_yelp_4_cores_3.l
condor_wait profiling_yelp_3_cores_1.l
condor_wait profiling_yelp_3_cores_2.l
condor_wait profiling_yelp_3_cores_3.l
condor_wait profiling_yelp_2_cores_1.l
condor_wait profiling_yelp_2_cores_2.l
condor_wait profiling_yelp_2_cores_3.l
condor_wait profiling_yelp_1_cores_1.l
condor_wait profiling_yelp_1_cores_2.l
condor_wait profiling_yelp_1_cores_3.l

for f in profiling_yelp_*.o; do mv $f $CURRENT_DIR/${f%.o}.txt; done;
rm *.sub *.e *.l
cd $CURRENT_DIR