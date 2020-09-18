/**
 * Example code for linking against astro-accelerate library.
 * 
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion.cpp -o test

Compile with:
g++ -std=c++11 -I/home/guest/abhinav/FRB_pipeline-1.7.10/include/ -L/home/guest/abhinav/FRB_pipeline-1.7.10/build/ -Wl,-rpath,/home/guest/abhinav/FRB_pipeline-1.7.10/build/ -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart filterbank_dedispersion.cpp -o test
*/

#include "aa_sigproc_input.hpp"
#include "aa_pipeline_wrapper_functions.hpp"

using namespace astroaccelerate;

int main() {
//  aa_sigproc_input       filterbank_datafile("/mnt/data/AstroAccelerate/filterbank/BenMeerKAT.fil");
//  aa_sigproc_input       filterbank_datafile("/data/jroy/data/HR_1850-48_ia_500_200_4096_4_1_8_14may2018.raw.fil");

// reading header information from config file
//  aa_sigproc_input       filterbank_datafile_config("/home/guest/abhinav/FRB_pipeline/input_files/config/FRB_DM100_163.84us_4K_7mP_4msW_10sig.header.config");

  aa_sigproc_input       filterbank_datafile_config("/home/guest/Rushikesh/FRB_pipeline-1.7.10.GMRT300-500/input_files/gpu.hdr");

//  aa_sigproc_input       filterbank_datafile("/data/jroy/data/FRB_DM1000_163.84us_4K_7mP_4msW_3sig.header.fil");
//  aa_sigproc_input       filterbank_datafile_config("/home/guest/abhinav/FRB_pipeline/input_files/config/HR_1850-48_ia_500_200_4096_4_1_8_14may2018.raw.config");

  aa_filterbank_metadata filterbank_metadata = filterbank_datafile_config.read_metadata();
//  aa_filterbank_metadata filterbank_metadata = filterbank_datafile.read_metadata();



//  aa_sigproc_input       filterbank_datafile("/data/jroy/data/FRB_DM1000_163.84us_4K_7mP_4msW_3sig.header.fil");
  aa_sigproc_input       filterbank_datafile("/data/jroy/data/FRB_samples/FRB_band3_DM1000_327.68us_4K_7mP_1msW_10sig.raw");

  filterbank_datafile.open();

  std::cout << "Test hai ye" << std::endl;

  aa_filterbank_metadata filterbank_metadata_config1 = filterbank_datafile_config.read_metadata_new(filterbank_metadata);

//  aa_filterbank_metadata filterbank_metadata_config1 = filterbank_datafile_config.read_metadata_new(filterbank_metadata);
  
//  printf("\n filterbank_dedispersion.cpp code ::  filterbank_metadata.strat:: is :: %d \n", filterbank_metadata.strat());
//  printf("\n filterbank_dedispersion.cpp code :: filterbank_metadata_config1.strat:: is :: %d \n", filterbank_metadata_config1.strat());


  int buf_count;
  unsigned long int f_pos;  
  
  for (buf_count=0; buf_count < 1; buf_count++)
  {
  if(!filterbank_datafile.read_new_signal(buf_count, filterbank_metadata, &f_pos)) {
    std::cout << "ERRORRR: Could not read telescope data." << std::endl;
    return 0;
  }

  std::vector<aa_ddtr_plan::dm> dm_ranges;
  aa_ddtr_plan::dm range1 = {0, 150, 0.1, 1, 1};
  aa_ddtr_plan::dm range2 = {150, 300, 0.2, 1, 1};
  aa_ddtr_plan::dm range3 = {300, 500, 0.25, 1, 1};
  aa_ddtr_plan::dm range4 = {500, 900, 0.4, 2, 2};
  aa_ddtr_plan::dm range5 = {900, 1200, 0.6, 4, 4};
  aa_ddtr_plan::dm range6 = {1200, 1500, 0.8, 4, 4};
  aa_ddtr_plan::dm range7 = {1500, 2000, 1.0, 4, 4};

/*
  aa_ddtr_plan::dm range3 = {300, 500, 0.25, 2, 2};
  aa_ddtr_plan::dm range4 = {500, 900, 0.4, 4, 4};
//  aa_ddtr_plan::dm range5 = {900, 1200, 0.6, 4, 4};
//  aa_ddtr_plan::dm range6 = {1200, 1500, 0.8, 8, 8};
//  aa_ddtr_plan::dm range7 = {1500, 2000, 1.0, 8, 8};
//  aa_ddtr_plan::dm range8 = {2000, 3000, 2.0, 8, 8};
*/


  dm_ranges.push_back(range1);
  dm_ranges.push_back(range2);
  dm_ranges.push_back(range3);
  dm_ranges.push_back(range4);
  dm_ranges.push_back(range5);
  dm_ranges.push_back(range6);
  dm_ranges.push_back(range7);
//  dm_ranges.push_back(range8);

  const aa_pipeline::pipeline_option pipeline_options = {aa_pipeline::component_option::zero_dm};

  printf("\n Buffer count is %d. Sending strat = %d \n", buf_count, filterbank_metadata.strat());
  dedisperse_telescope_data(filterbank_metadata, pipeline_options, dm_ranges, filterbank_datafile.input_buffer());  

  if (buf_count == 0) { 
  printf("\n filterbank_dedispersion.cpp code :: filterbank_metadata_config1.strat:: is :: %d \n", filterbank_metadata_config1.strat());
  filterbank_metadata = filterbank_metadata_config1;  
  }


  } // buf_count loop ends here----------------------------------------------------------
  std::cout << "NOTICE: Finished." << std::endl;


  return 0;

}
