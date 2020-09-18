/**
 * Example code for linking against astro-accelerate library.
 *
 * Compile with: g++ -std=c++11 -I/path/to/astro-accelerate/include/ -L/path/to/astro-accelerate/build -Wl,-rpath,/path/to/astro-accelerate/build -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion.cpp -o test

Compile with:
g++ -std=c++11 -I/home/guest/abhinav/FRB_pipeline-1.7.10/include/ -L/home/guest/abhinav/FRB_pipeline-1.7.10/build/ -Wl,-rpath,/home/guest/abhinav/FRB_pipeline-1.7.10/build/ -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion_and_analysis.cpp -o test_a

On TAPTI server, compile with:
g++ -std=c++11 -I/Data/rraj/abhinav/FRB_pipeline-1.7.10/include/ -L/Data/rraj/abhinav/FRB_pipeline-1.7.10/build/ -Wl,-rpath,/Data/rraj/abhinav/FRB_pipeline-1.7.10/build/ -lastroaccelerate -I/usr/local/cuda/include/ -L/usr/local/cuda/lib64 -Wl,-rpath,/usr/local/cuda/lib64 -I/usr/local/cuda-8.0/samples/common/inc/ -lcudart dedispersion_and_analysis.cpp -o test_a

Run the script using:
./astro-accelerate.sh ../input_files/GMRT.SPS.10m.txt /home/guest/abhinav/FRB_pipeline-1.7.10/output > log_file.txt

On TAPTI server, run the script using:
./astro-accelerate.sh ../input_files/GMRT.SPS.10m.txt /Data/rraj/abhinav/FRB_pipeline-1.7.10/output > log_file.txt


 */

#define _POSIX_SOURCE
#include <unistd.h>
#include <fstream> // for file-access
#include <string>
#include <vector>
#include <iostream>

#include "aa_ddtr_plan.hpp"
#include "aa_ddtr_strategy.hpp"
#include "aa_filterbank_metadata.hpp"
#include "aa_permitted_pipelines_2.hpp"
#include "aa_analysis_plan.hpp"
#include "aa_analysis_strategy.hpp"
#include "aa_log.hpp"
#include "aa_sigproc_input.hpp"
#include "aa_gpu_timer.hpp"
#include "frb_shm.h"   //Added  by Rushikesh
#include "externalLibraries.h"

#include <cstdlib>
#include <cstdio>
#include <array>

#include <wordexp.h>

#include <sys/types.h>
#include <sys/wait.h>

#include <sstream>

#include <sys/stat.h>  // to create directory

#include <algorithm>  // to remove new line charcter from string

//-----Added by Rushikesh----

using namespace std;

//----- End ----

using namespace astroaccelerate;

std::string timestamp(std::string source_fname);

//int main() {
int main(int argc,  char **argv) {

  unsigned long int f_pos;  // input file poistion variable

  int sysint, sysint1, sysint2, sysint3, sysint4, sysint5, sysint6;	// this is to store value returned by system command when script is called at the end of file.

// variables to read data from input file
  float sigma_cutoff_input, sigma_constant_input, max_boxcar_width_in_sec_input, periodicity_sigma_cutoff_input, periodicity_harmonics_input, max_dm_val;
  bool analysis_val, baselinenoise_val, debug_val;
  char file_name_input[128], config_path[128]; 
  baselinenoise_val = false; // yes, it is being used
  analysis_val = false; // not being used now
  debug_val = false; // not being used now

// buf count loop variables
  int buf_count = 0;
  int max_buf_count = 8;
  int new_file_flag = 1;

// vaiable to enable storing of each output global_peaks.dat file 
  int test_flg = 1;

  char cwd[256]; // to store current working directory

  std::string rm_cmnd, dir_cmnd, cat_cmnd;

// timing variables
    float time_tmp, time_first_read, time_buf_read, time_buf_process, time_conc_files_find_frb, time_total;
    aa_gpu_timer       timer_var, timer_var_total;

// timestamp of radec.dat
  std::string result_tmp, result_time_stamp;

// varibales that store the amount of data processed
   int nsamp_ovrlp_pt, nsamp_ovrlp_val;
   float tsamp_unique;


// variables to write the python file name
  std::stringstream stream;
  std::string s_ra, s_dec, s_mjd, s_bcount, py_fname, s_tsamp_unique, s_new_file_flag;

// finding the intial precision to be set later
  std::streamsize ss = std::cout.precision();
  std::cout << "Initial precision = " << ss << '\n';


  aa_ddtr_plan ddtr_plan;

// reading the input file for dm range and other input data::----------------------------------------------------
// extracting the input file name from the argument passed through the terminal to the astroaccelerate.sh bash script
  if (argc > 1) {
        printf("argv[1] = %s", argv[1]); printf("argv[2] = %s", argv[0]);
    } else {
        printf("No file name entered. Exiting...");
        return 0;
    }

// added a new function named read_metadata_input in aa_sigproc_input to read input file data
  aa_sigproc_input       filterbank_datafile_input(argv[1]);
  filterbank_datafile_input.read_metadata_input(&sigma_cutoff_input, &sigma_constant_input, &max_boxcar_width_in_sec_input, &periodicity_sigma_cutoff_input, &periodicity_harmonics_input, &ddtr_plan, &baselinenoise_val, &analysis_val, &debug_val, &file_name_input, &config_path);
  printf("\n ddtr_plan.user_dm(ddtr_plan.range()-1).high: %f \n", ddtr_plan.user_dm(ddtr_plan.range()-1).high);
  std::string file_name_input_str;

// input file read ----------------------------------------------------------------------------------------------
  printf("\n file_name_input :: %s \n", file_name_input);
  printf("\nbaselinenoise_val: is: %d \n", baselinenoise_val);

// assigning variables from the input file:
  const float sigma_cutoff = sigma_cutoff_input;
  const float sigma_constant = sigma_constant_input;
  const float max_boxcar_width_in_sec = max_boxcar_width_in_sec_input;

  const aa_analysis_plan::selectable_candidate_algorithm algo = aa_analysis_plan::selectable_candidate_algorithm::off; // not changed, being used as it was in the example file
// following constants added for FRB file (not being used now):
  const float periodicity_sigma_cutoff = periodicity_sigma_cutoff_input;
  const float periodicity_harmonics = periodicity_harmonics_input;

  

// input data file location and config file path ------------------------------------------------------------------------------------------------------------------------------

  std:: string fname_input(file_name_input);
  printf("\n fname_input::%s::test\n", fname_input.c_str());
  fname_input.erase(std::remove(fname_input.begin(), fname_input.end(), '\n'), fname_input.end());
  printf("\n fname_input::%s::test\n", fname_input.c_str());

  aa_sigproc_input       filterbank_datafile(fname_input);  

  filterbank_datafile.open();

  std:: string configpath_name(config_path);
  printf("\n configpath_name::%s::test\n", configpath_name.c_str());
  configpath_name.erase(std::remove(configpath_name.begin(), configpath_name.end(), '\n'), configpath_name.end());
  printf("\n configpath_name::%s::test\n", configpath_name.c_str());

  std::string source_fname, gpu_fname;
  source_fname = configpath_name + "radec.dat"; //RD: Feb 7 CHECK
  gpu_fname = configpath_name + "gpu.hdr";

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// reading metadata from gpu.hdr file
  aa_sigproc_input       filterbank_datafile_config(gpu_fname);
  aa_filterbank_metadata filterbank_metadata = filterbank_datafile_config.read_metadata(ddtr_plan.user_dm(ddtr_plan.range()-1).high, configpath_name);

// variables assigned after reading gu.hdr file
  const double tstart = filterbank_metadata.tstart();
  const double tsamp = filterbank_metadata.tsamp();
  const double nbits = filterbank_metadata.nbits();
  const double nsamples = filterbank_metadata.nsamples();
  const double fch1 = filterbank_metadata.fch1();
  const double foff = filterbank_metadata.foff();
  const double nchans = filterbank_metadata.nchans();
  
  aa_filterbank_metadata metadata(tstart, tsamp, nbits, nsamples, fch1, foff, nchans);

  printf("\n filterbank_metadata.ovrlp(): %f \n", filterbank_metadata.ovrlp());

// deciding how many time samples will be overlapped in each reading.
  if (filterbank_metadata.ovrlp() > 59){
      if (nsamples == 2*1024*1024) nsamp_ovrlp_val = 768*1024;	
      else if(nsamples == 1*1024*1024) nsamp_ovrlp_val = 384*1024;
      else if(nsamples == 512*1024) nsamp_ovrlp_val = 192*1024;
      else if(nsamples == 256*1024) nsamp_ovrlp_val = 96*1024;
      else if(nsamples == 128*1024) nsamp_ovrlp_val = 48*1024;
  }
  else{
      if (nsamples == 2*1024*1024) nsamp_ovrlp_val = 256*1024;	//RD: Updated for band 4 and 5.
      else if(nsamples == 1*1024*1024) nsamp_ovrlp_val = 128*1024;
      else if(nsamples == 512*1024) nsamp_ovrlp_val = 64*1024;
      else if(nsamples == 256*1024) nsamp_ovrlp_val = 32*1024;
      else if(nsamples == 128*1024) nsamp_ovrlp_val = 16*1024;
  }

  nsamp_ovrlp_pt = nsamples - nsamp_ovrlp_val;
  tsamp_unique = nsamp_ovrlp_pt*tsamp;
  printf("\n number of time samples to be processed only once (without overlap): %d \n", nsamp_ovrlp_pt);
  printf("\n amount of telescope time to be processed only once (without overlap): %f \n", tsamp_unique);
  

//  const size_t free_memory = 7000000000; // Free memory on the GPU in bytes
  const size_t free_memory = 12000000000; // Free memory on the GPU in bytes


  bool enable_analysis = true;       // The strategy will be optimised to run just dedispersion
  aa_ddtr_strategy ddtr_strategy(ddtr_plan, metadata, free_memory, enable_analysis); //RD : Look into more detail.

// read the timestamp of radec.dat
  result_tmp = timestamp(source_fname);
  printf("\n first check: timestamp of radec.dat is: %s \n", result_tmp.c_str());

   
  if(!(ddtr_strategy.ready())) {
    LOG(log_level::error, "ddtr_strategy not ready.");
    return 0;
  }
  
  aa_analysis_plan analysis_plan(ddtr_strategy, sigma_cutoff, sigma_constant, max_boxcar_width_in_sec, algo, true);
  aa_analysis_strategy analysis_strategy(analysis_plan); //RD : Look into more detail.

  if(!(analysis_strategy.ready())) {
    LOG(log_level::error, "analysis_strategy not ready.");
    return 0;
  }

// reading first buffer:

  timer_var_total.Start(); // starting timer for total time taken

  timer_var.Start();  // starting timer for first read

//---- File name for raw file for 1K@3ms ------Rushikesh edits
    stream << std::fixed << std::setprecision(6) << filterbank_metadata.src_raj();
    s_ra = stream.str();
    stream.str(std::string());
    stream.clear();

    stream << std::fixed << std::setprecision(6) << filterbank_metadata.src_dej();
    s_dec = stream.str();
    stream.str(std::string());
    stream.clear();

    stream << std::fixed << std::setprecision (ss) <<filterbank_metadata.tstart();
    s_mjd = stream.str();
    stream.str(std::string());
    stream.clear();

    system("pwd > ./pwd.txt");
    ifstream ifs("./pwd.txt");
    std::string line, path;
   
    while (getline(ifs, line))
         path += line;
    printf("Output file path :: %s \n", path.c_str());
    

    std::string FileNameRaw;
    std::string FilePath;
    FileNameRaw = "Rawfile_1K@3ms_" + s_mjd + "_" + s_ra + "_" + s_dec + ".raw";
    FilePath = path + "/" + FileNameRaw;
    printf("##### FilePath ::: %s", FilePath.c_str());

//------End------

    if(!filterbank_datafile.read_new_signal(buf_count, filterbank_metadata, &f_pos, new_file_flag, FilePath)) {
      std::cout << "ERROR: Could not read telescope data." << std::endl;
      return 0;
    }
  timer_var.Stop();

  time_first_read = timer_var.Elapsed() / 1000;
  printf("\n\n === First buffer read ===\n");
  //// RD: Check ////
  aa_permitted_pipelines_2<aa_pipeline::component_option::zero_dm, false> runner(ddtr_strategy, analysis_strategy, filterbank_datafile.input_buffer().data());


//  bool dump_to_disk = false;
  bool dump_to_disk = true;
  bool dump_to_user = true;
  std::vector<analysis_output> output;

/*
// making small edit to source.hdr to check if the timestep check is working well:
// only uncomment it to test weather the timestamp changes are being picked up!

    FILE * pFile_source_loc;

//    std::ofstream ofs;
//    ofs.open ("/home/guest/abhinav/FRB_pipeline-1.7.10/input_files/source.hdr");
    pFile_source_loc = fopen("/home/guest/abhinav/FRB_pipeline-1.7.10/input_files/source.hdr", "w");
    if(pFile_source_loc == NULL) {
      printf("\n ATTENTION:: source.hdr source cannot be opened! \n");
    }
   else {
   printf("\n File opened :: source.hdr for editing to simulate that the code will pick up changes in timestamp \n");
   fprintf(pFile_source_loc, "src_raj: 214400.0\n");
   fprintf(pFile_source_loc, "src_dej: -393300.0\n");
    }  // else condition ends -----------------
   if(!fclose(pFile_source_loc)) printf("\n source.hdr file successfully closed \n");
*/

//---------------------------------------------------------------------------------------------------------------------

  
  if(runner.setup()) {

  for (buf_count=0; true ; buf_count++)
  {

    printf("\n buf count:: :: %d \n", buf_count);

// read the timestamp of radec.dat
    if (buf_count==0){
    printf("\n source_fname:: :: %s \n", source_fname.c_str());  
    result_time_stamp = timestamp(source_fname);
    printf("\n buf_count: %d timestamp is: %s \n", buf_count, result_time_stamp.c_str());
    if(result_tmp != result_time_stamp)
        {    
         	printf("\n radec.dat has been modified. Time stamps do not match. Will read source information again. \n");
                aa_sigproc_input       filterbank_datafile_config(source_fname);
                aa_filterbank_metadata filterbank_metadata_source = filterbank_datafile_config.read_metadata_source(filterbank_metadata);
	        filterbank_metadata = filterbank_metadata_source;
	        printf("\n filterbank_metadata.src_raj(): %f \n", filterbank_metadata.src_raj());
                //----Rushikesh edits
	        new_file_flag = 1;
		printf(" radec.dat has been changed, creating a new filename\n");
		FileNameRaw = "Rawfile_1K@3ms_" + s_mjd + "_" + s_ra + "_" + s_dec + ".raw";
                FilePath = path + FileNameRaw;
    		printf("##### new FilePath ::: %s", FilePath.c_str());
		new_file_flag = 1;
		//-----end---
        }
    else	
        {
	        printf("\n radec.dat is still the same. Time stamps match \n"); 
		new_file_flag = 0;
        }
    }

// removing files from previous run and creating directories
// if test_flg is true, it will create a bcount$ directory for each buffer ($ = buf_count)    
    if (test_flg)
    {
    stream << buf_count;   
    s_bcount = stream.str();
    stream.str(std::string());
    stream.clear();

    rm_cmnd = "rm ./test/bcount" + s_bcount + "/global*";
    dir_cmnd = "test/bcount" + s_bcount;
    printf("\n buf count remove command: %s \n", rm_cmnd.c_str());
    printf("\n buf count mkdir command: %s \n", dir_cmnd.c_str());
    sysint2 = system("pwd");
    sysint2 = system("mkdir test/");
    mkdir(dir_cmnd.c_str(), ACCESSPERMS); 
    sysint2 = system(rm_cmnd.c_str());
    }

    sysint1 = system("rm -f analysed*"); 
    sysint2 = system("rm -f acc*");
    sysint3 = system("rm -f global*");
    sysint4 = system("rm -f fourier*");
    sysint5 = system("rm -f harmonic*");
    sysint6 = system("rm -f candidate*");
    sysint = system("rm -f peak*");
    printf("\n all files removed \n");


// read only if the buffer number is greater than 0 - because first buffer has been read above already:
    if(buf_count > 0) {
    printf("\n we are here with buf_count: %d \n", buf_count);

    timer_var.Start();
    filterbank_datafile.open();
    if(!filterbank_datafile.read_new_signal(buf_count, filterbank_metadata, &f_pos, new_file_flag, FilePath)) {
      std::cout << "ERROR: Could not read telescope data." << std::endl;
      return 0;
    }  

    aa_permitted_pipelines_2<aa_pipeline::component_option::zero_dm, false> runner(ddtr_strategy, analysis_strategy, filterbank_datafile.input_buffer().data());

    timer_var.Stop();
    time_buf_read = timer_var.Elapsed() / 1000;
    printf("\n\n === Next buffer read ===\n");
    } //RD: 'buf_count > 0' loop ends here

    timer_var.Start();

    while(runner.next(dump_to_disk, dump_to_user, output)) {
      LOG(log_level::notice, "Pipeline running over next chunk.");

      for(size_t i = 0; i < output.size(); i++) {
	std::cout << " dm_low " << output.at(i).dm_low << " dm_high " << output.at(i).dm_high << std::endl;
	for(size_t j = 0; j < output.at(i).pulses.size(); j++) {
	  analysis_pulse p = output.at(i).pulses.at(j);
	  std::cout << "dispersion measure " << p.dispersion_measure	 << " time " << p.time << " snr " << p.snr << " pulse_width " << p.pulse_width << std::endl;
	}
      }
    }

    timer_var.Stop();
    time_buf_process = timer_var.Elapsed() / 1000;
    printf("\n\n === Buffer processed ===\n");



// creating python output file name
    printf("\n RA: %f DEC: %f MJD: %f \n", filterbank_metadata.src_raj(), filterbank_metadata.src_dej(), filterbank_metadata.tstart());

    stream << new_file_flag;
    s_new_file_flag = stream.str();
    stream.str(std::string());
    stream.clear();

    stream << buf_count;
    s_bcount = stream.str();
    stream.str(std::string());
    stream.clear();

    stream << std::fixed << std::setprecision(6) << filterbank_metadata.src_raj();
    s_ra = stream.str();
    stream.str(std::string());
    stream.clear();

    stream << std::fixed << std::setprecision(6) << filterbank_metadata.src_dej();
    s_dec = stream.str();
    stream.str(std::string());
    stream.clear();

    stream << std::fixed << std::setprecision (ss) <<filterbank_metadata.tstart();
    s_mjd = stream.str();
    stream.str(std::string());
    stream.clear();

    stream << std::fixed << std::setprecision (ss) <<tsamp_unique;
    s_tsamp_unique = stream.str();
    stream.str(std::string());
    stream.clear();

    printf("\n After setting precision and string conversion: RA: %s DEC: %s MJD: %s buf_count: %s \n", s_ra.c_str(), s_dec.c_str(), s_mjd.c_str(), s_bcount.c_str());

    py_fname = "python "; // --bc 0 --ra 2.5 --dec -36 --mjd 55.5";
    py_fname += configpath_name;
    py_fname += "spsplotii_frb_detect.py global_peaks.dat";
    py_fname += " --mjd " + s_mjd + " --ra " + s_ra + " --dec " + s_dec + " --bc " + s_bcount + " --buf_t " + s_tsamp_unique + " --new_file_flag " + s_new_file_flag;
    printf("\n py_fname: %s \n", py_fname.c_str());

    // writing output files and running python script
    timer_var.Start();

     sysint1 = system("cat peak* > global_peaks.dat");
     sysint2 = system("cat analysed* > global_analysed_frb.dat");
     sysint3 = system("cat fourier-* > global_periods.dat");
     sysint4 = system("cat fourier_inter* > global_interbin.dat");
     sysint5 = system("cat harmo* > global_harmonics.dat");
     sysint6 = system("cat candidate* > global_candidates.dat");
     sysint = system(py_fname.c_str());

/*-------RD: To write python output file -------
     const char * ls_args[2] = { py_fname.c_str(), NULL} ;
     pid_t c_pid, pid;
     int status;

     c_pid = fork();

     if (c_pid == 0){
       /* CHILD */

/*       printf("Child: executing ls\n");
       execvp( ls_args[0], ls_args);
                                                                                                                                            
       perror("execve failed");
       }else if (c_pid > 0){
        /* PARENT */

/*        if( (pid = wait(&status)) < 0){
          perror("wait");
          _exit(1);
        }

        printf("Parent: finished\n");

      }else{
        perror("fork failed");
        _exit(1);
      }

//--------End------------------------*/

     printf("\nOutput files written\n");

// if test_flg is true, it will save every global_peaks.dat file from every buffer    
     if (test_flg)
     {
     cat_cmnd = "cat peak* > ./test/bcount" + s_bcount + "/global_peaks.dat";
     printf("\n Concatenation command:: %s \n", cat_cmnd.c_str());
     sysint2 = system(cat_cmnd.c_str());
     }



    timer_var.Stop();
    time_conc_files_find_frb = timer_var.Elapsed() / 1000;
    printf("\n\n === Output files created | Concatenated | Python script ran ===\n");

   if (buf_count == 0) 
   {
     printf("\n--------------------- Timing information after processing first buffer ---------------------\n");
     printf("\n Time to read first buffer: %f (GPU estimate) \n", time_first_read);
     printf("\n Time to process buffer number %d: %f (GPU estimate) \n", buf_count, time_buf_process);
     printf("\n Time to concatenate files and run python script for buffer number %d: %f (GPU estimate) \n", buf_count, time_conc_files_find_frb);
     printf("\n--------------------------------------------------------------------------------------------\n");
   }
   else
   {
     printf("\n--------------------- Timing information after processing new buffer ---------------------\n");
     printf("\n Time to read buffer number %d: %f (GPU estimate) \n", buf_count, time_buf_read);
     printf("\n Time to process buffer number %d: %f (GPU estimate) \n", buf_count, time_buf_process);
     printf("\n Time to concatenate files and run python script for buffer number %d: %f (GPU estimate) \n", buf_count, time_conc_files_find_frb);
     printf("\n-------------------------------------------------------------------------------------------\n");
   }


  }// buf_count loop ends here ----------------------------------------------------------------------------------

  }// setup condition ends here ------------------------------------------------------------------------------



  timer_var_total.Stop();
  time_total = timer_var_total.Elapsed() / 1000;

  printf("\n------------------------------- Final Timing Information  ---------------------------------\n");
  printf("\n Total number of buffers processed: %d \n", buf_count);
  printf("\n Total number of time samples processed uniquely: %d \n", buf_count*nsamp_ovrlp_pt);
  printf("\n Total telescope time processed: %f (GPU estimate) \n", buf_count*nsamp_ovrlp_pt*tsamp);
  printf("\n Total time taken: %f (GPU estimate) \n", time_total);
  printf("\n Total Real-time speedup factor: %f \n",  (buf_count*nsamp_ovrlp_pt*tsamp) / ( time_total ));
  printf("\n-------------------------------------------------------------------------------------------\n");

 
  LOG(log_level::notice, "Finished.");
  return 0;
} //runner()function ends here
std::string timestamp(std::string source_fname)
{
    std::string command_part = "stat -c '%y' ";
    command_part += source_fname;
// checking the last modified time of radec.dat ---------------------------------------------------
    std::string command = "stat -c '%y' ";
    command += source_fname;

//    std::string command("pwd");

    std::array<char, 128> buffer;
    std::string result, result_tmp;

    std::cout << "Inside timestamp function: Opening reading pipe" << std::endl;

    printf(" command ::::::: %s \n", command.c_str()); 
    FILE* pipe = popen(command.c_str(), "r"); 
    if (!pipe)
    {
        std::cerr << "Couldn't start command." << std::endl;
        return 0;
    }
    while (fgets(buffer.data(), 128, pipe) != NULL) {
        std::cout << "Reading timestamp" << std::endl;
        result += buffer.data();
        result_tmp += buffer.data();
    }
    auto returnCode = pclose(pipe);

    std::cout << result << std::endl;
    std::cout << returnCode << std::endl;
    return result;
// check complete ----------------------------------------------------------------------------------

}

