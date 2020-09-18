#include "aa_ddtr_strategy.hpp"

#include "aa_device_info.hpp"

#include <cstring>

#define MAX_OUTPUT_SIZE 8589934592

namespace astroaccelerate {
  /**
   * Trivial constructor for aa_ddtr_strategy which will result in a strategy that cannot be ready. 
   */
  aa_ddtr_strategy::aa_ddtr_strategy() : m_ready(false), m_strategy_already_calculated(false), m_configured_for_analysis(false), is_setup(false), m_maxshift(0), m_num_tchunks(0), m_total_ndms(0), m_max_dm(0.0), m_maxshift_high(0), m_max_ndms(0), m_power(0.0), m_enable_msd_baseline_noise(false) {
    
  }

  /**
   * Constructor for aa_ddtr_strategy that computes the strategy upon construction, and sets the ready state of the instance of the class.
   */
  aa_ddtr_strategy::aa_ddtr_strategy(const aa_ddtr_plan &plan, const aa_filterbank_metadata &metadata, const size_t &free_memory, const bool &enable_analysis) : m_ready(false), m_strategy_already_calculated(false), m_configured_for_analysis(enable_analysis), is_setup(false), m_metadata(metadata), m_maxshift(0), m_num_tchunks(0), m_total_ndms(0), m_max_dm(0.0), m_maxshift_high(0), m_max_ndms(0), m_power(plan.power()), m_enable_msd_baseline_noise(plan.enable_msd_baseline_noise()) {    
    strategy(plan, free_memory, enable_analysis);
  }

  /**
   * Implementation of the function that computes the strategy.
   */
  bool aa_ddtr_strategy::strategy(const aa_ddtr_plan &plan, const size_t &free_memory, const bool &enable_analysis) {
    /**
     * This method relies on defining points when nsamps is a multiple of
     * nchans - bin on the diagonal or a fraction of it.
     */

    printf("\n 1. in aa_ddtr_strategy.cpp :: m_metadata.strat(): %d \n", m_metadata.strat());

    if(m_strategy_already_calculated) {
//    if(m_strategy_already_calculated || m_metadata.strat()) { 						// also my addition -- abhinav

/*      m_ready = true;											// look right below
      m_strategy_already_calculated = true;								// my additions  -- abhinav
      printf("\n 2. in aa_ddtr_strategy.cpp :: m_metadata.strat(): %d \n", m_metadata.strat());		// look right above


      FILE * pFile;
      pFile = fopen ("myfile.txt","r");

      char strat_line[512], tempo_val[400];
      int range, scruch_value;
      float max_disp_delay, diagonal_dm;
//      const size_t range;

      while (fgets(strat_line, sizeof(strat_line), pFile)) {
        printf("\n **-*-** %.35s \n", strat_line);

        strcpy(tempo_val, &strat_line[17]);	
        if (strncmp(strat_line, "m_maxshift     : ", 17) == 0 )              // Time stamp of first sample (MJD) 
		{
		m_maxshift = strtol(tempo_val, NULL, 10);
		printf("\n m_masxshift from file: %s :: %d \n", tempo_val, m_maxshift);
		}
        if (strncmp(strat_line, "m_max_ndms     : ", 17) == 0 )              // Time stamp of first sample (MJD) 
		{
		m_max_ndms = strtol(tempo_val, NULL, 10);
		printf("\n m_max_ndms from file: %s :: %d \n", tempo_val, m_max_ndms);
		}
        if (strncmp(strat_line, "range          : ", 17) == 0 )              // Time stamp of first sample (MJD) 
		{
		range = strtol(tempo_val, NULL, 10);
		printf("\n range from file: %s :: %d \n", tempo_val, range);
		}
        if (strncmp(strat_line, "scrunch value  : ", 17) == 0 )              // Time stamp of first sample (MJD) 
		{
		scruch_value = strtol(tempo_val, NULL, 10);
		printf("\n scrunch value: %s :: %d \n", tempo_val, scruch_value);
		}

        if (strncmp(strat_line, "Max disp delay : ", 17) == 0 )              // Time stamp of first sample (MJD) 
		{
		max_disp_delay = strtof(tempo_val, NULL);
		printf("\n Max disp delay: %s :: %f \n", tempo_val, max_disp_delay);
		}
        if (strncmp(strat_line, "Diagonal DM    : ", 17) == 0 )              // Time stamp of first sample (MJD) 
		{
		diagonal_dm = strtof(tempo_val, NULL);
		printf("\n Diagonal DM: %s :: %f \n", tempo_val, diagonal_dm);
		}


        if (strncmp(strat_line, "m_num_tchunks  : ", 17) == 0 )              // Time stamp of first sample (MJD) 
		{
		m_num_tchunks = strtol(tempo_val, NULL, 10);
		printf("\n m_num_tchunks: %s :: %d \n", tempo_val, m_num_tchunks);
		}
        if (strncmp(strat_line, "dm limits ndm  : ", 17) == 0 )              // Time stamp of first sample (MJD) 
		{
//		printf("\n (l,h,s,i,o,n)  :: is > %s \n", tempo_val);
                float tmp_low, tmp_high, tmp_step;
                int tmp_inBin, tmp_outBin, tmp_ndms;
                char* ptr_line;
                for(size_t i = 0; i < range; i++) {
/*		  str_dm[i].low = strtof(tempo_val, &ptr_line);
		  str_dm[i].high = strtof(ptr_line, &ptr_line);
		  str_dm[i].step = strtof(ptr_line, &ptr_line);
		  str_dm[i].inBin = strtol(ptr_line, &ptr_line,10);
		  str_dm[i].outBin = strtol(ptr_line, &ptr_line,10);
		  m_ndms[i] = strtol(ptr_line, &ptr_line,10);
*/
//		  str_dm[i].low = strtof(tempo_val, NULL);
//	          fprintf (pFile, "%f %f %f %d %d %d\t",str_dm[i].low,str_dm[i].high,str_dm[i].step, str_dm[i].inBin,str_dm[i].outBin,m_ndms[i]);
/*		  printf("\n (l,h,s,i,o,n)  :  %f ::::::: %f \n", tmp_low, tmp_high);
		  printf("\n (l,h,s,i,o,n)  :  %f ::::::: %d \n", tmp_step, tmp_ndms);
 		  printf("\n (l,h,s,i,o,n)  :  %d ::::::: %d \n", tmp_inBin, tmp_outBin); 
*/
//                  }
//		}



//           fprintf (pFile, "m_num_tchunks  : %d \n", m_num_tchunks);
//           fprintf (pFile, "remainder      : %d \n", remainder);
/*	   fprintf (pFile, "(l,h,s,i,o,n)  : ");
           for(size_t i = 0; i < range; i++) {
           const aa_ddtr_plan::dm tmp = dm(i);
	   fprintf (pFile, "%f %f %f %d %d %d\t",tmp.low,tmp.high,tmp.step, tmp.inBin,tmp.outBin,ndms(i));
           }
           fprintf(pFile, "\n");
           fprintf (pFile, "t_procsd i j   : %d %d \n", t_processed().size(), t_processed().at(0).size());
           fprintf (pFile, "t_procsd arr   : ");
           for(size_t i = 0; i < t_processed().size(); i++) {
	     for(size_t j = 0; j < t_processed().at(i).size(); j++) {
               fprintf (pFile, "%d ", t_processed()[i][j]);
             }
           }
           fprintf(pFile, "\n");
*/

/*
      }

      fclose (pFile);

*/
      LOG(log_level::notice, "Strategy already calculated.");  // previous code ------------------------- 
      return true;
    }

    LOG(log_level::notice, "Calculating strategy."); //RD: What does it do?
	    
    //Part of the filterbank metadata
    const int nchans  = m_metadata.nchans();
    const int nsamp   = m_metadata.nsamples();
    const float fch1  = m_metadata.fch1();
    const float foff  = m_metadata.foff();
    const float tsamp = m_metadata.tsamp();

    printf("nchans :: %d\n",nchans);
    printf("nsamp :: %d\n",nsamp);
    printf("tsamp :: %f\n",tsamp);
    printf("fch1 :: %f\n",fch1);
    printf("foff :: %f\n",foff);

    
    if(!plan.range()) {
      //No user requested dm settings have been added, this is an invalid aa_ddtr_plan.
      return false;
    }
    
    //Plan requested DM settings
    const size_t range = plan.range();
    printf("\n range:: >> %d \n", range);
    m_ndms.resize(range);
    
    m_dmshifts.resize(nchans);
    
    //Strategy set DM settings
    str_dm.resize(range);

    const size_t gpu_memory = free_memory;
    printf("\n gpu_memory %d\n", gpu_memory);
    
    const double SPDT_fraction = 3.0/4.0; // 1.0 for MSD plane profile validation
    //Calculate maxshift, the number of dms for this bin and the highest value of dm to be calculated in this bin
    if (m_power != 2.0) {
      // Calculate time independent dm shifts
      for (int c = 0; c < nchans; c++) {
	m_dmshifts[c] = 4148.741601f * ( ( 1.0 / pow(( fch1 + ( foff * c ) ), m_power) ) - ( 1.0 / pow(fch1, m_power) ) );
      }
    }
    else {
      // Calculate time independent dm shifts
      for (int c = 0; c < nchans; c++) {
	m_dmshifts[c] = (float) ( 4148.741601f * ( ( 1.0 / pow((double) ( fch1 + ( foff * c ) ), m_power) ) - ( 1.0 / pow((double) fch1, m_power) ) ) );
      }
    }
    
    for(int i = 0; i < (int)range; i++)    {
      float n;
      modff(      ( ( (int) ( (  plan.user_dm(i).high - plan.user_dm(i).low ) / plan.user_dm(i).step ) + SDIVINDM ) / SDIVINDM )     , &n); // This calculates number of SDIVINDM blocks per DM range
      m_ndms[i] = (int) ( (int) n * SDIVINDM ); // This is number of DM trial per DM range. SDIVINDM is set to 21 in aa_params.hpp : RD
      if (m_max_ndms < m_ndms[i])
	m_max_ndms = m_ndms[i]; // looking for maximum number of DM trials for memory allocation
      m_total_ndms = m_total_ndms + m_ndms[i];
    }
    LOG(log_level::dev_debug, "mmMaximum number of dm trials in any of the range steps: " + std::to_string(m_max_ndms));
    
    str_dm[0].low = plan.user_dm(0).low;                        //
    str_dm[0].high = str_dm[0].low + ( m_ndms[0] * ( plan.user_dm(0).step ) );   // Redefines DM plan to suit GPU. //RD:Why needed?
    str_dm[0].step = plan.user_dm(0).step;                      //
    for (size_t i = 1; i < range; i++)    {
      str_dm[i].low = str_dm[i-1].high;
      str_dm[i].high = str_dm[i].low + m_ndms[i] * plan.user_dm(i).step;
      str_dm[i].step = plan.user_dm(i).step;
        
      if (plan.user_dm(i-1).inBin > 1) {
	m_maxshift = (int) ceil(( ( str_dm[i-1].low + str_dm[i-1].step * m_ndms[i - 1] ) * m_dmshifts[nchans - 1] ) / ( tsamp ));
	m_maxshift = (int) ceil((float) ( m_maxshift + ( (float) ( SDIVINT*2*SNUMREG ) ) ) / (float) plan.user_dm(i-1).inBin) / (float) ( SDIVINT*2*SNUMREG );
	m_maxshift = ( m_maxshift ) * ( SDIVINT*2*SNUMREG ) * plan.user_dm(i-1).inBin;
	if (( m_maxshift ) > m_maxshift_high)
	  m_maxshift_high = ( m_maxshift );
      }
    }//end of for loop.
    
    if (plan.user_dm(range-1).inBin > 1) { //it is calculated here
      m_maxshift = (int) ceil(( ( str_dm[range-1].low + str_dm[range-1].step * m_ndms[range - 1] ) * m_dmshifts[nchans - 1] ) / ( tsamp ));
      m_maxshift = (int) ceil((float) ( m_maxshift + ( (float) ( SDIVINT*2*SNUMREG ) ) ) / (float) plan.user_dm(range-1).inBin) / (float) ( SDIVINT*2*SNUMREG );
      m_maxshift = m_maxshift * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
      if (( m_maxshift ) > m_maxshift_high)
	m_maxshift_high = ( m_maxshift );
    }
    
    if (m_maxshift_high == 0)    {
      m_maxshift_high = (int) ceil(( ( str_dm[range-1].low + str_dm[range-1].step * ( m_ndms[range - 1] ) ) * m_dmshifts[nchans - 1] ) / tsamp);
    }
    m_max_dm = ceil(str_dm[range-1].high);
    
    m_maxshift = ( m_maxshift_high + ( SNUMREG * 2 * SDIVINT ) );
    LOG(log_level::dev_debug, "Range:\t"+std::to_string(range-1) + " MAXSHIFT:\t" + std::to_string(m_maxshift) + " Scrunch value:\t" + std::to_string(plan.user_dm(range-1).inBin));
    LOG(log_level::dev_debug, "Maximum dispersive delay:\t "+std::to_string(m_maxshift * tsamp)+ "(s).");

    if (m_maxshift >= nsamp)    {
      LOG(log_level::error, "Your maximum DM trial exceeds the number of samples you have. Reduce your maximum DM trial.");
      return false;
    }
    
    LOG(log_level::dev_debug, "Diagonal DM:\t" + std::to_string(( tsamp * nchans * 0.0001205 * powf(( fch1 + ( foff * ( nchans / 2 ) ) ), 3.0) ) / ( -foff * nchans )));
    if (m_maxshift >= nsamp)    {
      LOG(log_level::error, "Your maximum DM trial exceeds the number of samples you have. Reduce your maximum DM trial.");
      return false;
    }
    
    /* Four cases:
     * 1) nchans < m_max_ndms & nsamp fits in GPU RAM
     * 2) nchans > m_max_ndms & nsamp fits in GPU RAM
     * 3) nchans < m_max_ndms & nsamp does not fit in GPU RAM
     * 4) nchans > m_max_ndms & nsamp does not fit in GPU RAM
     */
    
    unsigned int max_tsamps;
    // Allocate memory to store the t_processed ranges:
    m_t_processed.resize(range);
    
    if (nchans < ( m_max_ndms )) {
      // This means that we can cornerturn into the allocated output buffer
      // without increasing the memory needed
      
      // Maximum number of samples we can fit in our GPU RAM is then given by:
      size_t SPDT_memory_requirements = (enable_analysis ? (sizeof(float)*(m_max_ndms)*SPDT_fraction) : 0 );
      max_tsamps = (unsigned int) ( (gpu_memory) / ( sizeof(unsigned short)*nchans + sizeof(float)*(m_max_ndms) + SPDT_memory_requirements )); // maximum number of timesamples we can fit into GPU memory
	size_t output_size_per_tchunk = ((size_t)max_tsamps)*sizeof(float)*((size_t)m_max_ndms);
	if (output_size_per_tchunk > MAX_OUTPUT_SIZE) max_tsamps = (unsigned int) ( (MAX_OUTPUT_SIZE) / (sizeof(float)*(m_max_ndms)));
      
      // Check that we dont have an out of range maxshift:
      if ((unsigned int)( m_maxshift ) > max_tsamps)    {
	LOG(log_level::error, "1. The selected GPU does not have enough memory for this number of dispersion trials. Reduce maximum dm or increase the size of dm step.");
	return false;
      }
      
      // Next check to see if nsamp fits in GPU RAM:
      if ((unsigned int)nsamp < max_tsamps)    {
	printf("We are in case 1");
	// Allocate memory to hold the values of nsamps to be processed
	unsigned long int local_t_processed = (unsigned long int) floor(( (double) ( nsamp - (m_maxshift) ) / (double) plan.user_dm(range-1).inBin ) / (double) ( SDIVINT*2*SNUMREG )); //number of timesamples per block
	local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
	for (size_t i = 0; i < range; i++)    {
	  m_t_processed[i].resize(1);
	  m_t_processed[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) plan.user_dm(i).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
	  m_t_processed[i][0] = m_t_processed[i][0] * ( SDIVINT*2*SNUMREG );
	}
	( m_num_tchunks ) = 1;
	LOG(log_level::dev_debug, "In 1");
      }
      else {
	printf("We are in case 3");
	// Work out how many time samples we can fit into ram
	int samp_block_size = max_tsamps - ( m_maxshift );
        
	// Work out how many blocks of time samples we need to complete the processing
	// upto nsamp-maxshift
	//int num_blocks = (int) floor(( (float) nsamp - ( *maxshift ) )) / ( (float) ( samp_block_size ) ) + 1;
        
	// Find the common integer amount of samples between all bins
	int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) plan.user_dm(range-1).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
	local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
        
	int num_blocks = (int) floor(( (float) nsamp - (float)( m_maxshift ) )) / ( (float) ( local_t_processed ) );
        
	// Work out the remaining fraction to be processed
	int remainder =  nsamp -  (num_blocks*local_t_processed ) - (m_maxshift) ;
	remainder = (int) floor((float) remainder / (float) plan.user_dm(range-1).inBin) / (float) ( SDIVINT*2*SNUMREG );
	remainder = remainder * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
        
	for (size_t i = 0; i < range; i++)    {
	  // Allocate memory to hold the values of nsamps to be processed
	  m_t_processed[i].resize(num_blocks + 1);
	  // Remember the last block holds less!
	  for (int j = 0; j < num_blocks; j++) {
	    m_t_processed[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) plan.user_dm(i).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
	    m_t_processed[i][j] = m_t_processed[i][j] * ( SDIVINT*2*SNUMREG );
	  }
	  // fractional bit
	  m_t_processed[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) plan.user_dm(i).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
	  m_t_processed[i][num_blocks] = m_t_processed[i][num_blocks] * ( SDIVINT*2*SNUMREG );
	}
	( m_num_tchunks ) = num_blocks + 1;
	LOG(log_level::dev_debug, "In 3");
	LOG(log_level::dev_debug, "num_blocks:\t" + std::to_string(num_blocks));
      }
    }
    else {
      // This means that we cannot cornerturn into the allocated output buffer . Case 2
      // without increasing the memory needed. Set the output buffer to be as large as the input buffer:
      
      // Maximum number of samples we can fit in our GPU RAM is then given by:
      size_t SPDT_memory_requirements = (enable_analysis ? (sizeof(float)*(m_max_ndms)*SPDT_fraction) : 0 );
      max_tsamps = (unsigned int) ( ( gpu_memory ) / ( nchans * ( sizeof(float) + sizeof(unsigned short) )+ SPDT_memory_requirements ));
	size_t output_size_per_tchunk = ((size_t)max_tsamps)*sizeof(float)*((size_t)nchans);
	if (output_size_per_tchunk > MAX_OUTPUT_SIZE) max_tsamps = (unsigned int)((MAX_OUTPUT_SIZE)/(sizeof(float)*nchans));

        printf("\n 2. aa_ddtr_Strategy.cpp:: max_tsamps:: %u m_maxshift:: %d m_max_ndms:: %d\n", max_tsamps, m_maxshift, m_max_ndms);
      
      // Check that we dont have an out of range maxshift:
      if (( m_maxshift ) > (int)max_tsamps) {
	LOG(log_level::error, "2. The selected GPU does not have enough memory for this number of dispersion trials. Reduce maximum dm or increase the size of dm step.");
	return false;
      }
      
      // Next check to see if nsamp fits in GPU RAM:
      if ((unsigned int)nsamp < max_tsamps) {
        printf("We are in case 2");
	// We have case 2)
	// Allocate memory to hold the values of nsamps to be processed
	int local_t_processed = (int) floor(( (float) ( nsamp - ( m_maxshift ) ) / (float) plan.user_dm(range-1).inBin ) / (float) ( SDIVINT*2*SNUMREG )); //RD : Why is it nsamp-m_maxshift. This is overlap.
	local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
	for (size_t i = 0; i < range; i++) {
	  m_t_processed[i].resize(1);
	  m_t_processed[i][0] = (int) floor(( (float) ( local_t_processed ) / (float) plan.user_dm(i).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
	  m_t_processed[i][0] = m_t_processed[i][0] * ( SDIVINT*2*SNUMREG );
	}
	( m_num_tchunks ) = 1;
	LOG(log_level::dev_debug, "In 2");
      }
      else {
	printf("We are in case 4");
	// Work out how many time samples we can fit into ram
	int samp_block_size = max_tsamps - ( m_maxshift );
        printf("We are in case 4: samp_block_size %d", samp_block_size);
        
	// Work out how many blocks of time samples we need to complete the processing
	// upto nsamp-maxshift
	//int num_blocks = (int) floor(( (float) nsamp - (float) ( *maxshift ) ) / ( (float) samp_block_size ));
        
	// Find the common integer amount of samples between all bins
	int local_t_processed = (int) floor(( (float) ( samp_block_size ) / (float) plan.user_dm(range-1).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
	local_t_processed = local_t_processed * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
        
	// samp_block_size was not used to calculate remainder instead there is local_t_processed which might be different
	int num_blocks = (int) floor(( (float) nsamp - (float) ( m_maxshift ) ) / ( (float) local_t_processed ));
        printf("No of blocks: %d", num_blocks);
        
	// Work out the remaining fraction to be processed
	int remainder = nsamp - ( num_blocks * local_t_processed ) - ( m_maxshift );
	remainder = (int) floor((float) remainder / (float) plan.user_dm(range-1).inBin) / (float) ( SDIVINT*2*SNUMREG );
	remainder = remainder * ( SDIVINT*2*SNUMREG ) * plan.user_dm(range-1).inBin;
        
	for (size_t i = 0; i < range; i++)    {
	  // Allocate memory to hold the values of nsamps to be processed
	  m_t_processed[i].resize(num_blocks + 1);
	  // Remember the last block holds less!
	  for (int j = 0; j < num_blocks; j++) {
	    m_t_processed[i][j] = (int) floor(( (float) ( local_t_processed ) / (float) plan.user_dm(i).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
	    m_t_processed[i][j] = m_t_processed[i][j] * ( SDIVINT*2*SNUMREG );
	  }
	  // fractional bit
	  m_t_processed[i][num_blocks] = (int) floor(( (float) ( remainder ) / (float) plan.user_dm(i).inBin ) / (float) ( SDIVINT*2*SNUMREG ));
	  m_t_processed[i][num_blocks] = m_t_processed[i][num_blocks] * ( SDIVINT*2*SNUMREG );
	}
	( m_num_tchunks ) = num_blocks + 1;
	LOG(log_level::dev_debug, "In 4");
      }
    }
    printf("\nNOTICE: Maxshift memory needed:\t%lu MB.", nchans * ( m_maxshift ) * sizeof(unsigned short) / 1024 / 1024);
    if (nchans < ( m_max_ndms ))    {
      printf("\nNOTICE: Output memory needed:\t%lu MB.", ( m_max_ndms ) * ( m_maxshift ) * sizeof(float) / 1024 / 1024);
    }
    else {
      printf("\nNOTICE: Output memory needed:\t%lu MB.\n", nchans * ( m_maxshift ) * sizeof(float) / 1024 / 1024);
    }

    // The memory that will be allocated on the GPU in function allocate_memory_gpu is given by gpu_inputsize + gpu_outputsize
    // gpu_inputsize
    aa_device_info& device_info = aa_device_info::instance();
    int time_samps = m_t_processed[0][0] + m_maxshift;
    device_info.request((size_t) time_samps * (size_t)nchans * sizeof(unsigned short));
    // gpu_outputsize depends on nchans
    if (nchans < m_max_ndms) {
      printf("\n nchans :# %d m_max_ndms :# %d \n", nchans, m_max_ndms);
      if(!device_info.request((size_t)time_samps * (size_t)m_max_ndms * sizeof(float))) {
	std::cout << "1. ERROR:  Could not request memory." << std::endl;
	return false;
      }
    }
    else {
      if(!device_info.request((size_t)time_samps * (size_t)nchans * sizeof(float))) {
	std::cout << "2. ERROR:  Could not request memory." << std::endl;
	return false;
      }
    }
    
    //Strategy does not change inBin, outBin.
    //Re-assign original inBin, outBin to the strategy.
    for(size_t i = 0; i < str_dm.size(); i++) {
      str_dm.at(i).inBin = plan.user_dm(i).inBin;
      str_dm.at(i).outBin = plan.user_dm(i).outBin;
    }

// additions by abhinav ---------------------writing strategy information to file so as to read it when metadata.strat() == 1 -----------------------------------------------
/*    if (m_metadata.strat() == 0) {

       FILE * pFile;
       pFile = fopen ("myfile.txt","w");

       if (pFile!=NULL)
       {
           fprintf (pFile, "m_maxshift     : %d \n", m_maxshift);
           fprintf (pFile, "m_max_ndms     : %d \n", m_max_ndms);
           fprintf (pFile, "range          : %d \n", (int) range);
           fprintf (pFile, "scrunch value  : %d \n", plan.user_dm(range-1).inBin);
           fprintf (pFile, "Max disp delay : %f \n", m_maxshift * tsamp);
           fprintf (pFile, "Diagonal DM    : %f \n", (tsamp * nchans * 0.0001205 * powf(( fch1 + ( foff * ( nchans / 2 ) ) ), 3.0) ) / ( -foff * nchans ));
           fprintf (pFile, "m_num_tchunks  : %d \n", m_num_tchunks);
//           fprintf (pFile, "remainder      : %d \n", remainder);
	   fprintf (pFile, "dm limits ndm  : ");
           for(size_t i = 0; i < range; i++) {
           const aa_ddtr_plan::dm tmp = dm(i);
	   fprintf (pFile, "%f %f %f %d %d %d ",tmp.low,tmp.high,tmp.step, tmp.inBin,tmp.outBin,ndms(i));
           }
           fprintf(pFile, "\n");
           fprintf (pFile, "t_procsd i j   : %d %d \n", t_processed().size(), t_processed().at(0).size());
           fprintf (pFile, "t_procsd arr   : ");
           for(size_t i = 0; i < t_processed().size(); i++) {
	     for(size_t j = 0; j < t_processed().at(i).size(); j++) {
               fprintf (pFile, "%d ", t_processed()[i][j]);
             }
           }
           fprintf(pFile, "\n");

           fprintf (pFile, "m_total_ndms   : %d \n", m_total_ndms);
           fprintf (pFile, "m_max_ndms     : %d \n", m_max_ndms);
           fprintf (pFile, "m_max_dm       : %f \n", m_max_dm);
           fprintf (pFile, "m_maxshift_high: %f \n", m_maxshift_high);


           fclose (pFile);
       }
 
    }
*/
// edits end ----------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    m_ready = true;
    m_strategy_already_calculated = true;
    return true;
  }

  /**
   * Returns the ready state of the instance of the class.
   */
  bool aa_ddtr_strategy::setup() {
    is_setup = true;
    if(is_setup && ready()) {
      return true;
    }
    else {
      return false;
    }
  }

} //namespace astroaccelerate
