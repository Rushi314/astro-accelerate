#include "aa_sigproc_input.hpp"
#include <cmath>

#include "aa_ddtr_plan.hpp" // added by abhinav

#include <string>
#include <wordexp.h>


namespace astroaccelerate {

  /** \brief Constructor for aa_sigproc_input. */
  aa_sigproc_input::aa_sigproc_input(const std::string &path) : header_is_read(false), data_is_read(false), file_path(path) {
    isopen = false;
  printf("\n constructor \n");
  }

  /** \brief Destructor for aa_sigproc_input. */
  aa_sigproc_input::~aa_sigproc_input() {
    if(isopen) {
      close();
    }
  printf("\n destructor \n");
  }

  /**
   * \brief Method to open the sigproc input file.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   */
  bool aa_sigproc_input::open() {
    fp = fopen(file_path.c_str(), "rb");
    if(fp == NULL) {
      return false;
    }
    printf("\n file opened: %s \n",file_path.c_str());
    isopen = true;
    return true;
  }

  /**
   * \brief Closes the input file.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   */
  bool aa_sigproc_input::close() {
    if(isopen) {
      if(fclose(fp) == -1) {
	return false;
      }
      else {
	isopen = false;
      }
    }
    return true;
  printf("\n file closed \n");
  }

  /**
   * \brief Method to read the metadata from the sigproc input file.
   * \returns An aa_filterbank_metadata object containing the metadata read from the sigproc input file. If the data could not be read, a trivial instance is returned.
   */
  aa_filterbank_metadata aa_sigproc_input::read_metadata() {
    if(!isopen) {
      if(!open()) {
	aa_filterbank_metadata empty;
        printf("\n returning empty metadata \n");
	return empty;
      }
    }

    printf("\n will read metadata \n");

    aa_filterbank_metadata metadata;
    get_file_data(metadata);
    header_is_read = true;
    return metadata;
  }


  aa_filterbank_metadata aa_sigproc_input::read_metadata(float max_dm_val) {
    if(!isopen) {
      if(!open()) {
	aa_filterbank_metadata empty;
        printf("\n returning empty metadata \n");
	return empty;
      }
    }

    printf("\n will read metadata \n");

    aa_filterbank_metadata metadata;
    get_file_data(metadata, max_dm_val);
    header_is_read = true;
    return metadata;
  }


  aa_filterbank_metadata aa_sigproc_input::read_metadata(float max_dm_val, std::string config_path) {
    if(!isopen) {
      if(!open()) {
	aa_filterbank_metadata empty;
        printf("\n returning empty metadata \n");
	return empty;
      }
    }

    printf("\n will read metadata \n");

    aa_filterbank_metadata metadata;
    get_file_data(metadata, max_dm_val, config_path);
    header_is_read = true;
    return metadata;
  }

  aa_filterbank_metadata aa_sigproc_input::read_metadata_new(aa_filterbank_metadata metadata) {
/*    if(!isopen) {
      if(!open()) {
	aa_filterbank_metadata empty;
        printf("\n returning empty metadata \n");
	return empty;
      }
    }
*/
    printf("\n will read new metadata \n");
    printf("\n filterbank_data.nchans(): >>> %d \n", metadata.nchans());

//    aa_filterbank_metadata metadata;
    get_file_data_new(metadata);
    header_is_read = true;
    return metadata;
  }

// following function updates source raj and dec in the metadata:
  aa_filterbank_metadata aa_sigproc_input::read_metadata_source(aa_filterbank_metadata metadata) {
    if(!isopen) {
      if(!open()) {
	aa_filterbank_metadata empty;
        printf("\n returning empty metadata from read_metadata_source \n");
	return empty;
      }
    }

    char line[128], tempo_val[64];
    float src_raj, src_dej; 

    printf("\n will update metadata with new source location \n");
//    printf("\n filterbank_data.nchans(): >>> %d \n", metadata.nchans());

//    aa_filterbank_metadata metadata;

//    printf("\n  source: metadata.nchans:: is :: %d \n", metadata.nchans());

//    printf("\n  m_meta.nchans:: is :: %d \n", m_meta.nchans());
//    printf("\n  metadata.nchans:: is :: %d \n", m_meta_new.nchans());

      while (fgets(line, sizeof(line), fp)) {
        printf("\n ** ** ** %.90s \n", line);

        strcpy(tempo_val, &line[9]);
        printf("\n temo_val is :: %s \n", tempo_val);	

	if (strncmp(line, "src_raj", 6) == 0 )              // Time stamp of first sample (MJD) 
		{
		src_raj = strtof(tempo_val, NULL);
		printf("\n raj: %f \n", src_raj);
		}
	else if (strncmp(line, "src_dej", 6) == 0 )              // Time stamp of first sample (MJD) 
		{
		src_dej = strtof(tempo_val, NULL);
		printf("\n dej: %f \n", src_dej);
		}
	}

//    nsamples = (int)((nsamples) ? nsamples : (nsamp) ? nsamp : 0);
    aa_filterbank_metadata meta(metadata.telescope_id(),
                                metadata.machine_id(),
                                metadata.data_type(),
                                "unknown",
                                "unknown",
                                metadata.barycentric(),
                                metadata.pulsarcentric(),
                                metadata.az_start(),
                                metadata.za_start(),
                                src_raj,
                                src_dej,
                                metadata.tstart(),
                                metadata.tsamp(),
                                metadata.nbits(),
                                metadata.nsamples(),
                                metadata.fch1(),
                                metadata.foff(),
                                0,
                                metadata.fchannel(),
                                0,
                                metadata.nchans(),
                                metadata.nifs(),
                                metadata.refdm(),
                                metadata.period(),
                                1,
                                metadata.ovrlp());

//    printf("\n  filterbank_data.nchans:: is :: %d \n", filterbank_data.nchans());
//    printf("\n  meta.nchans:: is :: %d \n", meta.nchans());

    metadata = meta;
    m_meta = meta;
    m_meta_new   = meta;
//    printf("\n  aa_sigproc_input :: metadata.nchans:: is :: %d \n", metadata.nchans());
//    printf("\n  aa_sigproc_input :: metadata.strat:: is :: %d \n", metadata.strat());
//    printf("\n  aa_sigproc_input :: metadata.nsamples:: is :: %d \n", metadata.nsamples());


    header_is_read = true;
    return metadata;
  }



  /** 
  * \brief Method to read in data from input file
  * \returns An ________ object containing the constants read from input data.
  */
  bool aa_sigproc_input::read_metadata_input(float *sigma_cutoff_input, float *sigma_constant_input, float *max_boxcar_width_in_sec_input, float *periodicity_sigma_cutoff_input, float *periodicity_harmonics_input, aa_ddtr_plan *ddtr_plan, bool *baselinenoise_val, bool *analysis_val, bool *debug_val, char (*file_name_input)[128], char (*config_path)[128]) {
    if(!isopen) {
      if(!open()) {
        printf("\n returning empty handed \n");
	return false;
      }
    }

    printf("\n will read input metadata \n");

	char line[128]; 
        aa_ddtr_plan d_p_input;
        float high_input, low_input, step_input;
	int inBin_input, outBin_input;
//	float *sig_cut_off, sigma_constant_input, max_boxcar_width_in_sec_input, periodicity_sigma_cutoff_input;
//	float *periodicity_harmonics_input;

	while (fgets(line, sizeof(line), fp)) {
//        printf("\n **-*-** %.35s \n", line);
//		printf("\n tempo_val %s \n", tempo_val);
		if (strncmp(line, "range", 5) == 0 )               // Frequency of channel 1 (MHz)
		{
		low_input = strtof(&line[9], NULL);
		high_input = strtof(&line[14], NULL);
		step_input = strtof(&line[19], NULL);
		inBin_input = strtol(&line[24], NULL, 10);
		outBin_input = strtol(&line[26], NULL, 10);
		printf("\n low : %f high : %f step : %f inBin : %d outBin : %d \n", low_input, high_input, step_input, inBin_input, outBin_input);
		ddtr_plan->add_dm(low_input, high_input, step_input, inBin_input, outBin_input);  // this use of -> is quite a find !! abhinav
		}
		else if (strncmp(line, "sigma_cutoff", 12) == 0 )               // Frequency of channel 1 (MHz)
		{
		*sigma_cutoff_input = strtof(&line[16], NULL);
		printf("\n sigma_cutoff :  ::: %f \n", *sigma_cutoff_input);
		}
		else if (strncmp(line, "sigma_constant", 14) == 0 )              // Time stamp of first sample (MJD) 
		{
		*sigma_constant_input = strtof(&line[16], NULL);
		printf("\n sigma_constant_input :  ::: %f \n", *sigma_constant_input);
		}
		else if (strncmp(line, "max_boxcar_width_in_sec", 23) == 0 )               // "Channel bandwidth"
		{
		*max_boxcar_width_in_sec_input = strtof(&line[24], NULL);
		printf("\n max_boxcar_width_in_sec : :::%f \n", *max_boxcar_width_in_sec_input);
		}
		else if (strncmp(line, "periodicity_sigma_cutoff", 24) == 0 )               // "Channel bandwidth"
		{
		*periodicity_sigma_cutoff_input = strtof(&line[25], NULL);
		printf("\n periodicity_sigma_cutoff_input: ::: %f \n", *periodicity_sigma_cutoff_input);
		}
		else if (strncmp(line, "periodicity_harmonics", 21) == 0 )               // "Channel bandwidth"
		{
		*periodicity_harmonics_input = strtof(&line[22], NULL);
		printf("\n periodicity_harmonics_input: ::: %f \n", *periodicity_harmonics_input);
		}
		else if (strncmp(line, "baselinenoise", 13) == 0 )               // "Channel bandwidth"
		{
		*baselinenoise_val = true;
		printf("\n baselinenoise_val: ::: %d \n", *baselinenoise_val);
		}
		else if (strncmp(line, "file", 4) == 0 )               // "Channel bandwidth"
		{


                strcpy(*file_name_input, &line[5]);
//                strcpy(tmp_line, &line[5]);

/*
	        FILE *fp_test = NULL;      // Path to the input data file

		char tmp_line[80];
		std::string tmp_fname, fpath_test;

                strcpy(tmp_line, &line[5]);
		*file_name_input_str = tmp_line;
		tmp_fname = tmp_line;

		fscanf(fp_test, "%s", tmp_fname);

//		*file_name_input = true;

	    // This command expands "~" to "home/username/"
	    wordexp_t expanded_string;
	    wordexp(tmp_line, &expanded_string, 0);
*/
//	    if (( fp_test = fopen(expanded_string.we_wordv[0], "rb") ) == NULL)
/*	    if (( fp_test = fopen("/data/jroy/data/FRB_DM1000_163.84us_4K_7mP_4msW_3sig.header.gpt", "rb") ) == NULL)
	      {
		printf("\n Invalid data file!\n");
		return false;
	      }
	    else {
	      //File is valid
	      //Close it again for later re-opening
		printf("\n file opened !!!!!!!!\n");
	      fpath_test = expanded_string.we_wordv[0];
	      fclose(fp_test);
	    }
	    wordfree(&expanded_string);
*/


		printf("\n line is: ::: %s :: file_name_input is :: %s :: \n", line, *file_name_input);
		}
		else if (strncmp(line, "config_file_path", 16) == 0 )               // "Channel bandwidth"
		{
                strcpy(*config_path, &line[17]);
//		*file_name_input = true;
		printf("\n config path line is: ::: %s \n", line);
		}

//	printf("\n %s ++--++-- \n", line); 
	}
    return true;
  }


 

  /**
   * \brief If the file is open, and the header has been read, and the data has not yet been read, then read the input data from the input file.
   * \details Reading the telescope input data can only be performed once, after which this method will always return false.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   * \warning The method will return true only once, that is the first time the data are read from the input successfully. At this point the input_buffer should be checked for data. 
   */
  bool aa_sigproc_input::read_signal() {
/*    printf("\n !isopen: %d !header_is_read: %d data_is_read: %d \n", !isopen, !header_is_read, data_is_read);
    if(!isopen || !header_is_read || data_is_read) {
      return false;
    }
*/
    printf("\n get recorded data to be called \n");  
    get_recorded_data(m_input_buffer);
    return true;
  }

  bool aa_sigproc_input::read_new_signal(int buf_count, const aa_filterbank_metadata &filterbank_data, unsigned long int *f_pos) {
/*    printf("\n !isopen: %d !header_is_read: %d data_is_read: %d \n", !isopen, !header_is_read, data_is_read);
    if(!isopen || !header_is_read || data_is_read) {
      return false;
    }
*/
    printf("\n get recorded data to be called \n");
//    aa_filterbank_metadata metadata;  
    if(get_new_recorded_data(m_input_buffer, buf_count, filterbank_data, &f_pos))
    {
      return true;
    }
    else
    {
      return false;
    }
  }


  /**
   * \brief Method to read the input data from the sigproc input file.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   */
  bool aa_sigproc_input::get_file_data(aa_filterbank_metadata &metadata) {
    double az_start = 0;
    double za_start = 0;
    double src_raj = 0;
    double src_dej = 0;
    double tstart = 0;
    double tsamp = 0;
    double refdm = 0;
    double period = 0;
    double fch1 = 0;
    double foff = 0;
    double fchannel = 0;
    
    int telescope_id = 0;
    int machine_id = 0;
    int data_type = 0;
    int barycentric = 0;
    int pulsarcentric = 0;
    int nbits = 0;
    int nsamples = 0;
    int nchans = 0;
    int nifs = 0;
    char *string = (char *) malloc(80 * sizeof(char));
// added by abhinav --------------------------------------------------
    char line[256], tempo_name[10], tempo_val[50], linetel[500];
    int count_num = 0;
    unsigned long int nsamp = 0;
    int dej1 = 0;
    int dej2 = 0;
    float dej3 = 0.0f;
    float dej, dej1f, dej2f;
    int raj1 = 0;
    int raj2 = 0;
    float raj3 = 0;
    float raj;
    int strat = 0;

    float sign_bw, bw_val, ch_band_width, time_res_tel, low_freq, high_freq;
    int int_val, ovrlp;


// edits end -------------------------------------------------------------    
    int nchar = 0;
    int nbytes = sizeof(int);
    
/*    while (1) {
      strcpy(string, "ERROR");
      if (fread(&nchar, sizeof(int), 1, fp) != 1) {
	fprintf(stderr, "\nError while reading file: block not equal to nchar\n");
	return false;
      }
      if (feof(fp)) {
	return false;
      }
        
      if (nchar > 1 && nchar < 80) {
	if (fread(string, nchar, 1, fp) != 1) {
	  fprintf(stderr, "\nError while reading file: size of nchar not equal to string\n");
	  return false;
	}
            
	string[nchar] = '\0';
	// For debugging only
	//printf("From .fil header:\t%d\t%s\n", nchar, string);
	nbytes += nchar;
            
	if (strcmp(string, "HEADER_END") == 0) {
	  break;
	}
            
	if(strcmp(string, "telescope_id") == 0) {
	  if (fread(&telescope_id, sizeof(telescope_id), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "machine_id") == 0) {
	  if (fread(&machine_id, sizeof(machine_id), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "data_type") == 0) {
	  if (fread(&data_type, sizeof(data_type), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "barycentric") == 0) {
	  printf("barycentric\n");
	  if (fread(&barycentric, sizeof(barycentric), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "pulsarcentric") == 0) {
	  if (fread(&pulsarcentric, sizeof(pulsarcentric), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "az_start") == 0) {
	  if (fread(&az_start, sizeof(az_start), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "za_start") == 0) {
	  if (fread(&za_start, sizeof(za_start), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "src_raj") == 0) {
	  if (fread(&src_raj, sizeof(src_raj), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if(strcmp(string, "src_dej") == 0) {
	  if (fread(&src_dej, sizeof(src_dej), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "tstart") == 0) {
	  if (fread(&tstart, sizeof(tstart), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "tsamp") == 0) {
	  printf("\n test test test\n");
	  if (fread(&tsamp, sizeof(tsamp), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "nbits") == 0) {
	  if (fread(&nbits, sizeof(nbits), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "nsamples") == 0) {
	  if (fread(&nsamples, sizeof(nsamples), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "fch1") == 0) {
	  if (fread(&fch1, sizeof(fch1), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "foff") == 0) {
	  if (fread(&foff, sizeof(foff), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "fchannel") == 0) {
	  if (fread(&fchannel, sizeof(fchannel), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "nchans") == 0) {
	  if (fread(&nchans, sizeof(nchans), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "nifs") == 0) {
	  if (fread(&nifs, sizeof(nifs), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "refdm") == 0) {
	  if (fread(&refdm, sizeof(refdm), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
	else if (strcmp(string, "period") == 0) {
	  if (fread(&period, sizeof(period), 1, fp) != 1) {
	    fprintf(stderr, "\nError while reading file\n");
	    return false;
	  }
	}
      }
    }

    printf("\n refdm: %f :: period: %f :: fchannel: %f :: az_start: %f :: za_start : %f :: barycentric: %d :: pulsarcentric: %d :: nchans: %d \n", refdm, period, fchannel, az_start, za_start, barycentric, pulsarcentric, nchans);
*/


/*      while (fgets(line, sizeof(line), fp)) {
//        printf("\n **-*-** %.35s \n", line);

        strcpy(tempo_val, &line[35]);	
/*        if (strncmp(line, "Time stamp of first sample (MJD) ", 33) == 0 )              // Time stamp of first sample (MJD) 
		{
		tstart = strtof(tempo_val, NULL);
//		printf("\n Time stamp of first sample: %f \n", tempo_val, tstart);
		}
*/
/*        if (strncmp(line, "Frequency of channel 1 (MHz)     ", 33) == 0 )               // Frequency of channel 1 (MHz)
		{
		fch1 = strtof(tempo_val, NULL);
//		printf("\n Frequency of channel 1 (MHz) >>>>>>>> %s >>> %f \n", tempo_val, fch1);
		}
	else if (strncmp(line, "Channel bandwidth      (MHz)     ", 33) == 0 )               // "Channel bandwidth"
		{
		foff = strtof(tempo_val, NULL);
//		printf("\n Channel bandwidth >>>>>>>> %s >>> %f \n", tempo_val, foff);
		}
	else if (strncmp(line, "Sample time (us)                 ", 33) == 0 )              // sample time
		{
		tsamp = strtof(tempo_val, NULL);
		tsamp = tsamp/1000000.0;
		printf("\n Sample time >>>>>>>> %s >>> %12.9f \n", tempo_val, tsamp);
		}

	else if (strncmp(line, "Number of channels               ", 33) == 0 )              // Time stamp of first sample (MJD) 
		{
		nchans = strtol(tempo_val, NULL, 10);
//		printf("\n Number of channelss (nchans): %d \n", nchans);
		}
/*	else if (strncmp(line, "Number of bits per sample        ", 33) == 0 )              // Time stamp of first sample (MJD) 
		{
		nbits = strtol(tempo_val, NULL, 10);
//		printf("\n Number of bits per sample (nbits): %d \n", nbits);
		}
	else if (strncmp(line, "Number of IFs                    ", 33) == 0 )              // Time stamp of first sample (MJD) 
		{
		nifs = strtol(tempo_val, NULL, 10);
//		printf("\n nifs: %d \n", nifs);
		}
	else if (strncmp(line, "Number of samples                ", 33) == 0 )              // Time stamp of first sample (MJD) 
		{
		nsamp = strtol(tempo_val, NULL, 10);
//		printf("\n <<<<<<<<<<<<nsamples>>>>>>>> %d \n", nsamp);
		}
	else if (strncmp(line, "Source RA (J2000)                ", 33) == 0 )              // Time stamp of first sample (MJD) 
		{
		raj1 = strtol(tempo_val, NULL, 10);
	 	strcpy(tempo_val, &tempo_val[3]);
                raj2 = strtol(tempo_val, NULL, 10);
	 	strcpy(tempo_val, &tempo_val[3]);
                raj3 = strtof(tempo_val, NULL);
		src_raj = raj1*10000+raj2*100+raj3;
		printf("\n raj1: %d :: raj2: %d :: raj3: %f :: raj: %f :: tempo_val: %s \n", raj1, raj2, raj3, src_raj, tempo_val);
		}
	else if (strncmp(line, "Source DEC (J2000)               ", 33) == 0 )              // Time stamp of first sample (MJD) 
		{
		dej1 = strtol(tempo_val, NULL, 10);
//		printf("\n dej1: %d :: tempo_val: %s \n", dej1, tempo_val);
	 	strcpy(tempo_val, &tempo_val[4]);
                dej2 = strtol(tempo_val, NULL, 10);
//		printf("\n dej2: %d :: tempo_val: %s \n", dej2, tempo_val);
	 	strcpy(tempo_val, &tempo_val[3]);
                dej3 = strtof(tempo_val, NULL);
//		printf("\n dej3: %f :: tempo_val: %s \n", dej3, tempo_val);
		dej1f = (float) dej1*10000.0;
//		printf("\n dej1f: %f \n", dej1f);
		dej2f = (float) dej2*100.0;
//		printf("\n dej2f: %f \n", dej2f);
		if (dej1 < 0) src_dej = dej1*10000.0 - dej2*100.0 - dej3;
		else src_dej = dej1*10000.0 + dej2*100.0 + dej3;
		printf("\n dej1: %d :: dej2: %d :: dej3: %f :: dej: %f :: tempo_val: %s \n", dej1, dej2, dej3, src_dej, tempo_val);
//		printf("\n <<<<<<<<<<<<dej1>>>>>>>> %d :: tempo_val: %s \n", dej1, tempo_val);
		}
*/
//    count_num = count_num+1;
//    printf("\n count is: %d \n", count_num);
//    if(count_num == 10) break;
//    } // while loop ends here ------------------------------------------------------------

//    open("abhinav/FRB_pipeline-1.7.10/input_files/time.hdr");          

    FILE * pFile;
    pFile = fopen("/home/guest/Rushikesh/FRB_pipeline-1.7.10.GMRT300-500/input_files/time.hdr", "r");
    if(pFile == NULL) {
      printf("\n ATTENTION:: time.hdr cannot be opened! \n");
    }
   else {
   printf("\n File opened :: pFILE \n");
   while (fgets(line, sizeof(line), pFile)) {
        printf("\n **-*-** %.90s \n", line);

        strcpy(tempo_val, &line[35]);
        printf("\n temo_val is :: %s \n", tempo_val);	
        if (strncmp(line, "Time stamp of first sample (MJD) ", 33) == 0 )              // Time stamp of first sample (MJD) 
		{
//		tstart = strtof(tempo_val, NULL);
		tstart = atof(tempo_val);
		printf("\n Time stamp of first sample:: %s :: tstart:: %12.9f \n", tempo_val, tstart);
		}
    }  // time.hdr - while loop ends here ------------------------------------------------
    }  // else condition ends ---------------------------------------------
//    fclose(pFile);

// Reading data from gpu.hdr file::

	while (fgets(line, sizeof(line), fp)) {
//        printf("\n **-*-** %.35s \n", line);
		strcpy(tempo_val, &line[14]);	
//		printf("\n line:: line:: %s \n", linetel);
//		printf("\n tempo_val %s \n", tempo_val);
		if (strncmp(line, "GPU_RF4", 7) == 0 )               // Frequency of channel 1 (MHz)
		{
		fch1 = strtof(tempo_val, NULL);
		printf("\n Frequency of channel 1 (MHz): %f \n", fch1);
		}
		else if (strncmp(line, "GPU_CHAN_MAX", 12) == 0 )              // Time stamp of first sample (MJD) 
		{
		nchans = strtol(tempo_val, NULL, 10);
		printf("\n Number of channels: %d \n", nchans);
		}
		else if (strncmp(line, "GPU_SB4", 7) == 0 )               // "Channel bandwidth"
		{
		sign_bw = (-1)*strtof(tempo_val, NULL);
		printf("\n Channel bandwidth sign: %f \n", sign_bw);
		}
		else if (strncmp(line, "GPU_ACQ_BW", 10) == 0 )               // "Channel bandwidth"
		{
		bw_val = strtof(tempo_val, NULL);
		printf("\n Channel bandwidth size: %f \n", bw_val);
		}
		else if (strncmp(line, "GPU_BM_INT", 10) == 0 )               // "Channel bandwidth"
		{
		int_val = strtol(tempo_val, NULL, 10);
		printf("\n Beam Integration, i.e. Number of FFT cycles: %d \n", int_val);
		}
//	printf("\n %s ++--++-- \n", line); 
	}
	time_res_tel = (int_val*(nchans)*2)/(2*bw_val);		// in microseconds
	time_res_tel = time_res_tel/1000000.0;				// in seconds 

	tsamp = time_res_tel;
	printf("\n ** Time resolution :: %12.9f \n",tsamp);
	ch_band_width = (sign_bw*bw_val)/(nchans);

	foff = ch_band_width;
	printf("\n ** Channel Bandwith (foff): %9.8f \n",foff);

	high_freq = fch1;
	fch1 = (fch1)+(0.5*(foff));
	printf("\n fch1 is : %9.8f \n", fch1);	
// -reading of gpu.hdr ends here -------------------------------------------------------------------

    free(string);

// hard coded values ::    
    nbits = 8;
    nifs = 1;
   
// following should be changed if input file is changed- will be given by the telescope later
//    src_raj = 214400.0;
//    src_dej = -393300.0;

    FILE * pFile_source;
    pFile_source = fopen("/home/guest/Rushikesh/FRB_pipeline-1.7.10.GMRT300-500/input_files/source.hdr", "r");
    if(pFile_source == NULL) {
      printf("\n ATTENTION:: source.hdr cannot be opened! \n");
    }
   else {
   printf("\n File opened :: source.hdr \n");
   while (fgets(line, sizeof(line), pFile_source)) {
        printf("\n ** ** ** %.90s \n", line);

        strcpy(tempo_val, &line[9]);
        printf("\n temo_val is :: %s \n", tempo_val);	

	if (strncmp(line, "src_raj", 6) == 0 )              // Time stamp of first sample (MJD) 
		{
		src_raj = strtof(tempo_val, NULL);
		printf("\n raj: %f :: tempo_val: %s \n", src_raj, tempo_val);
		}
	else if (strncmp(line, "src_dej", 6) == 0 )              // Time stamp of first sample (MJD) 
		{
		src_dej = strtof(tempo_val, NULL);
		printf("\n dej: %f :: tempo_val: %s \n", src_dej, tempo_val);
		}
	} // while loop ends here
	} // else condition ends here

   if(!fclose(pFile_source)) printf("\n source.hdr file successfully closed \n");

//--------------------------------------------------------------------edits end----------------------------------------
    
    // Check that we are working with one IF channel
    if (nifs != 1) {
      printf("ERROR: Astro-Accelerate can only work with one IF channel.\n");
      return false;
    }

/*
    FILE * pFilee;
    pFilee = fopen("/data/jroy/data/FRB_DM500_163.84us_4K_7mP_4msW_3sig.header.fil", "r");
    if(pFilee == NULL) {
      printf("\n ATTENTION:: pFile cannot be opened! \n");
    }
   else {
   printf("\n File opened :: pFILEE \n");
    
//    fpos_t file_loc;
//    fgetpos(fp, &file_loc); //Position of end of header
      fseek(pFilee, 399, SEEK_SET);
      printf("\n Setting pointer to 399.\n");
    
//    unsigned long int nsamp = 0;

    // Getting number of time-samples based on file size
    unsigned long int data_start = ftell(pFilee);
    if (fseek(pFilee, 0, SEEK_END) != 0) {
      printf("\nERROR!! Failed to seek to the end of data file\n");
      return false;
    }
    unsigned long int exp_total_data = ftell(pFilee);
    exp_total_data = exp_total_data - data_start;
    fseek(pFilee, data_start, SEEK_SET);

    if (( nbits ) == 32) {
      nsamp = exp_total_data/((nchans)*4);
    }
    else if (( nbits ) == 16) {
      nsamp = exp_total_data/((nchans)*2);
    }
    else if (( nbits ) == 8) {
      nsamp = exp_total_data/((nchans));
    }
    else if (( nbits ) == 4) {
      nsamp = 2*exp_total_data/((nchans));
    }
    else if (( nbits ) == 2) {
      nsamp = 4*exp_total_data/((nchans));
    }
    else if (( nbits ) == 1) {
      nsamp = 8*exp_total_data/((nchans));
    }
    else {
      printf("ERROR: Currently this code only runs with 1, 2, 4, 8, and 16 bit data\n");
      return false;
    }

    // Move the file pointer back to the end of the header
//    fsetpos(fp, &file_loc);
    fclose(pFilee);
    }
*/
//    nsamp = 1024*1024;

    float epsilon = 0.0000001;
// number of time samples to be read depends upon the value of time between each time sample is taken::
    if (fabs(tsamp - 0.00008192) < epsilon) nsamp = 2*1024*1024;	//sampling_time=81.92 us
    else if(fabs(tsamp - 0.00016384) < epsilon) nsamp = 1*1024*1024;    //sampling_time=163.84 us
    else if(fabs(tsamp - 0.00032768) < epsilon) nsamp = 512*1024;       //sampling_time=327.68 us
    else if(fabs(tsamp - 0.00065536) < epsilon) nsamp = 256*1024;       //sampling_time=655.36 us
    else if(fabs(tsamp - 0.00131072) < epsilon) nsamp = 128*1024;       //sampling_time=1310.72 us

    printf("\n nsamp is: %d \t tsamp: %12.9f \n", nsamp, tsamp);

// change the value hard coded below
    low_freq = high_freq + (sign_bw)*(bw_val);
    printf("\n low freq: %f high freq: %f \n", low_freq, high_freq);
    ovrlp = 4.15*1000.*2000*((1./(low_freq*low_freq))*(1./(high_freq*high_freq)));
    
    nsamples = (int)((nsamples) ? nsamples : (nsamp) ? nsamp : 0);
    aa_filterbank_metadata meta(telescope_id,
                                machine_id,
                                data_type,
                                "unknown",
                                "unknown",
                                barycentric,
                                pulsarcentric,
                                az_start,
                                za_start,
                                src_raj,
                                src_dej,
                                tstart,
                                tsamp,
                                nbits,
                                nsamples,
                                fch1,
                                foff,
                                0,
                                fchannel,
                                0,
                                nchans,
                                nifs,
                                refdm,
                                period,
                                strat,
                                ovrlp);

    metadata = meta;
    m_meta   = meta;

    return true;
  }


  /**
   * \brief Method to read the input data from the sigproc input file.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   */
  bool aa_sigproc_input::get_file_data(aa_filterbank_metadata &metadata, float max_dm_val) {
    double az_start = 0;
    double za_start = 0;
    double src_raj = 0;
    double src_dej = 0;
    double tstart = 0;
    double tsamp = 0;
    double refdm = 0;
    double period = 0;
    double fch1 = 0;
    double foff = 0;
    double fchannel = 0;
    
    int telescope_id = 0;
    int machine_id = 0;
    int data_type = 0;
    int barycentric = 0;
    int pulsarcentric = 0;
    int nbits = 0;
    int nsamples = 0;
    int nchans = 0;
    int nifs = 0;
    char *string = (char *) malloc(80 * sizeof(char));
// added by abhinav --------------------------------------------------
    char line[256], tempo_name[10], tempo_val[50], linetel[500];
    int count_num = 0;
    unsigned long int nsamp = 0;
    int dej1 = 0;
    int dej2 = 0;
    float dej3 = 0.0f;
    float dej, dej1f, dej2f;
    int raj1 = 0;
    int raj2 = 0;
    float raj3 = 0;
    float raj;
    int strat = 0;

    float sign_bw, bw_val, ch_band_width, time_res_tel, low_freq, high_freq, ovrlp;
    int int_val;


// edits end -------------------------------------------------------------    
    int nchar = 0;
    int nbytes = sizeof(int);
    

    FILE * pFile;
    pFile = fopen("/home/guest/Rushikesh/FRB_pipeline-1.7.10.GMRT300-500/input_files/time.hdr", "r");
    if(pFile == NULL) {
      printf("\n ATTENTION:: time.hdr cannot be opened! \n");
    }
   else {
   printf("\n File opened :: pFILE \n");
   while (fgets(line, sizeof(line), pFile)) {
        printf("\n **-*-** %.90s \n", line);

        strcpy(tempo_val, &line[35]);
        printf("\n temo_val is :: %s \n", tempo_val);	
        if (strncmp(line, "Time stamp of first sample (MJD) ", 33) == 0 )              // Time stamp of first sample (MJD) 
		{
//		tstart = strtof(tempo_val, NULL);
		tstart = atof(tempo_val);
		printf("\n Time stamp of first sample:: %s :: tstart:: %12.9f \n", tempo_val, tstart);
		}
    }  // time.hdr - while loop ends here ------------------------------------------------
    }  // else condition ends ---------------------------------------------
//    fclose(pFile);

// Reading data from gpu.hdr file::

	while (fgets(line, sizeof(line), fp)) {
//        printf("\n **-*-** %.35s \n", line);
		strcpy(tempo_val, &line[14]);	
//		printf("\n line:: line:: %s \n", linetel);
//		printf("\n tempo_val %s \n", tempo_val);
		if (strncmp(line, "GPU_RF4", 7) == 0 )               // Frequency of channel 1 (MHz)
		{
		fch1 = strtof(tempo_val, NULL);
		printf("\n Frequency of channel 1 (MHz): %f \n", fch1);
		}
		else if (strncmp(line, "GPU_CHAN_MAX", 12) == 0 )              // Time stamp of first sample (MJD) 
		{
		nchans = strtol(tempo_val, NULL, 10);
		printf("\n Number of channels: %d \n", nchans);
		}
		else if (strncmp(line, "GPU_SB4", 7) == 0 )               // "Channel bandwidth"
		{
		sign_bw = (-1)*strtof(tempo_val, NULL);
		printf("\n Channel bandwidth sign: %f \n", sign_bw);
		}
		else if (strncmp(line, "GPU_ACQ_BW", 10) == 0 )               // "Channel bandwidth"
		{
		bw_val = strtof(tempo_val, NULL);
		printf("\n Channel bandwidth size: %f \n", bw_val);
		}
		else if (strncmp(line, "GPU_BM_INT", 10) == 0 )               // "Channel bandwidth"
		{
		int_val = strtol(tempo_val, NULL, 10);
		printf("\n Beam Integration, i.e. Number of FFT cycles: %d \n", int_val);
		}
//	printf("\n %s ++--++-- \n", line); 
	}
	time_res_tel = (int_val*(nchans)*2)/(2*bw_val);		// in microseconds
	time_res_tel = time_res_tel/1000000.0;				// in seconds 

	tsamp = time_res_tel;
	printf("\n ** Time resolution :: %12.9f \n",tsamp);
	ch_band_width = (sign_bw*bw_val)/(nchans);

	foff = ch_band_width;
	printf("\n ** Channel Bandwith (foff): %9.8f \n",foff);

	high_freq = fch1;
	fch1 = (fch1)+(0.5*(foff));
	printf("\n fch1 is : %9.8f \n", fch1);	
// -reading of gpu.hdr ends here -------------------------------------------------------------------

    free(string);

// hard coded values ::    
    nbits = 8;
    nifs = 1;
   
// following should be changed if input file is changed- will be given by the telescope later
//    src_raj = 214400.0;
//    src_dej = -393300.0;

    FILE * pFile_source;
    pFile_source = fopen("/home/guest/Rushikesh/FRB_pipeline-1.7.10.GMRT300-500/input_files/source.hdr", "r");
    if(pFile_source == NULL) {
      printf("\n ATTENTION:: source.hdr cannot be opened! \n");
    }
   else {
   printf("\n File opened :: source.hdr \n");
   while (fgets(line, sizeof(line), pFile_source)) {
        printf("\n ** ** ** %.90s \n", line);

        strcpy(tempo_val, &line[9]);
        printf("\n temo_val is :: %s \n", tempo_val);	

	if (strncmp(line, "src_raj", 6) == 0 )              // Time stamp of first sample (MJD) 
		{
		src_raj = strtof(tempo_val, NULL);
		printf("\n raj: %f :: tempo_val: %s \n", src_raj, tempo_val);
		}
	else if (strncmp(line, "src_dej", 6) == 0 )              // Time stamp of first sample (MJD) 
		{
		src_dej = strtof(tempo_val, NULL);
		printf("\n dej: %f :: tempo_val: %s \n", src_dej, tempo_val);
		}
	} // while loop ends here
	} // else condition ends here

   if(!fclose(pFile_source)) printf("\n source.hdr file successfully closed \n");

//--------------------------------------------------------------------edits end----------------------------------------
    
    // Check that we are working with one IF channel
    if (nifs != 1) {
      printf("ERROR: Astro-Accelerate can only work with one IF channel.\n");
      return false;
    }

/*
    FILE * pFilee;
    pFilee = fopen("/data/jroy/data/FRB_DM500_163.84us_4K_7mP_4msW_3sig.header.fil", "r");
    if(pFilee == NULL) {
      printf("\n ATTENTION:: pFile cannot be opened! \n");
    }
   else {
   printf("\n File opened :: pFILEE \n");
    
//    fpos_t file_loc;
//    fgetpos(fp, &file_loc); //Position of end of header
      fseek(pFilee, 399, SEEK_SET);
      printf("\n Setting pointer to 399.\n");
    
//    unsigned long int nsamp = 0;

    // Getting number of time-samples based on file size
    unsigned long int data_start = ftell(pFilee);
    if (fseek(pFilee, 0, SEEK_END) != 0) {
      printf("\nERROR!! Failed to seek to the end of data file\n");
      return false;
    }
    unsigned long int exp_total_data = ftell(pFilee);
    exp_total_data = exp_total_data - data_start;
    fseek(pFilee, data_start, SEEK_SET);

    if (( nbits ) == 32) {
      nsamp = exp_total_data/((nchans)*4);
    }
    else if (( nbits ) == 16) {
      nsamp = exp_total_data/((nchans)*2);
    }
    else if (( nbits ) == 8) {
      nsamp = exp_total_data/((nchans));
    }
    else if (( nbits ) == 4) {
      nsamp = 2*exp_total_data/((nchans));
    }
    else if (( nbits ) == 2) {
      nsamp = 4*exp_total_data/((nchans));
    }
    else if (( nbits ) == 1) {
      nsamp = 8*exp_total_data/((nchans));
    }
    else {
      printf("ERROR: Currently this code only runs with 1, 2, 4, 8, and 16 bit data\n");
      return false;
    }

    // Move the file pointer back to the end of the header
//    fsetpos(fp, &file_loc);
    fclose(pFilee);
    }
*/
//    nsamp = 1024*1024;

    float epsilon = 0.0000001;
// number of time samples to be read depends upon the value of time between each time sample is taken::
    if (fabs(tsamp - 0.00008192) < epsilon) nsamp = 2*1024*1024;	//sampling_time=81.92 us
    else if(fabs(tsamp - 0.00016384) < epsilon) nsamp = 1*1024*1024;    //sampling_time=163.84 us
    else if(fabs(tsamp - 0.00032768) < epsilon) nsamp = 512*1024;       //sampling_time=327.68 us
    else if(fabs(tsamp - 0.00065536) < epsilon) nsamp = 256*1024;       //sampling_time=655.36 us
    else if(fabs(tsamp - 0.00131072) < epsilon) nsamp = 128*1024;       //sampling_time=1310.72 us 

    printf("\n nsamp is: %d \t tsamp: %12.9f \n", nsamp, tsamp);

// change the value hard coded below
    low_freq = high_freq + (sign_bw)*(bw_val);
    printf("\n low freq: %f high freq: %f  max_dm_val: %f \n", low_freq, high_freq, max_dm_val);
    ovrlp = 4.15*1000.*max_dm_val*((1./(low_freq*low_freq))-(1./(high_freq*high_freq)));
    printf("\n ovrlp:: %f \n", ovrlp);
    
    nsamples = (int)((nsamples) ? nsamples : (nsamp) ? nsamp : 0);
    aa_filterbank_metadata meta(telescope_id,
                                machine_id,
                                data_type,
                                "unknown",
                                "unknown",
                                barycentric,
                                pulsarcentric,
                                az_start,
                                za_start,
                                src_raj,
                                src_dej,
                                tstart,
                                tsamp,
                                nbits,
                                nsamples,
                                fch1,
                                foff,
                                0,
                                fchannel,
                                0,
                                nchans,
                                nifs,
                                refdm,
                                period,
                                strat,
                                ovrlp);

    metadata = meta;
    m_meta   = meta;

    return true;
  }



  /**
   * \brief Method to read the input data from the sigproc input file.
   * \returns A boolean flag to indicate whether the operation was successful (true) or not (false).
   */
  bool aa_sigproc_input::get_file_data(aa_filterbank_metadata &metadata, float max_dm_val, std::string config_path) {
    double az_start = 0;
    double za_start = 0;
    double src_raj = 0;
    double src_dej = 0;
    double tstart = 0;
    double tsamp = 0;
    double refdm = 0;
    double period = 0;
    double fch1 = 0;
    double foff = 0;
    double fchannel = 0;
    
    int telescope_id = 0;
    int machine_id = 0;
    int data_type = 0;
    int barycentric = 0;
    int pulsarcentric = 0;
    int nbits = 0;
    int nsamples = 0;
    int nchans = 0;
    int nifs = 0;
    char *string = (char *) malloc(80 * sizeof(char));
// added by abhinav --------------------------------------------------
    char line[256], tempo_name[10], tempo_val[50], linetel[500];
    int count_num = 0;
    unsigned long int nsamp = 0;
    int dej1 = 0;
    int dej2 = 0;
    float dej3 = 0.0f;
    float dej, dej1f, dej2f;
    int raj1 = 0;
    int raj2 = 0;
    float raj3 = 0;
    float raj;
    int strat = 0;

    float sign_bw, bw_val, ch_band_width, time_res_tel, low_freq, high_freq, ovrlp;
    int int_val;


// edits end -------------------------------------------------------------    
    int nchar = 0;
    int nbytes = sizeof(int);
    
    std::string time_file, source_file;
    time_file = config_path + "time.hdr";
    printf("\n time file::%s::test\n", time_file.c_str());

    FILE * pFile;
    pFile = fopen(time_file.c_str(), "r");
    if(pFile == NULL) {
      printf("\n ATTENTION:: time.hdr cannot be opened! \n");
    }
   else {
   printf("\n File opened :: pFILE \n");
   while (fgets(line, sizeof(line), pFile)) {
        printf("\n **-*-** %.90s \n", line);

        strcpy(tempo_val, &line[35]);
        printf("\n temo_val is :: %s \n", tempo_val);	
        if (strncmp(line, "Time stamp of first sample (MJD) ", 33) == 0 )              // Time stamp of first sample (MJD) 
		{
//		tstart = strtof(tempo_val, NULL);
		tstart = atof(tempo_val);
		printf("\n Time stamp of first sample:: %s :: tstart:: %12.9f \n", tempo_val, tstart);
		}
    }  // time.hdr - while loop ends here ------------------------------------------------
    }  // else condition ends ---------------------------------------------
//    fclose(pFile);

// Reading data from gpu.hdr file::

	while (fgets(line, sizeof(line), fp)) {
//        printf("\n **-*-** %.35s \n", line);
		strcpy(tempo_val, &line[14]);	
//		printf("\n line:: line:: %s \n", linetel);
//		printf("\n tempo_val %s \n", tempo_val);
		if (strncmp(line, "GPU_RF4", 7) == 0 )               // Frequency of channel 1 (MHz)
		{
		fch1 = strtof(tempo_val, NULL);
		printf("\n Frequency of channel 1 (MHz): %f \n", fch1);
		}
		else if (strncmp(line, "GPU_CHAN_MAX", 12) == 0 )              // Time stamp of first sample (MJD) 
		{
		nchans = strtol(tempo_val, NULL, 10);
		printf("\n Number of channels: %d \n", nchans);
		}
		else if (strncmp(line, "GPU_SB4", 7) == 0 )               // "Channel bandwidth"
		{
		sign_bw = (-1)*strtof(tempo_val, NULL);
		printf("\n Channel bandwidth sign: %f \n", sign_bw);
		}
		else if (strncmp(line, "GPU_ACQ_BW", 10) == 0 )               // "Channel bandwidth"
		{
		bw_val = strtof(tempo_val, NULL);
		printf("\n Channel bandwidth size: %f \n", bw_val);
		}
		else if (strncmp(line, "GPU_BM_INT", 10) == 0 )               // "Channel bandwidth"
		{
		int_val = strtol(tempo_val, NULL, 10);
		printf("\n Beam Integration, i.e. Number of FFT cycles: %d \n", int_val);
		}
//	printf("\n %s ++--++-- \n", line); 
	}
	time_res_tel = (int_val*(nchans)*2)/(2*bw_val);		// in microseconds
	time_res_tel = time_res_tel/1000000.0;				// in seconds 

	tsamp = time_res_tel;
	printf("\n ** Time resolution :: %12.9f \n",tsamp);
	ch_band_width = (sign_bw*bw_val)/(nchans);

	foff = ch_band_width;
	printf("\n ** Channel Bandwith (foff): %9.8f \n",foff);

	high_freq = fch1;
	fch1 = (fch1)+(0.5*(foff));
	printf("\n fch1 is : %9.8f \n", fch1);	
// -reading of gpu.hdr ends here -------------------------------------------------------------------

    free(string);

// hard coded values ::    
    nbits = 8;
    nifs = 1;
   
// following should be changed if input file is changed- will be given by the telescope later
//    src_raj = 214400.0;
//    src_dej = -393300.0;

    source_file = config_path + "source.hdr";
    printf("\n source file::%s::test\n", source_file.c_str());

    FILE * pFile_source;
    pFile_source = fopen(source_file.c_str(), "r");
    pFile_source = fopen("/home/guest/Rushikesh/FRB_pipeline-1.7.10.GMRT300-500/input_files/source.hdr", "r");
    if(pFile_source == NULL) {
      printf("\n ATTENTION:: source.hdr cannot be opened! \n");
    }
   else {
   printf("\n File opened :: source.hdr \n");
   while (fgets(line, sizeof(line), pFile_source)) {
        printf("\n ** ** ** %.90s \n", line);

        strcpy(tempo_val, &line[9]);
        printf("\n temo_val is :: %s \n", tempo_val);	

	if (strncmp(line, "src_raj", 6) == 0 )              // Time stamp of first sample (MJD) 
		{
		src_raj = strtof(tempo_val, NULL);
		printf("\n raj: %f :: tempo_val: %s \n", src_raj, tempo_val);
		}
	else if (strncmp(line, "src_dej", 6) == 0 )              // Time stamp of first sample (MJD) 
		{
		src_dej = strtof(tempo_val, NULL);
		printf("\n dej: %f :: tempo_val: %s \n", src_dej, tempo_val);
		}
	} // while loop ends here
	} // else condition ends here

   if(!fclose(pFile_source)) printf("\n source.hdr file successfully closed \n");

//--------------------------------------------------------------------edits end----------------------------------------
    
    // Check that we are working with one IF channel
    if (nifs != 1) {
      printf("ERROR: Astro-Accelerate can only work with one IF channel.\n");
      return false;
    }

/*
    FILE * pFilee;
    pFilee = fopen("/data/jroy/data/FRB_DM500_163.84us_4K_7mP_4msW_3sig.header.fil", "r");
    if(pFilee == NULL) {
      printf("\n ATTENTION:: pFile cannot be opened! \n");
    }
   else {
   printf("\n File opened :: pFILEE \n");
    
//    fpos_t file_loc;
//    fgetpos(fp, &file_loc); //Position of end of header
      fseek(pFilee, 399, SEEK_SET);
      printf("\n Setting pointer to 399.\n");
    
//    unsigned long int nsamp = 0;

    // Getting number of time-samples based on file size
    unsigned long int data_start = ftell(pFilee);
    if (fseek(pFilee, 0, SEEK_END) != 0) {
      printf("\nERROR!! Failed to seek to the end of data file\n");
      return false;
    }
    unsigned long int exp_total_data = ftell(pFilee);
    exp_total_data = exp_total_data - data_start;
    fseek(pFilee, data_start, SEEK_SET);

    if (( nbits ) == 32) {
      nsamp = exp_total_data/((nchans)*4);
    }
    else if (( nbits ) == 16) {
      nsamp = exp_total_data/((nchans)*2);
    }
    else if (( nbits ) == 8) {
      nsamp = exp_total_data/((nchans));
    }
    else if (( nbits ) == 4) {
      nsamp = 2*exp_total_data/((nchans));
    }
    else if (( nbits ) == 2) {
      nsamp = 4*exp_total_data/((nchans));
    }
    else if (( nbits ) == 1) {
      nsamp = 8*exp_total_data/((nchans));
    }
    else {
      printf("ERROR: Currently this code only runs with 1, 2, 4, 8, and 16 bit data\n");
      return false;
    }

    // Move the file pointer back to the end of the header
//    fsetpos(fp, &file_loc);
    fclose(pFilee);
    }
*/
//    nsamp = 1024*1024;

    float epsilon = 0.0000001;
// number of time samples to be read depends upon the value of time between each time sample is taken::
    if (fabs(tsamp - 0.00008192) < epsilon) nsamp = 2*1024*1024;	//sampling_time=81.92 us
    else if(fabs(tsamp - 0.00016384) < epsilon) nsamp = 1*1024*1024;    //sampling_time=163.84 us
    else if(fabs(tsamp - 0.00032768) < epsilon) nsamp = 512*1024;       //sampling_time=327.68 us
    else if(fabs(tsamp - 0.00065536) < epsilon) nsamp = 256*1024;       //sampling_time=655.36 us
    else if(fabs(tsamp - 0.00131072) < epsilon) nsamp = 128*1024;       //sampling_time=1310.72 us 

    printf("\n nsamp is: %d \t tsamp: %12.9f \n", nsamp, tsamp);

// change the value hard coded below
    low_freq = high_freq + (sign_bw)*(bw_val);
    printf("\n low freq: %f high freq: %f  max_dm_val: %f \n", low_freq, high_freq, max_dm_val);
    if (high_freq > low_freq)    ovrlp = 4.15*1000.*max_dm_val*((1./(low_freq*low_freq))-(1./(high_freq*high_freq)));
    else if (high_freq < low_freq)    ovrlp = 4.15*1000.*max_dm_val*((1./(high_freq*high_freq))-(1./(low_freq*low_freq)));
    printf("\n ovrlp:: %f \n", ovrlp);
    
    nsamples = (int)((nsamples) ? nsamples : (nsamp) ? nsamp : 0);
    aa_filterbank_metadata meta(telescope_id,
                                machine_id,
                                data_type,
                                "unknown",
                                "unknown",
                                barycentric,
                                pulsarcentric,
                                az_start,
                                za_start,
                                src_raj,
                                src_dej,
                                tstart,
                                tsamp,
                                nbits,
                                nsamples,
                                fch1,
                                foff,
                                0,
                                fchannel,
                                0,
                                nchans,
                                nifs,
                                refdm,
                                period,
                                strat,
                                ovrlp);

    metadata = meta;
    m_meta   = meta;

    return true;
  }


  bool aa_sigproc_input::get_file_data_new(aa_filterbank_metadata &metadata) {
/*    double az_start = 0;
    double za_start = 0;
    double src_raj = 0;
    double src_dej = 0;
    double tstart = 0;
    double tsamp = 0;
    double refdm = 0;
    double period = 0;
    double fch1 = 0;
    double foff = 0;
    double fchannel = 0;
    
    int telescope_id = 0;
    int machine_id = 0;
    int data_type = 0;
    int barycentric = 0;
    int pulsarcentric = 0;
    int nbits = 0;
    int nsamples = 0;
    int nchans = 0;
    int nifs = 0;
    char *string = (char *) malloc(80 * sizeof(char));
// added by abhinav --------------------------------------------------
    char line[256], tempo_name[10], tempo_val[50], linetel[500];
    int count_num = 0;
    unsigned long int nsamp = 0;
    int dej1 = 0;
    int dej2 = 0;
    float dej3 = 0.0f;
    float dej, dej1f, dej2f;
    int raj1 = 0;
    int raj2 = 0;
    float raj3 = 0;
    float raj;
    int strat = 1;
// edits end -------------------------------------------------------------    
    int nchar = 0;
    int nbytes = sizeof(int);
    
*/
//    printf("\n  filterbank_data.nchans:: is :: %d \n", filterbank_data.nchans());
    printf("\n  metadata.nchans:: is :: %d \n", metadata.nchans());

//    printf("\n  m_meta.nchans:: is :: %d \n", m_meta.nchans());
//    printf("\n  metadata.nchans:: is :: %d \n", m_meta_new.nchans());


//    nsamples = (int)((nsamples) ? nsamples : (nsamp) ? nsamp : 0);
    aa_filterbank_metadata meta(metadata.telescope_id(),
                                metadata.machine_id(),
                                metadata.data_type(),
                                "unknown",
                                "unknown",
                                metadata.barycentric(),
                                metadata.pulsarcentric(),
                                metadata.az_start(),
                                metadata.za_start(),
                                metadata.src_raj(),
                                metadata.src_dej(),
                                metadata.tstart(),
                                metadata.tsamp(),
                                metadata.nbits(),
                                metadata.nsamples(),
                                metadata.fch1(),
                                metadata.foff(),
                                0,
                                metadata.fchannel(),
                                0,
                                metadata.nchans(),
                                metadata.nifs(),
                                metadata.refdm(),
                                metadata.period(),
                                1,
                                metadata.ovrlp());

//    printf("\n  filterbank_data.nchans:: is :: %d \n", filterbank_data.nchans());
//    printf("\n  meta.nchans:: is :: %d \n", meta.nchans());

    metadata = meta;
    m_meta = meta;
    m_meta_new   = meta;
    printf("\n  aa_sigproc_input :: metadata.nchans:: is :: %d \n", metadata.nchans());
    printf("\n  aa_sigproc_input :: metadata.strat:: is :: %d \n", metadata.strat());
    printf("\n  aa_sigproc_input :: metadata.nsamples:: is :: %d \n", metadata.nsamples());


    return true;
  }

/*
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

  bool aa_sigproc_input::get_input_file_data(aa_filterbank_metadata &metadata_input) {
    


    return true;
  }

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
*/

  template <typename T>
  bool aa_sigproc_input::get_recorded_data(std::vector<T> &input_buffer) {
  printf("\n Warning: get_recorded_data function is deprecated! Use get_recorded_data_new instead! \n");

    const size_t inputsize = (size_t)m_meta.nsamples() * (size_t)m_meta.nchans();
    input_buffer.resize(inputsize);
    int c;
    
    unsigned long int total_data;
    const int nchans = m_meta.nchans();
    const int nbits = m_meta.nbits();

/*    const size_t inputsize = 15000928256;
    input_buffer.resize(inputsize);
    int c;
    
    unsigned long int total_data;
    const int nchans = 4096;
    const int nbits = 8;
    const int nsamples = 3662336;
*/

//    aa_filterbank_metadata metadata;
//    printf("\nnbits:: %lu \n", &metadata.nbits);

//     printf("\n size is:: %zu \n", inputsize);
//     printf("\n nbits: %ld\n", nbits);

    //{{{ Load in the raw data from the input file and transpose
    if (nbits == 32) {
      // Allocate a tempory buffer to store a line of frequency data
      float *temp_buffer = (float *) malloc(nchans * sizeof(float));
        
      // Allocate floats to hold the mimimum and maximum value in the input data
      float max = -10000000000;
      float min = 10000000000;
        
      // Allocate a variable to hold the file pointer position
      fpos_t file_loc;
        
      // Set the file pointer position
      fgetpos(fp, &file_loc);
        
      // Find the minimum and maximum values in the input file.
      while (!feof(fp)) {
	if (fread(temp_buffer, sizeof(float), nchans, fp) != (size_t)nchans)
	  break;
	for (c = 0; c < nchans; c++) {
	  if(temp_buffer[c] > max) max = temp_buffer[c];
	  if(temp_buffer[c] < min) min = temp_buffer[c];
	}
      }
        
      // Calculate the bin size in a distribution of unsigned shorts for the input data.
      float bin = (max - min) / 65535.0f;
        
      printf("\n Conversion Parameters: %f\t%f\t%f", min, max, bin);
        
      // Move the file pointer back to the end of the header
      fsetpos(fp, &file_loc);
        
      // Read in the data, scale to an unsigned short range, transpose it and store it in the input buffer
      total_data = 0;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(float), nchans, fp) != (size_t)nchans)
	  break;
	for (c = 0; c < nchans; c++) {
	  ( input_buffer )[c + total_data * ( nchans )] = (unsigned short) ((temp_buffer[c]-min)/bin);
	}
	total_data++;
      }
      free(temp_buffer);
    }
    else if (nbits == 16) {
      // Allocate a tempory buffer to store a line of frequency data
      unsigned short *temp_buffer = (unsigned short *) malloc(nchans * sizeof(unsigned short));
        
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned short), nchans, fp) != (size_t)nchans)
	  break;
	for (c = 0; c < nchans; c++) {
	  ( input_buffer )[c + total_data * ( nchans )] = (unsigned short) temp_buffer[c];
	}
	total_data++;
      }
      free(temp_buffer);
    }
    else if (nbits == 8) {

      printf("\n We are here with nbits: %d\n",nbits);

      fseek(fp, 399, SEEK_SET);
      printf("\n Setting pointer to 399.\n");
        
      // Allocate a tempory buffer to store a line of frequency data
      unsigned char *temp_buffer = (unsigned char *) malloc(nchans * sizeof(unsigned char));
        
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;

      while (!feof(fp)) {
//      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
//	if(total_data % (int)(nsamples*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, nsamples, (int)(100.0*((float)total_data/(float)nsamples)));
	if (fread(temp_buffer, sizeof(unsigned char), nchans, fp) != (size_t)nchans) {
	  break;
	}
	for (c = 0; c < nchans; c++) {
	  ( input_buffer )[c + total_data * ( nchans )] = (unsigned short) temp_buffer[c];
//          if(c < 700 && total_data == 0) printf("\n %hu \n", ( input_buffer )[c + total_data * ( nchans )]);
	}
	total_data++;
/*        if (total_data == nsamples) 
         {
           break;
         }
*/      }
      free(temp_buffer);
    }
    else if (nbits == 4) {
      // Allocate a temporary buffer to store a line of frequency data
      // each byte stores 2 frequency data
      int nb_bytes = nchans/2;
      unsigned char *temp_buffer = (unsigned char *) malloc(nb_bytes * sizeof(unsigned char));
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      // 00001111
      char mask = 0x0f;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
	  break;
	for (c = 0; c < nb_bytes; c++) {
	  // (n >> a) & ( (1 << a) - 1) -> right shift by 'a' bits, then keep the last 'b' bits
	  // Here, we right shift 4 bits and keep the last 4 bits
	  ( input_buffer )[ (c*2) + total_data * ( nchans )]     = (unsigned short)( (temp_buffer[c] >> 4) & mask );
	  // n & ( (1 << a ) - 1)
	  // Here we keep the last 4 bits
	  ( input_buffer )[ (c*2) + 1 + total_data * ( nchans )] = (unsigned short)( temp_buffer[c] & mask );
	}
	total_data++;
      }
      free(temp_buffer);
    }
    else if (nbits == 2) {
      // Allocate a temporary buffer to store a line of frequency data
      // each byte stores 4 frequency data
      int nb_bytes = nchans/4;
      unsigned char *temp_buffer = (unsigned char *) malloc(nb_bytes * sizeof(unsigned char));
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      // 00001111
      char mask = 0x03;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
	  break;
	for (c = 0; c < nb_bytes; c++) {
	  ( input_buffer )[ (c*4) + total_data * ( nchans )]     = (unsigned short)( (temp_buffer[c] >> 6) & mask );
	  ( input_buffer )[ (c*4) + 1 + total_data * ( nchans )] = (unsigned short)( (temp_buffer[c] >> 4) & mask );
	  ( input_buffer )[ (c*4) + 2 + total_data * ( nchans )] = (unsigned short)( (temp_buffer[c] >> 2) & mask );
	  ( input_buffer )[ (c*4) + 3 + total_data * ( nchans )] = (unsigned short)( temp_buffer[c] & mask );
	}
	total_data++;
      }
      free(temp_buffer);
    }
    else if (nbits == 1) {
      // each byte stores 8 frequency data
      int nb_bytes = nchans/8;
        
      // Allocate a temporary buffer to store a line of frequency data
      unsigned char *temp_buffer = (unsigned char *) malloc(nb_bytes * sizeof(unsigned char));
        
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
        
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
	  break;
            
	for (c = 0; c < nb_bytes; c++) {
	  for(int i=0; i<8; i++) {
	    unsigned char mask =  1 << i;
	    unsigned char masked_char = temp_buffer[c] & mask;
	    unsigned char bit = masked_char >> i;
	    ( input_buffer )[ (c*8) + i + total_data * ( nchans )] = (unsigned short)bit;
	  }
	}
	total_data++;
      }
      free(temp_buffer);
    }
    else {
      printf("ERROR: Invalid number of bits in input data.\n");
      return false;
    }

    data_is_read = true;
    return true;
  }


  template <typename T>
  bool aa_sigproc_input::get_new_recorded_data(std::vector<T> &input_buffer, int buf_count, const aa_filterbank_metadata filterbank_data, unsigned long int **fp_pos) {
  printf("\nwe are here in get_new_recorded_data\n");
/*
    const size_t inputsize = (size_t)m_meta.nsamples() * (size_t)m_meta.nchans();
    input_buffer.resize(inputsize);
    int c;
    
    unsigned long int total_data;
    const int nchans = m_meta.nchans();
    const int nbits = m_meta.nbits();
*/

    printf("\n get_new_recorded_data::filterbank_data::nchans:: %d \n", filterbank_data.nchans());

    int c;
    
    unsigned long int total_data;
//    const int nchans = 4096;
//    const int nbits = 8;
//    const int nsamples = 3662336;
    const int nchans = filterbank_data.nchans();
    const int nbits = filterbank_data.nbits();
    const int nsamples = filterbank_data.nsamples();
    const int nsamp = nsamples;
    const int foff = filterbank_data.foff();
    double tsamp = filterbank_data.tsamp();
//  ----- Added by Rushikesh -----
    static short int data[IABeamBufferSize]; //IABeamBufferSize is the size of one_block_SM
    static short int data1[IABeamBufferSize];
    unsigned int curRec = 0, curBlock;
    int fp_PA;
    int i, j, iteration = 0;

//    const size_t inputsize = 15000928256;
//    const size_t inputsize = 34393296896;
//    const size_t inputsize = nsamp*nchans;
    const size_t inputsize = (size_t)filterbank_data.nsamples() * (size_t)filterbank_data.nchans();
    input_buffer.resize(inputsize);


    printf("\n get_new_recorded_data::::::nchans:: %d \n", nchans);
    printf("\n get_new_recorded_data::::::nbits:: %d \n", nbits);
    printf("\n get_new_recorded_data::::::nsamp:: %ld \n", nsamp);
    printf("\n get_new_recorded_data::::::inputsize:: %zu \n", inputsize);

    int data_block_size = 8192;   // reading 8192 time samples at a time to reduce file I/O operations

    static int nsamp_ovrlp_val;
    int t_samples_per_block;
    int nsamp_ovrlp_pt, block_num; 

// deciding how many time samples will be overlapped in each reading.
    float epsilon = 1.0;
    if (buf_count == 0)
    {
    if (filterbank_data.ovrlp() >= 59 ) 
    {
    printf("\n ovralp is greater than 59 seconds. buf_count is %d \n", buf_count); 
    if (nsamp == 2*1024*1024) nsamp_ovrlp_val = 768*1024;	
    else if(nsamp == 1*1024*1024) nsamp_ovrlp_val = 384*1024;
    else if(nsamp == 512*1024) nsamp_ovrlp_val = 192*1024;
    else if(nsamp == 256*1024) nsamp_ovrlp_val = 96*1024;
    else if(nsamp == 128*1024) nsamp_ovrlp_val = 48*1024;
    }
    else
    {
    printf("\n ovralp is less than 59 seconds. buf_count is %d \n", buf_count);
    if (nsamp == 2*1024*1024) nsamp_ovrlp_val = 256*1024;	
    else if(nsamp == 1*1024*1024) nsamp_ovrlp_val = 128*1024;
    else if(nsamp == 512*1024) nsamp_ovrlp_val = 64*1024;
    else if(nsamp == 256*1024) nsamp_ovrlp_val = 32*1024;
    else if(nsamp == 128*1024) nsamp_ovrlp_val = 16*1024;
    }
    }  // buf_count condition ends
    printf("\n buf count::nsamp_ovrlp_val:: %d \n", nsamp_ovrlp_val);

//    nsamp_ovrlp_pt = (nsamp - nsamp_ovrlp_val)/data_block_size;
    nsamp_ovrlp_pt = nsamp - nsamp_ovrlp_val;
    block_num = nsamp/data_block_size;              // Number of blocks in each part of buf_count loop
//    t_samples_per_block = data_block_size*nchans; 	// Number of time samples in each block = block number * 4096 channels

    printf("\n Number of time samples: %d \n", nsamp);
    printf("\n Time samples overlap between subsequent readings: %d\n", nsamp_ovrlp_val);
    printf("\n Overlapped data beigns after first %d samples \n",nsamp_ovrlp_pt);
    printf("\n Number of blocks %d to be read at a time \n",block_num);
//    printf("\n Maximum number of iterations for reading the file: %d \n", total_iter);


    //{{{ Load in the raw data from the input file and transpose
    if (nbits == 32) {
      // Allocate a tempory buffer to store a line of frequency data
      float *temp_buffer = (float *) malloc(nchans * sizeof(float));
        
      // Allocate floats to hold the mimimum and maximum value in the input data
      float max = -10000000000;
      float min = 10000000000;        
        
      // Find the minimum and maximum values in the input file.
      while (!feof(fp)) {
	if (fread(temp_buffer, sizeof(float), nchans, fp) != (size_t)nchans)
	  break;
	for (c = 0; c < nchans; c++) {
	  if(temp_buffer[c] > max) max = temp_buffer[c];
	  if(temp_buffer[c] < min) min = temp_buffer[c];
	}
      }
        
      // Calculate the bin size in a distribution of unsigned shorts for the input data.
      float bin = (max - min) / 65535.0f;
        
      printf("\n Conversion Parameters: %f\t%f\t%f", min, max, bin);
        
        
      // Read in the data, scale to an unsigned short range, transpose it and store it in the input buffer
      total_data = 0;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(float), nchans, fp) != (size_t)nchans)
	  break;
	for (c = 0; c < nchans; c++) {
	  ( input_buffer )[c + total_data * ( nchans )] = (unsigned short) ((temp_buffer[c]-min)/bin);
	}
	total_data++;
      }
      free(temp_buffer);
    }
    else if (nbits == 16) {
      // Allocate a tempory buffer to store a line of frequency data
      unsigned short *temp_buffer = (unsigned short *) malloc(nchans * sizeof(unsigned short));
        
      // Read in the data, transpose it and store it in the input buffer
      total_data = 0;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned short), nchans, fp) != (size_t)nchans)
	  break;
	for (c = 0; c < nchans; c++) {
	  ( input_buffer )[c + total_data * ( nchans )] = (unsigned short) temp_buffer[c];
	}
	total_data++;
      }
      free(temp_buffer);
    }
    else if (nbits == 8) {

      printf("\n We are here with nbits: %d\n",nbits);

      if (buf_count ==0){
      fseek(fp, 0, SEEK_SET); // use this with GPT files
      printf("\n Setting pointer to start of file !!\n");
      }
      else {
      fseek(fp, **fp_pos, SEEK_SET);  // use this with GPT files
      printf("\n Setting pointer to poisiton: %lu.\n", **fp_pos);
      printf("\n Setting pointer after %ld time samples.\n", buf_count*nsamp_ovrlp_pt);
      }
  
      const size_t one_block_SM = (size_t)filterbank_data.nsamples() * (size_t)filterbank_data.nchans() * 0.125; //RD: Each block is 21.474 sec
      const size_t telescope_time_one_block_SM = nsamp*tsamp/8; 
      printf("\n one_block_SM: %d\n", one_block_SM);
      printf("\n telescope_time_one_block_SM: %f\n", telescope_time_one_block_SM);//RD: Should always be 21.474636648 sec
      printf("\n filterbank_data.nsamples(): %d\n", filterbank_data.nsamples()); 
      printf("\n nsamp  :    %d\n", nsamp);
      

      unsigned char *temp_buffer = (unsigned char *) malloc(one_block_SM * sizeof(unsigned char));
     
     total_data = 0;
     
      while (!feof(fp)) {
        //if (fread(temp_buffer, sizeof(unsigned char), one_block_SM, fp) != (size_t)one_block_SM) 

        if (memcpy(temp_buffer, dataBuf->temp_buffer+curBlock*one_block_SM, one_block_SM)){
          printf("\n Number of records read: %d\n", total_data);
	  break;
	}
       
	  for (c = 0; c < one_block_SM; c++) {
	    ( input_buffer )[c + total_data * ( one_block_SM )] = (unsigned short) temp_buffer[c];
            if(c < 3 && total_data == 0) printf("\n #!# %hu \n", ( input_buffer )[c + total_data * ( one_block_SM )]); 
	  }
	
	
	total_data++;
	printf("\n Number of data blocks read: %d\n", total_data);
	
	if (total_data == 8)
	{
	  if (filterbank_data.ovrlp() >= 59) **fp_pos = 5*(buf_count+1)*one_block_SM;	
          else if(filterbank_data.ovrlp() < 59) **fp_pos = 7*(buf_count+1)*one_block_SM;             
        printf("\n Setting file pointer after reading data of %ld data blocks. File position: %lu\n", total_data, **fp_pos);
        printf("\nBreaking after reading %d blocks\n", total_data);
	break;
	} 
      } 
     
      if (total_data != 8)
	{
	printf("\n Not enough data to process \n");
        data_is_read = false;
        return false;
	}

      free(temp_buffer);
//      free(temp_buffer_2);


    }   /*############### nbits == 8 ends here ##################*/
    else if (nbits == 4) {
      int nb_bytes = nchans/2;
      unsigned char *temp_buffer = (unsigned char *) malloc(nb_bytes * sizeof(unsigned char));
      total_data = 0;
      // 00001111
      char mask = 0x0f;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
	  break;
	for (c = 0; c < nb_bytes; c++) {
	  ( input_buffer )[ (c*2) + total_data * ( nchans )]     = (unsigned short)( (temp_buffer[c] >> 4) & mask );
	  ( input_buffer )[ (c*2) + 1 + total_data * ( nchans )] = (unsigned short)( temp_buffer[c] & mask );
	}
	total_data++;
      }
      free(temp_buffer);
    }
    else if (nbits == 2) {
      int nb_bytes = nchans/4;
      unsigned char *temp_buffer = (unsigned char *) malloc(nb_bytes * sizeof(unsigned char));
      total_data = 0;
      // 00001111
      char mask = 0x03;
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
	  break;
	for (c = 0; c < nb_bytes; c++) {
	  ( input_buffer )[ (c*4) + total_data * ( nchans )]     = (unsigned short)( (temp_buffer[c] >> 6) & mask );
	  ( input_buffer )[ (c*4) + 1 + total_data * ( nchans )] = (unsigned short)( (temp_buffer[c] >> 4) & mask );
	  ( input_buffer )[ (c*4) + 2 + total_data * ( nchans )] = (unsigned short)( (temp_buffer[c] >> 2) & mask );
	  ( input_buffer )[ (c*4) + 3 + total_data * ( nchans )] = (unsigned short)( temp_buffer[c] & mask );
	}
	total_data++;
      }
      free(temp_buffer);
    }
    else if (nbits == 1) {
      int nb_bytes = nchans/8;
      unsigned char *temp_buffer = (unsigned char *) malloc(nb_bytes * sizeof(unsigned char));
      total_data = 0;
        
      while (!feof(fp)) {
	if(total_data % (int)(m_meta.nsamples()*0.1) == 0) printf("Reading record %10lu / %10d \t (%2d%% complete)\n", total_data, m_meta.nsamples(), (int)(100.0*((float)total_data/(float)m_meta.nsamples())));
	if (fread(temp_buffer, sizeof(unsigned char), nb_bytes, fp) != (size_t)nb_bytes)
	  break;
            
	for (c = 0; c < nb_bytes; c++) {
	  for(int i=0; i<8; i++) {
	    unsigned char mask =  1 << i;
	    unsigned char masked_char = temp_buffer[c] & mask;
	    unsigned char bit = masked_char >> i;
	    ( input_buffer )[ (c*8) + i + total_data * ( nchans )] = (unsigned short)bit;
	  }
	}
	total_data++;
      }
      free(temp_buffer);
    }
    else {
      printf("ERROR: Invalid number of bits in input data.\n");
      return false;
    }

    data_is_read = true;
    return true;
  }

} //namespace astroaccelerate

void initialise () {

  shmHId = shmget( DasHeaderKey, sizeof( DataHeader ), SHM_RDONLY ); /*shmget() returns an identifier for the shared memory segment
                                 associated with DasHeaderKey. Here, in place of data header we can use the 'one_block_sm'*/
  shmBId = shmget( DasBufferKey, sizeof( DataBufferIA ), SHM_RDONLY );
  printf(shmHID, shmBID)

  if( shmHId < 0 || shmBId < 0 ) {
    fprintf(stderr, "Error in attaching shared memory..\n");
    exit(-1);
  }

  dataHdr = (DataHeader *) shmat( shmHId, 0, 0 ); /*shmat(). void *shmat(int shmid ,void *shmaddr ,int shmflg) is used to attach shared the memory to the dataHdr; where, shmid is shared memory id. shmaddr specifies specific address to use but we should set it to zero and OS will automatically choose the address.*/
  dataBuf = (DataBufferIA *) shmat( shmBId, 0, 0 ); /*shmat attaches a shared memory segment with identifies shmBID to dataBuf*/
}
