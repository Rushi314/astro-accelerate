#ifndef ASTRO_ACCELERATE_AA_SIGPROC_INPUT_HPP
#define ASTRO_ACCELERATE_AA_SIGPROC_INPUT_HPP

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "aa_input.hpp"

#include "aa_ddtr_plan.hpp" // added by abhinav

#include <string>

namespace astroaccelerate {

  /**
   * \class aa_sigproc_input aa_sigproc_input.hpp "include/aa_sigproc_input.hpp"
   * \brief The aa_sigproc_input class is used to parse a sigproc (.fil) file.
   * \details The aa_sigproc_input class is used to parse a sigproc (.fil) file.
   * \details Implementation details provided in source src/aa_sigproc_input.cpp.
   * \details Data description from "SIGPROC-v3.7 (Pulsar) Signal Processing Programs"
   * \details Source: http://sigproc.sourceforge.net/sigproc.pdf
   * \author Cees Carels.
   * \date 5 November 2018.
   */
  class aa_sigproc_input : public aa_input {
  public:

    /** \brief Constructor for aa_sigproc_input, requires a file path. */
    aa_sigproc_input(const std::string &path);

    /** \brief Destructor for aa_sigproc_input. */
    virtual ~aa_sigproc_input();

    /** \brief Method to open the input file. */
    bool open();

    /** \brief Method to close the input file. */
    bool close();

    /**
     * \brief Method to read only the metadata from the input file.
     * \returns an aa_filterbank_metadata object with the metadata from the input file. */
//    aa_filterbank_metadata read_metadata();
    aa_filterbank_metadata read_metadata();
    aa_filterbank_metadata read_metadata(float max_dm_val);
    aa_filterbank_metadata read_metadata(float max_dm_val, std::string configpath_name);
    aa_filterbank_metadata read_metadata_new(aa_filterbank_metadata metadata);
    aa_filterbank_metadata read_metadata_source(aa_filterbank_metadata metadata);
    bool read_metadata_input(float *sigma_cutoff_input, float *sigma_constant_input, float *max_boxcar_width_in_sec_input, float *periodicity_sigma_cutoff_input, float *periodicity_harmonics_input, aa_ddtr_plan *ddtr_plan, bool *baselinenoise_val, bool *analysis_val, bool *debug_val, char (*file_name_input)[128], char (*config_path)[128]);


    /** \brief Method to read only the input data from the input file. */
    bool read_signal();
    bool read_new_signal(int buf_count, const aa_filterbank_metadata &filterbank_data, unsigned long int *f_pos);

    /** \returns The input data from the telescope. */
    const std::vector<unsigned short>& input_buffer() const {
      return m_input_buffer;
    }

    /** \returns A modifiable reference to the input data from the telescope. */
    std::vector<unsigned short>& input_buffer_modifiable() {
      return m_input_buffer;
    }

    /** \brief Method to check if the input data from the input file has already been read.
     * \returns A boolean flag to indicate if the data has been read (true) or not (false). */
    bool did_read_signal() const {
      return data_is_read;
    }
    

  private:
    
    bool get_file_data(aa_filterbank_metadata &metadata); /** Reads the metadata from the filterbank input file. */
    bool get_file_data(aa_filterbank_metadata &metadata, float max_dm_val); /** Reads the metadata from the filterbank input file. */
    bool get_file_data(aa_filterbank_metadata &metadata, float max_dm_val, std::string config_path); /** Reads the metadata from the filterbank input file. */
    aa_filterbank_metadata m_meta; /** Stores the metadata associated with the input file. */

    bool get_file_data_new(aa_filterbank_metadata &metadata); /** Reads the metadata from the filterbank input file. */
    aa_filterbank_metadata m_meta_new; /** Stores the metadata associated with the input file. */

    bool get_input_file_data(aa_filterbank_metadata &metadata_input); /** Reads the metadata from the filterbank input file. */
    aa_filterbank_metadata m_meta_input; /** Stores the metadata associated with the input file. */

    std::vector<unsigned short> m_input_buffer; /** Stores the data in the sigproc file. */
    
    template <typename T>
    bool get_recorded_data(std::vector<T> &input_buffer);

    template <typename T>
    bool get_new_recorded_data(std::vector<T> &input_buffer, int buf_count, const aa_filterbank_metadata filterbank_data, unsigned long int **f_pos);
    
    bool header_is_read; /** Flag to indicate whether the input file header (metadata) has been read. */
    bool data_is_read; /** Flag to indicate whether the input file data has been read. */

    
    std::string file_path; /** Stores the file path to the input file. */
    FILE *fp; /** File pointer to the input file. */
    
  };

} // namespace astroaccelerate

#endif // ASTRO_ACCELERATE_AA_SIGPROC_INPUT_HPP
