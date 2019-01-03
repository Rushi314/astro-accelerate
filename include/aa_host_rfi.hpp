#ifndef ASTRO_ACCELERATE_AA_HOST_RFI_HPP
#define ASTRO_ACCELERATE_AA_HOST_RFI_HPP

namespace astroaccelerate {

  /** \brief Function that performs RFI mitigation. */
  void rfi(int nsamp, int nchans, unsigned short **input_buffer);

} // namespace astroaccelerate
  
#endif // ASTRO_ACCELERATE_AA_HOST_RFI_HPP

