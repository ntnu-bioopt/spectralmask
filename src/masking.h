//==============================================================================
// Copyright 2015 Asgeir Bjorgan, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//==============================================================================

#ifndef MASKING_H_DEFINED
#define MASKING_H_DEFINED


/**
 * Masking parameters. Reference spectra and so on.  
 **/
typedef struct{
	/// Number of reference spectra 
	int num_masking_spectra;
	/// Number of bands 
	int num_bands;
	/// Original reference spectra, as input into the initializator 
	float **orig_spectra;
	/// Spectra that are updated with new information as the image is segmented as skin 
	float **updated_spectra;
	/// Number of samples used in updated_spectra 
	long *num_samples_in_spectra;
	/// Threshold values for SAM 
	float *sam_thresh;
	/// Start band for SAM calculations 
	int start_band_ind;
	/// End band for SAM calculations 
	int end_band_ind;
} masking_t;

/**
 * Masking error values.
 **/
enum masking_err_t{
	MASKING_NO_ERR = 0, 
	MASKING_LIBRARY_ERR = -1
};

/**
 * Internal datatype for controlling segmentations. Samples along the first table direction, ref. spectrum number along the next.
 **/
typedef bool** mask_thresh_t;

/** 
 * Allocate mask_thresh_t object.
 * \param mask_param Masking parameters
 * \param num_samples Number of samples in image
 **/
mask_thresh_t masking_allocate_thresh(const masking_t *mask_param, int num_samples);

/**
 * Free mask_thresh_t object. Everything will be freed, no need to do anything from outside. 
 **/
void masking_free_thresh(mask_thresh_t *mask_thresh_t, int num_samples);

/**
 * Check whether specified pixel belongs according to the thresholded SAM values. 
 * \param mask_param Masking parameters
 * \param threshed Thresholded values obtained from masking_thresh()
 * \param sample Sample coordinate along image
 * \return true if pixel belongs to the segmented image
 **/
bool masking_pixel_belongs(const masking_t *mask_param, mask_thresh_t threshed, int sample); 

/**
 * Data type for initialization. Chooses which directory to read reference spectra from. 
 **/
enum masking_input_data_type_t{REFLECTANCE_MASKING, TRANSMITTANCE_MASKING};

/**
 * Initialize masking parameters. 
 * \param num_wlens Number of bands in image to segment
 * \param wlens Wavelengths
 * \param masking_type Whether transmittance or reflectance masking
 * \param mask_param Output masking parameters
 **/
masking_err_t masking_init(int num_wlens, float *wlens, masking_input_data_type_t masking_type, masking_t *mask_param);

/** 
 * Do masking thresholding according to parameter specifications and update reference spectra according to segmented parts. 
 * \param mask_param Masking parameters
 * \param num_samples Number of samples in image
 * \param line_data Input hyperspectral data
 * \param ret_thresh Return segmented values.
 **/
void masking_thresh(masking_t *mask_param, int num_samples, float *line_data, mask_thresh_t *ret_thresh);

/**
 * Free memory associated with masking parameters. 
 **/
void masking_free(masking_t *mask_param);

#warning "NB! Fixed masking spectra directories have been defined. Please check (line approx. 94, src/masking.h)."
#define REFLECTANCE_MASKING_SPECTRA_DIRECTORY "/home/hyspex/IACOBUS/processing_libraries/masking/reflectance_spectra/"
#define TRANSMITTANCE_MASKING_SPECTRA_DIRECTORY "/home/hyspex/IACOBUS/processing_libraries/masking/transmittance_spectra/"


#endif
