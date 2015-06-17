//==============================================================================
// Copyright 2015 Asgeir Bjorgan, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//==============================================================================

#ifndef SPECTRAL_H_DEFINED
#define SPECTRAL_H_DEFINED


/**
 * Defines behavior for each chromophore/spectrum. 
 **/
typedef struct{
	/* number of values in value array */
	int num_values;
	/* wavelength corresponding to first index */
	float start_wlen; 
	/* wavelength increment throughout the value array */
	float step_wlen; 
	/* contained values */
	float *values; 
} spectrum_t;

/**
 * Returned error values from calls to spectral_ functions.
 **/
enum spectral_err_t{
	SPECTRAL_NO_ERR = 0, 
	SPECTRAL_FILE_NOT_FOUND = -1,
	SPECTRAL_NOT_VALID = -2,
	SPECTRAL_DIRECTORY_NOT_FOUND = -3,
	SPECTRAL_DIRECTORY_FILE_ERROR = -4
};

/**
 * Read spectral information from file. 
 *
 * \param filename Input filename
 * \param spectrum Output spectrum
 * \return Error value, SPECTRAL_NO_ERR on success
 **/
spectral_err_t spectral_read_file(const char *filename, spectrum_t *spectrum);

/**
 * Free internal data associated with spectrum object. Will still have to deallocate spectrum outside if it was allocated on the heap. 
 *
 * \param spectrum Spectrum to free.
 * \return Error value, SPECTRAL_NO_ERR on success
 **/
spectral_err_t spectral_free(spectrum_t *spectrum);

/**
 * Copy data contained in one spectrum object to another.
 *
 * \param destination Destination for data copy
 * \param source Source for data copy
 **/
void spectral_copy(spectrum_t *destination, const spectrum_t *source);

/**
 * Get value of spectrum at specified wavelength, using linear interpolation. O(1)
 *
 * \param spectrum Spectral data
 * \param wlen Wavelength
 * \param ret_val Return value
 * \return Error value, SPECTRAL_NO_ERR on success
 **/
spectral_err_t spectral_get_value(const spectrum_t *spectrum, float wlen, float *ret_val);

/**
 * Get array of values at specified wavelengths, using linear interpolation. 
 *
 * \param spec Spectral data
 * \param start_wlen Start wavelength for data array
 * \param step_wlen Wavelength step difference
 * \param num_wlens Number of wavelengths (i.e. size of array)
 * \param res Array of size num_wlens, in which data is returned
 * \return Error value, SPECTRAL_NO_ERR on success
 **/
spectral_err_t spectral_get_values_array(const spectrum_t *spec, float start_wlen, float step_wlen, int num_wlens, float *res);

/**
 * Defines the library of all available chromophores/spectra and values.
 **/
typedef struct{
	/* Number of spectra */
	int num_spectra;
	/* Array over the available spectra */
	spectrum_t *spectra;
} spectral_library_t;

/**
 * Construct spectral chromophore library from files contained in specified directory.
 *
 * \param directory Input directory. Should end with the directory separator unique to the operating system in use (e.g. '/' for Leegnux and '\' for Windows (probably) 
 * \param library Output spectral library 
 * \return Error value, SPECTRAL_NO_ERR on success
 **/
spectral_err_t spectral_construct_library_from_directory(const char *directory, spectral_library_t *library);

/**
 * Construct spectral chromophore library from list of filenames. 
 *
 * \param num_files Number of filenames
 * \param filenames Array of filename strings
 * \param library Output spectral library
 * \return Error value, SPECTRAL_NO_ERR on sucess
 **/
spectral_err_t spectral_construct_library_from_files(int num_files, char **filenames, spectral_library_t *library);

/** 
 * Free internal memory associated with spectral library. 
 * \param library Library data to free
 **/
void spectral_free_library(spectral_library_t *library);


#endif
