//==============================================================================
// Copyright 2015 Asgeir Bjorgan, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//==============================================================================

#include "masking.h"
#include "spectral.h"
#include <cmath>
#include <iostream>
using namespace std;

#define SAM_THRESH_DEFAULT 0.3
#define SAM_THRESH_TRANSMITTANCE 0.10
masking_err_t masking_init(int num_wlens, float *wlens, masking_input_data_type_t masking_type, masking_t *mask_param){
	spectral_library_t library;
	
	//read from specified library directory according to masking type
	float sam_thresh = SAM_THRESH_DEFAULT;
	spectral_err_t retval = SPECTRAL_NO_ERR;
	switch (masking_type){
		case REFLECTANCE_MASKING:
			retval = spectral_construct_library_from_directory(REFLECTANCE_MASKING_SPECTRA_DIRECTORY, &library);
			if (retval != SPECTRAL_NO_ERR) {
				return MASKING_REFLECTANCE_LIBRARY_ERR;
			}
		break;

		case TRANSMITTANCE_MASKING:
			retval = spectral_construct_library_from_directory(TRANSMITTANCE_MASKING_SPECTRA_DIRECTORY, &library);
			if (retval != SPECTRAL_NO_ERR) {
				return MASKING_TRANSMITTANCE_LIBRARY_ERR;
			}
			sam_thresh = SAM_THRESH_TRANSMITTANCE;
		break;
	}

	//generate masking spectra from the spectral library
	mask_param->num_masking_spectra = library.num_spectra;
	mask_param->num_bands = num_wlens;
	mask_param->orig_spectra = new float*[library.num_spectra];
	mask_param->updated_spectra = new float*[library.num_spectra];
	mask_param->sam_thresh = new float[library.num_spectra]();
	mask_param->start_band_ind = 0;
	mask_param->end_band_ind = num_wlens - 1;
	mask_param->num_samples_in_spectra = new long[library.num_spectra]();

	for (int i=0; i < mask_param->num_masking_spectra; i++){
		mask_param->orig_spectra[i] = new float[num_wlens]();
		mask_param->updated_spectra[i] = new float[num_wlens]();
		for (int j=0; j < num_wlens; j++){
			spectral_get_value(&(library.spectra[i]), wlens[j], &(mask_param->orig_spectra[i][j]));
			spectral_get_value(&(library.spectra[i]), wlens[j], &(mask_param->updated_spectra[i][j]));
		}
		mask_param->sam_thresh[i] = sam_thresh;
	}

	spectral_free_library(&library);
	return MASKING_NO_ERR;
}

void masking_free(masking_t *mask_param){
	for (int i=0; i < mask_param->num_masking_spectra; i++){
		delete [] mask_param->orig_spectra[i];
		delete [] mask_param->updated_spectra[i];
	}
	delete [] mask_param->orig_spectra;
	delete [] mask_param->updated_spectra;
	delete [] mask_param->sam_thresh;
	delete [] mask_param->num_samples_in_spectra;
}

void masking_thresh(masking_t *mask_param, int num_samples, float *line_data, mask_thresh_t *ret_thresh){
	//calculate norms of reference spectra
	float *ref_norms_orig = new float[mask_param->num_masking_spectra]();
	float *ref_norms_updated = new float[mask_param->num_masking_spectra]();
	for (int i=0; i < mask_param->num_masking_spectra; i++){
		for (int j=mask_param->start_band_ind; j <= mask_param->end_band_ind; j++){
			ref_norms_orig[i] += mask_param->orig_spectra[i][j]*mask_param->orig_spectra[i][j];
			ref_norms_updated[i] += mask_param->updated_spectra[i][j]*mask_param->updated_spectra[i][j];
		}
		ref_norms_orig[i] = sqrt(ref_norms_orig[i]);
		ref_norms_updated[i] = sqrt(ref_norms_updated[i]);
	}
	
	for (int j=0; j < num_samples; j++){
		//get pixel band values, calculate norm of pixel spectrum
		float pixel_norm = 0;
		float *pixel_vals = new float[mask_param->num_bands];
		for (int i=mask_param->start_band_ind; i <= mask_param->end_band_ind; i++){
			pixel_vals[i] = line_data[i*num_samples + j];
			pixel_norm += pixel_vals[i]*pixel_vals[i];
		}
		pixel_norm = sqrt(pixel_norm);

		//calculate sam values against all available spectra
		for (int k=0; k < mask_param->num_masking_spectra; k++){
			float samval_orig = 0;
			float samval_updated = 0;
			for (int i=mask_param->start_band_ind; i <= mask_param->end_band_ind; i++){
				samval_orig += pixel_vals[i]*mask_param->orig_spectra[k][i];
				samval_updated += pixel_vals[i]*mask_param->updated_spectra[k][i];
			}
			samval_orig /= pixel_norm*ref_norms_orig[k];
			samval_orig = acos(samval_orig);
			samval_updated /= pixel_norm*ref_norms_updated[k];
			samval_updated = acos(samval_updated);

			//compare against thresholds, save to return array in separate slots
			bool pixel_belong = (samval_orig < mask_param->sam_thresh[k]) || (samval_updated < mask_param->sam_thresh[k]);
			(*ret_thresh)[j][k] = pixel_belong;

			//update the updated spectra with new information if above threshold
			if (pixel_belong){
				long n = mask_param->num_samples_in_spectra[k];
				n++;
				ref_norms_updated[k] = 0;
				for (int i=mask_param->start_band_ind; i <= mask_param->end_band_ind; i++){
					double delta = pixel_vals[i] - mask_param->updated_spectra[k][i];

					//update reference spectrum
					mask_param->updated_spectra[k][i] += delta/(n*1.0);

					//update norm of reference spectrum
					ref_norms_updated[k] += mask_param->updated_spectra[k][i]*mask_param->updated_spectra[k][i];
				}
				ref_norms_updated[k] = sqrt(ref_norms_updated[k]);
				mask_param->num_samples_in_spectra[k] = n;
			}
		}
		delete [] pixel_vals;
	}
	delete [] ref_norms_orig;
	delete [] ref_norms_updated;
}

mask_thresh_t masking_allocate_thresh(const masking_t *mask_param, int num_samples){
	bool **ret_val = new bool*[num_samples];
	for (int i=0; i < num_samples; i++){
		ret_val[i] = new bool[mask_param->num_masking_spectra];
	}
	return ret_val;
}


void masking_free_thresh(mask_thresh_t *mask_thresh, int num_samples){
	for (int i=0; i < num_samples; i++){
		delete [] (*mask_thresh)[i];
	}
	delete [] (*mask_thresh);
}

bool masking_pixel_belongs(const masking_t *mask_param, mask_thresh_t threshed, int sample){
	bool belongs = false;
	for (int i=0; i < mask_param->num_masking_spectra; i++){
		belongs = belongs || threshed[sample][i];
	}
	return belongs;
}

const char *masking_error_message(masking_err_t errcode)
{
	switch (errcode) {
		case MASKING_NO_ERR:
			return "No error.";
		case MASKING_REFLECTANCE_LIBRARY_ERR:
			return "Error in constructing spectral library from masking spectra in " REFLECTANCE_MASKING_SPECTRA_DIRECTORY;
		case MASKING_TRANSMITTANCE_LIBRARY_ERR:
			return "Error in constructing spectral library from masking spectra in " TRANSMITTANCE_MASKING_SPECTRA_DIRECTORY;
	}
}
