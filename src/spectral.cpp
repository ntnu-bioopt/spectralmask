//==============================================================================
// Copyright 2015 Asgeir Bjorgan, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//==============================================================================

#include <fstream>
#include "spectral.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstring>
#include <sstream>
using namespace std;


spectral_err_t spectral_read_file(const char *filename, spectrum_t *spectrum){
	spectrum->num_values = 0;
	spectrum->values = NULL;
	spectrum->start_wlen = 0;
	spectrum->step_wlen = 0;

	//read from file
	ifstream file;
	file.open(filename);

	if (file.fail() || (!file.good())){
		return SPECTRAL_FILE_NOT_FOUND;
	}

	string line;
	vector<float> wlens;
	vector<float> vals;
	bool first_line = true;
	while (true){
		getline(file, line);
		if ((file.fail() || (!file.good())) && first_line){
			return SPECTRAL_NOT_VALID;
		}
		istringstream stream(line);

		float val1, val2;
		stream >> val1;
		stream >> val2;
		if (file.eof()){
			break;
		}

		if (val1 > 0){
			wlens.push_back(val1);
			vals.push_back(val2);
		}

		first_line = false;
	}
	file.close();

	if (wlens.size() == 0){
		return SPECTRAL_FILE_NOT_FOUND;
	}

	if (wlens.size() != vals.size()){
		return SPECTRAL_NOT_VALID;
	}

	//find minimum step wavelength
	float min_step = wlens[0];
	float min_wlen = wlens[0];
	float max_wlen = wlens[wlens.size()-1];
	for (int i=1; i < wlens.size(); i++){
		float step = wlens[i] - wlens[i-1];
		if ((step < min_step) && (step > 0)){
			min_step = step;
		}
	}
	float step_wlen = min_step;

	//create value array at regular wavelength steps where values falling between two wavelengths are interpolated linearly
	float wlen = min_wlen;
	int wlen_upper_ind = 0;
	vector<float> values_interpolated;
	while (wlen < max_wlen){
		if (wlen >= wlens[wlen_upper_ind]){
			wlen_upper_ind++;
		}
		
		int lower_ind = wlen_upper_ind-1;
		int upper_ind = wlen_upper_ind;

		float lower_wlen = wlens[lower_ind];
		float upper_wlen = wlens[upper_ind];

		float val = vals[lower_ind] + (wlen - lower_wlen)/(upper_wlen - lower_wlen)*(vals[upper_ind] - vals[lower_ind]);
		values_interpolated.push_back(val);

		wlen += step_wlen;
	}

	if (values_interpolated.size() == 0){
		return SPECTRAL_FILE_NOT_FOUND;
	}

	spectrum->values = new float[values_interpolated.size()]();
	for (int i=0; i < values_interpolated.size(); i++){
		spectrum->values[i] = values_interpolated[i];
	}

	spectrum->start_wlen = min_wlen;
	spectrum->step_wlen = step_wlen;
	spectrum->num_values = values_interpolated.size();

	return SPECTRAL_NO_ERR;
}

void spectral_copy(spectrum_t *destination, const spectrum_t *source){
	destination->num_values = source->num_values;
	destination->start_wlen = source->start_wlen;
	destination->step_wlen = source->step_wlen;
	destination->values = new float[destination->num_values];
	memcpy(destination->values, source->values, sizeof(float)*destination->num_values);
}

spectral_err_t spectral_free(spectrum_t *spectrum){
	delete [] spectrum->values;
	return SPECTRAL_NO_ERR;
}

spectral_err_t spectral_get_value(const spectrum_t *spectrum, float wlen, float *ret_value){
	if (spectrum->num_values == 0){
		return SPECTRAL_NOT_VALID;
	}

	//find lower and upper indices in value array
	int lower_ind = floor((wlen - spectrum->start_wlen)/spectrum->step_wlen);
	int upper_ind = lower_ind+1;
	if (lower_ind < 0){
		*ret_value = spectrum->values[0];
		return SPECTRAL_NO_ERR;
	} else if ((lower_ind >= spectrum->num_values) || (upper_ind >= spectrum->num_values)){
		*ret_value = spectrum->values[spectrum->num_values - 1];
	} else {
		//interpolate linearly
		float lower_wlen = spectrum->start_wlen + spectrum->step_wlen*lower_ind;
		float upper_wlen = spectrum->start_wlen + spectrum->step_wlen*upper_ind;
		*ret_value = spectrum->values[lower_ind] + (wlen - lower_wlen)/(upper_wlen - lower_wlen)*(spectrum->values[upper_ind] - spectrum->values[lower_ind]);
	}
	return SPECTRAL_NO_ERR;
}

spectral_err_t spectral_get_values_array(const spectrum_t *spec, float start_wlen, float step_wlen, int num_wlens, float *res){
	for (int i=0; i < num_wlens; i++){
		spectral_get_value(spec, start_wlen + i*step_wlen, &(res[i]));
	}
	return SPECTRAL_NO_ERR;
}


#ifdef _WIN32
#include <windows.h>
#else 
#include <sys/types.h>
#include <dirent.h>
#endif

/** 
 * Get list of files from a directory. 
 *
 * \param directory Specified directory
 * \param filenames Output filenames
 * \return Error value, SPECTRAL_NO_ERR on success
 **/
spectral_err_t spectral_get_files_in_directory(const char *directory, vector<string> *filenames){
	filenames->clear();

	//loop through directory 
	#ifdef _WIN32
	string win_dir = string(directory) + "/*";
	WIN32_FIND_DATA f;
	HANDLE h = FindFirstFile(win_dir, &f);
	if(h != INVALID_HANDLE_VALUE){
		do {
			filenames->push_back(f.cFileName);
		} while(FindNextFile(h, &f));
	} else {
		return SPECTRAL_DIRECTORY_NOT_FOUND;
	}

	#else
	DIR *dir = opendir(directory);
	if(dir){
		struct dirent *ent;
		while((ent = readdir(dir)) != NULL){
			filenames->push_back(ent->d_name);
		}
	} else {
		return SPECTRAL_DIRECTORY_NOT_FOUND;
	}
	closedir(dir);
	#endif

	return SPECTRAL_NO_ERR;
}

spectral_err_t spectral_construct_library_from_directory(const char *directory, spectral_library_t *library){
	//get list of filenames from directory
	vector<string> filenames;
	spectral_err_t retval = spectral_get_files_in_directory(directory, &filenames);
	if (retval != SPECTRAL_NO_ERR){
		return retval;
	}

	//convert to char* array
	char** c_filenames = new char*[filenames.size()];
	for (int i=0; i < filenames.size(); i++){
		string abspath = string(directory) + filenames[i];
		c_filenames[i] = new char[abspath.length()+1];
		memcpy(c_filenames[i], abspath.c_str(), sizeof(char)*(abspath.length()+1));
	}

	//create spectral library
	retval = spectral_construct_library_from_files(filenames.size(), c_filenames, library);

	for (int i=0; i < filenames.size(); i++){
		delete [] c_filenames[i];
	}
	delete [] c_filenames;
	return retval;
}

spectral_err_t spectral_construct_library_from_files(int num_files, char **filenames, spectral_library_t *library){
	vector<spectrum_t*> spectra;
	int num_valid_files = 0;
	for (int i=0; i < num_files; i++){
		spectrum_t *spectrum = new spectrum_t;
		spectral_err_t retval = spectral_read_file(filenames[i], spectrum);
		if (retval == SPECTRAL_NO_ERR){
			spectra.push_back(spectrum);
			num_valid_files++;
		} else {
			spectral_free(spectrum);
			delete spectrum;
		}
	}

	if (num_valid_files == 0){
		library->spectra = NULL;
		library->num_spectra = 0;
		return SPECTRAL_DIRECTORY_FILE_ERROR;
	}

	library->spectra = new spectrum_t[num_valid_files];
	for (int i=0; i < num_valid_files; i++){
		spectral_copy(&(library->spectra[i]), spectra[i]);
		spectral_free(spectra[i]);
		delete spectra[i];
	}
	library->num_spectra = num_valid_files;

	return SPECTRAL_NO_ERR;
}

void spectral_free_library(spectral_library_t *library){
	for (int i=0; i < library->num_spectra; i++){
		spectral_free(&(library->spectra[i]));
	}
	delete [] library->spectra;
}
