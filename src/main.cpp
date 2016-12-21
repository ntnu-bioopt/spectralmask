//==============================================================================
// Copyright 2015 Asgeir Bjorgan, Norwegian University of Science and Technology
// Distributed under the MIT License.
// (See accompanying file LICENSE or copy at
// http://opensource.org/licenses/MIT)
//==============================================================================

#include "readimage.h"
#include "masking.h"
#include "spectral.h"
#include <iostream>
#include <sys/time.h>
using namespace std;

int main(int argc, char *argv[]){
	if (argc < 2) {
		fprintf(stderr, "Usage: %s hyperspectral_filename.\n", argv[0]);
		exit(1);
	}

	char* filename = argv[1];

	//read hyperspectral image header
	HyspexHeader header;
	hyperspectral_read_header(filename, &header);
	int start_band = 0;
	int end_band = header.bands;
	float *wlens = new float[end_band - start_band];
	for (int i=start_band; i < end_band; i++){
		wlens[i - start_band] = header.wlens[i];
	}

	masking_t mask_param;
	masking_err_t errcode = masking_init(end_band - start_band, wlens, REFLECTANCE_MASKING, &mask_param);
	if (errcode != MASKING_NO_ERR) {
		fprintf(stderr, "Error in initializing masking parameters: %s\n", masking_error_message(errcode));
		exit(1);
	}

	mask_thresh_t thresh_val = masking_allocate_thresh(&mask_param, header.samples);

	//read image and mask
	for (int i=0; i < header.lines; i++){
		ImageSubset subset;
		subset.startSamp = 0;
		subset.endSamp = header.samples;
		subset.startLine = i;
		subset.endLine = i+1;
		subset.startBand = 0;
		subset.endBand = header.bands;
	
		//read line
		float *line = new float[header.samples*(end_band - start_band)];
		hyperspectral_read_image(filename, &header, subset, line);
		
		//mask line
		masking_thresh(&mask_param, header.samples, line, &thresh_val);
		for (int i=0; i < header.samples; i++){
			cout << masking_pixel_belongs(&mask_param, thresh_val, i) << " ";
		}
		
		cout << endl;
		delete [] line;
	}
	masking_free_thresh(&thresh_val, header.samples);
	masking_free(&mask_param);
	delete [] wlens;
}	
