libspectralmask
===============

Skin segmentation from hyperspectral images based on SAM (Spectral Angle
Mapper).  The segmentation is based on reference spectra located in
REFLECTANCE_MASKING_SPECTRA_DIRECTORY or
TRANSMITTANCE_MASKING_SPECTRA_DIRECTORY, hard-coded in src/masking.h. These datafiles
should be in plain-text ASCII format with wavelengths in nanometers along the
first column and reflectance or transmittance values along the second column. 

If you find this software useful, we would be grateful if you could cite the following paper: 

A. Bjorgan and L. L. Randeberg, "Towards real-time medical diagnostics using hyperspectral imaging technology", Proc. SPIE 9537 (in press). 
