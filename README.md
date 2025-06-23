# BlockMetrology [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15430176.svg)](https://doi.org/10.5281/zenodo.15430176)

This repository is part of the work 
"Enabling data-driven design of block copolymer self-assembly"
by Chiara Magosso, Irdi Murataj, Michele Perego, Gabriele Seguini, Debra J. Audus, Gianluca Milano,
and Federico Ferrarese Lupi, Scientific Data (2025), DOI: 10.1038/s41597-025-05379-w.

The repository contains 3 files:
- The file named Metadata_GUI.ipynb can be run on Google Colab and allows you to insert the process parameters in the Scanning Electron Microscope (SEM) images previously uploaded to a Google Drive folder.
- The file named sem_image_automated_preprocessing_for_ADAblok.py allows you to automatically analyze SEM images in order to extract morphological parameters. The file can be run locally and also requires the installation of ImageJ and modified ADAblock as per the following instructions.
- The file named sem_image_metadata_reader.py can be run locally and is used to generate a txt file with a copy of the metadata present in the SEM images. 

The ImageJ The software can be downloaded at https://imagej.net/ - We recommend version ImageJ 1.51w
The modified script ADAblock.ijm can be found at https://github.com/ChiaraMagosso/ADAblock/blob/ijMacro_updates_only/ADAblock.ijm
ADAblock was originally published by Murphy, J. N., Harris, K. D. & Buriak, J. M. Automated Defect and Correlation Length Analysis of Block Copolymer Thin Film Nanopatterns. PLOS ONE 10, e0133088 (2015) - https://doi.org/10.1371/journal.pone.0133088. 
To set up ADAblock, follow "Set up ImageJ" in "S1 Instructions - Use of ADAblock" of the above article. To understand the outputs of ADAblock, refer to "Sorting Through The Data - List of output files" of "S1 Instructions - Use of ADAblock" of the above article.

Before running the python files locally make sure you have installed the list of libraries with version reported in the requirements.txt file.
To run the files correctly, you need to change the path where indicated (commented: "#modify accordingly") according to where it is installed on your computer.
