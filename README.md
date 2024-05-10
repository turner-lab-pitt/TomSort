# TomSort
	About TomSort
	  TomSort is an automated spike sorter utilizing the super paramagnetic clustering (SPC) algorithm. 
 	We leveraged the SPC code implemented by the Quiroga lab at the University of Leicester (see the reference for further details). 
	In the Turner lab at the University of Pittsburgh, we use TomSort for spike sorting, then convert the data into the nex (NeuroExplorer, Plexon) file format. 
 	Subsequently, we curate the results using OfflineSorter (Plexon), and evaluate the quality of unit isolation. 
	This repository contains the codes for TomSort and conversion script.

	Data structure
	  The data file needs two variables: extracellular data names ad "signal" or "hp_cont" and the sampling rate (Hz) names as "samplerate"

	MatLab requirements
	  Parallel Computing Toolbox:	
	    https://jp.mathworks.com/products/parallel-computing.html	
	  extrema.m function:	
	    https://jp.mathworks.com/matlabcentral/fileexchange/12275-extrema-m-extrema2-m	
	  wave_clus (SPC):	
	    https://github.com/csn-le/wave_clus
	  Code to read and write NeuroExplorer data files:
	    https://www.neuroexplorer.com/downloadspage/

	Reference
	  Chaure, F. J. and Rey, H. G. and Quian Quiroga, R. "A novel and fully automatic spike sorting implementation with variable number of features", Journal of Neurophysiology, vol. 120-4, pg. 1859-1871, 2018

