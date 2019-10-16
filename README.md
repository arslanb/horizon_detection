This is the MATLAB software package for horizon detection on a single image paper [Xu et al. CVPR 2013](Xu_A_Minimum_Error_2013_CVPR_paper.pdf).
Use the following citation:

	Xu, Yiliang, Sangmin Oh, and Anthony Hoogs. 
	"A minimum error vanishing point detection approach for uncalibrated monocular images of man-made environments." 
	In Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition, pp. 1376-1383. 2013.

BibTex

	@InProceedings{Xu_2013_CVPR,
		author = {Xu, Yiliang and Oh, Sangmin and Hoogs, Anthony},
		title = {A Minimum Error Vanishing Point Detection Approach for Uncalibrated Monocular Images of Man-Made Environments},
		booktitle = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
		month = {June},
		year = {2013}
	}

Copyright 2014 by Kitware, Inc. All Rights Reserved. 
Please refer to LICENSE.TXT for licensing information, or contact General Counsel, Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065

*Software Structure*

horizon_detection.m contains the main function. 
	- extract_lsd.m: calls Line Segment Detector (LSD) C++ executable "LSD_VC.exe" to compute line segments in the input image. LSD_VC.exe read in an image and output a text file containing all detected line segments information. The Matlab script reads in this text file. This module can be substituted by other line segment detection algorithm.
	- vp_probability_EM.m: compute vanishing points given line segments that are detected.
	- horizon.m: compute horizon given vanishing points that are detected.


*Dependency*

	- Line Segment Detector (LSD) C++ executable "LSD_VC.exe" is included in this package. The source files can be found at http://www.ipol.im/pub/art/2012/gjmr-lsd/.
	- MATLAB basic package
	- MATLAB image processing toolbox (only required for visualizing results).

