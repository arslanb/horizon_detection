This is the MATLAB software package for horizon detection on a single image. 

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
