Pore Tension Calculation
========================

What's this about
-----------------

The initial procedure for calculating vesicle pore edge tension from electroporation 
images was written by Thomas Portet while his stay in Max-Planck-Institute of 
Colloids and Inerfaces between April and June 2009.

Initial procedure uses a mix of ImageJ and MATLAB to do the job. 
While the first one being FLOSS is acceptable, MATLAB is an expensive 
proprietary tool, and I had re-implemented the MATLAB part in 
free and cross-platform Python plus some packages.


Image acquisition and processing
--------------------------------

Image acquisition and storage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*This section is copy-pasted with minor editions from T.Portet's manual.*

Images must be transferred from the camera in tiff format. 
Acquisition must be performed with a linear timescale. Images relative to an 
experiment (i.e. to the application of one pulse on one vesicle) must be saved 
as an image stack in tiff format. The user must ensure that the cathode-facing 
hemisphere of the vesicle is situated on the right-hand side of the images. 
Converting images to a stack and rotating them if necessary can be performed with ImageJ.

Image processing with ImageJ for membrane detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*This section is copy-pasted with some editions from T.Portet's manual.*

You first have to open this stack with ImageJ, define a rectangular selection 
enclosing the vesicle, and crop the image. You have to make sure that for 
each slice of the stack, the center of the image lies inside the vesicle,
and that the membrane is not porated neither above nor below the center of the 
image.

As a first step in this membrane detection stage, perform an image background 
correction by using the background subtraction rolling-ball algorithm 
[Sternberg(1983)], with ball radius set to 50 pixels4 (Process, 
Subtract Background, Sliding Paraboloid, 50 pixels). A ball radius of 
50 pixels gave good results with processed images, other values can also 
be tried if problems. Then locate the membrane position by using a common Sobel 
edge detector to highlight sharp changes in intensity (Process, Find Edges). 
The last step is the image binarization (Process, Binary, Make Binary, Calculate
Threshold for each Image, Black Background) where the threshold value is 
calculated using the classical Isodata algorithm [Ridler and Calvard(1978)]. 
Of course all these three operations should be performed on the entire stack. 

There is a macro for performing automatically these three procedures 
in the file macros_processing.txt. For using this macro, first install it in ImageJ 
(Plugin, Macros, Install, choose file macros processing.txt), then simply use 
the shortcut p. (This macro named process image will also appear in the menu
Plugins, Macros.)

**Update:** The whole toolbar for poration image preparation is now available.
Just copy "Poration Toolset.txt" to "[ImageJ folder]/macros/toolsets/" folder, 
restart ImageJ, click on right-most ">>" button on the ImageJ toolbar and choose
"Poration Toolset". This will fill the right part of the toolbar with buttons 
and menues for most used tools and functions when preparing images of 
poration experiments to help you minimize deep ImageJ submenu navigation. 
It also includes button providing the same macros T.Portet has included before,
so you do not need to install it anymore. Since the icons on tools might 
be somewht cryptic, pay attention to the ImageJ tooltips for these buttons, 
they show what the buttons mean.

In some cases, we have to remove residual white pixels appearing inside the 
vesicle, because they could lead to errors in the pore radius measurements. 
Such pixels are problematic if they are located on the vertical line passing 
through the center of the image and if they are not connected to the membrane. 
You could just draw a white rectangle which fills the black part between these 
residuals pixels and the membrane. This solution is often easier to apply than 
the one consisting in removing the residual pixels, and can be easily applied 
to the whole stack.

User manual for Python implementation
-------------------------------------

Prerequisites (version tested with is in parentheses):

- Python (2.7.2)
- NumPy (1.6.1)
- SciPy (0.10)
- PIL Python Image Library (1.1.7)
- MatPlotLib (1.1.0)

Put pores.py and fitpore.py somewhere in the PATH, or adjust the PATH to include them.
Run pores.py -h and fitpore.py -h to see available command line options and explanations.

pores.py
~~~~~~~~

::

	usage: pores.py [-h] [--skip SKIP] [--test TEST] filename

	Extract pore positions and radii.

	positional arguments:
	  filename     Name of the multipage b/w tiff file

	optional arguments:
	  -h, --help   show this help message and exit
	  --skip SKIP  Analyse only every SKIPth frame (default is every frame)
	  --test TEST  Run the procedure TEST times and report the minimal of them

	Outputs TSV text file named imagename_skipXXX.txt, with 6 columns: pore radius
	in pixels, 4 columns for x and y coordinates of pore edges, corresponding
	frame number.
  
fitpore.py
~~~~~~~~~~

::

	usage: fitpore.py [-h] [-s S] [-e E] [-r R] [-f F] [-v V] filename

	Get pore tension in 2 steps.

	positional arguments:
	  filename    Name of the input file with pore radii and frame numbers

	optional arguments:
	  -h, --help  show this help message and exit
	  -s S        Start of the linear region (frame number)
	  -e E        End of the linear region (frame number)
	  -r R        Initial radius of the vesicle in microns.
	  -f F        Speed of image acquisition in frames per second
	  -v V        Viscosity of the bulk solution in Pa*s (defaults 1e-3 Pa*s for water)

	First run with only a filename as input and remember the boundaries of the
	linear stage. Then run again supplying all arguments to get the linear region
	fitted and pore tension displayed.


Technical details
-----------------

Pore finding procedure
~~~~~~~~~~~~~~~~~~~~~~

The implementation in Python very closely follows to the MATLAB one, except 
using a fast library for cluster detection instead of brute-force high-level 
code in MATLAB. Result is **33-fold increase in speed**, while the difference between 
two implementations for all but few frames of 5 test images supplied with 
MATLAB code is close to zero, and even in those few the difference is 
in the order of half pixel.

Also with Python implementation it is possible to count the number of frames 
in the TIFF file programmaticaly (although at the cost of some relatively short time), 
so this parameter is no longer needed.

Below is the insight on workings of algorithm in respect to a single frame 
of single image file. The image is supposed to be rotated so that 
the horizontal midsection always goes through the pore with the pore located 
on the right side (if the pore is present that is). 

#. Find center of the image
#. Blacken the left half of the image.
#. Find the innermost intersection points between vesicle and vertical midsection
   (these most likely are on those squares put onto images as described in the 
   manual for MATLAB code).
#. Find the continuous clusters those innermost points belong to. Continuous means 
   that every point of the cluster has at least one nearest neighbour 
   in any of 8 directions.
#. Find indices (i.e. coordinates) of all non-zero elements of the clusters found.
#. Find the (signed) angles between the center of image, 
   positive x-direction (right) and the nonzero points of the clusters.
#. For nonzero elements in upper-right quadrant take element and its position 
   with the minimal angle.
#. For nonzero elements in lower-right quadrant take element and its position 
   with the maximal angle.
#. Find distance between these two points, filtering out possible overlapping cases.

MATLAB files
~~~~~~~~~~~~

Here is my idea of what those Matlab files are specifically for:

- affiche.m - displays sets of R**2 * ln(r) lines for user to visually determine 
  the linear regime boundaries (english: display)
- chargement.m - loads data from txt file user creates from MS Excel file;
  aslo stores names of corresponding image files (english: load)
- fit_lineaire.m - make linear fit of data (self-explanatory)
- pentes.m - makes series of linear fits and extracts tension values from them 
  (english: slopes)
- trous.m - performs image analysis to find pore radius (english: holes)