========================
Pore Tension Calculation
========================
:Info: http:/sites.google.com/shchelokovskyy
:Author: Pavlo Shchelokovskyy <shchelokovskyy@gmail.com>
:Description: use Docutils or http://rst2a.com to convert to HTML, PDF or other formats.

What's this about
-----------------

The initial procedure for calculating vesicle pore edge tension from electroporation 
images was written by Thomas Portet while his stay in Max-Planck-Institute of 
Colloids and Inerfaces between April and June 2009. 
Some parts of this document are heavily based on the document "howto.pdf"
written by him as a documentation for his project.

Initial procedure uses a mix of ImageJ and MATLAB to do the job. 
While the first one being FLOSS is acceptable, MATLAB is an expensive 
proprietary tool, and I had re-implemented the MATLAB part in 
free and cross-platform Python plus some packages.


Image acquisition and processing
--------------------------------

Image acquisition and storage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Images must be transferred from the camera in tiff format. 
Acquisition must be performed with a linear timescale. Images related to an 
experiment (i.e. to the application of one pulse on one vesicle) must be saved 
as an image stack in tiff format. The user must ensure that the cathode-facing 
hemisphere (wher the pore appears) of the vesicle is situated on the right-hand side of the images. 
Converting images to a stack and rotating them if necessary can be performed with ImageJ.

Image processing with ImageJ for membrane detection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You first have to open this stack with ImageJ, define a rectangular selection 
enclosing the vesicle, and crop the image. You have to make sure that for 
each slice of the stack, the center of the image lies inside the vesicle,
and that the membrane is not porated neither above nor below the center of the 
image.

As a first step in this membrane detection stage, perform an image background 
correction by using the background subtraction rolling-ball algorithm 
[Sternberg1983]_, with ball radius set to 50 pixels (Process, 
Subtract Background, Sliding Paraboloid, 50 pixels). A ball radius of 
50 pixels gave good results with processed images, other values can also 
be tried if problems. Then locate the membrane position by using a common Sobel 
edge detector to highlight sharp changes in intensity (Process, Find Edges). 
The last step is the image binarization (Process, Binary, Make Binary, Calculate
Threshold for each Image, Black Background) where the threshold value is 
calculated using the classical Isodata algorithm [Ridler_Calvard1978]_. 
Of course all these three operations should be performed on the entire stack. 

The whole toolbar for poration image preparation is available to simplify 
the image processing with ImageJ.
Just copy "Poration Toolset.txt" to "[ImageJ folder]/macros/toolsets/" folder, 
restart ImageJ, click on right-most ">>" button on the ImageJ toolbar and choose
"Poration Toolset". This will fill the right part of the toolbar with buttons 
and menus for most used tools and functions when preparing images of 
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

.. [Sternberg1983] Sternberg, S., 1983. 
   Biomedical Image Processing. 
   *IEEE Comput.* **16**:22-34.

.. [Ridler_Calvard1978] Ridler, T. W., and S. Calvard, 1978. 
   Picture Thresholding Using an Iterative Selection Method. 
   *IEEE Trans. Syst. Man. Cybern.* **8**:630-632.

User manual for Python implementation
-------------------------------------

Prerequisites (version developed and tested with is given in parentheses):

- Python (2.7.2) - programming language interpreter and its standard library
- NumPy (1.6.1) - provides fast numerical arrays
- SciPy (0.10) - provides rich set of scientific functions based on NumPy
- PIL Python Image Library (1.1.7) - for image loading
- MatPlotLib (1.1.0) - for plotting
- wxPython (2.8.12.1) - for Graphical User Interface


Run pores.py -h to see available command line options and explanations.

pores.py
~~~~~~~~

::

    usage: pores.py [-h] [--file FILE] [--skip SKIP] [--test TEST] [--nogui]
    
    Extract pore positions and radii.
    
    optional arguments:
      -h, --help            show this help message and exit
      --file FILE, -f FILE  Name of the multipage b/w tiff file to process
      --skip SKIP, -s SKIP  Analyse only every SKIPth frame (default is every
                            frame)
      --test TEST, -t TEST  Run the procedure TEST times and report the minimal of
                            them
      --nogui, -n           Run in no GUI mode suitable for batch processing.
    
    with --nogui outputs TSV text file named imagename_skipXXX.txt, with 6
    columns: pore radius in pixels, 4 columns for x and y coordinates of pore
    edges, corresponding frame number.

When you launch pores.py, you will be presented with a window containing 
a toolbar, an empty plot with a plot toolbar and two horizontal sliders 
under it, a parameters panel on the right side and a status bar at the bottom. 
Also the a standard console window will open in background, showing errors and 
warnings when they happen.

Toolbar has the following buttons:

- Open image file - opens and immediately analyses the image with settings 
  found on parameters panel.
- Open pores data file - opens a text data file previousely saved by pores.py 
  or wxpores.py
- Save pores data file - saves a text data file for later reference/analysis, 
  coompatible with files saved by pores.py and old MATLAB procedure.
- Show found pores - opens an extra window where position of pore edges found
  is visualized overlayed on image being analysed.

Parameters panel can set following parameters:

- Take every - take every n-th frame for calculations. 1 means all frames 
  are analysed, 2 means every other frame and so on.
- Radius - radius of the vesicle after poration in micrometers. Currently clipped 
  at 200 micrometers, which is a quite unrealistically huge vesicle.
- Speed - speed of image acquisition in frames per second. Currently clipped 
  at 50000 fps, which is beyond the speed of the fast camera available in our lab.
- Viscosity - viscosity of bulk media in mPa*s = 0.1 Pois. Currently clipped at 
  2000 mPa*s, which is way above viscosity even of pure glycerol 
  (1.2 Pa*s at room temperature).
- Autozoom - if enabled the plot will automatically zoom to the region 
  defined by two sliders.

Status bar shows toolbar items hints and coordinates of the cursor when over the plot.

With the plot toolbar (below the plot) you can pan and zoom the plot, 
revert to original pan and zoom settings and save the image in variety of formats, 
both vector and raster.

To analyse an image prepared as described in section on image processing, 
open it by pressing "Open Image" button and choosing the image. The image will be
immediately analysed (no visual clues for now, so it might look like the program 
hangs), taking parameters as set on Parameters panel. You could also open 
a previousely saved data file. 

In both cases you will be presented with the plot of ln(Rp) vs time. 
First, adjust the parameters on the Parameters panel to the desired values. 
Than, using two sliders below the plot, define the linear region (it will be 
visualized as two vertical dashed lines on the plot). If the linear fitting of 
this region succeeds, the plot will also show the fitted line, and the plot title 
will be adjusted to display the calculated value of edge tensions (in picoNewtons) 
and its standard error (derived from the fit), the frame interval where the fit 
was performed and values of other material parameters used for fit. If the value 
of edge tension is displayed as "nan" (i.e. not a number), it means that the 
fitting has failed due to presence of pore radius zero somewhere in the defined 
region of fit. Such data also produce warnings in the background console, 
something like

::

    C:\pores.py:572: RuntimeWarning: divide by zero encountered in log
      lnr = np.log(self.data[0])
    C:\Python27\lib\site-packages\numpy\lib\function_base.py:1989: RuntimeWarning: invalid value encountered in subtract
      X -= X.mean(axis=1-axis)[tup]
    C:\Python27\lib\site-packages\scipy\stats\stats.py:2810: RuntimeWarning: invalid value encountered in absolute
      prob = distributions.t.sf(np.abs(t),df)*2


Just ignore this. However if there are other messages not of this type, 
it may be a bug. In this case contact me and I will try to investigate.

*Note:* in this program the first frame of the multi-page TIFF image is numbered 1, 
as done by ImageJ. This can be different from frame numbers by image acquisition 
software or other analysis tools, as they can assign number 0 to the first frame. 
Keep that in mind when doing frame-by-frame comparisons or searching for a particular 
frame with other tools.

If you are interested in how well pore detection algorithm had performed, press 
"Show found pores" button on the toolbar. If you have opened the experiment 
from the text file you will be prompted to open an image file corresponding to it. 
In any way, you will see a window where you can browse through the multi-page TIFF 
with the help of the slider, and the found pore will be shown as a line 
joining the edges of the pore established by the algorithm (or a single dot 
in the center of the image if no pore was found). The title of the image 
will show you the current frame number and found pore radius in pixels.


Technical details
-----------------

Pore finding procedure
~~~~~~~~~~~~~~~~~~~~~~

The implementation in Python very closely follows to the MATLAB one, except 
using a fast library for cluster detection instead of brute-force high-level 
code in MATLAB. Result is **33-fold increase in speed** as measured on several 
test images, while the difference between two implementations for all but few 
frames of 5 test images supplied with MATLAB code is close to zero, 
and even in those few the difference is in the order of half pixel.

Also with Python implementation it is possible to count the number of frames 
in the TIFF file programmaticaly (although at the cost of some relatively short time), 
so this parameter is no longer needed.

Below is the insight on workings of algorithm in respect to a single frame 
of single image file. The image is supposed to be rotated with the pore located 
on the right side (see section on image processing). 

#. Find center of the image
#. Blacken the left half of the image.
#. Find the innermost intersection points between vesicle and vertical midsection
   (these most likely are on those squares put onto images as described in the 
   section on image processing).
#. Find continuous clusters those innermost points belong to. Continuous means 
   that every point of the cluster has at least one nearest neighbour 
   in any of 8 directions.
#. Find the (signed) angles between the center of image, 
   positive x-direction (right) and each point of the clusters.
#. For nonzero elements in upper-right quadrant take element and its position 
   with the minimal angle.
#. For nonzero elements in lower-right quadrant take element and its position 
   with the maximal angle.
#. Find distance between these two points, filtering out possible overlapping cases.

