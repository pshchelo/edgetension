Converting Pore Tension calc to Python
======================================

The initial procedure was written by Thomas Portet while his stay in 
Max-Planck-Institute of colloids and inerfaces.

Initial procedure uses a mix of ImageJ, Matlab (and MS Excel?) 
to do the job. While the first one being FLOSS is OK, the other two, 
specifically Matlab, are expensive proprietary tools, and the idea is 
to implement at least the Matlab part in free and cross-platform Python 
plus some packages.

The main complication is that most comments and names of functions and 
variables are in French or use abbreviations presumably stemming from French.

Matlab files
------------
Here is my initial guess of Matlab files present:

- affiche.m - displays sets of R**2*ln(r) lines for user to visually determine 
  the linear regime boundaries (english: display)
- chargement.m - loads data from txt file user creates from MS Excel file
  aslo stores names of corresponding image files (english: load)
- fit_lineaire.m - make linear fit of data (self-explanatory)
- pentes.m - makes a lot of linear fits and extracts tension value from them 
  (english: slopes)
- trous.m - performs image analysis to find pore radius (english: holes)

Pore finding procedure
----------------------
This is the insight on workings of algorythm in trous.m in respect to a single 
frame of single image file

The image is supposed to be rotated so that the horizontal midsection 
always goes through the pore with the pore located on the right side 
(if the pore is present that is).

# Find center of the image
# Blacken the left half of the image - that leaves either two separated 
  arc-like clusters or a single big hemi-circular one.
# Find indices (i.e. coordinates) of all non-zero elements
# for all nonzero elements find an angle between the element, center of image 
  and horizontal right (+) direction
# for nonzero elements in upper-right quadrant take element and its position 
  with the minimal angle
# for nonzero elements in lower-right quadrant take element and its position 
  with the maximal angle
# find distance between these two points. if it is (one or zero?) - 
  there is no pore, otherwise it is a pore diameter