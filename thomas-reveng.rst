Converting Pore Tension calc to Python
======================================

The initial procedure was written by Thomas Portet while his stay in 
Max-Planck-Institute of colloids and inerfaces.

Initial procedure uses complicated mix of ImageJ, Matlab (and MS Excel?) 
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
- chargement.m - loads data from txt file, stores names of corresponding 
  image files (english: load)
- fit_lineaire.m - make linear fit of data
- pentes.m - makes a lot of linear fits and extracts temsion value from them 
  (english: slopes)
- trous.m - performs image analysis to find pore radius (english: holes)

