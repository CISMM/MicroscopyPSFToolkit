Microscopy PSF Toolkit

What it is
----------

This project provides analytical image models of microscope
point-spread functions. It is set up as an ITK external module and
extends the capabilities of ITK.

Who's responsible for this?
---------------------------

Cory Quammen
cquammen@cs.unc.edu
Department of Computer Science
UNC Chapel Hill

Funding provided by:

ARRA - Adding Deconvolution Algorithms to ITK
Award No: A10-1515-002
Sponsor Award No: HHSN276201000581P
Awarding Agency: National Library of Medicine

Acknowledgements
----------------

The Microscopy PSF Toolkit uses the open-source COSMOS code from
Dr. Crysanthe Preza's Computational Imaging Research Laboratory at The
University of Memphis (http://cirl.memphis.edu/cosmos.php).

License
------- 

The Microscopy PSF Toolkit is release under the terms of the GNU
Public License 2.

Build Instructions
------------------

The Microscopy PSF Toolkit is an external module for ITK. To use it, follow
these steps.

1. Clone ITK

 git clone git::/itk.org/ITK.git

2. Clone the Microscopy PSF Toolkit

 cd ITK/Modules/External
 git clone git://github.com/cquammen/MicroscopyPSFToolkit.git

3. Configure ITK

 cd ../../..
 mkdir ITK-build
 cd ITK-build
 ccmake ../ITK

The classes itkHaeberleCOSMOSPointSpreadFunctionImageSource and
itkGibsonLanniCOSMOSPointSpreadFunctionImageSource will now be
available in your build of ITK.
 
