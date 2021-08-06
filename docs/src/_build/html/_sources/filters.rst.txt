=====================
The filters module.
=====================
This module contains the code to generate a selection of 2-dimensional filters:

* Gaussian - specified by :math:`\sigma`, the standard deviation of the Gaussian.
* Spectral Cutoff - specified by the wavenumber. A rough equivalence with
  the Gaussian filter is :math:`wavenumber = \pi/(2\sigma)`. Hence :math:`wavelength=4\sigma`.
* Running mean - specified by the width in grid points. A rough equivalence with
  the Gaussian filter is 
  :math:`width = int(\sigma /dx \times \pi \times 2.0/3.0)+1`, where
  :math:`dx` is the grid spacing.

.. topic:: New at 0.3

    #. The filters.filter_2D class has been replaced with :py:class:`filters.Filter`. This now accepts an optional argument ndim when creating a Filter instance. This may be 1 or 2 and defaults to 2. The use_ave option is no longer supported.


Detailed Module Contents
------------------------
The entire module is documented below.

.. automodule:: subfilter.filters
   :member-order: bysource
   :members:
   :undoc-members:

