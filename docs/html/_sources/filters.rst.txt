=====================
The filters module.
=====================
This module contains the code to generate a selection of 2-dimensional filters:

* Gaussian - specified by :math:`\sigma`, the standard deviation of the Gaussian. Note that, in MONC, we have found that the Smagorinsky sub-filter model corresponds roughly with a Gaussian filter with :math:`\sigma = \Delta`. 
* Spectral Cutoff - specified by the wavenumber. A rough equivalence with
  the Gaussian filter is :math:`wavenumber = \pi/(2\sigma)`. Hence :math:`wavelength=4\sigma`.
* Spectral Cylindrical Cutoff - specified by the wavenumber. A rough equivalence with
  the Gaussian filter is :math:`wavenumber = \pi/(2\sigma)`. Hence :math:`wavelength=4\sigma`.
* 2D version of the 1-2-1 filter. Note: if ``options['FFT_type']`` is set to ``'DIRECT'``, this is calculated directly, not using FFTs.
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

