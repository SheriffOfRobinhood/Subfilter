=====================
The subfilter module.
=====================
This module implements overall filtering.

The main functions are:

* :py:func:`subfilter.subfilter.filter_variable_list`.
* :py:func:`subfilter.subfilter.quadratic_subfilter`.

Before use, two files have to be created, using:

* :py:func:`subfilter.subfilter.setup_derived_data_file`.
* :py:func:`subfilter.subfilter.setup_filtered_data_file`.

.. topic:: New at 0.3
    
    #. The :py:func:`subfilter.subfilter.filter_variable_pair_list` function outputs filtered pairs :math:`\phi,\psi` inder the name "s(:math:`\phi,\psi`)_on_grid" where "grid" will be "u", "v", "w" or "p".


Detailed Module Contents
------------------------
The entire module is documented below.

.. automodule:: subfilter.subfilter
   :member-order: bysource
   :members:
   