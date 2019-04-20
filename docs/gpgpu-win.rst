GPGPU on Windows 10
===================

.. caution::

	We assume core turboWAVE **and** tools have already been installed according to the documentation, with no steps omitted.

Driver
------

Update to the latest graphics drivers, following the guidance on the vendor's website.

Compile with OpenCL
-------------------

#. Edit :samp:`{twroot}\\core\\source\\win.make`
#. Change ``win.make`` to support OpenCL compilation (TBD)
#. From the developer/intel prompt navigate to :samp:`{twroot}`:samp:`\\core\\source`
#. :samp:`nmake /F win.make`
