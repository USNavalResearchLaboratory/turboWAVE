Input File: PIC
===============

Field Solvers
-------------

Field solver modules have priority 3 in the update sequence.

.. py:function:: new electrostatic field solver { directives }

	Installs a general electrostatic field solver

	:param block directives: The following directives are supported:

		.. py:function:: elliptical solver = slv

		 	:param enum slv: can be ``facr``, ``iterative``, or ``eigenmode``

		.. py:function:: poisson boundary condition coord = ( bc1 , bc2 )

			:param enum coord: can be ``x``, ``y``, ``z``
			:param enum bc1: boundary condition on lower side, can be ``open``, ``dirichlet``, ``neumann``.
			:param enum bc2: boundary condition on lower side, can be ``open``, ``dirichlet``, ``neumann``.

			At present ``open`` only works for ``facr`` or ``eigenmode`` elliptical solvers.

		.. py:function:: external potential = ( val1 , val2 )

			In case of dirichlet boundary condition fix the lower and upper potential at the z boundaries at the given values.

			:param float val1: potential at the lower z boundary
			:param float val2: potential at the upper z boundary

.. _coulomb-solver:
.. py:function:: new coulomb electromagnetic module { directives }

 	Installs a Coulomb gauge electromagnetic field solver.  Similar to WAVE field solver, but assumes continuity of sources (no divergence cleaning structure).  Cartesian coordinates only.

	:param block directives: The following directives are supported:

		.. py:function:: elliptical solver = slv

			:param enum slv: can be ``facr``, ``iterative``, or ``eigenmode``

		.. py:function:: poisson boundary condition z = ( bc1 , bc2 )

			:param enum bc1: boundary condition on lower side, can be ``open``, ``dirichlet``, ``neumann``.
			:param enum bc2: boundary condition on lower side, can be ``open``, ``dirichlet``, ``neumann``.

			At present ``open`` only works for ``facr`` or ``eigenmode`` elliptical solvers.

		.. py:function::	dipole center = (x,y,z)

			Reference point for dipole moment diagnostic

		.. py:function::	gamma beam = g

		 	Reference gamma for initializing beam potentials


.. py:function:: new curvilinear coulomb module { directives }

	Create a WAVE type EM solver that works in arbitrary coordinates.  Directives are the same as :ref:`coulomb electromagnetic module <coulomb-solver>`.  Elliptical solver cannot be :samp:`facr`.

.. _direct-solver:
.. py:function:: new direct electromagnetic module { directives }

	Create an EM module that advances Maxwell's curl equations directly, relying on continuity of sources to preserve divergence conditions.  The elliptical solver is only used for initialization.  Cartesian only.

	:param block directives: The following directives are supported:

		.. py:function:: elliptical solver = slv

			:param enum slv: can be ``facr``, ``iterative``, or ``eigenmode``

		.. py:function::	dipole center = (x,y,z)

			Reference point for dipole moment diagnostic

		.. py:function:: layer thickness = L

			:param int L: number of cells (in a single strip) occupied by absorbing layers.  If moving window is in use, layers are not added to the z boundaries.

		.. py:function:: layers = ( x0,x1,y0,y1,z0,z1 )

			Allows for control of layers at each individual boundary wall.

			:param int x0: number of cells in a single strip occupied by absorbing layers adjacent to the lower boundary in the x direction.  If 0 there are no PML media at this boundary.  Other 5 parameters are analogous.

		.. py:function:: reflection coefficient = R

		 	:param float R: Desired fraction of AMPLITUDE reflected.  If actual reflection is larger than requested, try increasing the number of layers.

.. py:function:: new curvilinear direct module { directives }

	Same as :ref:`direct electromagnetic module <direct-solver>` except for arbitrary coordinate system. Elliptical solver cannot be ``facr``.

.. py:function:: new pgc laser module { directives }

	Create an enveloped field solver suitable for use with ponderomotive guiding center simulations.

	:param block directives: The following directives are supported:

		.. py:function:: carrier frequency = f

			:param float f: base frequency ratio for the laser radiation

		.. py:function::	polarization = p

			:param enum p: can be ``linear``, ``circular``, or ``radial``

		.. py:function:: propagator = prop

			:param enum prop: can be ``spectral`` or ``adi``.

		.. py:function:: file = fname

			:param str fname: Name of file to use for reading initial profile.  The file is expected to be in DataViewer binary format (.dvdat file).  The scalar data is :math:`a^2`.

		.. py:function:: byte reversal

			If directive is present, bytes are reversed when reading binary profile data


Particle Species
----------------

Particle species can be used in electromagnetic PIC or as Bohmian particles in :doc:`bak-quantum`.

.. py:function:: new species name { directives }

	:param str name: name given to the species
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`boundaries <boundaries>`, :ref:`ionization`

		.. py:function:: mass = m0

			:param float m0: mass of the particle, default = 1.0

		.. py:function:: charge = q0

			:param float q0: charge of the particle, default = -1.0

		.. py:function:: particles per cell = ( Nx , Ny , Nz ) when density = n0

			Lays out particles on a subgrid of dimension :math:`N_x \times N_y \times N_z` within a cell.  The particles are weighted so that the density in the cell is ``n0``.  If particle weights are variable, the density specification is ignored (but still required), and the requested profile density is achieved in every cell.

		.. py:function:: minimum density = nmin

			:param float nmin: suppress creation of particles with density less than this

		.. py:function:: emission temperature = ( Tx , Ty , Tz )

			Thermal momentum of particles re-emitted from the boundaries

		.. py:function:: mobile = tst

			:param bool tst: set to false to hold this species immobile (defaults to true)

		.. py:function:: accelerate to pz in dt

		 	:param float pz: desired momentum of particle after acceleration
			:param float dt: time over which to accelerate particle

		.. py:function:: radiation damping = tst

			:param bool tst: set to true to apply radiation damping to the particles (default = false)

Nonlinear Optics
----------------

Bound particles treated as anharmonic oscillators can be used in the electromagnetic PIC environment.

.. py:function:: new bound name { directives }

	:param str name: name given to the bound species
	:param block directives: The following directives are supported:

		.. py:function:: mass = m0

			:param float m0: mass of the particle, default = 1.0

		.. py:function:: charge = q0

			:param float q0: charge of the particle, default = -1.0

		.. py:function:: basis = ( u1,u2,u3,v1,v2,v3,w1,w2,w3)

			Defines the :math:`{\bf u}`, :math:`{\bf v}`, and :math:`{\bf w}` unit vectors which define the principal axes of the crystal.

		.. py:function:: resonance = ( w1 , w2 , w3 )

			:param float w1: resonant frequency along u axis
			:param float w2: resonant frequency along v axis
			:param float w3: resonant frequency along w axis

		.. py:function:: damping = ( d1 , d2 , d3 )

			:param float d1: damping frequency along u axis
			:param float d2: damping frequency along v axis
			:param float d3: damping frequency along w axis

		.. py:function:: strength = ( f1 , f2 , f3 )

			:param float f1: oscillator strength along u axis
			:param float f2: oscillator strength along v axis
			:param float f3: oscillator strength along w axis

		.. py:function:: a1 = ( a11 , a12 , a13 , a14 , a15 , a16)

			First row of the second order anharmonic tensor

		.. py:function:: a2 = ( a21 , a22 , a23 , a24 , a25 , a26)

			Second row of the second order anharmonic tensor

		.. py:function:: a3 = ( a31 , a32 , a33 , a34 , a35 , a36)

			Third row of the second order anharmonic tensor

		.. py:function:: b = b0

			:param float b0: cubic anharmonic coefficient

		.. py:function:: d = d0

			:param float d0: quintic anharmonic coefficient

		.. py:function:: phi = q1

			:param float q1: Rotation about z in degrees.  Initial orientation has principal axes aligned with standard basis.  This rotation happens before the theta rotation.

		.. py:function:: theta = q2

			:param float q2: Rotation about y in degrees.  Initial orientation has principal axes aligned with standard basis.  This rotation happens after the phi rotation.
