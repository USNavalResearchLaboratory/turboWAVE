Input File: PIC
===============

Smoothing
---------

Any field solver or particle species module can accept the following smoothing directives.  If applied to a field solver, sources are smoothed prior to advancing the fields, and the physics is affected.  If applied to a particle species, only the single species density diagnostic is affected (note that the sources may be smoothed, even if the density diagnostic is not).  The original WAVE smoother used ``smoothing=(4,4,4)`` and ``compensation=(1,1,1)``.

.. py:function:: smoothing = (smx,smy,smz)

	Used to perform smoothing passes (0.25, 0.5, 0.25) on the source functions.

	:param int smx: number of smoothing passes in x-direction.
	:param int smy: number of smoothing passes in y-direction.
	:param int smz: number of smoothing passes in z-direction.

.. py:function:: compensation = (cnx,cny,cnz)

	Apply compensation passes (-1.25 , 3.5 , -1.25) after smoother.

	:param int cnx: number of compensation passes in x-direction.
	:param int cny: number of compensation passes in y-direction.
	:param int cnz: number of compensation passes in z-direction.

Field Solvers
-------------

.. py:function:: new electrostatic [field solver] [<name>] { <directives> }

	Installs a general electrostatic field solver.  Inhomogeneous boundary conditions are imposed using
	:ref:`Conductor <conductor>` objects.  Conductor objects may fill either ghost cells or interior cells.
	The topology of conductor objects is unrestricted.  However, elliptical solvers have varying restrictions
	on the complexity of the topology.  An elliptic solver tool is required.

	:param block directives: The following directives are supported:

		Installable tools: :ref:`elliptic`

.. _coulomb-solver:
.. py:function:: new coulomb [field solver] [<name>] { <directives> }

	Installs a Coulomb gauge electromagnetic field solver.  Similar to WAVE field solver, but assumes continuity of sources (no divergence cleaning structure).  Cartesian coordinates only.  An elliptic solver tool is required.  At present only the ``facr`` elliptic solver is recommended.

	:param block directives: The following directives are supported:

		Installable tools: :ref:`elliptic`, :ref:`radiation`

		.. py:function::	dipole center = (x,y,z)

			Reference point for dipole moment diagnostic

		.. py:function::	gamma beam = g

		 	Reference gamma for initializing beam potentials


.. _direct-solver:
.. py:function:: new direct [field solver] [<name>] { <directives> }

	Create an EM module that advances Maxwell's curl equations directly, relying on continuity of sources to preserve divergence conditions.  An elliptical solver is used for initialization.  Cartesian only.

	:param block directives: The following directives are supported:

		Installable tools: :ref:`elliptic`, :ref:`radiation`

		.. py:function::	dipole center = (x,y,z)

			Reference point for dipole moment diagnostic

		.. py:function:: layers = ( x0,x1,y0,y1,z0,z1 )

			thickness of absorbing layers.

			:param int x0: number of cells in a single strip occupied by absorbing layers adjacent to the lower boundary in the x direction.  If 0 there are no PML media at this boundary.  Other 5 parameters are analogous.

		.. py:function:: reflection coefficient = Rx0 Rx1 Ry0 Ry1 Rz0 Rz1

		 	:param float Rx0: Desired fraction of AMPLITUDE reflected from lower x boundary.  If actual reflection is larger than requested, try increasing the number of layers. Other 5 parameters are analogous.

.. py:function:: new curvilinear direct [field solver] [<name>] { <directives> }

	Same as :ref:`direct electromagnetic module <direct-solver>` except for arbitrary coordinate system. Elliptical solver should be ``eigenmode``.

.. py:function:: new pgc [laser solver] [<name>] { <directives> }

	Create an enveloped field solver suitable for use with ponderomotive guiding center simulations.

	:param block directives: The following directives are supported:

		Installable tools: :ref:`propagator`, :ref:`radiation`

		.. py:function:: carrier frequency = f

			:param float f: base frequency ratio for the laser radiation

		.. py:function::	polarization = p

			:param enum p: can be ``linear``, ``circular``, or ``radial``


Particle Species
----------------

Particle species can be used in electromagnetic PIC or as Bohmian particles in :doc:`bak-quantum`.

.. py:function:: new species [<name>] { <directives> }

	:param str name: name given to the species
	:param block directives: The following directives are supported:

		Shared directives: see :ref:`boundaries <boundaries>`

		Installable tools: :ref:`ionization`

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

.. py:function:: new bound [<name>] { <directives> }

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

			:param float q1: Rotation about z in radians.  Initial orientation has principal axes aligned with standard basis.  This rotation happens before the theta rotation.

		.. py:function:: theta = q2

			:param float q2: Rotation about y in radians.  Initial orientation has principal axes aligned with standard basis.  This rotation happens after the phi rotation.
