Quantum Optics Modules
======================

The quantum optics modules solve the four basic quantum mechanical wave equations under time dependent excitation by an arbitrary electromagnetic field.  The primary application is to solve for the response of an atom to a strong laser field.  Admittedly, "quantum optics" is strictly a misnomer, since the radiation field is treated classically.

Wave Equations
--------------

The wave equations that are treated are the Schroedinger equation, the Pauli equation, the Klein-Gordon equation, and the Dirac equation.  These can be presented in a :math:`2\times 2` table which categorizes the spin and the energy-momentum (dispersion) relation

.. csv-table:: Table I. Wave Equations.

	"", spin = 0, spin = :math:`\hbar/2`
	:math:`\hbar\omega = (\hbar k)^2/2m`, "Schroedinger", "Pauli"
	:math:`(\hbar\omega)^2 = (mc^2)^2 + (c\hbar k)^2`, "Klein-Gordon", "Dirac"

In the case of the relativistic equations, it is sometimes useful to consider cylindrical atoms which allow for 2D geometry.  This should not be confused with the use of cylindrical geometry, which can be used to simulate a spherical atom in the non-relativistic case, in the dipole approximation.  In contrast, relativistic cylindrical atoms are simulated in 2D Cartesian coordinates, with the ignorable coordinate lying along the atom's symmetry axis.

Quantum Units
-------------

For the non-relativistic wave equations, atomic units are used for all input file parameters and output file data.  For the relativistic wave equations, natural units are used for all input file parameters and output file data.  Tables II and III give the unit conversions.  To convert an atomic or natural unit to Gaussian or SI units, multiply by the formula in the cgs or mks column, respectively.

.. csv-table:: Table II. Atomic Units to Gaussian or SI Units.
	:header: "Quantity", "Unit", "cgs", "mks"

	"Charge", "electronic charge", :math:`|e|`, :math:`|e|`
	"Mass", "electronic mass", :math:`m`, :math:`m`
	"Time", "atomic time", :math:`t_a = \hbar/m(\alpha c)^2`, :math:`t_a = \hbar/m(\alpha c)^2`
	"Space", "Bohr Radius", :math:`a_0 = \hbar/m(\alpha c)`, :math:`a_0 = \hbar/m(\alpha c)`
	"Velocity", "",  :math:`\alpha c`, :math:`\alpha c`
	"Energy", "Hartree",  :math:`Ha = m(\alpha c)^2`, :math:`Ha = m(\alpha c)^2`
	"Momentum", "",  :math:`m\alpha c`, :math:`m\alpha c`
	"Scalar Potential", "",  :math:`Ha/|e|`, :math:`Ha/|e|`
	"Vector Potential", "",  :math:`Ha/\alpha|e|`, :math:`Ha/\alpha|e|c`
	"Electric Field", "", :math:`Ha/a_0|e|`, :math:`Ha/a_0|e|`
	"Magnetic Field", "", :math:`Ha/\alpha a_0|e|`, :math:`Ha/\alpha a_0|e|c`

It should be mentioned there is another way to normalize the vector potential. In the scheme of table II a factor of :math:`\alpha` is eliminated from the Schroedinger equation, compared with the other possible scheme.  The way to understand this is to look at the Schroedinger equation written in various units:

	**cgs** --- :math:`i\hbar\partial\psi/\partial t = \left[\left(i\hbar\nabla + Q{\bf A}/c\right)^2/2M + Q\Phi\right]\psi`

	**mks** --- :math:`i\hbar\partial\psi/\partial t = \left[\left(i\hbar\nabla + Q{\bf A}\right)^2/2M + Q\Phi\right]\psi`

	**TW** --- :math:`i\partial\psi/\partial t = \left[\left(i\nabla + Q{\bf A}\right)^2/2M + Q\Phi\right]\psi`

Here, Q and M are the charge and mass of the simulated particle, written in capital letters to distinguish them from the electronic charge and mass.  By including the factor of :math:`\alpha` in the normalization of the vector potential, it disappears from the TW Schroedinger equation.  If this had not been done, then the TW equation would resemble more the cgs equation, with :math:`\alpha` playing the part of the speed of light.

.. csv-table:: Table III. Natural Units to Gaussian or SI Units.
	:header: "Quantity", "Unit", "cgs", "mks"

	"Charge", "Expanded charge", :math:`q = \alpha^{-1/2}|e|`, :math:`q = \alpha^{-1/2}|e|`
	"Mass", "electronic mass", :math:`m`, :math:`m`
	"Time", "QED time", :math:`\hbar/mc^2`, :math:`\hbar/mc^2`
	"Space", "Reduced Compton Wavelength", :math:`\Lambda = \hbar/mc`, :math:`\Lambda = \hbar/mc`
	"Velocity", "Speed of Light",  :math:`c`, :math:`c`
	"Energy", "",  :math:`mc^2`, :math:`mc^2`
	"Momentum", "",  :math:`mc`, :math:`mc`
	"Scalar Potential", "",  :math:`mc^2/q`, :math:`mc^2/q`
	"Vector Potential", "",  :math:`mc^2/q`, :math:`mc/q`
	"Electric Field", "", :math:`mc^2/\Lambda q`, :math:`mc^2/\Lambda q`
	"Magnetic Field", "", :math:`mc^2/\Lambda q`, :math:`mc/\Lambda q`

To illustrate the effect of natural units on a wave equation, the Klein-Gordon equation is presented in various units:

	**cgs** --- :math:`\left[\left(i\hbar\partial/\partial t - Q\Phi\right)^2 - \left(i\hbar c\nabla + Q{\bf A}\right)^2 - (Mc^2)^2\right]\psi = 0`

	**mks** --- :math:`\left[\left(i\hbar\partial/\partial t - Q\Phi\right)^2 - \left(i\hbar c\nabla + Qc{\bf A}\right)^2 - (Mc^2)^2\right]\psi = 0`

	**TW** --- :math:`\left[\left(i\partial/\partial t - Q\Phi\right)^2 - \left(i\nabla + Q{\bf A}\right)^2 - M^2\right]\psi = 0`

The question of how to normalize the vector potential does not come up because in natural units :math:`c=1`, whereas in atomic units :math:`c=\alpha`.

.. note::

	The Schwinger field is not the unit of electric field, because of the use of the expanded charge in the unit of electric field.  In natural units the Schwinger field is :math:`q/|e| = \alpha^{-1/2}`.

Quantum Numbers
----------------

In order to initialize most quantum optics simulations, it is necessary to create one or more stationary states.  Each stationary state is associated with a given static field.  Typical static fields include the Coulomb field of the nucleus, an approximation of the Coulomb field called a soft-core potential, or a nuclear field superposed with a uniform magnetic field.

The quantum numbers identifying a given stationary state differ depending on the wave equation that is selected.  In order to bring as much uniformity as possible to the situation, we adopt a nomenclature that is based on the most general of the four equations, the Dirac equation.  The four quantum numbers that are used to characterize any turboWAVE quantum state are :math:`n_r, j, l, m`.  The values they take in various situations are listed in table IV.  There is also an auxiliary number, :math:`\kappa`, that can be helpful in formulating the state.

.. csv-table:: Table IV. Quantum Numbers.
	:header: "Theory", "Symbol", "Interpretation", "Values", "Comments"
	:delim: ;

	"All"; :math:`n_r`; Radial Number; :math:`[0,1,2,...]`; cannot be 0 for :math:`\kappa>0`
	"3D Spin 0";j;Total Angular Momentum; :math:`[0,1,2,...]`; :math:`l=j`
	"3D Spin 0";l;Orbital Angular Momentum; :math:`[0,1,2,...]`; :math:`l=j`
	"3D Spin 0";m;Angular Momentum Projection; :math:`[-l,-l+1,...,l]`; ""
	"3D Spin 1/2";j;Total Angular Momentum; :math:`[1/2,3/2,5/2,...]`; ""
	"3D Spin 1/2";l;Parity Number; :math:`[j+1/2,j-1/2]`; :math:`P=(-1)^l`
	"3D Spin 1/2"; m; Angular Momentum Projection; :math:`[-j,-j+1,...,j]`; ""
	"3D Spin 1/2"; :math:`\kappa`; auxiliary number; :math:`2(l-j)(j+0.5)`; :math:`\omega = \omega(n_r,\kappa)`
	"2D Spin 0";j;Total Angular Momentum; integers; :math:`m=l=j`
	"2D Spin 0";l;Orbital Angular Momentum; integers; :math:`m=l=j`
	"2D Spin 0";m;Angular Momentum Projection; integers; :math:`m=l=j`
	"2D Spin 1/2";j;Total Angular Momentum; half integers; :math:`m=j`
	"2D Spin 1/2";l;Parity Number; :math:`[m+1/2,m-1/2]`; :math:`P=(-1)^l`
	"2D Spin 1/2";m;Angular Momentum Projection; half integers; :math:`m=j`
	"2D Spin 1/2"; :math:`\kappa`; auxiliary number; :math:`2(l-m)m`; :math:`\omega = \omega(n_r,\kappa)`

.. csv-table:: Table V. Quantum Number Variable Mapping
	:header: "Symbol", "Code Variable"

	:math:`n_r`, :samp:`nr`
	:math:`j`, :samp:`Jam`
	:math:`l`, :samp:`Lam`
	:math:`m`, :samp:`jzam`
