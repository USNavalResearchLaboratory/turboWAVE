Electromagnetic Wave Modes
==========================

Electromagnetic wave modes are are monochromatic, continuous wave solutions of Maxwell's equations.  TurboWAVE makes use of such solutions to inject radiation into the simulation.  Since the analytical expressions are typically approximate in some way, turboWAVE numerically forces the initial or boundary fields into compliance with Maxwell's divergence equations.  The total solution is then guaranteed to satisfy Maxwell's equations, to within a discretization error.

The initial or boundary fields are always derived from a vector potential in the Coulomb gauge, regardless of the field solver that is used.  It follows that :math:`\nabla\cdot {\bf B} = 0` is automatically satisfied (divergence of a curl is identically zero).  During initialization, an elliptical solver is used to guarantee :math:`\nabla\cdot {\bf A} = 0`.  The elliptical solver is also used to generate electrostatic fields in case there is a charge separation in the initial condition.

Coordinates
-----------

A point in Cartesian coordinates is denoted

:math:`(x,y,z)`

A point in cylindrical coordinates is denoted

:math:`(\varrho,\varphi,z)`

A point in spherical coordinates is denoted

:math:`(r,\varphi,\theta)`

Note that the azimuthal variable is always :math:`\varphi`.  The polar angle is :math:`\theta`.

Phase and Envelope Functions
----------------------------

All the modes considered here are in the form

:math:`A({\bf r})T(\tau)e^{i\psi}`

where *A* is some component of the vector potential (not necessarily Cartesian), *T* is a pulse envelope function, and :math:`\psi` is a phase function that may have a complex spatial dependence.  The argument of the pulse function, :math:`\tau`, can also have a complex spatial dependence.  The envelope and phase function arguments are typically closely related.  To express this fact, the phase is always put in the form

:math:`\psi = \psi_0 + \psi' - \omega\tau - \chi\tau^2`

Here :math:`\omega` is the center frequency of the wave, :math:`\chi` is the ``chirp`` parameter, and :math:`\psi_0` is a reference phase set by the user independently for each wave in the system.  The anomalous phase term :math:`\psi'` allows for arbitrary modifications to the relationship between the phase and amplitude contours.

In general, all functions may depend on the refractive index, :math:`\eta`.

The pulse envelope factor *T* is put into the monochromatic solution by hand, and will generally introduce an error, except in the eikonal limit. Since the divergence error is corrected numerically, the factor *T* introduces no error into the overall solution, in cases where a Maxwell field solver is used to advance the fields in time.

.. tip::

	The origin of space and time are set by the user through the :ref:`wave <wave-obj>` parameters ``focus position``, ``delay`` and ``risetime``.  The coordinate system can also be rotated arbitrarily.

.. note::

	Normalized units are used in the following.

Plane Wave
----------

For a plane wave we have the constant amplitude

:math:`A({\bf r}) = a_0`

and the pulse function argument

:math:`\tau = t - \eta z`

with no anomalous phase.

Bessel Beams
------------------

The Bessel beam modes are exact Maxwell solutions that are separable in cylindrical coordinates.  As of this writing, only the lowest order is supported:

:math:`A({\bf r}) = a_0 J_0(\varrho/r_0)`

The pulse and phase variables are the same as for a plane wave.

Hermite Gaussian
------------------

The Hermite-Gaussian modes are paraxial beam modes that are separable in Cartesian coordinates.  The amplitude and pulse function argument are in the forms

:math:`A({\bf r}) = a_0 A_m(x,z)A_n(y,z)`

:math:`\tau = t - \eta z + \tau_m(x,z) + \tau_n(y,z)`

where *m* and *n* are the mode numbers along the *x* and *y* axes.  The anomalous phase term is twice the Guoy phase

:math:`\psi'({\bf r}) = 2\Gamma_m(z) + 2\Gamma_n(z)`

The factor of two reverses the sign of an equal and opposite Guoy phase that is put in the envelope function to make it subluminal in the confocal region.

For the *x* direction mode factors:

:math:`A_m = \sqrt{\frac{1}{2^m m!}} \sqrt{\frac{x_0}{x_m}} H_m\left(\frac{\sqrt{2}x}{x_m}\right)e^{-x^2/x_m^2}`

:math:`\tau_m = \Gamma_m/\omega - \frac{\eta z}{2}\frac{x^2}{z^2+z_x^2}`

where

:math:`z_x = \frac{1}{2} \eta \omega x_0^2`

:math:`x_m = x_0\sqrt{1 + \frac{z^2}{z_x^2}}`

:math:`\Gamma_m = -\left(m+\frac{1}{2}\right)\tan^{-1}\frac{z}{z_x}`

Here the inputs are the index of refraction :math:`\eta`, the waist radius :math:`x_0`, and the peak vector potential :math:`a_0`.  The factors for the *y* direction are exactly analogous.

Laguerre Gaussian
------------------

The Laguerre-Gaussian modes are paraxial beam modes that are separable in cylindrical coordinates.  The amplitude factor and pulse function argument are in the forms

:math:`A({\bf r}) = a_0 A_{nm}(\varrho,z)`

:math:`\tau = t - \eta z + \tau_{nm}(\varrho,z)`

where *n* and *m* are the radial and azimuthal mode numbers.  The anomalous phase is

:math:`\psi'({\bf r}) = 2\Gamma_{nm}(z) - m\varphi`

As in the Hermite case the factor of two multiplying the Guoy phase compensates for a reverse Guoy phase in the pulse function argument.  The radial factors are:

:math:`A_{nm} = a_0 \sqrt{\frac{n!}{(n+m)!}} \left(\frac{\sqrt{2}\varrho}{r_m}\right)^m \frac{r_0}{r_m} L_{nm}\left(\frac{2\varrho^2}{r_m^2}\right)e^{-\varrho^2/r_m^2}`

:math:`\tau_{nm} = \Gamma_{nm}/\omega - \frac{\eta z}{2}\frac{\varrho^2}{z^2+z_R^2}`

where

:math:`z_R = \frac{1}{2} \eta \omega r_0^2`

:math:`r_m = r_0\sqrt{1 + \frac{z^2}{z_R^2}}`

:math:`\Gamma_{nm} = -\left(2n+m+1\right)\tan^{-1}\frac{z}{z_R}`

Here the inputs are the index of refraction :math:`\eta`, the waist radius :math:`r_0`, and the peak vector potential :math:`a_0`.

Multipole Fields
------------------
The multipole fields are exact Maxwell solutions that are separable in spherical coordinates.  They are also the states of definite angular momentum of the photon.  As of this writing, only the lowest order magnetic multipole (dipole) is supported.  Multipole radiation is most compactly expressed as a standing wave.  For a magnetic multipole of any order we have

:math:`{\bf A} \propto j_l(kr) {\bf \Phi}_{lm}(\theta,\varphi) e^{-i\omega t}`

where :math:`{\bf \Phi}_{lm}` is the appropriate vector spherical harmonic function, and :math:`j_l` is the spherical Bessel function.  The spherical Bessel function is real valued and carries the rapidly varying spatial dependence.

For a magnetic dipole field, it is straightforward to decompose the Bessel function into incoming and outgoing waves, which allows the field to be put in the standard form above.  The two solutions are:

:math:`{\bf A}_\pm = \frac{a_0}{2} \frac{\sin\theta}{j_1(x_1)} R_\pm(\omega_0,r) T(t\mp r) e^{-i\omega_0(t\mp r)} {\bf e}_\varphi`

Here :math:`x_1` is the coordinate of the first maximum of the spherical Bessel function, and

:math:`R_\pm = \pm\frac{1}{i\omega^2 r^2} - \frac{1}{\omega r}`
