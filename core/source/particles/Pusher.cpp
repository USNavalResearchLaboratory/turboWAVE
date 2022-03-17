#include "meta_base.h"
#include "computeTool.h"
#include "bundle.h"
#include "pusher.h"
#include "slicer.h"
#include "tiler.h"
#include "mover.h"

void BundlePusherBoris::Push()
{
	impulse(u,F);
	rotation1(t,u,F);
	rotation2(s,t);
	rotation3(s,t,vel,u);
	impulse(u,F);
	velocity(vel,u);
	translate(x,vel);
	load_j4(J,number,vel);
}

void BundlePusherUnitary::Lambda()
{
	// Left multiply the spinor by Lambda (time translation operator)
	copy_spinor(zf,zi);

	copy_spinor(z,zi);
	set_psi24(a,F);
	scalar_mul(a,z);
	add_spinor(zf,z);

	copy_spinor(z,zi);
	left_mul_sig1(z);
	set_psi_1(a,F);
	scalar_mul(a,z);
	add_spinor(zf,z);

	copy_spinor(z,zi);
	left_mul_sig2(z);
	set_psi_2(a,F);
	scalar_mul(a,z);
	add_spinor(zf,z);

	copy_spinor(z,zi);
	left_mul_sig3(z);
	set_psi_3(a,F);
	scalar_mul(a,z);
	add_spinor(zf,z);
}

void BundlePusherUnitary::Push()
{
	// We rewrite L*z*L^dag as (L*(L*z)^dag)^dag
	to_spinor(u,zi);
	estimate_ds(ds,F,u); // ds gets buried in F; at present must go after to_spinor
	Lambda();
	copy_spinor(zi,zf);
	dagger(zi);
	Lambda();
	dagger(zf);
	to_vector(u,zf,a);

	velocity(vel,u);
	translate(x,vel);
	load_j4(J,number,vel);
}

void BundlePusherPGC::Push()
{
	avg_gam_1(avgGam,vel,u,F,las);
	impulse(u,F,las,avgGam);
	rotation1(t,F,avgGam);
	rotation2(s,t);
	rotation3(s,t,vel,u);
	impulse(u,F,las,avgGam);
	avg_gam_2(avgGam,u,las);
	velocity(vel,u,avgGam);
	translate(x,vel);
	load_j4(J,number,vel);
	load_chi(chi,number,avgGam);
}

void BundlePusherBohmian::Push()
{
	// estimate vel(n+1/2) using u(n-1) and J(n), and update u(n-1) to u(n)
	bohm_velocity(vel,u,J);
	// update x(n) to x(n+1) using vel(n+1/2)
	translate(x,vel);
}
