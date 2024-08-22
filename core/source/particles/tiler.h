struct BundleTilerEM2D : BundleSlicerEM
{
	float F_tile[3][3][6];
	float J_tile[5][5][4];

	BundleTilerEM2D(Mover *owner) : ParticleBundle(owner), BundleSlicerEM(owner) {}
	void LoadFTile();
	void ResetJTile();
	void StoreJTile();
	/// @brief load field tensor from tile for N particles
	/// @param F this is shorthand for F*q*dt/2/m, where F is a field amplitude
	/// @param w0 weight factors for N particles
	/// @param l0 wall weights for N particles
	/// @param qmdth (q/m)*(dt/2)
	void GatherF(float F[6][N],const float w0[3][3][N],const float l0[3][3][N],const float qmdth);
	/// @brief scatter current for N particles into tile
	/// @param J amps to deposit, the tile will receive coulombs in cell and amps in wall
	/// @param w0 weight factors for N particles at start of step
	/// @param w1 weight factors for N particles at end of step
	/// @param cellMask used to mask out particles that need scalar processing
	/// @param dti inverse time step
	void ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti);
};

struct BundleTilerEM3D : BundleSlicerEM
{
	float F_tile[3][3][3][6];
	float J_tile[5][5][5][4];

	BundleTilerEM3D(Mover *owner) : ParticleBundle(owner), BundleSlicerEM(owner) {}
	void LoadFTile();
	void ResetJTile();
	void StoreJTile();
	/// @brief load field tensor from tile for N particles
	/// @param F this is shorthand for F*q*dt/2/m, where F is a field amplitude
	/// @param w0 weight factors for N particles
	/// @param l0 wall weights for N particles
	/// @param qmdth (q/m)*(dt/2)
	void GatherF(float F[6][N],const float w0[3][3][N],const float l0[3][3][N],const float qmdth);
	/// @brief scatter current for N particles into tile
	/// @param J amps to deposit, the tile will receive coulombs in cell and amps in wall
	/// @param w0 weight factors for N particles at start of step
	/// @param w1 weight factors for N particles at end of step
	/// @param cellMask used to mask out particles that need scalar processing
	/// @param dti inverse time step
	void ScatterJ4(const float J[4][N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N],const float& dti);
};

struct BundleTilerPGC2D : BundleSlicerPGC
{
	float las_tile[3][3][8];
	float chi_tile[5][5];

	BundleTilerPGC2D(Mover *owner) : ParticleBundle(owner), BundleSlicerPGC(owner)  {}
	void LoadLaserTile();
	void ResetChiTile();
	void StoreChiTile();
	/// @brief load laser envelope data from tile for N particles
	/// @param las [q2m2dth*grad(a^2(n)), q2m2dth*grad(a^2(n+1/2)), q2m2h*a^2(n), q2m2h*a^2(n+1/2)]
	/// @param w0 weight factors for N particles
	/// @param q2m2dth (q^2/m^2)*(dt/2), n.b. we are usually updating u = p/m
	/// @param q2m2h (q^2/m^2)*(1/2)
	void GatherLaser(float las[8][N],const float w0[3][3][N],const float q2m2dth,const float q2m2h);
	/// @brief scatter current susceptibility for N particles into tile
	/// @param chi current susceptibility defined by -q0*q0*number/m0/avgGam
	/// @param w0 weight factors for N particles at start of step
	/// @param w1 weight factors for N particles at end of step
	/// @param cellMask used to mask out particles that need scalar processing
	void ScatterChi(const float chi[N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N]);
};

struct BundleTilerPGC3D : BundleSlicerPGC
{
	float las_tile[3][3][3][8];
	float chi_tile[5][5][5];

	BundleTilerPGC3D(Mover *owner) : ParticleBundle(owner), BundleSlicerPGC(owner)  {}
	void LoadLaserTile();
	void ResetChiTile();
	void StoreChiTile();
	/// @brief load laser envelope data from tile for N particles
	/// @param las [q2m2dth*grad(a^2(n)), q2m2dth*grad(a^2(n+1/2)), q2m2h*a^2(n), q2m2h*a^2(n+1/2)]
	/// @param w0 weight factors for N particles
	/// @param q2m2dth (q^2/m^2)*(dt/2), n.b. we are usually updating u = p/m
	/// @param q2m2h (q^2/m^2)*(1/2)
	void GatherLaser(float las[8][N],const float w0[3][3][N],const float q2m2dth,const float q2m2h);
	/// @brief scatter current susceptibility for N particles into tile
	/// @param chi current susceptibility defined by -q0*q0*number/m0/avgGam
	/// @param w0 weight factors for N particles at start of step
	/// @param w1 weight factors for N particles at end of step
	/// @param cellMask used to mask out particles that need scalar processing
	void ScatterChi(const float chi[N],const float w0[3][3][N],const float w1[3][3][N],const float cellMask[N]);
};

struct BundleTilerBohmian2D : BundleSlicerBohmian
{
	float tile[3][3][4];
	BundleTilerBohmian2D(Mover *owner) : ParticleBundle(owner), BundleSlicerBohmian(owner)  {}
	void LoadTile();
	void GatherJ4(float J[4][N],const float w0[3][3][N]);
};

struct BundleTilerBohmian3D : BundleSlicerBohmian
{
	float tile[3][3][3][4];
	BundleTilerBohmian3D(Mover *owner) : ParticleBundle(owner), BundleSlicerBohmian(owner)  {}
	void LoadTile();
	void GatherJ4(float J[4][N],const float w0[3][3][N]);
};

//////////
inline void BundleTilerEM2D::LoadFTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int k=0;k<3;k++)
			for (tw::Int s=0;s<6;s++)
				F_tile[i][k][s] = BundleSlicerEM::Fx(ijk0[1]-1+i,0,ijk0[3]-1+k,s);
}
inline void BundleTilerEM2D::ResetJTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int k=0;k<5;k++)
			for (tw::Int s=0;s<4;s++)
				J_tile[i][k][s] = 0.0f;
}
inline void BundleTilerEM2D::StoreJTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int k=0;k<5;k++)
			for (tw::Int s=0;s<4;s++)
				BundleSlicerEM::Jx(ijk0[1]-2+i,0,ijk0[3]-2+k,s) += J_tile[i][k][s];
}

//////////
inline void BundleTilerEM3D::LoadFTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
				for (tw::Int s=0;s<6;s++)
					F_tile[i][j][k][s] = BundleSlicerEM::Fx(ijk0[1]-1+i,ijk0[2]-1+j,ijk0[3]-1+k,s);
}
inline void BundleTilerEM3D::ResetJTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
				for (tw::Int s=0;s<4;s++)
					J_tile[i][j][k][s] = 0.0f;
}
inline void BundleTilerEM3D::StoreJTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
				for (tw::Int s=0;s<4;s++)
					BundleSlicerEM::Jx(ijk0[1]-2+i,ijk0[2]-2+j,ijk0[3]-2+k,s) += J_tile[i][j][k][s];
}

//////////
inline void BundleTilerPGC2D::LoadLaserTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int k=0;k<3;k++)
			for (tw::Int s=0;s<8;s++)
				las_tile[i][k][s] = BundleSlicerPGC::lasx(ijk0[1]-1+i,0,ijk0[3]-1+k,s);
}
inline void BundleTilerPGC2D::ResetChiTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int k=0;k<5;k++)
			chi_tile[i][k] = 0.0f;
}
inline void BundleTilerPGC2D::StoreChiTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int k=0;k<5;k++)
			BundleSlicerPGC::chix(ijk0[1]-2+i,0,ijk0[3]-2+k,0) += chi_tile[i][k];
}

//////////
inline void BundleTilerPGC3D::LoadLaserTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
				for (tw::Int s=0;s<8;s++)
					las_tile[i][j][k][s] = BundleSlicerPGC::lasx(ijk0[1]-1+i,ijk0[2]-1+j,ijk0[3]-1+k,s);
}
inline void BundleTilerPGC3D::ResetChiTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
				chi_tile[i][j][k] = 0.0f;
}
inline void BundleTilerPGC3D::StoreChiTile()
{
	for (tw::Int i=0;i<5;i++)
		for (tw::Int j=0;j<5;j++)
			for (tw::Int k=0;k<5;k++)
				BundleSlicerPGC::chix(ijk0[1]-2+i,ijk0[2]-2+j,ijk0[3]-2+k,0) += chi_tile[i][j][k];
}

//////////
inline void BundleTilerBohmian2D::LoadTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int k=0;k<3;k++)
			for (tw::Int s=0;s<4;s++)
				tile[i][k][s] = BundleSlicerBohmian::Jx(ijk0[1]-1+i,0,ijk0[3]-1+k,s);
}

//////////
inline void BundleTilerBohmian3D::LoadTile()
{
	for (tw::Int i=0;i<3;i++)
		for (tw::Int j=0;j<3;j++)
			for (tw::Int k=0;k<3;k++)
				for (tw::Int s=0;s<4;s++)
					tile[i][j][k][s] = BundleSlicerBohmian::Jx(ijk0[1]-1+i,ijk0[2]-1+j,ijk0[3]-1+k,s);
}
