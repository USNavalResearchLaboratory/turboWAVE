struct Mover;

struct ParticleBundle
{
	// The purpose of the particle bundle is to allow for efficient vectorization of the pusher
	// Several inlined functions are needed due to OpenMP inability to work with member variables
	// Evidently this is related to uninstantiated variables being tagged as aligned
	// The workaround is to pass members back in as arguments
	Mover *owner;
	tw::Float q0,m0,k[4];
	tw::Int num,cell0,ijk0[4];
	static const tw::Int N = tw::max_bundle_size;
	static const tw::Int AB = tw::vec_align_bytes;
	// arrays storing particle state
	alignas(AB) tw::Int cell[N];
	alignas(AB) float x[4][N];
	alignas(AB) tw::Float u[4][N];
	alignas(AB) float number[N];
	// temporary arrays used during calculation
	alignas(AB) float J[4][N];
	alignas(AB) float domainMask[N];
	alignas(AB) float cellMask[N];
	alignas(AB) tw::Int ijk[4][N];
	alignas(AB) tw::Float vel[4][N];
	alignas(AB) float w0[3][3][N];
	alignas(AB) float w1[3][3][N];
	alignas(AB) float l0[3][3][N];

	std::valarray<Particle*> refs;

	ParticleBundle(Mover *owner);
	void PadBundle();
	void Reset();
	bool Complete(const Particle& par);
	void Append(Particle& par);
	void CopyBack();
	void PrepareGather();
	void PrepareScatter();
	void translate(float x[4][N],tw::Float vel[4][N],tw::Float dt);
	void set_cell_mask(float cellMask[N],tw::Int cell0,tw::Int cell[N]);
	void get_cell_displ(tw::Int dc[3],tw::Int n);
	void load_j4(float J[4][N],const float number[N],const tw::Float vel[4][N]);
	void ZeroArray(float q[][N],tw::Int s1,tw::Int s2);
	void ZeroArray(tw::Float q[][N],tw::Int s1,tw::Int s2);
};

inline void ParticleBundle::PadBundle()
{
	for (int ax=0;ax<4;ax++)
		for (int i=num;i<N;i++)
		{
			x[ax][i] = 0.0;
			u[ax][i] = 0.0;
		}
	for (int i=num;i<N;i++)
	{
		u[0][i] = m0==0.0 ? 1.0 : m0; // u needs to be a good 4-vector
		u[3][i] = m0==0.0 ? 1.0 : 0.0;
		cell[i] = cell[0];
		number[i] = 0.0;
	}
}
inline void ParticleBundle::Reset()
{
	num = 0;
}
inline bool ParticleBundle::Complete(const Particle& par)
{
	return (num && par.q.cell!=cell[0]) || num==N;
}
inline void ParticleBundle::Append(Particle& par)
{
	refs[num] = &par;
	cell[num] = par.q.cell;
	x[0][num] = par.q.x[0];
	x[1][num] = par.q.x[1];
	x[2][num] = par.q.x[2];
	x[3][num] = par.q.x[3];
	u[0][num] = par.p[0]*(m0+tw::tiny)/(sqr(m0)+tw::tiny);
	u[1][num] = par.p[1]*(m0+tw::tiny)/(sqr(m0)+tw::tiny);
	u[2][num] = par.p[2]*(m0+tw::tiny)/(sqr(m0)+tw::tiny);
	u[3][num] = par.p[3]*(m0+tw::tiny)/(sqr(m0)+tw::tiny);
	number[num] = par.number;
	num++;
}
inline void ParticleBundle::translate(float x[4][N],tw::Float vel[4][N],tw::Float dt)
{
	for (int c=0;c<4;c++)
		#pragma omp simd aligned(x,vel:AB)
		for (int i=0;i<N;i++)
			x[c][i] += vel[c][i]*dt*k[c];
}
inline void ParticleBundle::set_cell_mask(float cellMask[N],tw::Int cell0,tw::Int cell[N])
{
	for (int i=0;i<N;i++)
		cellMask[i] = float(cell[i]==cell0);
}
inline void ParticleBundle::get_cell_displ(tw::Int dc[3],tw::Int n)
{
	dc[0] = ijk[1][n] - ijk0[1];
	dc[1] = ijk[2][n] - ijk0[2];
	dc[2] = ijk[3][n] - ijk0[3];
}
inline void ParticleBundle::load_j4(float J[4][N],const float number[N],const tw::Float vel[4][N])
{
	#pragma omp simd aligned(J,number:AB)
	for (int i=0;i<N;i++)
		J[0][i] = q0*number[i];
	#pragma omp simd aligned(J,number,vel:AB)
	for (int i=0;i<N;i++)
		J[1][i] = q0*number[i]*vel[1][i]*k[1];
	#pragma omp simd aligned(J,number,vel:AB)
	for (int i=0;i<N;i++)
		J[2][i] = q0*number[i]*vel[2][i]*k[2];
	#pragma omp simd aligned(J,number,vel:AB)
	for (int i=0;i<N;i++)
		J[3][i] = q0*number[i]*vel[3][i]*k[3];
}
inline void ParticleBundle::ZeroArray(float q[][N],tw::Int s1,tw::Int s2)
{
	for (tw::Int s=s1;s<=s2;s++)
		for (tw::Int par=0;par<N;par++)
			q[s][par] = 0.0f;
}
inline void ParticleBundle::ZeroArray(tw::Float q[][N],tw::Int s1,tw::Int s2)
{
	for (tw::Int s=s1;s<=s2;s++)
		for (tw::Int par=0;par<N;par++)
			q[s][par] = 0.0f;
}
