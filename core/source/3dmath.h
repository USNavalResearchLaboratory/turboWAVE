namespace tw
{

	struct vec3
	{
		tw::Float x,y,z;

		vec3() noexcept
		{
			x = 0.0;
			y = 0.0;
			z = 0.0;
		}

		vec3(const tw::Float& x,const tw::Float& y,const tw::Float& z) noexcept
		{
			this->x = x;
			this->y = y;
			this->z = z;
		}

		vec3(const vec3& v) noexcept
		{
			x = v.x;
			y = v.y;
			z = v.z;
		}

		vec3(const tw::Float& a) noexcept
		{
			x = a;
			y = a;
			z = a;
		}


		vec3(const tw::Float *a) noexcept
		{
			x = a[0];
			y = a[1];
			z = a[2];
		}

		friend tw::Float Magnitude(const vec3& v)
		{
			return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
		}

		friend tw::Float Norm(const vec3& v)
		{
			return v.x*v.x + v.y*v.y + v.z*v.z;
		}

		friend tw::Float Normalize(vec3& v)
		{
			tw::Float mag;
			mag = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
			v.x /= mag;
			v.y /= mag;
			v.z /= mag;
			return mag;
		}

		const tw::Float& operator [] (const tw::Int& i) const
		{
			return (&x)[i];
		}

		tw::Float& operator [] (const tw::Int& i)
		{
			return (&x)[i];
		}

		vec3& operator - ()
		{
			x = -x;
			y = -y;
			z = -z;
			return *this;
		}

		vec3& operator = (const vec3& v) noexcept
		{
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		}

		vec3& operator = (const tw::Float& a) noexcept
		{
			x = a;
			y = a;
			z = a;
			return *this;
		}

		vec3& operator += (const vec3& v)
		{
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}

		vec3& operator -= (const vec3& v)
		{
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}

		vec3& operator *= (const tw::Float& a)
		{
			x *= a;
			y *= a;
			z *= a;
			return *this;
		}

		vec3& operator *= (const vec3& v)
		{
			x *= v.x;
			y *= v.y;
			z *= v.z;
			return *this;
		}

		vec3& operator /= (const tw::Float& a)
		{
			x /= a;
			y /= a;
			z /= a;
			return *this;
		}

		vec3& operator /= (const vec3& v)
		{
			x /= v.x;
			y /= v.y;
			z /= v.z;
			return *this;
		}

		friend vec3 operator + (const vec3& lhs,const vec3& rhs)
		{
			return vec3(lhs.x + rhs.x,lhs.y + rhs.y,lhs.z + rhs.z);
		}

		friend vec3 operator - (const vec3& lhs,const vec3& rhs)
		{
			return vec3(lhs.x - rhs.x,lhs.y - rhs.y,lhs.z - rhs.z);
		}

		friend tw::Float operator ^ (const vec3& lhs,const vec3& rhs)		// Dot Product
		{
			return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
		}

		friend vec3 operator | (const vec3& lhs,const vec3& rhs)		// Cross Product
		{
			return vec3(lhs.y*rhs.z - lhs.z*rhs.y,
						- lhs.x*rhs.z + lhs.z*rhs.x,
						lhs.x*rhs.y - lhs.y*rhs.x);
		}

		friend vec3 operator * (const vec3& lhs,const vec3& rhs)		// Direct Product
		{
			return vec3(lhs.x*rhs.x,lhs.y*rhs.y,lhs.z*rhs.z);
		}

		friend vec3 operator / (const vec3& lhs,const vec3& rhs)		// Direct Divide
		{
			return vec3(lhs.x/rhs.x,lhs.y/rhs.y,lhs.z/rhs.z);
		}

		friend vec3 operator * (const vec3& v,const tw::Float& a)	// Scaling
		{
			return vec3(v.x*a, v.y*a, v.z*a);
		}

		friend vec3 operator * (const tw::Float& a,const vec3& v)
		{
			return vec3(v.x*a, v.y*a, v.z*a);
		}

		friend vec3 operator / (const vec3& v,const tw::Float& a)
		{
			return vec3(v.x/a, v.y/a, v.z/a);
		}

		void RotateZ(tw::Float angle)
		{
			tw::Float st,ct,xOld;
			st = sin(angle);
			ct = cos(angle);
			xOld = x;
			x = xOld*ct - y*st;
			y = xOld*st + y*ct;
		}

		void RotateY(tw::Float angle)
		{
			tw::Float st,ct,zOld;
			st = sin(angle);
			ct = cos(angle);
			zOld = z;
			z = zOld*ct - x*st;
			x = zOld*st + x*ct;
		}

		void RotateX(tw::Float angle)
		{
			tw::Float st,ct,yOld;
			st = sin(angle);
			ct = cos(angle);
			yOld = y;
			y = yOld*ct - z*st;
			z = yOld*st + z*ct;
		}

		friend std::ostream& operator << (std::ostream& os, const vec3& v)
		{
			os << '(' << v.x << ',' << v.y << ',' << v.z << ')';
			return os;
		}
	};

	struct cvec3
	{
		tw::Complex x,y,z;

		cvec3()
		{
			x = 0.0;
			y = 0.0;
			z = 0.0;
		}

		cvec3(const tw::Complex& x,const tw::Complex& y,const tw::Complex& z)
		{
			this->x = x;
			this->y = y;
			this->z = z;
		}

		cvec3(const tw::Float& x,const tw::Float& y,const tw::Float& z)
		{
			this->x = x;
			this->y = y;
			this->z = z;
		}

		cvec3(const cvec3& v)
		{
			x = v.x;
			y = v.y;
			z = v.z;
		}

		cvec3(const tw::Complex& a)
		{
			x = a;
			y = a;
			z = a;
		}

		cvec3(const tw::Float& a)
		{
			x = a;
			y = a;
			z = a;
		}

		friend tw::Float Magnitude(const cvec3& v)
		{
			return sqrt(norm(v.x) + norm(v.y) + norm(v.z));
		}

		friend tw::Float Norm(const cvec3& v)
		{
			return norm(v.x) + norm(v.y) + norm(v.z);
		}

		friend tw::Float Normalize(cvec3& v)
		{
			tw::Float mag;
			mag = sqrt(norm(v.x) + norm(v.y) + norm(v.z));
			v.x /= mag;
			v.y /= mag;
			v.z /= mag;
			return mag;
		}

		cvec3& operator - ()
		{
			x = -x;
			y = -y;
			z = -z;
			return *this;
		}

		cvec3& operator = (const cvec3& v)
		{
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		}

		cvec3& operator = (const tw::Float& a)
		{
			x = a;
			y = a;
			z = a;
			return *this;
		}

		cvec3& operator += (const cvec3& v)
		{
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}

		cvec3& operator -= (const cvec3& v)
		{
			x -= v.x;
			y -= v.y;
			z -= v.z;
			return *this;
		}

		cvec3& operator *= (const tw::Float& a)
		{
			x *= a;
			y *= a;
			z *= a;
			return *this;
		}

		cvec3& operator *= (const cvec3& v)
		{
			x *= v.x;
			y *= v.y;
			z *= v.z;
			return *this;
		}

		cvec3& operator /= (const tw::Float& a)
		{
			x /= a;
			y /= a;
			z /= a;
			return *this;
		}

		friend cvec3 operator + (const cvec3& lhs,const cvec3& rhs)
		{
			return cvec3(lhs.x + rhs.x,lhs.y + rhs.y,lhs.z + rhs.z);
		}

		friend cvec3 operator - (const cvec3& lhs,const cvec3& rhs)
		{
			return cvec3(lhs.x - rhs.x,lhs.y - rhs.y,lhs.z - rhs.z);
		}

		friend tw::Complex operator ^ (const cvec3& lhs,const cvec3& rhs)		// Dot Product
		{
			return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
		}

		friend cvec3 operator | (const cvec3& lhs,const cvec3& rhs)		// Cross Product
		{
			return cvec3(lhs.y*rhs.z - lhs.z*rhs.y,
						- lhs.x*rhs.z + lhs.z*rhs.x,
						lhs.x*rhs.y - lhs.y*rhs.x);
		}

		friend cvec3 operator * (const cvec3& lhs,const cvec3& rhs)		// Direct Product
		{
			return cvec3(lhs.x*rhs.x,lhs.y*rhs.y,lhs.z*rhs.z);
		}

		friend cvec3 operator / (const cvec3& lhs,const cvec3& rhs)		// Direct Divide
		{
			return cvec3(lhs.x/rhs.x,lhs.y/rhs.y,lhs.z/rhs.z);
		}

		friend cvec3 operator * (const cvec3& v,const tw::Float& a)	// Scaling
		{
			return cvec3(v.x*a, v.y*a, v.z*a);
		}

		friend cvec3 operator * (const tw::Float& a,const cvec3& v)
		{
			return cvec3(v.x*a, v.y*a, v.z*a);
		}

		friend cvec3 operator / (const cvec3& v,const tw::Float& a)
		{
			return cvec3(v.x/a, v.y/a, v.z/a);
		}

		void RotateZ(tw::Float angle)
		{
			tw::Complex xOld;
			tw::Float st,ct;
			st = sin(angle);
			ct = cos(angle);
			xOld = x;
			x = xOld*ct - y*st;
			y = xOld*st + y*ct;
		}

		void RotateY(tw::Float angle)
		{
			tw::Complex zOld;
			tw::Float st,ct;
			st = sin(angle);
			ct = cos(angle);
			zOld = z;
			z = zOld*ct - x*st;
			x = zOld*st + x*ct;
		}

		void RotateX(tw::Float angle)
		{
			tw::Complex yOld;
			tw::Float st,ct;
			st = sin(angle);
			ct = cos(angle);
			yOld = y;
			y = yOld*ct - z*st;
			z = yOld*st + z*ct;
		}

		friend std::ostream& operator << (std::ostream& os, const cvec3& v)
		{
			os << '(' << v.x << ',' << v.y << ',' << v.z << ')';
			return os;
		}
	};

	struct basis
	{
		vec3 u,v,w;

		basis()
		{
		}

		void SetWithEulerAngles(tw::Float alpha,tw::Float beta,tw::Float gamma)
		{
			u = vec3(1,0,0); v = vec3(0,1,0); w = vec3(0,0,1);
			u.RotateZ(gamma); v.RotateZ(gamma); w.RotateZ(gamma);
			u.RotateX(beta); v.RotateX(beta); w.RotateX(beta);
			u.RotateZ(alpha); v.RotateZ(alpha); w.RotateZ(alpha);
		}

		void ExpressInBasis(vec3* r) const
		{
			// we assume it's orthogonal
			vec3 ans;
			ans.x = u.x*r->x + u.y*r->y + u.z*r->z;
			ans.y = v.x*r->x + v.y*r->y + v.z*r->z;
			ans.z = w.x*r->x + w.y*r->y + w.z*r->z;
			*r = ans;
		}

		void ExpressInBasis(cvec3* r) const
		{
			// we assume it's orthogonal
			cvec3 ans;
			ans.x = u.x*r->x + u.y*r->y + u.z*r->z;
			ans.y = v.x*r->x + v.y*r->y + v.z*r->z;
			ans.z = w.x*r->x + w.y*r->y + w.z*r->z;
			*r = ans;
		}

		void ExpressInStdBasis(vec3* r) const
		{
			vec3 ans;
			ans.x = u.x*r->x + v.x*r->y + w.x*r->z;
			ans.y = u.y*r->x + v.y*r->y + w.y*r->z;
			ans.z = u.z*r->x + v.z*r->y + w.z*r->z;
			*r = ans;
		}

		void ExpressInStdBasis(cvec3* r) const
		{
			cvec3 ans;
			ans.x = u.x*r->x + v.x*r->y + w.x*r->z;
			ans.y = u.y*r->x + v.y*r->y + w.y*r->z;
			ans.z = u.z*r->x + v.z*r->y + w.z*r->z;
			*r = ans;
		}

		friend void Normalize(basis& b)
		{
			Normalize(b.u);
			Normalize(b.v);
			Normalize(b.w);
		}
	};

	struct vec4
	{
		// relativistic 4-vector

		tw::Float array[4];

		vec4()
		{
			array[0] = 0.0;
			array[1] = 0.0;
			array[2] = 0.0;
			array[3] = 0.0;
		}

		vec4(const tw::Float& t,const vec3& r)
		{
			array[0] = t;
			array[1] = r.x;
			array[2] = r.y;
			array[3] = r.z;
		}

		vec4(const tw::Float& a,const tw::Float& b,const tw::Float& c,const tw::Float& d)
		{
			array[0] = a;
			array[1] = b;
			array[2] = c;
			array[3] = d;
		}

		vec4(const vec4& v)
		{
			array[0] = v.array[0];
			array[1] = v.array[1];
			array[2] = v.array[2];
			array[3] = v.array[3];
		}

		vec4(const tw::Float& a)
		{
			array[0] = a;
			array[1] = a;
			array[2] = a;
			array[3] = a;
		}

		vec4(tw::Float *a)
		{
			array[0] = a[0];
			array[1] = a[1];
			array[2] = a[2];
			array[3] = a[3];
		}

		tw::vec3 spatial()
		{
			return tw::vec3(array[1],array[2],array[3]);
		}

		friend tw::Float Inner(const vec4& v1,const vec4& v2)
		{
			// Be sure to raise one of the two vectors to form the Minkowski product.
			// Otherwise we have the Euclidean inner product.
			return v1.array[0]*v2.array[0] + v1.array[1]*v2.array[1] + v1.array[2]*v2.array[2] + v1.array[3]*v2.array[3];
		}

		tw::vec4 raise_pmmm()
		{
			return tw::vec4(array[0],-array[1],-array[2],-array[3]);
		}

		tw::vec4 raise_mppp()
		{
			return tw::vec4(-array[0],array[1],array[2],array[3]);
		}

		void zBoost(const tw::Float& g,const tw::Float& sgn)
		{
			tw::vec4 v(*this);
			tw::vec4 L0(g,0.0,0.0,sgn*sqrt(g*g-1.0));
			tw::vec4 L3(sgn*sqrt(g*g-1.0),0.0,0.0,g);
			// Here L is a pure matrix, no need to raise the vector.
			array[0] = Inner(L0,v);
			array[3] = Inner(L3,v);
		}

		void Boost(const tw::vec4& gb)
		{
			// Boost in an arbitrary direction with relativistic 4-velocity gb (gamma*beta)
			tw::vec4 v(*this);
			tw::Float b2 = gb[1]*gb[1] + gb[2]*gb[2] + gb[3]*gb[3];
			auto Lij = [&] (tw::Int i,tw::Int j) { return gb[i]*gb[j]*(gb[0]-1)/(tw::small_pos+b2); };
			tw::vec4 L0(gb[0],gb[1],gb[2],gb[3]);
			tw::vec4 L1(gb[1],1+Lij(1,1),Lij(1,2),Lij(1,3));
			tw::vec4 L2(gb[2],Lij(2,1),1+Lij(2,2),Lij(2,3));
			tw::vec4 L3(gb[3],Lij(3,1),Lij(3,2),1+Lij(3,3));
			// Here L is a pure matrix, no need to raise the vector.
			array[0] = Inner(L0,v);
			array[1] = Inner(L1,v);
			array[2] = Inner(L2,v);
			array[3] = Inner(L3,v);
		}

		const tw::Float& operator [] (const tw::Int& i) const
		{
			return array[i];
		}

		tw::Float& operator [] (const tw::Int& i)
		{
			return array[i];
		}

		vec4& operator - ()
		{
			array[0] = -array[0];
			array[1] = -array[1];
			array[2] = -array[2];
			array[3] = -array[3];
			return *this;
		}

		vec4& operator = (const vec4& v)
		{
			array[0] = v.array[0];
			array[1] = v.array[1];
			array[2] = v.array[2];
			array[3] = v.array[3];
			return *this;
		}

		vec4& operator = (const tw::Float& a)
		{
			array[0] = a;
			array[1] = a;
			array[2] = a;
			array[3] = a;
			return *this;
		}

		vec4& operator += (const vec4& v)
		{
			array[0] += v.array[0];
			array[1] += v.array[1];
			array[2] += v.array[2];
			array[3] += v.array[3];
			return *this;
		}

		vec4& operator -= (const vec4& v)
		{
			array[0] -= v.array[0];
			array[1] -= v.array[1];
			array[2] -= v.array[2];
			array[3] -= v.array[3];
			return *this;
		}

		vec4& operator *= (const tw::Float& a)
		{
			array[0] *= a;
			array[1] *= a;
			array[2] *= a;
			array[3] *= a;
			return *this;
		}

		vec4& operator *= (const vec4& v)
		{
			array[0] *= v.array[0];
			array[1] *= v.array[1];
			array[2] *= v.array[2];
			array[3] *= v.array[3];
			return *this;
		}

		vec4& operator /= (const tw::Float& a)
		{
			array[0] /= a;
			array[1] /= a;
			array[2] /= a;
			array[3] /= a;
			return *this;
		}

		friend vec4 operator + (const vec4& lhs,const vec4& rhs)
		{
			return vec4(lhs.array[0] + rhs.array[0],lhs.array[1] + rhs.array[1],lhs.array[2] + rhs.array[2],lhs.array[3] + rhs.array[3]);
		}

		friend vec4 operator - (const vec4& lhs,const vec4& rhs)
		{
			return vec4(lhs.array[0] - rhs.array[0],lhs.array[1] - rhs.array[1],lhs.array[2] - rhs.array[2],lhs.array[3] - rhs.array[3]);
		}

		friend tw::Float operator ^ (const vec4& lhs,const vec4& rhs)		// Spatial Dot Product
		{
			return lhs.array[1]*rhs.array[1] + lhs.array[2]*rhs.array[2] + lhs.array[3]*rhs.array[3];
		}

		friend vec3 operator | (const vec4& lhs,const vec4& rhs)		// Spatial Cross Product
		{
			return vec3(lhs.array[2]*rhs.array[3] - lhs.array[3]*rhs.array[2],
						- lhs.array[1]*rhs.array[3] + lhs.array[3]*rhs.array[1],
						lhs.array[1]*rhs.array[2] - lhs.array[2]*rhs.array[1]);
		}

		friend vec4 operator * (const vec4& lhs,const vec4& rhs)		// Direct Product
		{
			return vec4(lhs.array[0]*rhs.array[0],lhs.array[1]*rhs.array[1],lhs.array[2]*rhs.array[2],lhs.array[3]*rhs.array[3]);
		}

		friend vec4 operator / (const vec4& lhs,const vec4& rhs)		// Direct Divide
		{
			return vec4(lhs.array[0]/rhs.array[0],lhs.array[1]/rhs.array[1],lhs.array[2]/rhs.array[2],lhs.array[3]/rhs.array[3]);
		}

		friend vec4 operator * (const vec4& v,const tw::Float& a)	// Scaling
		{
			return vec4(v.array[0]*a, v.array[1]*a, v.array[2]*a, v.array[3]*a);
		}

		friend vec4 operator * (const tw::Float& a,const vec4& v)
		{
			return vec4(v.array[0]*a, v.array[1]*a, v.array[2]*a, v.array[3]*a);
		}

		friend vec4 operator / (const vec4& v,const tw::Float& a)
		{
			return vec4(v.array[0]/a, v.array[1]/a, v.array[2]/a, v.array[3]/a);
		}

		friend std::ostream& operator << (std::ostream& os, const vec4& v)
		{
			tw::Int i;
			os << '(';
			for (i=0;i<3;i++)
				os << v.array[i] << ',';
			os << v.array[3] << ')';
			return os;
		}
	};

	struct spinor
	{
		tw::Complex array[2];

		spinor()
		{
			array[0] = 1.0;
			array[1] = 0.0;
		}

		spinor(const tw::Complex& a,const tw::Complex& b)
		{
			array[0] = a;
			array[1] = b;
		}

		spinor(const spinor& v)
		{
			array[0] = v.array[0];
			array[1] = v.array[1];
		}

		spinor(const tw::vec3& n,const tw::Float& sgn)
		{
			// Construct a helicity state
			// n = unit vector in direction of motion
			// sgn = sign determines sign of the helicity
			// Work out spherical coordinates and half-angles
			tw::Float ctheta = n.z;
			tw::Float stheta = sqrt(1.0-n.z*n.z);
			tw::Float cphi = n.x/stheta;
			tw::Float sphi = sqrt(1.0 - cphi*cphi);
			tw::Float ctheta2 = sqrt(0.5*(1+ctheta));
			tw::Float stheta2 = sqrt(0.5*(1-ctheta));
			tw::Float cphi2 = sqrt(0.5*(1+cphi));
			tw::Float sphi2 = sqrt(0.5*(1-cphi));
			if (sgn>0.0)
			{
				array[0] = ctheta2 * (cphi2-ii*sphi2);
				array[1] = stheta2 * (cphi2+ii*sphi2);
			}
			else
			{
				array[0] = -stheta2 * (cphi2-ii*sphi2);
				array[1] = ctheta2 * (cphi2+ii*sphi2);
			}
		}

		const tw::Complex& operator [] (const tw::Int& i) const
		{
			return array[i];
		}

		tw::Complex& operator [] (const tw::Int& i)
		{
			return array[i];
		}

		spinor& operator - ()
		{
			array[0] = -array[0];
			array[1] = -array[1];
			return *this;
		}

		spinor& operator = (const spinor& v)
		{
			array[0] = v.array[0];
			array[1] = v.array[1];
			return *this;
		}

		spinor& operator += (const spinor& v)
		{
			array[0] += v.array[0];
			array[1] += v.array[1];
			return *this;
		}

		spinor& operator -= (const spinor& v)
		{
			array[0] -= v.array[0];
			array[1] -= v.array[1];
			return *this;
		}

		spinor& operator *= (const tw::Complex& a)
		{
			array[0] *= a;
			array[1] *= a;
			return *this;
		}

		spinor& operator *= (const spinor& v)
		{
			array[0] *= v.array[0];
			array[1] *= v.array[1];
			return *this;
		}

		spinor& operator /= (const tw::Complex& a)
		{
			array[0] /= a;
			array[1] /= a;
			return *this;
		}

		friend spinor operator + (const spinor& lhs,const spinor& rhs)
		{
			return spinor(lhs.array[0] + rhs.array[0],lhs.array[1] + rhs.array[1]);
		}

		friend spinor operator - (const spinor& lhs,const spinor& rhs)
		{
			return spinor(lhs.array[0] - rhs.array[0],lhs.array[1] - rhs.array[1]);
		}

		friend spinor operator * (const spinor& lhs,const spinor& rhs)		// Direct Product
		{
			return spinor(lhs.array[0]*rhs.array[0],lhs.array[1]*rhs.array[1]);
		}

		friend spinor operator / (const spinor& lhs,const spinor& rhs)		// Direct Divide
		{
			return spinor(lhs.array[0]/rhs.array[0],lhs.array[1]/rhs.array[1]);
		}

		friend spinor operator * (const spinor& v,const tw::Complex& a) // phase factor
		{
			return spinor(v.array[0]*a, v.array[1]*a);
		}

		friend spinor operator * (const tw::Complex& a,const spinor& v)
		{
			return spinor(v.array[0]*a, v.array[1]*a);
		}

		friend spinor operator / (const spinor& v,const tw::Complex& a)
		{
			return spinor(v.array[0]/a, v.array[1]/a);
		}

		friend std::ostream& operator << (std::ostream& os, const spinor& v)
		{
			os << '(' << v.array[0] << ',' << v.array[1] << ')';
			return os;
		}

		void TransformAxis(tw::Int ax)
		{
			// Transform spinor components to a coordinate axis
			const tw::Complex pauli_mat[3][2][2] = {{{0.0,1.0},{1.0,0.0}}, {{0.0,-ii},{ii,0.0}}, {{1.0,0.0},{0.0,-1.0}}};
			const tw::Complex temp = array[0];
			array[0] = pauli_mat[ax][0][0]*temp + pauli_mat[ax][0][1]*array[1];
			array[1] = pauli_mat[ax][1][0]*temp + pauli_mat[ax][1][1]*array[1];
		}

		void TransformAxis(const vec3& n)
		{
			// Transform spinor components to an arbitrary axis
			spinor temp,ans(0.0,0.0);
			for (tw::Int i=1;i<=3;i++)
			{
				temp = *this;
				temp.TransformAxis(i);
				ans += n[i-1]*temp;
			}
			*this = ans;
		}
	};

	struct bispinor
	{
		tw::Complex array[4];

		bispinor()
		{
			array[0] = 1.0;
			array[1] = 0.0;
			array[2] = 0.0;
			array[3] = 0.0;
		}

		bispinor(const tw::Complex& a,const tw::Complex& b,const tw::Complex& c,const tw::Complex& d)
		{
			array[0] = a;
			array[1] = b;
			array[2] = c;
			array[3] = d;
		}

		bispinor(const bispinor& v)
		{
			array[0] = v.array[0];
			array[1] = v.array[1];
			array[2] = v.array[2];
			array[3] = v.array[3];
		}

		bispinor(const tw::spinor& w)
		{
			array[0] = w[0];
			array[1] = w[1];
			array[2] = 0.0;
			array[3] = 0.0;
		}

		void Boost(const tw::vec4& p)
		{
			tw::spinor phi(array[0],array[1]),chi(array[2],array[3]);
			tw::spinor nsigmaphi(phi),nsigmachi(chi);
			tw::vec3 vel = tw::vec3(-p[1],-p[2],-p[3])/p[0];
			tw::Float vmag = Magnitude(vel);
			tw::vec3 n = vel/vmag;
			tw::Float rot = std::atanh(vmag);
			nsigmaphi.TransformAxis(n);
			nsigmachi.TransformAxis(n);
			phi = std::cosh(0.5*rot)*phi - std::sinh(0.5*rot)*nsigmachi;
			chi = std::cosh(0.5*rot)*chi - std::sinh(0.5*rot)*nsigmaphi;
			array[0] = phi[0];
			array[1] = phi[1];
			array[2] = chi[0];
			array[3] = chi[1];
		}

		void ChargeConjugate()
		{
			tw::spinor phi(std::conj(array[0]),std::conj(array[1]));
			tw::spinor chi(std::conj(array[2]),std::conj(array[3]));
			phi.TransformAxis(2);
			chi.TransformAxis(2);
			array[0] = chi[0];
			array[1] = chi[1];
			array[2] = -phi[0];
			array[3] = -phi[1];
		}

		const tw::Complex& operator [] (const tw::Int& i) const
		{
			return array[i];
		}

		tw::Complex& operator [] (const tw::Int& i)
		{
			return array[i];
		}

		bispinor& operator - ()
		{
			array[0] = -array[0];
			array[1] = -array[1];
			array[2] = -array[2];
			array[3] = -array[3];
			return *this;
		}

		bispinor& operator = (const bispinor& v)
		{
			array[0] = v.array[0];
			array[1] = v.array[1];
			array[2] = v.array[2];
			array[3] = v.array[3];
			return *this;
		}

		bispinor& operator += (const bispinor& v)
		{
			array[0] += v.array[0];
			array[1] += v.array[1];
			array[2] += v.array[2];
			array[3] += v.array[3];
			return *this;
		}

		bispinor& operator -= (const bispinor& v)
		{
			array[0] -= v.array[0];
			array[1] -= v.array[1];
			array[2] -= v.array[2];
			array[3] -= v.array[3];
			return *this;
		}

		bispinor& operator *= (const tw::Complex& a)
		{
			array[0] *= a;
			array[1] *= a;
			array[2] *= a;
			array[3] *= a;
			return *this;
		}

		bispinor& operator *= (const bispinor& v)
		{
			array[0] *= v.array[0];
			array[1] *= v.array[1];
			array[2] *= v.array[2];
			array[3] *= v.array[3];
			return *this;
		}

		bispinor& operator /= (const tw::Complex& a)
		{
			array[0] /= a;
			array[1] /= a;
			array[2] /= a;
			array[3] /= a;
			return *this;
		}

		friend bispinor operator + (const bispinor& lhs,const bispinor& rhs)
		{
			return bispinor(lhs.array[0] + rhs.array[0],lhs.array[1] + rhs.array[1],lhs.array[2] + rhs.array[2],lhs.array[3] + rhs.array[3]);
		}

		friend bispinor operator - (const bispinor& lhs,const bispinor& rhs)
		{
			return bispinor(lhs.array[0] - rhs.array[0],lhs.array[1] - rhs.array[1],lhs.array[2] - rhs.array[2],lhs.array[3] - rhs.array[3]);
		}

		friend bispinor operator * (const bispinor& lhs,const bispinor& rhs)		// Direct Product
		{
			return bispinor(lhs.array[0]*rhs.array[0],lhs.array[1]*rhs.array[1],lhs.array[2]*rhs.array[2],lhs.array[3]*rhs.array[3]);
		}

		friend bispinor operator / (const bispinor& lhs,const bispinor& rhs)		// Direct Divide
		{
			return bispinor(lhs.array[0]/rhs.array[0],lhs.array[1]/rhs.array[1],lhs.array[2]/rhs.array[2],lhs.array[3]/rhs.array[3]);
		}

		friend bispinor operator * (const bispinor& v,const tw::Complex& a)
		{
			return bispinor(v.array[0]*a, v.array[1]*a, v.array[2]*a, v.array[3]*a);
		}

		friend bispinor operator * (const tw::Complex& a,const bispinor& v)
		{
			return bispinor(v.array[0]*a, v.array[1]*a, v.array[2]*a, v.array[3]*a);
		}

		friend bispinor operator / (const bispinor& v,const tw::Complex& a)
		{
			return bispinor(v.array[0]/a, v.array[1]/a, v.array[2]/a, v.array[3]/a);
		}

		friend std::ostream& operator << (std::ostream& os, const bispinor& v)
		{
			os << '(' << v.array[0] << ',' << v.array[1] << ',' << v.array[2] << ',' << v.array[3] << ')';
			return os;
		}

	};

}
