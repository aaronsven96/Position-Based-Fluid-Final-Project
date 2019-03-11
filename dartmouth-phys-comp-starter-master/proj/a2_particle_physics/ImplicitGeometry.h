#ifndef __ImplicitGeometry_h__
#define __ImplicitGeometry_h__

template<int d> class ImplicitGeometry
{using VectorD=Vector<real,d>;
public:
	virtual real Phi(const VectorD& pos) const {return 0.;}
	virtual VectorD Normal(const VectorD& pos) const {return VectorD::Zero();}
};

template<int d> class Bowl : public ImplicitGeometry<d>
{using VectorD=Vector<real,d>;
public:
	VectorD center;
	real radius;
	Bowl(VectorD _center=VectorD::Zero(),real _radius=1.):center(_center),radius(_radius){}
	// box - min of (distance to each line) - has direction (within box is neg, outside is pos, boundary is zero)
	virtual real Phi(const VectorD& pos) const {return radius-(pos-center).norm();}
	virtual VectorD Normal(const VectorD& pos) const {return (center-pos).normalized();}
};

template<int d> class ParticleSphere<d>: public ImplicitGeometry<d>
{using VectorD=Vector<real,d>;
public:
	VectorD center;
	real radius;
	Sphere(VectorD _center=VectorD::Zero(),real _radius=1.):center(_center),radius(_radius){}
	virtual real Phi(const VectorD& pos) const {return (pos-center).norm()-radius;}
	virtual VectorD Normal(const VectorD& pos) const {return (pos-center).normalized();}

	virtual VectorD<> getParticleArray(double rad){



		double arc = asin(rad/radius)*2 * radius;
		real pi=3.1415927;
		int n=floor((2*pi*radius)/arc);
        
		VectorD boundpos<n>;


        for(int i=0;i<n;i++){
            real theta=2.pi(real)i/(real)n;
            VectorD pos=VectorD::Zero();
            pos[0]=particles.X(0)[0]+rcos(theta);
            pos[1]=particles.X(0)[1]+rsin(theta);
            boundpos.push_back(pos);}

            return boundpos;
	}
};

// template<int d> class Cube : public ImplicitGeometry<d>
// {using VectorD=Vector<real,d>;
// public:
// 	VectorD corners;
// 	real side_length;
// 	Cube (VectorD _center=VectorD::Zero(),real _radius=1.):center(_center),radius(_radius){}
// 	// box - min of (distance to each line) - has direction (within box is neg, outside is pos, boundary is zero)
// 	virtual real Phi(const VectorD& pos) const {return radius-(pos-center).norm();}
// 	virtual VectorD Normal(const VectorD& pos) const {return (center-pos).normalized();}
// };

template<int d> class Sphere : public ImplicitGeometry<d>
{using VectorD=Vector<real,d>;
public:
	VectorD center;
	real radius;
	Sphere(VectorD _center=VectorD::Zero(),real _radius=1.):center(_center),radius(_radius){}
	virtual real Phi(const VectorD& pos) const {return (pos-center).norm()-radius;}
	virtual VectorD Normal(const VectorD& pos) const {return (pos-center).normalized();}
};

#endif
