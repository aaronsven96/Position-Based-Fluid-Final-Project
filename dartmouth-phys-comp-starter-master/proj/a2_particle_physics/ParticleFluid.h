//#####################################################################
// Particle Fluid (SPH)
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#ifndef __ParticleFluid_h__
#define __ParticleFluid_h__
#include "Common.h"
#include <cstdlib>
#include "Particles.h"
#include "ImplicitGeometry.h"
#include "Math.h"

// offline - record positions of particles 

//////////////////////////////////////////////////////////////////////////
////Kernel function
template<int d> class Kernel
{using VectorD=Vector<real,d>;
public:
	////precomputed coefs;
	real h;
	real coef_Wspiky;
	real coef_h2;
	real coef_dWspiky;
	real coef_Wvis;
	real coef_d2Wvis;
	real coef_1poly6;
	real coef_1poly4;
	real coef_1poly2;
	real coef_9poly6;
	real pi=3.1415927;

	void Precompute_Coefs(real _h)
	{
		h=_h;
		coef_1poly6 = pow(h, 6);
		coef_h2 = pow(h, 2);
		coef_Wspiky=15.0/(pi*coef_1poly6);
		coef_dWspiky=-45.0/(pi*coef_1poly6);
		coef_Wvis=2*pi*pow(h,3);
		coef_d2Wvis=45.0/(pi*coef_1poly6);
		coef_1poly4 = 3 * pow(h, 4);
		coef_1poly2 = 3 * pow(h, 2);
		coef_9poly6 = 315/(pi*64* pow(h, 9));

	}
	////Kernel Poly6
	real Poly6(const VectorD& xji) {
		real r=xji.norm();
		if (r >= 0 && r <= h) { return (coef_9poly6 * pow(coef_h2 -pow(r, 2), 3)); }//return (coef_1poly6-coef_1poly4*pow(r,2)+coef_1poly2*pow(r,4)-pow(r,6)); }
		else { return 0; }
	}

	real Poly6_scorr(real h, real coeff) {
		real r=coeff*h;
		if (r >= 0 && r <= h) { return (coef_9poly6 * pow(pow(h, 2)-pow(r, 2), 3)); }//return (coef_1poly6-coef_1poly4*pow(r,2)+coef_1poly2*pow(r,4)-pow(r,6)); }
		else { return 0; }
	}

	////Kernel Spiky
	real Wspiky(const VectorD& xji)
	{
		real r=xji.norm();
		if(r>=0&&r<=h){return coef_Wspiky *pow(h-r,3);}
		else{return 0;}
	}

	real Wspiky_scorr(real h, real coeff)
	{
		real r=coeff*h;
		if(r>=0&&r<=h){return coef_Wspiky *pow(h-r,3);}
		else{return 0;}
	}

	VectorD gradientWspiky(const VectorD& v){
		real r=v.norm();
		if(r<= h&&r>0){return coef_dWspiky *pow(h-r,2)*v/r;}
		else{return VectorD::Zero();}
	}

	////Kernel viscosity
	real Wvis(const VectorD& xji){
		real r=xji.norm();
		if(r>=0&&r<=h){return 15.0/(2*pi*pow(h,3))*((-pow(r,3)/(2*pow(h,3))+r*r/(h*h)+h/(2*r)-1));}
		else{return 0;}
	}
	real laplacianWvis(const VectorD& v){
		real r=v.norm();
		if(r<=h&&r>0){return 45.0/(pi*pow(h,6))*(h-r);}
		else{return 0;}
	}
};

//////////////////////////////////////////////////////////////////////////
////Spatial hashing


template<int d> class SpatialHashing
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;
public:
	real dx=1.;	////grid cell size
	Hashtable<VectorDi,Array<int> > voxels;

	void Update_Voxels(const Array<VectorD>& points)
	{Clear_Voxels();for(int i=0;i<(int)points.size();i++)Add_Point(i,points[i]);}

	void Clear_Voxels(){voxels.clear();}

	bool Add_Point(const int point_idx,const VectorD& point_pos)
	{
		VectorDi cell=Cell_Coord(point_pos);
		auto iter=voxels.find(cell);
		if(iter==voxels.end())iter=voxels.insert(std::make_pair(cell,Array<int>())).first;
		Array<int>& bucket=iter->second;
		bucket.push_back(point_idx);
		return true;
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): find all the neighboring particles within the "kernel_radius" around "pos" and record their indices in "nbs", the position of the particles are given in "points"
	////You need to traverse all the 3^d neighboring cells in the background grid around the cell occupied by "pos", and then check the distance between each particle in each neighboring cell and the given "pos"
	////Use the helper function Cell_Coord to get the cell coordinates for a given "pos"
	////Use the helper function Nb_R to get the cell coordinates of the ith neighboring cell around the cell "coord"
	bool Find_Nbs(const VectorD& pos,const Array<VectorD>& points,const real kernel_radius,/*returned result*/Array<int>& nbs,Array<int>& nbs_b,Array<bool>& isboundary) const
	{
		/* Your implementation start */
		VectorDi cur= Cell_Coord(pos);
		for (int i = 0; i < pow(3, d);i++) {
			VectorDi curCell;
			curCell = Nb_R(cur, i);
			auto iter = voxels.find(curCell);
			if (iter != voxels.end()) {
				Array<int> negbors = iter->second;
				for (int j = 0; j < negbors.size(); j++) {
					int idx = negbors[j];
					double val = (pos-points[idx]).norm();
					if ((points[idx]-pos ).norm() < kernel_radius) {
						if(isboundary[idx]){
							nbs_b.push_back(idx);
						}
						else{ 
							nbs.push_back(idx); 
						}
						
					}
				}
			}
		}
		/* Your implementation end */
		return nbs.size()>0;
	}


protected:	////Helper functions
	VectorDi Cell_Coord(const VectorD& pos) const
	{VectorD coord_with_frac=(pos)/dx;return coord_with_frac.template cast<int>();}
	Vector2i Nb_R(const Vector2i& coord,const int index) const
	{assert(index>=0&&index<9);int i=index/3;int j=index%3;return coord+Vector2i(-1+i,-1+j);}
	Vector3i Nb_R(const Vector3i& coord,const int index) const
	{assert(index>=0&&index<27);int i=index/9;int m=index%9;int j=m/3;int k=m%3;return coord+Vector3i(-1+i,-1+j,-1+k);}
};

//////////////////////////////////////////////////////////////////////////
////Particle fluid simulator
template<int d> class ParticleFluid
{using VectorD=Vector<real,d>;
public:
	Particles<d> particles;
	Array<Array<int> > neighbors;
	Array<Array<int> > neighbors_boundary;
	Array<bool> is_boundary;
	
	Array<VectorD> last_positions;			////temp positions in solver]]
	Array<VectorD> delta_positions;			////change in positions in solver
	Array<real> lambda_i;					////array for lambda values
	Array<VectorD> weight_i;						////weight of the particles
	
	//Boundary Particle Volume
	Array<real> b_volume;

	//User Interaction
	Array<VectorD> flow_direction;
	bool is_flow = false;
	Array<ImplicitGeometry<d>* > flow_objects;

	SpatialHashing<d> spatial_hashing;
	
	Kernel<d> kernel;
	int solver_iterations = 15;				////solver iterations
	real kernel_radius=(real).8;			////kernel radius
	real pressure_density_coef=(real)1e1;	////pressure-density-relation coefficient, used in Update_Pressure()
	real density_0=(real)20.;				////rest density, used in Update_Pressure()
	real viscosity_coef=(real)5;		////viscosity coefficient, used in Update_Viscosity_Force()
	real kd=(real)1e2;						////stiffness for environmental collision response
	VectorD g=VectorD::Unit(1)*(real)-1.;	////gravity
	
	////realx coef
	
	real relax = 1;

	//// Max Vel
	real velocity_max = 100;

	////values used in scorr 
	real k = 0.01;
	real n = 1;
	real coeff = 0.1;
	
	real denom_coef;
	real XSPH_c = 0.01;
	
	////Environment objects
	Array<ImplicitGeometry<d>* > env_objects;

	virtual void Initialize()
	{
		if (d == 2) {
			density_0 = 10;
		}

		kernel.Precompute_Coefs(kernel_radius);
		denom_coef =kernel.Wspiky_scorr(kernel_radius, coeff);
	}

	virtual void Update_Neighbors()
	{
		
		spatial_hashing.Clear_Voxels();
		spatial_hashing.Update_Voxels(particles.XRef());
		
		is_boundary.resize(particles.Size());
		neighbors_boundary.resize(particles.Size());
		neighbors.resize(particles.Size());
		for(int i=0;i<particles.Size();i++){
			Array<int> nbs;
			Array<int> nbs_b;
			spatial_hashing.Find_Nbs(particles.X(i),particles.XRef(),kernel_radius,nbs,nbs_b, is_boundary);
			neighbors[i]=nbs;
			neighbors_boundary[i] = nbs_b;
		}
	}

	/*virtual void Update_Voracity_Force() {
		for (int i = 0; i < particles.Size(); i++) {
			particles.F(i)= VectorD::Zero();
			VectorD n = VectorD::Zero();
			for (int idx : neighbors[i]) {
				n += (weight_i[idx]- weight_i[i]).dot(kernel.gradientWspiky(particles.X(i) - particles.X(idx)));
			}
			particles.F(i) = 1 * ((n / n.norm()).cross(weight_i[i]));
		}
	}*/

	/*virtual void Update_Voracity_Weight()
	{
		for(int i=0; i<particles.Size(); i++) {
			weight_i[i] = VectorD::Zero();
			for(int idx : neighbors[i]){
				weight_i[i] += (particles.V(i) - particles.V(idx)).cross(kernel.gradientWspiky(particles.X(i) - particles.X(idx)));
			}
		}
	}*/

	virtual void Update_Viscosity() {
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				VectorD gradV = VectorD::Zero();
				for (int idx : neighbors[i]) {
					gradV += (particles.V(idx) - particles.V(i))*kernel.Wspiky(particles.X(i) - particles.X(idx));
				}
				particles.V(i) = particles.V(i) + XSPH_c * gradV;
			}
		}
	}

	virtual void Update_Flow_Velocity() {
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				for (int j = 0; j < flow_objects.size(); j++) {
					real phi = flow_objects[j]->Phi(particles.X(i));
					if (phi < particles.R(i)) {
						particles.V(i) += flow_direction[j];
					}
				}
			}
		}
	}
	virtual void Update_Pressure()
	{
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				particles.P(i) = (particles.D(i) - density_0) *pressure_density_coef;
			}
		}
	}
	virtual void Update_Boundary_Pressure_Force() {
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				for (int idx : neighbors_boundary[i]) {
					particles.F(i) -= b_volume[idx] * (particles.P(i) / pow(particles.D(i), 2))*
						kernel.gradientWspiky(particles.X(i) - particles.X(idx));
				}
			}
		}
	}
	virtual void Update_Boundary_Friction_Force() {
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				for (int idx : neighbors_boundary[i]) {
					VectorD force= b_volume[idx] * (particles.P(i) / pow(particles.D(i), 2))*
						kernel.gradientWspiky(particles.X(i) - particles.X(idx));
					particles.F(i) -= 100 *force;
					//std::cout << "Force" << force << "\n";
					/*std::cout << "Bvol" << b_volume[idx] << "\n";
					std::cout << "Pressure" << particles.P(i) << "\n";
					std::cout << "Density" << particles.D(i) << "\n";*/
				}
			}
		}
	}

	virtual void Add_Boundary_Density() {
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				for (int idx : neighbors_boundary[i]) {
					particles.D(i) += b_volume[idx] * kernel.Wspiky(particles.X(i) - particles.X(idx));
				}
			}
		}
	}

	virtual void Update_Boundary_Volume(){
		for (int i = 0; i < particles.Size(); i++) {
			if (is_boundary[i]) {
				b_volume[i] = 0;
				for (int idx : neighbors_boundary[i]) {
					b_volume[i] += kernel.Wspiky(particles.X(i)-particles.X(idx));
				}
				b_volume[i] = (1 / b_volume[i])*density_0;
			}
		}
	}

	virtual void Advance(const real dt)
	{
		b_volume.resize(particles.Size());
		last_positions.resize(particles.Size());
		delta_positions.resize(particles.Size());
		lambda_i.resize(particles.Size());
		weight_i.resize(particles.Size());
		
		Update_Neighbors();
		
		for(int i=0;i<particles.Size();i++){
			particles.F(i)=VectorD::Zero();}


		Update_Body_Force();
		Update_Boundary_Collision_Force();

		//Particle Boundary Collision
		Update_Density();
		Update_Boundary_Volume();
		Add_Boundary_Density();
		Update_Pressure();
		Update_Boundary_Pressure_Force();
		Update_Boundary_Friction_Force();


		//predict postitions
		for (int i = 0; i < particles.Size(); i++) {
			particles.V(i) += particles.F(i) *dt;
		}
		
		if(is_flow){ Update_Flow_Velocity(); }
		
		Check_Boundary_Conditions();

		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
			for (int j = 0; j < d; j++) {
				last_positions[i][j] = particles.X(i)[j];
			}
			//if(i==1)std::cout << "Before change x:" << last_positions[i][1] << " Before chnage x:" << particles.X(i)[1] << "\n";
				particles.X(i) = last_positions[i] + particles.V(i)*dt;
			}
			//if(i==1)std::cout << "before x:" << last_positions[i][1] << " After x:" << particles.X(i)[1] << "\n";
		}

		Update_Neighbors();

		//Constraint Solver 
		for (int i = 0; i < 10; i++) {
			//std::cout << i;
			Update_Density();
			Update_Lambda();
			Update_Postion_Change();
			Update_Temp_Position();
		}
		/*solver_iterations*/
		
		//update position and Velocity
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				//std::cout << "x:"<<particles.X(i)[0] << " y:" << particles.X(i)[1]
				particles.V(i) = (particles.X(i) - last_positions[i]) / dt;
				//Update_Voracity_Weight();
				//Update_Voracity_Force();
			}
		}
		Update_Viscosity();
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): update the density (particles.D(i)) of each particle based on the kernel function (Wspiky)
	
	void Update_Postion_Change() {
		
		// delta_q.resize(particles.X(i);
		//real denom = kernel.Poly6_scorr(kernel_radius); 
		// delta_q(2) = 0.3*kernel_radius;
		real s_corr = 0;
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				delta_positions[i] = VectorD::Zero();
				for (int idx : neighbors[i]) {
					//std::cout << delta_positions[i];
					//if(particles.X(i))

					//s_corr = 0;
					s_corr = -k * pow((kernel.Wspiky(particles.X(i) - particles.X(idx)) / denom_coef), n);
					// if (neighbors[i].size() <12) {
					// 	// std::cout<<neighbors[i].size()<<"\n";
					// real fract = (kernel.Wspiky(particles.X(i) - particles.X(idx)) / denom_coef);
					// s_corr = -k * pow(fract, n);
					// // if (fract>1){
					// // 	// std::cout<<fract;
					// // 	s_corr = -k*fract;
					// // 	}


					// // //std::cout << s_corr << "\n";
					// }
					//s_corr = 0;
					// std::cout << s_corr << "\n";
					//std::cout << lambda_i[i] + lambda_i[idx] << "\n";

					delta_positions[i] += (lambda_i[i] + lambda_i[idx] + s_corr) * kernel.gradientWspiky(particles.X(i) - particles.X(idx));
				}
				delta_positions[i] = delta_positions[i] / density_0;
				// std::cout << "x:" << delta_positions[i][0] << " y:" << delta_positions[i][1];
			}
		}
	}
	void Update_Lambda() {
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				real Ci = (particles.D(i) / density_0) - 1;
				real sum = 0;//VectorD::Zero();
				for (int idx2 : neighbors[i]) {
					if (i != idx2) {
						sum += pow((kernel.gradientWspiky(particles.X(i) - particles.X(idx2)) / density_0).norm(), 2);
					}
				}
				lambda_i[i] = -1 * (Ci / (sum * 2 + relax));//(sum/density_0));
				//std::cout << "lambda:" << lambda_i[i] << "  ";
			}
		}
	}

	void Update_Temp_Position() {
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				particles.X(i) = particles.X(i) + delta_positions[i];
			}
		}
	}

	void Update_Density()
	{
		/* Your implementation start */
		
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				particles.D(i) = 0;
				for (int idx : neighbors[i]) {
					particles.D(i) += kernel.Poly6(particles.X(idx) - particles.X(i));
				}
			}
		}
		/* Your implementation end */
	}

	

	void Update_Body_Force()
	{
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				particles.F(i) += particles.D(i)*g;
			}
		}
	}
	void Check_Boundary_Conditions() {
		for (int i = 0; i < particles.Size(); i++) {
			VectorD velocity = particles.V(i);
			if (velocity.norm() > velocity_max) {
				particles.V(i) = (velocity / velocity.norm()) * velocity_max;
			}
		}
	}

	void Update_Boundary_Collision_Force()
	{
		for (int i = 0; i < particles.Size(); i++) {
			if (!is_boundary[i]) {
				for (int j = 0; j < env_objects.size(); j++) {
					real phi = env_objects[j]->Phi(particles.X(i));
					if (phi < particles.R(i)) {
						//std::cout<<phi<<"\n";
						VectorD normal = env_objects[j]->Normal(particles.X(i));
						particles.F(i) += normal * kd*(particles.R(i) - phi)*particles.D(i);
					}
				}
			}
		}
	}
};

#endif
