//#####################################################################
// Particle Fluid Driver
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#ifndef __ParticleFluidDriver_h__
#define __ParticleFluidDriver_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "ParticleFluid.h"

// download bunny mesh, put particles in the mesh, and make radius large enough - create another particle system
// treat in particle particle collision detection

template<int d> class ParticleFluidDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
	real dt=.02;
	ParticleFluid<d> fluid;
	Array<OpenGLSolidCircle*> opengl_circles;
	Array<OpenGLSphere*> opengl_spheres;								////spheres

	Bowl<d>* bowl=nullptr;
	Sphere<d>* sphere= nullptr;
	OpenGLSphere* moving_sphere;

	//interaction Booleans
	bool collision_cirle = false;
	bool water_flow = false;
	bool set_flow = false;

	//FlowCircles
	VectorD center_new;
	real radius_new;
	int click = 1;

	//Stop double clicking
	bool up_click = true;

public:
	virtual void Initialize()
	{
		////driver initialization, initialize simulation data
		real dx=.35;int nx=10;int ny=10;int nz=10;
		for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++){
				for(int k=0;k<nz;k++){
					VectorD pos;pos[0]=(real)i*dx-1.;pos[1]=(real)j*dx+3.;pos[2]=k*dx-1.;
					Add_Particle(pos);}}}

		bowl=new Bowl<d>(VectorD::Unit(1)*8,8);
		sphere = new Sphere<d>(VectorD::Unit(1), 2);
		fluid.env_objects.push_back(bowl);
		fluid.env_objects.push_back(sphere);

		fluid.Initialize();

		////viewer initialization, initialize visualization data
		OpenGLViewer::Initialize();
	}

	////synchronize simulation data to visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		if(bowl){
			auto opengl_circle=Add_Interactive_Object<OpenGLCircle>();
			opengl_circle->n=64;
			opengl_circle->pos=V3(bowl->center);
			opengl_circle->radius=bowl->radius;
			opengl_circle->color=OpenGLColor(1.f,.6f,.2f);
			opengl_circle->line_width=4.f;
			opengl_circle->Set_Data_Refreshed();
			opengl_circle->Initialize();
		}

		if(sphere){
			moving_sphere=Add_Interactive_Object<OpenGLSphere>();
			Set_Color(moving_sphere,OpenGLColor(.0+.05,1.,.0+.05,1.));
			moving_sphere->pos=V3(sphere->center);
			moving_sphere->radius=sphere->radius;
			moving_sphere->Set_Data_Refreshed();
			moving_sphere->Initialize();			
		}

		for(int i=0;i<fluid.particles.Size();i++){
			Add_Sphere(i);}

		// for(int i=0;i<fluid.particles.Size();i++){
		// 	Add_Solid_Circle(i);}
	}

	void Sync_Simulation_And_Visualization_Data()
	{	
		moving_sphere->pos=V3(sphere->center);
		moving_sphere->Set_Data_Refreshed();
		
		// for(int i=0;i<fluid.particles.Size();i++){
		// 	auto opengl_circle=opengl_circles[i];
		// 	opengl_circle->pos=V3(fluid.particles.X(i));
		// 	opengl_circle->Set_Data_Refreshed();}
		
		for(int i=0;i<fluid.particles.Size();i++){
			auto opengl_sphere=opengl_spheres[i];
			opengl_sphere->pos=V3(fluid.particles.X(i));
			opengl_sphere->Set_Data_Refreshed();}
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		fluid.Advance(dt);
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
	}

	virtual void Run()
	{
		OpenGLViewer::Run();
	}

	virtual void Initialize_Common_Callback_Keys()
	{
		OpenGLViewer::Initialize_Common_Callback_Keys();
		Bind_Callback_Key('a', &Keyboard_Event_A_Func, "press A");

		Bind_Callback_Key('w', &Keyboard_Event_W_Func, "press W");

		Bind_Callback_Key('s', &Keyboard_Event_S_Func, "press S");

		Bind_Callback_Key('d', &Keyboard_Event_D_Func, "press D");

		Bind_Callback_Key('q', &Keyboard_Event_Q_Func, "press Q");

		Bind_Callback_Key('e', &Keyboard_Event_E_Func, "press E");
	}

	virtual void Keyboard_Event_Q()
	{
		std::cout << "Q pressed" << std::endl;
		if (fluid.XSPH_c < .2)fluid.XSPH_c += .01;
		std::cout << fluid.XSPH_c << std::endl;

	}

	virtual void Keyboard_Event_E()
	{
		std::cout << "E pressed" << std::endl;
		if (fluid.XSPH_c > 0)fluid.XSPH_c -= .01;
		std::cout << fluid.XSPH_c << std::endl;

	}

	virtual void Keyboard_Event_A()
	{
		std::cout << "A pressed" << std::endl;
		collision_cirle = !collision_cirle;
	}

	virtual void Keyboard_Event_W()
	{
		std::cout << "W pressed" << std::endl;
		water_flow = !water_flow;
	}
	virtual void Keyboard_Event_S()
	{
		std::cout << "S pressed" << std::endl;
		set_flow = !set_flow;
		fluid.is_flow = !fluid.is_flow;
	}
	virtual void Keyboard_Event_D()
	{
		std::cout << "D pressed" << std::endl;
		if (fluid.flow_objects.size() > 0) {
			fluid.flow_objects.pop_back();
			fluid.flow_direction.pop_back();
		}
	}

	Define_Function_Object(ParticleFluidDriver, Keyboard_Event_S);
	Define_Function_Object(ParticleFluidDriver, Keyboard_Event_A);
	Define_Function_Object(ParticleFluidDriver, Keyboard_Event_W);
	Define_Function_Object(ParticleFluidDriver, Keyboard_Event_D);
	Define_Function_Object(ParticleFluidDriver, Keyboard_Event_Q);
	Define_Function_Object(ParticleFluidDriver, Keyboard_Event_E);



	virtual bool Mouse_Click(int left, int right, int mid, int x, int y, int w, int h)
	{
		up_click = !up_click;
		if (left != 1 || !set_flow || up_click) { return false; }

		Vector3f win_pos = opengl_window->Project(Vector3f::Zero());
		Vector3f pos = opengl_window->Unproject(Vector3f((float)x, (float)y, win_pos[2]));
		VectorD p_pos; for (int i = 0; i < d; i++)p_pos[i] = (real)pos[i];

		switch (click) {
		case 1: {
			std::cout << "Pick Center \n";
			center_new = p_pos;
			click = 2;
		}break;
		case 2: {
			std::cout << "Set Radius \n";
			radius_new = (center_new - p_pos).norm();
			click = 3;
		}break;
		case 3: {
			std::cout << "Pick Direction \n";
			fluid.flow_direction.push_back((p_pos - center_new) / 5);
			Sphere<d>* sphere_new = new Sphere<d>(center_new, radius_new);
			fluid.flow_objects.push_back(sphere_new);
			click = 1;
		}break;
		}
	}
	virtual bool Mouse_Drag(int x, int y, int w, int h) {
		if (!water_flow && !collision_cirle) {
			return false;
		}

		Vector3f win_pos = opengl_window->Project(Vector3f::Zero());
		Vector3f pos = opengl_window->Unproject(Vector3f((float)x, (float)y, win_pos[2]));
		VectorD p_pos; for (int i = 0; i < d; i++)p_pos[i] = (real)pos[i];

		if (collision_cirle) {
			sphere->center = p_pos;
			return true;
		}
		if (water_flow) {
			real r = (rand() % 5) - 2.5;
			Add_Particle(p_pos);
			Add_Solid_Circle(fluid.particles.Size() - 1);
			VectorD vel;
			for (int i = 0; i < d; i++) {
				vel[i] = p_pos[i];
			}
			vel[1] += 5;
			vel[0] += r;

			fluid.particles.V(fluid.particles.Size() - 1) = vel;
			return true;
		}
	}

protected:
	void Add_Particle(VectorD pos,real m=1.)
	{
		int i=fluid.particles.Add_Element();	////return the last element's index
		fluid.particles.X(i)=pos;
		fluid.particles.V(i)=VectorD::Zero();
		fluid.particles.R(i)=.1;
		fluid.particles.M(i)=m;
		fluid.particles.D(i)=1.;
	}

	void Add_Solid_Circle(const int i)
	{
		OpenGLColor c;
		for(int i=0;i<3;i++){
			c.rgba[i]=static_cast<float>(rand()%1000)/1000.f;}
		auto opengl_circle=Add_Interactive_Object<OpenGLSolidCircle>();
		opengl_circles.push_back(opengl_circle);
		opengl_circle->pos=V3(fluid.particles.X(i));
		opengl_circle->radius=fluid.particles.R(i);
		opengl_circle->color=c;
		opengl_circle->Set_Data_Refreshed();
		opengl_circle->Initialize();	
	}

	void Add_Sphere(const int i)
	{	
		OpenGLColor c;
		for(int i=0;i<3;i++){
			c.rgba[i]=static_cast<float>(rand()%1000)/1000.f;}
		OpenGLSphere* opengl_sphere=Add_Interactive_Object<OpenGLSphere>();
		opengl_sphere->pos=V3(fluid.particles.X(i));
		opengl_sphere->radius=fluid.particles.R(i);
		// Set_Color(opengl_sphere,OpenGLColor(.0+.05*(real)i,1.,.0+.05*(real)i,1.));
		opengl_sphere->color=c;
		opengl_sphere->Set_Data_Refreshed();
		opengl_sphere->Initialize();
		opengl_spheres.push_back(opengl_sphere);

	}

	////Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2){return Vector3(v2[0],v2[1],.0);}
	Vector3 V3(const Vector3& v3){return v3;}
};
#endif
