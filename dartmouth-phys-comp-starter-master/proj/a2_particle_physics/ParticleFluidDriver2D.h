//#####################################################################
// Particle Fluid Driver
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#ifndef __ParticleFluidDriver2D_h__
#define __ParticleFluidDriver2D_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "Particles.h"
#include "OpenGLMesh.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLParticles.h"

// download bunny mesh, put particles in the mesh, and make radius large enough - create another particle system
// treat in particle particle collision detection

template<int d> class ParticleFluidDriver2D : public Driver, public OpenGLViewer
{
	using VectorD = Vector<real, d>; using VectorDi = Vector<int, d>; using Base = Driver;
	real dt = .02;

	//interaction Booleans
	bool collision_cirle = false;
	bool water_flow = false;
	bool set_flow = false;
	bool color = false;
	
	//FlowCircles
	VectorD center_new;
	real radius_new;
	int click=1;

	//Stop double clicking
	bool up_click = true;

	ParticleFluid<d> fluid;

	Array<OpenGLSolidCircle*> opengl_circles;

	Bowl<d>* bowl = nullptr;
	Sphere<d>* sphere = nullptr;

public:
	virtual void Initialize()
	{
		////driver initialization, initialize simulation data
		real dx = .35; int nx = 50; int ny = 50;
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				VectorD pos; pos[0] = (real)i*dx - 5.; pos[1] = (real)j*dx + 3.; pos[2] = .1*pos[1];
				Add_Particle(pos);
			}
		}

		bowl = new Bowl<d>(VectorD::Unit(1) * 8, 16);
		sphere = new Sphere<d>(VectorD::Unit(1), 2);
		sphere->center[0] += 1000;
		fluid.env_objects.push_back(bowl);
		fluid.env_objects.push_back(sphere);

		fluid.Initialize();

		////viewer initialization, initialize visualization data
		OpenGLViewer::Initialize();
	}

	////synchronize simulation data to visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		if (bowl) {
			auto opengl_circle = Add_Interactive_Object<OpenGLCircle>();
			opengl_circle->n = 64;
			opengl_circle->pos = V3(bowl->center);
			opengl_circle->radius = bowl->radius;
			opengl_circle->color = OpenGLColor(1.f, .6f, .2f);
			opengl_circle->line_width = 4.f;
			opengl_circle->Set_Data_Refreshed();
			opengl_circle->Initialize();
		}

		for (int i = 0; i < fluid.particles.Size(); i++) {
			Add_Solid_Circle(i);
		}
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		for (int i = 0; i < fluid.particles.Size(); i++) {
			auto opengl_circle = opengl_circles[i];
			opengl_circle->pos = V3(fluid.particles.X(i));
			if (color) {
				float green = fluid.particles.V(i).norm() / 5;
				OpenGLColor c((green - 1) / 3, green, 2 - green, 1.f);
				opengl_circle->color = c;
			}
			else {
				OpenGLColor c(.5f, .5f,.5f, 1);
				opengl_circle->color = c;
			}
			
			opengl_circle->Set_Data_Refreshed();
		}
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

	////User interaction
	virtual void Initialize_Common_Callback_Keys()
	{
		OpenGLViewer::Initialize_Common_Callback_Keys();
		Bind_Callback_Key('a', &Keyboard_Event_A_Func, "Move Ball");

		Bind_Callback_Key('w', &Keyboard_Event_W_Func, "Add water");

		Bind_Callback_Key('s', &Keyboard_Event_S_Func, "Add flow");

		Bind_Callback_Key('d', &Keyboard_Event_D_Func, "Del flow");

		Bind_Callback_Key('q', &Keyboard_Event_Q_Func, "Add vicosity");

		Bind_Callback_Key('e', &Keyboard_Event_E_Func, "lower viscosity");

		Bind_Callback_Key('c', &Keyboard_Event_C_Func, "Add color");

		Bind_Callback_Key('z', &Keyboard_Event_Z_Func, "stop scorr");
	}

	virtual void Keyboard_Event_Z()
	{
		std::cout << "Z pressed" << std::endl;
		fluid.do_s_corr = !fluid.do_s_corr;
	}

	virtual void Keyboard_Event_C()
	{
		std::cout << "C pressed" << std::endl;
		color = !color;
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
		if(fluid.XSPH_c>0)fluid.XSPH_c -= .01;
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
		set_flow= !set_flow;
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

	Define_Function_Object(ParticleFluidDriver2D, Keyboard_Event_Z);
	Define_Function_Object(ParticleFluidDriver2D, Keyboard_Event_C);
	Define_Function_Object(ParticleFluidDriver2D, Keyboard_Event_S);
	Define_Function_Object(ParticleFluidDriver2D, Keyboard_Event_A);
	Define_Function_Object(ParticleFluidDriver2D, Keyboard_Event_W);
	Define_Function_Object(ParticleFluidDriver2D, Keyboard_Event_D);
	Define_Function_Object(ParticleFluidDriver2D, Keyboard_Event_Q);
	Define_Function_Object(ParticleFluidDriver2D, Keyboard_Event_E);

	

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
				fluid.flow_direction.push_back((p_pos - center_new)/5);
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
			real r = (rand() % 5)-2.5;
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
	void Add_Particle(VectorD pos, real m = 1.)
	{
		int i = fluid.particles.Add_Element();	////return the last element's index
		fluid.particles.X(i) = pos;
		fluid.particles.V(i) = VectorD::Zero();
		fluid.particles.R(i) = .1;
		fluid.particles.M(i) = m;
		fluid.particles.D(i) = 1.;
	}

	void Add_Solid_Circle(const int i)
	{
		OpenGLColor c(0.5f, 0.5f, 0.5f, 1.f);
		auto opengl_circle = Add_Interactive_Object<OpenGLSolidCircle>();
		opengl_circles.push_back(opengl_circle);
		opengl_circle->pos = V3(fluid.particles.X(i));
		opengl_circle->radius = fluid.particles.R(i);
		opengl_circle->color = c;
		opengl_circle->Set_Data_Refreshed();
		opengl_circle->Initialize();
	}


	////Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2) { return Vector3(v2[0], v2[1], .0); }
	Vector3 V3(const Vector3& v3) { return v3; }
};
#endif