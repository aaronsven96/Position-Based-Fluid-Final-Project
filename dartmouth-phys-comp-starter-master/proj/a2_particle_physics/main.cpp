//#####################################################################
// Main
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#include <iostream>
#include "ParticleFluidDriver.h"
#include "ParticleFluidDriver2D.h"
#include "ParticleFluidDriver2DExtended.h"

#ifndef __Main_cpp__
#define __Main_cpp__

int main(int argc,char* argv[])
{

	int driver=2;
	
	switch(driver){
	case 1:{
		ParticleFluidDriver2D<2> driver;
		driver.Initialize();
		driver.Run();	
	}break;
	case 2:{
		ParticleFluidDriver<3> driver;
		driver.Initialize();
		driver.Run();		
	}break;
	case 3: {
		ParticleFluidDriver2DExtended<2> driver;
		driver.Initialize();
		driver.Run();
	}break;
	}
}

#endif
