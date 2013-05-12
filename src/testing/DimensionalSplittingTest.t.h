#include <cxxtest/TestSuite.h>

#include "blocks/SWE_Block.hh"
#include "blocks/SWE_DimensionalSplitting.hh"
#include "testing/SWE1D/WavePropagation.h"
#include "testing/testing_scenario.hh"

#define private public


class DimensionalSplittingTest : public CxxTest::TestSuite {
public:
    static const int row  = 0;
    static const float dt = 0.01;
    static const float accuracy = 1.0E-6;
    
    static const int nx = 200;
    static const int ny = 1;
    static const float dx = 1.f;
    static const float dy = 1.f;
    
    
    void testCompareNetUpdates() {
        SWE_DimensionalSplitting dimsplit_solver(nx, ny, dx, dy);
        
        SWE_TestingScenario testingScenario;
        dimsplit_solver.initScenario(0.f, 0.f, testingScenario);

        float h[nx];
        float hu[nx];
        float b[nx];
        
        for(int i=0; i<nx; i++) {
            h[i]  = testingScenario.getWaterHeight(row, i);
            hu[i] = h[i] * testingScenario.getVeloc_u(row, i);
            b[i]  = testingScenario.getBathymetry(row, i);
        }
        
        WavePropagation solver1d(h, hu, b, nx, dx);
        
        // calculate one timestep with the dimsplit solver
	    dimsplit_solver.computeNumericalFluxes();
	    dimsplit_solver.updateUnknowns(dt);
	    
	    // calculate one timestep with the 1D solver
	    solver1d.computeNumericalFluxes();
	    solver1d.updateUnknowns(dt);
	    
	    // read the data from the dimension-splitting solver
	    Float2D dimsplit_height = dimsplit_solver.getWaterHeight();
	    Float2D dimsplit_momentum = dimsplit_solver.getDischarge_hu();
	    
	    // get the heights momentums as one-dimensional arrays
	    Float1D heights = dimsplit_height.getRowProxy(row);
	    Float1D momentums = dimsplit_momentum.getRowProxy(row);
	    
	    // compare the results
	    for(int i=0; i<nx; i++) {
	        TS_ASSERT_DELTA(heights[i], h[i], accuracy);
	        TS_ASSERT_DELTA(momentums[i], hu[i], accuracy);
	    }
    }
    
    void testHorizontalSteadyState() {
        SWE_DimensionalSplitting dimsplit_solver(nx, ny, dx, dy);
        
        SWE_TestingScenario testingScenario;
        dimsplit_solver.initScenario(0.f, 0.f, testingScenario);

        // calculate one timestep with the dimsplit solver
	    dimsplit_solver.computeNumericalFluxes();
	    dimsplit_solver.updateUnknowns(dt);
	    
	    // read the data from the dimension-splitting solver
	    Float2D dimsplit_height = dimsplit_solver.getWaterHeight();
	    Float2D dimsplit_momentum = dimsplit_solver.getDischarge_hu();
	    
	    // get the heights momentums as one-dimensional arrays
	    Float1D heights = dimsplit_height.getRowProxy(row);
	    Float1D momentums = dimsplit_momentum.getRowProxy(row);
	    
	    // compare the results
	    for(int i=0; i<nx; i++) {
	        TS_ASSERT_DELTA(heights[i], 0.f, accuracy);
	        TS_ASSERT_DELTA(momentums[i], 0.f, accuracy);
	    }
	}
	
};

