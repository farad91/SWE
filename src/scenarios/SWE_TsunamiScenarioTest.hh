#include <cxxtest/TestSuite.h>

#include "SWE_TsunamiScenario.hh"
#define private public

class SWE_TsunamiScenarioTest : public CxxTest::TestSuite {
public:
    void testgetBoundaryPos(void){
        SWE_TsunamiScenario test;
        test.readNetCDF("test.nc");
        BoundaryEdge edge = BND_RIGHT;
        float r = test.getBoundaryPos(BND_RIGHT);
        TS_ASSERT_DELTA(r, 700 , 0.0f);
        float a = test.getBoundaryPos(BND_LEFT);
        TS_ASSERT_DELTA(a, -200, 0.0f);
        float b = test.getBoundaryPos(BND_TOP);
        TS_ASSERT_DELTA(b, -120 , 0.0f);
        float g = test.getBoundaryPos(BND_BOTTOM);
        TS_ASSERT_DELTA(g, -80 , 0.0f);   
    };
    void testgetBathymetry(void) {
        SWE_TsunamiScenario test;
        test.readNetCDF("test.nc");
        float r = test.getBathymetry(700.f,-80.f);
        TS_ASSERT_DELTA(r, 49 , 0.0f);
        r = test.getBathymetry(700.f,-81.5f);
        TS_ASSERT_DELTA(r, 49 , 0.0f); 
    };
};
