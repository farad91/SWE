#include <cxxtest/TestSuite.h>

#include <SWE_TsunamiScenario.hh>
#define private public

class SWE_TsunamiScenarioTest : public CxxTest::TestSuite {
public:
    void testgetBoundaryPos(){
        SWE_Scenario test("test.nc");
        float r = test.getBoundaryPos(BND_RIGHT);
        TS_ASSERT_DELTA(r, 90 , 0.0f);    
    }; 
};
