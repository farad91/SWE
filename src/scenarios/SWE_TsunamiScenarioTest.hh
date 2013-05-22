#include <cxxtest/TestSuite.h>
#define private public

#include "SWE_TsunamiScenario.hh"


class SWE_TsunamiScenarioTest : public CxxTest::TestSuite {
public:
    void testgetBoundaryPos(void){
        SWE_TsunamiScenario test;
        test.readNetCDF("test.nc","d_test.nc");
        float r = test.getBoundaryPos(BND_RIGHT);
        TS_ASSERT_DELTA(r, 700 , 0.0f);
        float l = test.getBoundaryPos(BND_LEFT);
        TS_ASSERT_DELTA(l, -200, 0.0f);
        float t = test.getBoundaryPos(BND_TOP);
        TS_ASSERT_DELTA(t, -80 , 0.0f);
        float b = test.getBoundaryPos(BND_BOTTOM);
        TS_ASSERT_DELTA(b, -120 , 0.0f);   
    };
    void testgetBathymetry(void) {
        SWE_TsunamiScenario test;
        test.readNetCDF("test.nc","d_test.nc");
        float r = test.getBathymetry(700.f,-80.f);
        TS_ASSERT_DELTA(r, 49 , 0.0f);
        r = test.getBathymetry(749.9f,-84.9f);
        TS_ASSERT_DELTA(r, 49 , 0.0f);
        r = test.getBathymetry(650.1f,-75.1f);
        TS_ASSERT_DELTA(r, 49 , 0.0f); 
    };
    void testcornners(void){
    SWE_TsunamiScenario test;
        test.readNetCDF("test.nc","d_test.nc");
        float r = test.getBoundaryPos(BND_RIGHT);
        float l = test.getBoundaryPos(BND_LEFT);
        float t = test.getBoundaryPos(BND_TOP);
        float b = test.getBoundaryPos(BND_BOTTOM);
        float bath = test.getPurBathymetry(l,b);
        TS_ASSERT_DELTA(bath,0.f, 0.0f);
        bath = test.getPurBathymetry(l,t);
        TS_ASSERT_DELTA(bath,4.f, 0.0f);
        bath = test.getPurBathymetry(r,b);
        TS_ASSERT_DELTA(bath,45, 0.0f);
        bath = test.getPurBathymetry(r,t);
        TS_ASSERT_DELTA(bath,49.f, 0.0f);
    };
    
    void testpossibleScenario(void){
        SWE_TsunamiScenario test;
        test.readNetCDF("test.nc","d_test.nc");
        float r = test.getBoundaryPos(BND_RIGHT);
        float l = test.getBoundaryPos(BND_LEFT);
        float t = test.getBoundaryPos(BND_TOP);
        float b = test.getBoundaryPos(BND_BOTTOM);
        for(int x = 0; x < 10; x++){
            float in_x = l + ((r-l)/9)*x;
            int a = 0;
            if ( in_x > -115.f && in_x < 95.f)
                a = 1; 
            for(int y = 0; y < 5; y++){
                float in_y = b + ((t-b)/4)*y;
                int c = 0;
                if( in_y > -103.5f && in_y < -98.5f)
                    c = 1;
                float bath = test.getBathymetry(in_x,in_y);
                TS_ASSERT_DELTA(bath,((x*5)+y)+(a*c), 0.0f);
            }
        }
    };
    
    void testgetPurBathymetry(void){
        SWE_TsunamiScenario test;
        test.readNetCDF("test.nc","d_test.nc");
        float r = test.getBoundaryPos(BND_RIGHT);
        float l = test.getBoundaryPos(BND_LEFT);
        float t = test.getBoundaryPos(BND_TOP);
        float b = test.getBoundaryPos(BND_BOTTOM);
        for(int x = 0; x < 10; x++){
            float in_x = l + ((r-l)/9)*x;
            for(int y = 0; y < 5; y++){
                float in_y = b + ((t-b)/4)*y;
                float bath = test.getPurBathymetry(in_x,in_y);
                TS_ASSERT_DELTA(bath,((x*5)+y), 0.0f);
            }
        }
    };
    
    void testgetWaterHeight(void){
        SWE_TsunamiScenario test;
        test.readNetCDF("test.nc","d_test.nc");
        float r = test.getBoundaryPos(BND_RIGHT);
        float l = test.getBoundaryPos(BND_LEFT);
        float t = test.getBoundaryPos(BND_TOP);
        float b = test.getBoundaryPos(BND_BOTTOM);
        for(int x = 0; x < 10; x++){
            float in_x = l + ((r-l)/9)*x;
            for(int y = 0; y < 5; y++){
                float in_y = b + ((t-b)/4)*y;
                float water = test.getWaterHeight(in_x,in_y);
                TS_ASSERT_DELTA(water,-((x*5)+y), 0.0f);
            }
        }
    };
    
    void testgetVelco(void){
        SWE_TsunamiScenario test;
        test.readNetCDF("test.nc","d_test.nc");
        float r = test.getBoundaryPos(BND_RIGHT);
        float l = test.getBoundaryPos(BND_LEFT);
        float t = test.getBoundaryPos(BND_TOP);
        float b = test.getBoundaryPos(BND_BOTTOM);
        for(int x = 0; x < 10; x++){
            float in_x = l + ((r-l)/9)*x;
            for(int y = 0; y < 5; y++){
                float in_y = b + ((t-b)/4)*y;
                float hu = test.getVeloc_u(in_x,in_y);
                float hv = test.getVeloc_v(in_x,in_y);
                TS_ASSERT_DELTA(hv,0, 0.0f);
                TS_ASSERT_DELTA(hu,0, 0.0f);
            }
        }
    };
};
