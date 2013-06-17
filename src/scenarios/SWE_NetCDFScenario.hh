
#ifndef __SWE_NETCDFSCENARIO_H
#define __SWE_NETCDFSCENARIO_H
#include "SWE_Scenario.hh"

class SWE_NetCDFScenario : public SWE_Scenario {

 public:

    virtual float getWaterHeight(float x, float y) { return 10.0f; };
    virtual float getVeloc_u(float x, float y) { return 0.0f; };
    virtual float getVeloc_v(float x, float y) { return 0.0f; };
    virtual float getBathymetry(float x, float y) { return 0.0f; };
    virtual int readNetCDF(const char *file_bathy, const char *file_displ) { return 0;} ;
    virtual float waterHeightAtRest() { return 10.0f; };
    virtual float getDynamicBathymetry(float x, float y, float time) { return 0.0f; };
    virtual float getEruptionDuration() { return 0.0f; };
    virtual float getEruptionResolution() { return 0.0f; };
    
    virtual float endSimulation() { return 0.1f; };
    
    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };
    virtual float getBoundaryPos(BoundaryEdge edge) {
       if (edge==BND_LEFT || edge==BND_BOTTOM)
          return 0.0f;
       else
          return 1.0f; 
    };
    
    virtual ~SWE_NetCDFScenario() {};

};


#endif
