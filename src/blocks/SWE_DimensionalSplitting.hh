/**
 * @author Raphael Dümig <duemig@in.tum.de>
 */

#ifndef DIMSPLIT_H
#define DIMSPLIT_H

#include "blocks/SWE_Block.hh"
#include "tools/help.hh"
#include "solvers/FWave.hpp"

using namespace std;


class SWE_DimensionalSplitting: public SWE_Block {
private:
    solver::FWave<float> fWaveSolver;
    
    Float2D hNetUpdatesLeft;
    Float2D hNetUpdatesRight;
    Float2D hNetUpdatesAbove;
    Float2D hNetUpdatesBelow;
    Float2D huNetUpdatesLeft;
    Float2D huNetUpdatesRight;
    Float2D hvNetUpdatesAbove;
    Float2D hvNetUpdatesBelow;
public:
    // the constructor and destructor
    SWE_DimensionalSplitting(int   l_nx, int   l_ny,
                             float l_dx, float l_dy);
    virtual ~SWE_DimensionalSplitting();
    
    void simulateTimestep(float dt);
    float simulate(float tStart, float tEnd);
    void computeNumericalFluxes();
    void updateUnknowns(float dt);
    
};

#endif
