/**
 * @author Raphael Dümig <duemig@in.tum.de>
 */

#include "blocks/SWE_DimensionalSplitting.hh"

#include <limits>
#include <cmath>
#include <cassert>

using namespace std;


SWE_DimensionalSplitting::SWE_DimensionalSplitting(int   l_nx, int   l_ny,
                                                   float l_dx, float l_dy)
  : SWE_Block(l_nx, l_ny, l_dx, l_dy),
    // initialize the 2-dimensional arrays for the net-updates
    // nx, ny (= l_nx, l_ny): defined in super-class SWE_Block
    hNetUpdatesLeft  (nx+1, ny),
    hNetUpdatesRight (nx+1, ny),
    hNetUpdatesAbove (nx, ny+1),
    hNetUpdatesBelow (nx, ny+1),
    huNetUpdatesLeft (nx+1, ny),
    huNetUpdatesRight(nx+1, ny),
    hvNetUpdatesAbove(nx, ny+1),
    hvNetUpdatesBelow(nx, ny+1)
{
    // initialize the maximum timestep
    // maxTimestep = std::numeric_limits<float>::max();
}


SWE_DimensionalSplitting::~SWE_DimensionalSplitting()
{
}


void SWE_DimensionalSplitting::simulateTimestep(float dt)
{
    // update the net-updates
    computeNumericalFluxes();
    // apply the updates
    updateUnknowns(dt);
    
    return;
}


float SWE_DimensionalSplitting::simulate(float tStart, float tEnd)
{
    float time = tStart;
    float dt;
    
    while(time < tEnd) {
        // get the ghost cells content from another block
        setGhostLayer();
        
        // simulate a timestep with dynamic length
        // (depends on the maximum wave speed)
        computeNumericalFluxes();
        
        dt = maxTimestep;
        updateUnknowns(dt);
        
        // update the time
        time += dt;
    }
    
    return time;
}

/**
 * calculate the net-updates for the current state of the fluid,
 * that can be applied later by #updateUnknowns
 */
void SWE_DimensionalSplitting::computeNumericalFluxes()
{
    float maxWaveSpeed = 0.f;
    
    // execute the f-wave solver: first horizontally
    for(int i=1; i<nx+2; i++)
        for(int j=1; j<ny+1; j++)
        {
            float maxEdgeSpeed;
            
            fWaveSolver.computeNetUpdates(h[i][j-1],  h[i][j],
                                          hu[i][j-1], hu[i][j],
                                          b[i][j-1],  b[i][j],
                                          hNetUpdatesLeft[i-1][j-1],  hNetUpdatesRight[i-1][j-1],
                                          huNetUpdatesLeft[i-1][j-1], huNetUpdatesRight[i-1][j-1],
                                          maxEdgeSpeed );
            
            // increase maxWaveSpeed if necessary
            // this is needed for calculating the maximum timestep
            maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
        }
    
    // execute the f-wave solver: now vertically
    for(int i=1; i<nx+1; i++)
        for(int j=1; j<ny+2; j++)
        {
            float maxEdgeSpeed;
            
            fWaveSolver.computeNetUpdates(h[i][j-1],  h[i][j],
                                          hv[i][j-1], hv[i][j],
                                          b[i][j-1],  b[i][j],
                                          hNetUpdatesBelow[i-1][j-1],  hNetUpdatesAbove[i-1][j-1],
                                          hvNetUpdatesBelow[i-1][j-1], hvNetUpdatesAbove[i-1][j-1],
                                          maxEdgeSpeed );
            
            // increase maxWaveSpeed if necessary
            // - the same as in the loop directly above
            maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
        }
    
    // calculate the maximum timestep that can be simulated at once from the maximum wave speed
    if(maxWaveSpeed > 0.001)
        maxTimestep = std::min(dx, dy) / maxWaveSpeed;
    else
        maxTimestep = std::numeric_limits<float>::max();
    
    assert(maxTimestep >= 0.0);
    return;
}

/**
 * apply the net-updates calculated with #computeNumericalFluxes
 */
void SWE_DimensionalSplitting::updateUnknowns(float dt)
{
    // apply the net-updates
    
    for(int i=1; i<nx+1; i++)
        for(int j=1; j<ny+1; j++) {
            // update the height of the water columns
            h[i][j] -=   dt/dx * (hNetUpdatesRight[i-1][j-1] + hNetUpdatesLeft[i][j-1])
                       + dt/dy * (hNetUpdatesAbove[i-1][j-1] + hNetUpdatesBelow[i-1][j]);
            // we cannot handle cells that got dry
            assert(h[i][j] > 0.0);
            
            // update the horizontal (hu) and vertical (hv) momentum
            hu[i][j] -= dt/dx * (huNetUpdatesRight[i-1][j-1] + huNetUpdatesLeft[i][j-1]);
            hv[i][j] -= dt/dy * (hvNetUpdatesAbove[i-1][j-1] + hvNetUpdatesBelow[i-1][j]);
        }
    
    return;
}
