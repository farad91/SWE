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
    hNetUpdatesLeft  (nx+2, ny+2),
    hNetUpdatesRight (nx+2, ny+2),
    hNetUpdatesAbove (nx+2, ny+2),
    hNetUpdatesBelow (nx+2, ny+2),
    huNetUpdatesLeft (nx+2, ny+2),
    huNetUpdatesRight(nx+2, ny+2),
    hvNetUpdatesAbove(nx+2, ny+2),
    hvNetUpdatesBelow(nx+2, ny+2)
{
    // initialize the maximum timestep
    // maxTimestep = std::numeric_limits<float>::max();
}


SWE_DimensionalSplitting::~SWE_DimensionalSplitting()
{
}


void SWE_DimensionalSplitting::simulateTimestep(float dt)
{
    float maxWaveSpeed = 0.f;
    float maxEdgeSpeed = 0.f;
    // execute the f-wave solver: first horizontally
    maxEdgeSpeed = computeHorizontalFluxes();
    maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
    updateHorizontal(dt);
    // execute the f-wave solver: now vertically
    maxEdgeSpeed = computeVerticalFluxes();
    maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
    updateVertical(dt);
    // calculate the maximum timestep that can be simulated at once from the maximum wave speed
    //if(maxWaveSpeed > 0.0000000000001)
        maxTimestep = std::min(dx, dy) / maxWaveSpeed;
    //else
      //  maxTimestep = std::numeric_limits<float>::max();
    
    assert(maxTimestep > 0.01);    
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
        dt = maxTimestep;
        simulateTimestep(dt);
        //computeNumericalFluxes();
        //updateUnknowns(dt);
        
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
    float maxEdgeSpeed = 0.f;
    // execute the f-wave solver: first horizontally
    maxEdgeSpeed = computeHorizontalFluxes();
    maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
        if(maxWaveSpeed > 0.001)
        maxTimestep = std::min(dx, dy) / maxWaveSpeed;
    else
        maxTimestep = std::numeric_limits<float>::max();
    assert(maxTimestep >= 0.0);
    updateHorizontal(maxTimestep);
    // execute the f-wave solver: now vertically
    maxEdgeSpeed = computeVerticalFluxes();
    maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
    // calculate the maximum timestep that can be simulated at once from the maximum wave speed

    

    return;
}
/**
 * Function to compute the horizontal Fluxes 
 *
 * This private function computes the horizontal Fluxes for each Vertical Row including the Ghostcells and returns the maximum Edge speed.
 *
 *@return maximum Edge speed 
 */
float SWE_DimensionalSplitting::computeHorizontalFluxes(){
    float edgeSpeed = 0.f;
    float maxEdgeSpeed = 0.f;
    for(int y=0; y<ny+2; y++){
        for(int x=0; x<nx+1; x++){
            fWaveSolver.computeNetUpdates(h[x][y],  h[x+1][y],
                                          hu[x][y], hu[x+1][y],
                                          b[x][y],  b[x+1][y],
                                          hNetUpdatesLeft[x][y],  hNetUpdatesRight[x][y],
                                          huNetUpdatesLeft[x][y], huNetUpdatesRight[x][y],
                                          edgeSpeed);
            maxEdgeSpeed = std::max(edgeSpeed, maxEdgeSpeed);
        }
    }
    return maxEdgeSpeed;
}
/**
 * Function to compute the vertical Fluxes
 *
 * This private function computes the vertical Fluxes for each horizontal Row excluding the Ghostcells and returns the maximum Edge speed.
 *
 *@return maximum Edge speed  
 */
float SWE_DimensionalSplitting::computeVerticalFluxes(){
    float edgeSpeed = 0.f;
    float maxEdgeSpeed = 0.f;
    for(int x=0; x<nx+2; x++){
        for(int y=0; y<ny+1; y++){
            fWaveSolver.computeNetUpdates(h[x][y],  h[x][y+1],
                                          hv[x][y], hv[x][y+1],
                                          b[x][y],  b[x][y+1],
                                          hNetUpdatesBelow[x][y],  hNetUpdatesAbove[x][y],
                                          hvNetUpdatesBelow[x][y], hvNetUpdatesAbove[x][y],
                                          edgeSpeed);
            maxEdgeSpeed = std::max(edgeSpeed, maxEdgeSpeed);
        }
    }
    return maxEdgeSpeed;
}

/**
 * apply the net-updates calculated with #computeNumericalFluxes
 */
void SWE_DimensionalSplitting::updateUnknowns(float dt)
{
    // apply the net-updates
    //updateHorizontal(dt);
    updateVertical(dt);
    return;
}
/**
 * Updates the cells for the horizontal floating Water
 */
void SWE_DimensionalSplitting::updateHorizontal(float dt){
    for(int y = 0; y<ny+2; y++){
        for(int x = 1; x<nx+1; x++){
            h[x][y] -= dt/dx * (hNetUpdatesRight[x-1][y] + hNetUpdatesLeft[x][y]);
            hu[x][y] -= dt/dx * (huNetUpdatesRight[x-1][y] + huNetUpdatesLeft[x][y]);
            assert(h[x][y] > 0.0 || b[x][y] > 0.0);
        }
    }    
}
/**
 * Updates the cells for the vertical floating Water
 */
void SWE_DimensionalSplitting::updateVertical(float dt){
    for(int x = 1; x<nx+1; x++){
        for(int y = 1; y<ny+1; y++){
            h[x][y] -= dt/dy * (hNetUpdatesAbove[x][y-1] + hNetUpdatesBelow[x][y]);
            hv[x][y] -= dt/dy * (hvNetUpdatesAbove[x][y-1] + hvNetUpdatesBelow[x][y]);
            assert(h[x][y] > 0.0 || b[x][y] > 0.0);
        }
    }    
}
