/**
 * @author Raphael Dümig <duemig@in.tum.de>
 * @author Thomas Blocher <blocher@in.tum.de>
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

 /**Simulates a timestep of dt and doesn't check if this time is out of the conditones
 */
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

    
    assert(maxTimestep > 0.01);    
    return;
}
 /** This funktion calculates and applays all changes for one Timestep 
 */
void SWE_DimensionalSplitting::runTimestep()
{
    float maxWaveSpeed = 0.f;
    float maxEdgeSpeed = 0.f;
    // execute the f-wave solver: first horizontally
    maxEdgeSpeed = computeHorizontalFluxes();
    maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
    maxTimestep = 0.4f * std::min(dx, dy) / maxWaveSpeed;
    updateHorizontal(maxTimestep);
    // execute the f-wave solver: now vertically
    maxEdgeSpeed = computeVerticalFluxes();
    maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
    updateVertical(maxTimestep);
    assert(maxTimestep > 0.0001);
    // assert just controlls y because it just can take y updates with invalide values
    assert(maxTimestep <= 0.5f * dy / maxWaveSpeed);    
    return;
}
 /** This funktion calculates and applays all changes for one Timestep with a maximum Stepwith
 * @param tmax Maximum calculation Time 
 */
void SWE_DimensionalSplitting::runTimestep(float tmax)
{
    float maxWaveSpeed = 0.f;
    float maxEdgeSpeed = 0.f;
    // execute the f-wave solver: first horizontally
    maxEdgeSpeed = computeHorizontalFluxes();
    maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
    maxTimestep = 0.4f * std::min(dx, dy) / maxWaveSpeed;
    maxTimestep = std::min(maxTimestep,tmax);
    updateHorizontal(maxTimestep);
    // execute the f-wave solver: now vertically
    maxEdgeSpeed = computeVerticalFluxes();
    maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
    updateVertical(maxTimestep);
    assert(maxTimestep > 0.0001);
    // assert just controlls y because it just can take y updates with invalide values
    assert(maxTimestep <= 0.5f * dy / maxWaveSpeed);    
    return;
}

/**This methode runs a simulation for the time intervall from tStart to tEnd
 *
 * @param tStart Start Time of the intervall
 * @param tEnd End Time of the intervall  
 */
float SWE_DimensionalSplitting::simulate(float tStart, float tEnd)
{
    float time = tStart;
    
    while(time < tEnd) {        
        // simulate a timestep with dynamic length
        // (depends on the maximum wave speed)
        runTimestep();
        
        // update the time
        time += maxTimestep;
    }
    
    return time;
}

/**
 * calculate the net-updates for the current state of the fluid, 
 * that can be applied later by #updateUnknowns 
 *
 * Important if you change maxTimestep between this function an #updateUnkowns the accurancy of the calculation is not given any more.
 */
void SWE_DimensionalSplitting::computeNumericalFluxes()
{
    float maxWaveSpeed = 0.f;
    float maxEdgeSpeed = 0.f;
    
    // execute the f-wave solver: first horizontally
    maxEdgeSpeed = computeHorizontalFluxes();
    
    // calculate the maximum timestep that can be simulated at once from the maximum wave speed
    maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
    if(maxWaveSpeed > 0.001)
        maxTimestep = 0.4f * std::min(dx, dy) / maxWaveSpeed;
    else
        maxTimestep = std::numeric_limits<float>::max();
    maxEdgeSpeed = computeVerticalFluxes(maxTimestep);
    assert(maxTimestep >= 0.0);
    assert(maxTimestep <= 0.5f * dy / maxEdgeSpeed);
    

    

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
    #pragma omp parallel
    float maxEdgeSpeed = 0.f;
    float edgeSpeed = 0.f;
    float maxInnerEdgeSpeed = 0.f;
    float u[nx+2];
    #pragma omp for 
    for(int y=0; y<ny+2; y++){
        u[0] = hu[0][y] / h[0][y];
        for(int x=0; x<nx+1; x++){
        u[x+1] = hu[x+1][y] / h[x+1][y];
            fWaveSolver.computeNetUpdates(h[x][y],  h[x+1][y],
                                          hu[x][y], hu[x+1][y],
                                          b[x][y],  b[x+1][y],
                                          u[x], u[x+1],
                                          hNetUpdatesLeft[x][y],  hNetUpdatesRight[x][y],
                                          huNetUpdatesLeft[x][y], huNetUpdatesRight[x][y],
                                          edgeSpeed);
            maxInnerEdgeSpeed = std::max(edgeSpeed, maxInnerEdgeSpeed);
        }
        maxEdgeSpeed = std::max(maxInnerEdgeSpeed, maxEdgeSpeed);
    }
    return maxEdgeSpeed;
}
/**
 * Function to compute the vertical Fluxes in dependency of the results form @computeHorizontalFluxes 
 *
 * This private function computes the vertical Fluxes for each horizontal Row excluding the Ghostcells and returns the maximum Edge speed.
 * @param dt MaxTimestep from @computeHorizontalFluxes
 *
 * @return maximum Edge speed  
 */
float SWE_DimensionalSplitting::computeVerticalFluxes(float dt){
    float maxEdgeSpeed = 0.f;
    #pragma omp parallel for schedule(dynamic,16)
    for(int x=1; x<nx+1; x++){
        float edgeSpeed = 0.f;
        float maxInnerEdgeSpeed = 0.f;
        float u[ny+2];
        u[0] = hv[x][0] / h[x][0] ;
        for(int y=0; y<ny+1; y++){
            u[y+1] = hv[x][y+1] / h[x][y+1];
            fWaveSolver.computeNetUpdates(h[x][y] -(dt/dx * (hNetUpdatesRight[x-1][y] + hNetUpdatesLeft[x][y])),
                                          h[x][y+1] -(dt/dx*(hNetUpdatesRight[x-1][y+1] + hNetUpdatesLeft[x][y+1])),
                                          hv[x][y], hv[x][y+1],
                                          b[x][y],  b[x][y+1],
                                          u[y],u[y+1],
                                          hNetUpdatesBelow[x][y],  hNetUpdatesAbove[x][y],
                                          hvNetUpdatesBelow[x][y], hvNetUpdatesAbove[x][y],
                                          edgeSpeed);
            maxInnerEdgeSpeed = std::max(edgeSpeed, maxInnerEdgeSpeed);
        }
        maxEdgeSpeed = std::max(maxInnerEdgeSpeed, maxEdgeSpeed);
    }
    return maxEdgeSpeed;
}
/**
 * Function to compute the vertical Fluxes 
 *
 * This private function computes the vertical Fluxes for each horizontal Row excluding the Ghostcells and returns the maximum Edge speed.
 *
 * @return maximum Edge speed  
 */
float SWE_DimensionalSplitting::computeVerticalFluxes(){
    #pragma omp parallel 
    float maxEdgeSpeed = 0.f;
    float edgeSpeed = 0.f;
    float maxInnerEdgeSpeed = 0.f;
    float u[ny+2];
    #pragma omp for schedule(dynamic,16)
    for(int x=1; x<nx+1; x++){
        u[0] = hv[x][0] / h[x][0];
        for(int y=0; y<ny+1; y++){
            u[y+1] = hv[x][y+1] / h[x][y+1];
            fWaveSolver.computeNetUpdates(h[x][y], h[x][y+1],
                                          hv[x][y], hv[x][y+1],
                                          b[x][y],  b[x][y+1],
                                          u[y],u[y+1],
                                          hNetUpdatesBelow[x][y],  hNetUpdatesAbove[x][y],
                                          hvNetUpdatesBelow[x][y], hvNetUpdatesAbove[x][y],
                                          edgeSpeed);
            maxInnerEdgeSpeed = std::max(edgeSpeed, maxInnerEdgeSpeed);
        }
        maxEdgeSpeed = std::max(maxInnerEdgeSpeed, maxEdgeSpeed);
    }
    return maxEdgeSpeed;
}

/**
 * apply the net-updates calculated with #computeNumericalFluxes
 * @param dt MaxTimeStep 
 */
void SWE_DimensionalSplitting::updateUnknowns(float dt)
{
    // apply the net-updates horizontal
    updateHorizontal(dt);
    // execute the f-wave solver: now vertically
    updateVertical(dt);
    return;
}
/**
 * Updates the cells for the horizontal floating Water
 * @param dt MaxTimeStep 
 */
void SWE_DimensionalSplitting::updateHorizontal(float dt){
#pragma omp parallel for schedule(dynamic,16)
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
 * @param dt MaxTimeStep 
 */
void SWE_DimensionalSplitting::updateVertical(float dt){
#pragma omp parallel for schedule(dynamic,16)
    for(int x = 1; x<nx+1; x++){
        for(int y = 1; y<ny+1; y++){
            h[x][y] -= dt/dy * (hNetUpdatesAbove[x][y-1] + hNetUpdatesBelow[x][y]);
            hv[x][y] -= dt/dy * (hvNetUpdatesAbove[x][y-1] + hvNetUpdatesBelow[x][y]);
            assert(h[x][y] > 0.0 || b[x][y] > 0.0);
        }
    }    
}
/** This funktion Updates the Bathymetry data in the SWE_Block
 *  
 * @param &scenario an reference to the scenario with the Data for the Bathymetry
 * @param time the timestemp for the requested Bathymetry   
 */
void SWE_DimensionalSplitting::updateBathymetry(SWE_Scenario &scenario, float time){
    // regulate the updat on the Displacement so the Dataflow is minimized 
    int xStart = std::max(getXpos(scenario.getBoundaryPosDispl(BND_LEFT)) - 1,1);
    int xEnd = std::min(getXpos(scenario.getBoundaryPosDispl(BND_RIGHT)) + 1,nx);
    int yStart = std::max(getYpos(scenario.getBoundaryPosDispl(BND_BOTTOM)) - 1,1);
    int yEnd = std::min(getYpos(scenario.getBoundaryPosDispl(BND_TOP)) + 1,ny);
    for(int i=xStart; i<=xEnd; i++) {
        for(int j=yStart; j<=yEnd; j++) {
          //make sure dry datapoint doesn't get wet
          if(b[i][j] < 0.)
            b[i][j] = scenario.getDynamicBathymetry( offsetX + (i-0.5f)*dx,
                                                       offsetY + (j-0.5f)*dy, time);
        }
    }
}
/**This function returns to an x-coordinate the nearest Value for x-DataPoints
 * @param x requested x-Coordinate
 * @return Position of the x-Coordinate in the Block arrays
 */
int SWE_DimensionalSplitting::getXpos(float x){
    int result = (int) ((x - offsetX)/dx);
    if(result < 0)
        result = 0;
    if(result > nx)
        result = nx;
    return result;
}
/**This function returns to an y-coordinate the nearest Value for y-DataPoints
 * @param y requested y-Coordinate
 * @return Position of the y-Coordinate in the Block arrays
 */
int SWE_DimensionalSplitting::getYpos(float y){
    int result = (int) ((y - offsetY)/dy);
    if(result < 0)
        result = 0;
    if(result > ny)
        result = ny;
    return result;
}
