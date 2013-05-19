/**
 * @file
 * This file is part of SWE.
 *
 * @author Raphael Dümig
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * TODO
 */

#ifndef __SWE_ARTIFICIAL_TSUNAMI_SCENARIO_H
#define __SWE_ARTIFICIAL_TSUNAMI_SCENARIO_H

#include <cmath>

#define PI 3.1415926535897932384626433832795

// all values in this file can be interpreted as meters

class SWE_ArtificialTsunamiScenario : public SWE_Scenario {
    
public:
    
    float getWaterHeight(float x, float y) { 
        return (-1) * getOrigBathymetry(x, y);
    };
    
    float getBathymetry(float x, float y) {
        // get the original bathymetry
        float bath = getOrigBathymetry(x, y);
        
        // the rest of this function will determine, whether we need to
        // add a displacement to the original value
        
        // coordinates of the center
        float xc, yc;
        
        // calculate the center coordinates
        xc = (getBoundaryPos(BND_LEFT) + getBoundaryPos(BND_RIGHT)) / 2;
        yc = (getBoundaryPos(BND_BOTTOM) + getBoundaryPos(BND_TOP)) / 2;
        
        // add the displacement to the bathymetry, if we are close to the center
        // if the simulated area (max 500 meters distance in both directions)
        if(std::abs(x - xc) <= 500.f && std::abs(y - yc) <= 500.f)
            bath += getDisplacement(x - xc, y - yc);
        
        return bath;
    };
    
    BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };
    
    float getBoundaryPos(BoundaryEdge i_edge) {
        if ( i_edge == BND_LEFT || i_edge == BND_BOTTOM)
            return     0.f;
        else
            return 10000.f;
    };
    
    float endSimulation() { return 180.f; };
    
private:
    /**
     * getDisplacement implements the function d(x,y) from the paper of the
     * third assignment
     * 
     * @param x x-coordinate: has to be in the range (-500, 500)
     * @param y y-coordinate: has to be in the range (-500, 500)
     * 
     * @result the displacement
     */
    float getDisplacement(float x, float y) {
        float dx = std::sin((x/500.f + 1.f) * PI);
        
        float c = y/500.f;
        float dy = (-1) * c * c + 1.f;
        
        return 5.f * dx * dy;
    };
    
    /**
     * get the bathymetry before the displacement has occurred
     * 
     * @param x x-coordinate
     * @param y y-coordinate
     * 
     * @result the bathymetry before the displacement
     */
    float getOrigBathymetry(float x, float y) {
        return -100.f;
    }
};


#endif
