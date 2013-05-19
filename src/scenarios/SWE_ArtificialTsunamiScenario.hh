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

// all values in this file can be interpreted as meters

class SWE_ArtificialTsunamiScenario {
    
public:
    
    float getWaterHeight(float x, float y) { 
        float height = waterHeightAtRest();
        
        // coordinates of the center
        float xc, yc;
        
        // calculate the center coordinates
        xc = (getBoundaryPos(BND_LEFT) + getBoundaryPos(BND_RIGHT)) / 2;
        yc = (getBoundaryPos(BND_BOTTOM) + getBoundaryPos(BND_TOP)) / 2;
        
        // add the displacement to the water height, if we are close to the center
        // (max 500 meters distance)
        if(std::abs(x - xc <= 500.f) && std::abs(y - yc <= 500.f))
            height += getDisplacement(x - xc, y - yc);
        
        return height;
    };
    
    float getBathymetry(float x, float y) { return -100.0f; };
    
    virtual float waterHeightAtRest() { return (-1) * getBathymetry(x, y); };
    
    virtual float endSimulation() { return 0.1f; };
    
    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };
    
    virtual float getBoundaryPos(BoundaryEdge edge) {
        // values are meters
        if (edge==BND_LEFT || edge==BND_BOTTOM)
            // left and bottom boundary
            return 0.0f;
        else
            // right and top boundary
            return 10000.0f; 
    };
    
    virtual ~SWE_Scenario() {};

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
        float dx = std::sin((x/500.f + 1.f) * pi);
        
        float c = y/500.f;
        float dy = - c * c + 1.f;
        
        return 5.f * dx * dy;
    };
    
};


#endif
