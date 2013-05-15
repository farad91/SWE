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


class SWE_ArtificialTsunamiScenario {
    
public:
    
    float getWaterHeight(float x, float y) { return (-1) * getBathymetry(x, y); };
    float getBathymetry(float x, float y) { return 0.0f; };
    
    virtual float waterHeightAtRest() { return 10.0f; };
    
    virtual float endSimulation() { return 0.1f; };
    
    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };
    virtual float getBoundaryPos(BoundaryEdge edge) {
        if (edge==BND_LEFT || edge==BND_BOTTOM)
            return 0.0f;
        else
            return 1.0f; 
    };
    
    virtual ~SWE_Scenario() {};
    
};


#endif
