/**
 * @file
 * This file is part of SWE.
 *
 * @author Thomas Blocher, Raphael Dümig
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

#include "scenarios/SWE_Scenario.hh"

#include <netcdf.h>

#include <iostream>
#include <cstdlib>

using namespace std;


class SWE_TsunamiScenario : public SWE_Scenario {
    
public:
    
    
    float getWaterHeight(float x, float y) { 
        int err_val;
        int index[3];
        float result;
        
        toGridCoordinates(x, y, &index[0], &index[1]);
        index[2] = 0;
        
        if(err_val = nc_get_var1_float(nc_id, h_id, index, &result))
            cerr <<  nc_strerror(err_val) << endl;
        
        return result;
    };
    
    
    float getBathymetry(float x, float y) {
        int err_val;
        int index[3];
        float result;
        
        toGridCoordinates(x, y, &index[0], &index[1]);
        index[2] = 0;
        
        if(err_val = nc_get_var1_float(nc_id, b_id, index, &result))
            cerr <<  nc_strerror(err_val) << endl;
        
        return result;
    };
    
    
    float getVeloc_u(float x, float y) {
        int err_val;
        int index[3];
        float h, hu;
        
        toGridCoordinates(x, y, &index[0], &index[1]);
        index[2] = 0;
        
        // get the value of 'h'
        if(err_val = nc_get_var1_float(nc_id, h_id, index, &h))
            cerr <<  nc_strerror(err_val) << endl;
        
        // get the value of 'hu'
        if(err_val = nc_get_var1_float(nc_id, hu_id, index, &hu))
            cerr <<  nc_strerror(err_val) << endl;
        
        return hu / h;
    };
    
    
    float getVeloc_v(float x, float y) {
        int err_val;
        int index[3];
        float h, hv;
        
        toGridCoordinates(x, y, &index[2], &index[1]);
        index[0] = 0;
        
        // get the value of 'h'
        if(err_val = nc_get_var1_float(nc_id, h_id, index, &h))
            cerr <<  nc_strerror(err_val) << endl;
        
        // get the value of 'hu'
        if(err_val = nc_get_var1_float(nc_id, hv_id, index, &hv))
            cerr <<  nc_strerror(err_val) << endl;
        
        return hv / h;
    };
    
    // get boundary position
    float getBoundaryPos(BoundaryEdge i_edge) {
        if( i_edge == BND_RIGHT )
            return x_size;
        else if( i_edge == BND_TOP )
            return y_size;
        else
            // left or bottom boundary
            return 0.f;
    };
    
    BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };
    
private:
    // file id
    int nc_id;
    
    // variable ids
    int x_id, y_id, h_id, hu_id, hv_id, b_id;
    
    // float x0, y0;
    float x_size, y_size;
    
    
    /**
     * readNetCDF will initialize the ids of the nc file and the ids of all
     * the variables which are being used
     * 
     * @param filename the name of the nc-file to be opened
     * @return 0 if successful, else the error value of the netcdf-library
     */
    int readNetCDF(char *filename) {
        // error values will be stored in this variable
        int err_val;
        
        // number of dimensions of the variables 'x' and 'y'
        int dim_x, dim_y;
        
        
        // open .nc file
        if(err_val = nc_open(filename, NC_NOWRITE, &nc_id)) {
            cerr << nc_strerror(err_val) << endl;
            return err_val;
        }
        
        // get the variables IDs 
        if(err_val = nc_inq_varid(nc_id, "x",  &x_id ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(nc_id, "y",  &y_id ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(nc_id, "h",  &h_id ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(nc_id, "hu", &hu_id))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(nc_id, "hv", &hv_id))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(nc_id, "b",  &b_id ))
            cerr << nc_strerror(err_val) << endl;
        
        // get the number of dimensions contained by 'x' and 'y'
        if(err_val = nc_inq_dimid(nc_id, "x", &dim_x))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_dimid(nc_id, "y", &dim_y))
            cerr << nc_strerror(err_val) << endl;
        
        // get length x, y;
        nc_inq_dimlen(nc_id, dim_x, &x_size);
        nc_inq_dimlen(nc_id, dim_y, &y_size);
        
        // TODO: read the grid values of the x- and y-axis
        // nc_get_vara_int(nc_id, x_id, start[], count[], x_coords);
        // nc_get_vara_int(nc_id, y_id, start[], count[], y_coords);
        
        return 0;
    };
    
    
    void toGridCoordinates(float x_in, float y_in, int* x_out, int* y_out) {
        // TODO
    };
};
