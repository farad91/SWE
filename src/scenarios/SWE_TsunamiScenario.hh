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
    
    /**
     * load a scenario from a netCDF file
     * 
     * @param file the netCDF file to load
     */
    SWE_TsunamiScenario(char *file) {
        readNetCDF(file);
    };
    
    // destructor
    SWE_TsunamiScenario~() {
        nc_close(nc_id);
    };
    
    float getWaterHeight(float x, float y) { return -getBathymetry(x,y); };
    
    
    float getBathymetry(float x, float y) {
        int err_val;
        int index[2];
        float result;
        
        toGridCoordinates(x, y, &index[0], &index[1]);

        if(err_val = nc_get_var1_float(nc_id, z_id, index, &result))
            cerr <<  nc_strerror(err_val) << endl;
        
        return result;
    };
    
    
    float getVeloc_u(float x, float y) { return 0.0f; };
    
    
    float getVeloc_v(float x, float y) { return 0.0f; };
    
    // get boundary position TODO
    float getBoundaryPos(BoundaryEdge i_edge) {
        float result;
        if( i_edge == BND_RIGHT )
            if(err_val = nc_get_var1_float(nc_id, x_id, x_size, &result))
            cerr <<  nc_strerror(err_val) << endl;
        else if( i_edge == BND_TOP )
            if(err_val = nc_get_var1_float(nc_id, y_id, y_size, &result))
            cerr <<  nc_strerror(err_val) << endl;
        else if( i_edge == BND_BOTTOM )
            if(err_val = nc_get_var1_float(nc_id, y_id, 0, &result))
            cerr <<  nc_strerror(err_val) << endl;
        else
            if(err_val = nc_get_var1_float(nc_id, x_id, 0, &result))
            cerr <<  nc_strerror(err_val) << endl;
        return result;
    };
    
    BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };
    
    
private:
    // file id
    int nc_id;
    
    // variable ids
    int x_id, y_id, z_id;
    
    // float x0, y0;
    float x_size, y_size;
    // float delta x[n] x[n+1];
    float x_start, x_delta, y_start, y_delta;
    
    
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
        if(err_val = nc_inq_varid(nc_id, "z",  &z_id ))
            cerr << nc_strerror(err_val) << endl;
        
        // get the number of dimensions contained by 'x' and 'y'
        if(err_val = nc_inq_dimid(nc_id, "x", &dim_x))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_dimid(nc_id, "y", &dim_y))
            cerr << nc_strerror(err_val) << endl;
        
        // get length x, y;
        nc_inq_dimlen(nc_id, dim_x, &x_size);
        nc_inq_dimlen(nc_id, dim_y, &y_size);
        
        //get additional informations about the resulution of the x and y axes
        if(err_val = nc_get_var1_float(nc_id, y_id, 0, &y_start))
            cerr <<  nc_strerror(err_val) << endl;
        if(err_val = nc_get_var1_float(nc_id, y_id, 1, &y_delta))
            cerr <<  nc_strerror(err_val) << endl;
        y_delta -= y_start;
        if(err_val = nc_get_var1_float(nc_id, x_id, 0, &x_start))
            cerr <<  nc_strerror(err_val) << endl;
        if(err_val = nc_get_var1_float(nc_id, x_id, 1, &x_delta))
            cerr <<  nc_strerror(err_val) << endl;
        x_delta -= x_start;
        return 0;
    };
    
    
    void toGridCoordinates(float x_in, float y_in, int* x_out, int* y_out) {
        // TODO i think the calculation of getboundarypos is wrong but if this is ment to be right this is quit like this
        *y_out = (int) ((in_y-y_start)/y_delta)+0.5f; 
        *x_out = (int) ((x_in-x_start)/x_delta)+0.5f;
    };
};
