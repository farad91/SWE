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

#include "SWE_NetCDFScenario.hh"

#include <netcdf.h>

#include <iostream>
#include <cstdlib>

using namespace std;


class SWE_NetCDFCheckpointScenario : public SWE_NetCDFScenario {
    
public:
    
    /**
     * load a scenario from a netCDF file
     * 
     * @param file the netCDF file to load
     */
    SWE_NetCDFCheckpointScenario() {
    };
    
    // destructor
    ~SWE_NetCDFCheckpointScenario() {
        nc_close(nc_id);
    };
    
    /**
     * readNetCDF will initialize the ids of the nc file and the ids of all
     * the variables which are being used
     * 
     * @param data_file the name of the nc-file to be opened
     * @param CPFile filename of the checkpoint file
     * @return 0 if successfull, else the error value of the netcdf-library
     */ 
    int readNetCDF(const char *data_file, const char *CPFile) {
        // error values will be stored in this variable
        int err_val;
        
        // number of dimensions of the variables 'x','y' and Time
        int dim_x, dim_y, dim_time;
        
        
        // open .nc file
        if(err_val = nc_open(data_file, NC_WRITE, &nc_id)) {
            cerr << nc_strerror(err_val) << endl;
            return err_val;
        }
        if(err_val = nc_open(CPFile, NC_NOWRITE, &cp_id)) {
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
        if(err_val = nc_inq_varid(nc_id, "time", &time_id ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(cp_id, "EndTime", &EndTime_id ))
            cerr << nc_strerror(err_val) << endl;
        
        // get the number of dimensions contained by 'x' and 'y'
        if(err_val = nc_inq_dimid(nc_id, "x", &dim_x))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_dimid(nc_id, "y", &dim_y))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_dimid(nc_id, "time", &dim_time))
            cerr << nc_strerror(err_val) << endl;
            
        //get additional informations about the resulution of the x and y axes
        size_t u0 = 0;
        size_t u1 = 1;
        if(err_val = nc_get_var1_float(nc_id, y_id, &u0, &y_start))
            cerr <<  nc_strerror(err_val) << endl;
        if(err_val = nc_get_var1_float(nc_id, y_id, &u1, &y_delta))
            cerr <<  nc_strerror(err_val) << endl;
        y_delta -= y_start;
        if(err_val = nc_get_var1_float(nc_id, x_id, &u0, &x_start))
            cerr <<  nc_strerror(err_val) << endl;
        if(err_val = nc_get_var1_float(nc_id, x_id, &u1, &x_delta))
            cerr <<  nc_strerror(err_val) << endl;
        x_delta -= x_start;
        
        // get length x, y;
        nc_inq_dimlen(nc_id, dim_x, &x_size);
        nc_inq_dimlen(nc_id, dim_y, &y_size);
        nc_inq_dimlen(nc_id, dim_time, &CP_Number);
        CP_Number -= 1;   
        return 0;
    };
    
    float getWaterHeight(float x, float y) { 
        int err_val;
        size_t index[3];
        float result;
        
        toGridCoordinates(x, y, &index[1], &index[2]);
        index[0] = CP_Number;
        
        err_val = nc_get_var1_float(nc_id, h_id, index, &result);
        if( err_val )
            cerr <<  nc_strerror(err_val) << endl;
        
        return result;
    };
    
    
    float getBathymetry(float x, float y) {
        int err_val;
        size_t index[3];
        float result;
        
        toGridCoordinates(x, y, &index[1], &index[2]);
        index[0] = CP_Number;
        err_val = nc_get_var1_float(nc_id, b_id, index, &result);
        if( err_val )
            cerr <<  nc_strerror(err_val) << endl;
        
        return result;
    };
    
    
    float getVeloc_u(float x, float y) {
        int err_val;
        size_t index[3];
        float h, hu;
        
        toGridCoordinates(x, y, &index[1], &index[2]);
        index[0] = CP_Number;
        
        // get the value of 'h'
        err_val = nc_get_var1_float(nc_id, h_id, index, &h);
        if( err_val )
            cerr <<  nc_strerror(err_val) << endl;
        
        // get the value of 'hu'
        err_val = nc_get_var1_float(nc_id, hu_id, index, &hu);
        if( err_val )
            cerr <<  nc_strerror(err_val) << endl;
        
        return hu/(h*h);
    };
    
    
    float getVeloc_v(float x, float y) {
        int err_val;
        size_t index[3];
        float h, hv;
        
        toGridCoordinates(x, y, &index[2], &index[1]);
        index[0] = CP_Number;
        
        // get the value of 'h'
        err_val = nc_get_var1_float(nc_id, h_id, index, &h);
        if( err_val )
            cerr <<  nc_strerror(err_val) << endl;
        
        // get the value of 'hu'
        err_val = nc_get_var1_float(nc_id, hv_id, index, &hv);
        if( err_val )
            cerr <<  nc_strerror(err_val) << endl;
        
        return hv/(h*h);
    };
    
    float getTime() {
        int err_val;
        float time;
        
        err_val = nc_get_var1_float(nc_id, time_id, &CP_Number, &time);
        if( err_val )
            cerr << nc_strerror(err_val) << endl;
        
        return time;
    }
    
    // get boundary position  
    float getBoundaryPos(BoundaryEdge i_edge) {
        int err_val;
        float ret;
        size_t x1 = x_size-1;
        size_t y1 = y_size-1;
        size_t u0 = 0;
        if( i_edge == BND_RIGHT ) {
            if(err_val = nc_get_var1_float(nc_id, x_id, &x1 , &ret))
                cerr << nc_strerror(err_val) << endl;
            return ret+x_delta/2;
        }
        else if( i_edge == BND_TOP ) {
            if(err_val = nc_get_var1_float(nc_id, y_id, &y1, &ret))
                cerr << nc_strerror(err_val) << endl;
            return ret+y_delta/2;
        }
        if( i_edge == BND_LEFT ) {
            if(err_val = nc_get_var1_float(nc_id, x_id, &u0 , &ret))
                cerr << nc_strerror(err_val) << endl;
            return ret-x_delta/2;
        }
        else if( i_edge == BND_BOTTOM ) {
            if(err_val = nc_get_var1_float(nc_id, y_id, &u0, &ret))
                cerr << nc_strerror(err_val) << endl;
            return ret-y_delta/2;
        }
        return 0;
        
    };
    
    float endSimulation() {
        int err_val;
        float ret;
        
        err_val = nc_get_var1_float(cp_id, EndTime_id, 0, &ret);
        if( err_val )
            cerr << nc_strerror(err_val) << endl;
        
        return ret;    
    }
    
    //TODO
    BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };
    
    
private:
    // Check point Number
    size_t CP_Number;
    // file id
    int nc_id, cp_id;
    
    // variable ids
    int x_id, y_id, h_id, hu_id, hv_id, b_id, time_id, Bound_id, EndTime_id;
    
    // needed to make tho fitting of the requested coordinates to the grid
    float x_start, x_delta, y_start, y_delta;
    
    // float x0, y0;
    size_t x_size, y_size;
    
    /**
     * toGridCoordinates calculates the indices for accessing the x and y coordinates of the data file
     * 
     * @param x_coord x-coordinate
     * @param y_coord y-coordinate
     * @param x_index place to write the index for the x-axis
     * @param y_index place to write the index for the y-axis
     * @return true if the coordinate is inside the boundaries of the data, false otherwise
     */
bool toGridCoordinates(float x_coord, float y_coord,
                       size_t* x_index, size_t* y_index) {
        
        bool coordinate_valid = true;
        
        *y_index = (size_t) ( ( (y_coord - y_start) / y_delta ) + 0.5f ); 
        *x_index = (size_t) ( ( (x_coord - x_start) / x_delta ) + 0.5f );
        
        
        // check if any of the indices is outside of the boundaries
        // and set output to the nearest legal value
        if(*x_index < 0){
            *x_index = 0;
            coordinate_valid = false;
        }
        if(*x_index >= x_size){
            *x_index = x_size-1;
            coordinate_valid = false;
        }
        if(*y_index < 0){
            *y_index = 0;
            coordinate_valid = false;
        }
        if(*y_index >= y_size){
            *y_index = y_size-1;
            coordinate_valid = false;
        }
        
        
        return coordinate_valid;
    };
};
