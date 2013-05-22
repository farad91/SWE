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

#include "SWE_Scenario.hh"

#include <netcdf.h>

#include <iostream>
#include <cstdlib>
#include <cassert>

using namespace std;


class SWE_TsunamiScenario : public SWE_Scenario {
    
public:
    
    /**
     * load a scenario from a netCDF file
     * 
     * @param file the netCDF file to load
     */
    SWE_TsunamiScenario() {
    };
    
    // destructor
    ~SWE_TsunamiScenario() {
        nc_close(nc_id);
    };
    
    float getWaterHeight(float x, float y) { return -getPurBathymetry(x,y); };
    
    
    float getBathymetry(float x, float y) {
        int err_val;
        float result = 0.f;
        size_t index[2];        
        if(d_toGridCoordinates(x,y,&index[0], &index[1])){
            if(err_val = nc_get_var1_float(d_nc_id, d_z_id, index, &result))
                cerr <<  nc_strerror(err_val) << endl;
        }
        return result+getPurBathymetry(x,y);
    };
    
    
    float getVeloc_u(float x, float y) { return 0.0f; };
    
    
    float getVeloc_v(float x, float y) { return 0.0f; };
    
    // get boundary position TODO
    float getBoundaryPos(BoundaryEdge i_edge) {
        int err_val;
        float result;
        size_t u0 = 0;
        size_t x = x_size-1;
        size_t y = y_size-1;
        if( i_edge == BND_RIGHT ){
            if(err_val = nc_get_var1_float(nc_id, x_id, &x, &result))
                cerr <<  nc_strerror(err_val) << endl;
                result += (x_delta*0.5f);
            }
        else if( i_edge == BND_TOP ){
            if(err_val = nc_get_var1_float(nc_id, y_id, &y, &result))
                cerr <<  nc_strerror(err_val) << endl;
                result += (y_delta*0.5f);
                }
        else if( i_edge == BND_BOTTOM ){
            if(err_val = nc_get_var1_float(nc_id, y_id, &u0, &result))
                cerr <<  nc_strerror(err_val) << endl;
                result -= (y_delta*0.5f);
            }
        else if( i_edge == BND_LEFT ){
            if(err_val = nc_get_var1_float(nc_id, x_id, &u0, &result))
                cerr <<  nc_strerror(err_val) << endl;
                result -= (x_delta*0.5f);
            }
        return result;
    };
    
    BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };

       /**
     * readNetCDF will initialize the ids of the nc file and the ids of all
     * the variables which are being used
     * 
     * @param filename the name of the nc-file to be opened
     * @return 0 if successful, else the error value of the netcdf-library
     */
    int readNetCDF(const char *filename, const char *d_filename) {
        // error values will be stored in this variable
        int err_val;
        
        // number of dimensions of the variables 'x' and 'y'
        int dim_x, dim_y, d_dim_x, d_dim_y;
        
        
        // open .nc file
        if(err_val = nc_open(filename, NC_NOWRITE, &nc_id)) {
            cerr << nc_strerror(err_val) << endl;
            return err_val;
        }
        if(err_val = nc_open(d_filename, NC_NOWRITE, &d_nc_id)) {
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
        
        if(err_val = nc_inq_varid(d_nc_id, "x",  &d_x_id ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(d_nc_id, "y",  &d_y_id ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(d_nc_id, "z",  &d_z_id ))
            cerr << nc_strerror(err_val) << endl;
        
        // get the number of dimensions contained by 'x' and 'y'
        if(err_val = nc_inq_dimid(nc_id, "x", &dim_x))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_dimid(nc_id, "y", &dim_y))
            cerr << nc_strerror(err_val) << endl;
        
        if(err_val = nc_inq_dimid(d_nc_id, "x", &d_dim_x))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_dimid(d_nc_id, "y", &d_dim_y))
            cerr << nc_strerror(err_val) << endl;
        
        // get length x, y;
        nc_inq_dimlen(nc_id, dim_x, &x_size);
        nc_inq_dimlen(nc_id, dim_y, &y_size);
        
        nc_inq_dimlen(d_nc_id, d_dim_x, &d_x_size);
        nc_inq_dimlen(d_nc_id, d_dim_y, &d_y_size);
        
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
        
        if(err_val = nc_get_var1_float(d_nc_id, d_y_id, &u0, &d_y_start))
            cerr <<  nc_strerror(err_val) << endl;
        if(err_val = nc_get_var1_float(d_nc_id, d_y_id, &u1, &d_y_delta))
            cerr <<  nc_strerror(err_val) << endl;
        d_y_delta -= d_y_start;
        if(err_val = nc_get_var1_float(d_nc_id, d_x_id, &u0, &d_x_start))
            cerr <<  nc_strerror(err_val) << endl;
        if(err_val = nc_get_var1_float(d_nc_id, d_x_id, &u1, &d_x_delta))
            cerr <<  nc_strerror(err_val) << endl;
        d_x_delta -= d_x_start;
        return 0;
    }; 
    
private:
    // all variables are also declared for the displacement data and marked wit d_*    
    
    // file id
    int nc_id, d_nc_id;
    
    // variable ids
    int x_id, y_id, z_id, d_x_id, d_y_id, d_z_id ;
    
    // float x0, y0;
    size_t x_size, y_size, d_x_size, d_y_size;
    
    // float delta x[n] x[n+1];
    float x_start, x_delta, y_start, y_delta, d_x_start, d_x_delta, d_y_start, d_y_delta;
    
    // This function makes a Fitting of the x an y coordinates to the grid of the basic Bathymetry Data 
    void toGridCoordinates(float x_in, float y_in, size_t* x_out, size_t* y_out) {
        *y_out = (size_t) (((y_in-y_start)/y_delta)+0.5f); 
        *x_out = (size_t) (((x_in-x_start)/x_delta)+0.5f);
        if(*x_out<0)
            *x_out=0;
        if(*y_out < 0)
            *y_out = 0;
        if(*x_out >= x_size)
            *x_out = x_size-1;
        if(*y_out >=y_size)
            *y_out = y_size-1;
        assert(*y_out <  y_size && *x_out < x_size);
    };
    
    // This funktion checks if displacemant is relevant for the cordinates and if so it makes an fitting for them to the grid of the Displacement
    // @return bool  ( displacement avaible for this position;
    bool d_toGridCoordinates(float x_in, float y_in, size_t* x_out, size_t* y_out) {
        if((x_in > d_x_start-(d_x_delta*0.5f)) && (x_in < d_x_start + (d_x_delta*(d_x_size-0.5))) && 
            (y_in > d_y_start-(d_y_delta*0.5f)) && (y_in < d_y_start + (d_y_delta*(d_y_size-0.5)))){ 
            *y_out = (size_t) (((y_in-d_y_start)/d_y_delta)+0.5f); 
            *x_out = (size_t) (((x_in-d_x_start)/d_x_delta)+0.5f);
            return true;
        }
        else 
            return false;
    };
    
    //retuns Bathymetry data without displacement
    float getPurBathymetry(float x, float y) {
        int err_val;
        size_t index[2];
        float result;
                
        toGridCoordinates(x, y, &index[0], &index[1]);

        if(err_val = nc_get_var1_float(nc_id, z_id, index, &result))
            cerr <<  nc_strerror(err_val) << endl;
        
        return result;
    };
};
