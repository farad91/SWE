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
#include <cassert>

using namespace std;

typedef enum { BATHYMETRY, DISPLACEMENT } DataSource;


class SWE_TsunamiScenario : public SWE_NetCDFScenario {
    
public:
    
    /**
     * minimum elevation of the landmass and minimum depth of the water
     * This is needed for simulating the coastlines realistically.
     */
    const float bath_min_zero_offset;
    
    
    SWE_TsunamiScenario() : bath_min_zero_offset(20.f)
    {};
    
    // destructor
    ~SWE_TsunamiScenario() {
        nc_close(ncid_bathy);
        nc_close(ncid_displ);
    };
    
    
    float getWaterHeight(float x, float y) {
        float bath = getOriginalBathymetry(x,y);
        
        // apply the minimum height
        if( bath > bath_min_zero_offset )
            return bath_min_zero_offset;
        else
            return (-1) * bath;
    };
    
    
    float getBathymetry(float x, float y) {
        int   err_val;
        float displacement = 0.f;
        float result;
        
        size_t index[2];        
        
        // try to get the indices for the displacement data
        if( !toGridCoordinates(DISPLACEMENT, x ,y, &index[0], &index[1]) ) {
            // we are outside the area where the displacement is relevant
            result = getOriginalBathymetry(x,y);
        }
        else {
            // get the displacement
            err_val = nc_get_var1_float(ncid_displ, z_id_displ, index, &displacement);
            if( err_val )
                cerr << nc_strerror(err_val) << endl;
            
            result = getOriginalBathymetry(x,y) + displacement;
        }
        
        // apply the minimum elevation or depth
        if( result > (-1) * bath_min_zero_offset && result < 0)
            result = (-1) * bath_min_zero_offset;
        else if( result < bath_min_zero_offset && result >= 0)
            result = bath_min_zero_offset;
        
        return result;
    };
    
    
    float getVeloc_u(float x, float y) { return 0.0f; };
    
    
    float getVeloc_v(float x, float y) { return 0.0f; };
    
    // get boundary position TODO
    float getBoundaryPos(BoundaryEdge i_edge) {
        //! error values from the netcdf calls
        int err_val;
        
        float result;
        
        size_t u0 = 0;
        size_t x = x_size_bathy-1;
        size_t y = y_size_bathy-1;
        
        if( i_edge == BND_RIGHT ) {
            err_val = nc_get_var1_float(ncid_bathy, x_id_bathy, &x, &result);
            if( err_val )
                cerr <<  nc_strerror(err_val) << endl;
            
            result += (x_delta_bathy*0.5f);
        }
        else if( i_edge == BND_TOP ) {
            err_val = nc_get_var1_float(ncid_bathy, y_id_bathy, &y, &result);
            if( err_val )
                cerr <<  nc_strerror(err_val) << endl;
            
            result += (y_delta_bathy*0.5f);
        }
        else if( i_edge == BND_BOTTOM ) {
            err_val = nc_get_var1_float(ncid_bathy, y_id_bathy, &u0, &result);  
            if( err_val )
                cerr <<  nc_strerror(err_val) << endl;
            
            result -= (y_delta_bathy*0.5f);
        }
        else if( i_edge == BND_LEFT ) {
            err_val = nc_get_var1_float(ncid_bathy, x_id_bathy, &u0, &result);
            if( err_val )
                cerr <<  nc_strerror(err_val) << endl;
            
            result -= (x_delta_bathy*0.5f);
        }
        
        return result;
    };
    
    BoundaryType getBoundaryType(BoundaryEdge edge) { return WALL; };

    /**
     * readNetCDF will initialize the ids of the nc file and the ids of all
     * the variables which are being used
     * 
     * @param file_bathy the name of the nc-file containing the bathymetry
     * @param file_displ the name of the nc-file containing the displacements
     * @return 0 if successful, else the error value of the netcdf-library
     */
    int readNetCDF(const char *file_bathy, const char *file_displ) {
        // error values will be stored in this variable
        int err_val;
        
        // number of dimensions of the variables 'x' and 'y'
        int dim_x_bathy, dim_y_bathy, dim_x_displ, dim_y_displ;
        
        
        // open the bathymetry file
        if(err_val = nc_open(file_bathy, NC_NOWRITE, &ncid_bathy)) {
            cerr << nc_strerror(err_val) << endl;
            return err_val;
        }
        
        // open the displacements file
        if(err_val = nc_open(file_displ, NC_NOWRITE, &ncid_displ)) {
            cerr << nc_strerror(err_val) << endl;
            return err_val;
        }
        
        // get the variables IDs 
        if(err_val = nc_inq_varid(ncid_bathy, "x", &x_id_bathy ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(ncid_bathy, "y", &y_id_bathy ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(ncid_bathy, "z", &z_id_bathy ))
            cerr << nc_strerror(err_val) << endl;
        
        if(err_val = nc_inq_varid(ncid_displ, "x", &x_id_displ ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(ncid_displ, "y", &y_id_displ ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_varid(ncid_displ, "z", &z_id_displ ))
            cerr << nc_strerror(err_val) << endl;
        
        // get the number of dimensions contained by 'x' and 'y'
        if(err_val = nc_inq_dimid(ncid_bathy, "x", &dim_x_bathy))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_dimid(ncid_bathy, "y", &dim_y_bathy))
            cerr << nc_strerror(err_val) << endl;
        
        if(err_val = nc_inq_dimid(ncid_displ, "x", &dim_x_displ))
            cerr << nc_strerror(err_val) << endl;
        if(err_val = nc_inq_dimid(ncid_displ, "y", &dim_y_displ))
            cerr << nc_strerror(err_val) << endl;
        
        // get length x, y;
        nc_inq_dimlen(ncid_bathy, dim_x_bathy, &x_size_bathy);
        nc_inq_dimlen(ncid_bathy, dim_y_bathy, &y_size_bathy);
        
        nc_inq_dimlen(ncid_displ, dim_x_displ, &x_size_displ);
        nc_inq_dimlen(ncid_displ, dim_y_displ, &y_size_displ);
        
        
        // get the distance of the first two values on the x- and y-axis and
        // calculate the resolution of the axes from them
        // (does only work for equidistant values in the grid)
        size_t u0 = 0;
        size_t u1 = 1;
        
        float x_start_bathy, y_start_bathy;
        float x_next_bathy,  y_next_bathy;
        
        float x_start_displ, y_start_displ;
        float x_next_displ,  y_next_displ;
        
        // bathymetry file
        if(err_val = nc_get_var1_float(ncid_bathy, y_id_bathy, &u0, &y_start_bathy))
            cerr <<  nc_strerror(err_val) << endl;
        if(err_val = nc_get_var1_float(ncid_bathy, y_id_bathy, &u1, &y_next_bathy))
            cerr <<  nc_strerror(err_val) << endl;
        
        if(err_val = nc_get_var1_float(ncid_bathy, x_id_bathy, &u0, &x_start_bathy))
            cerr <<  nc_strerror(err_val) << endl;
        if(err_val = nc_get_var1_float(ncid_bathy, x_id_bathy, &u1, &x_next_bathy))
            cerr <<  nc_strerror(err_val) << endl;
        
        y_delta_bathy = y_next_bathy - y_start_bathy;
        x_delta_bathy = x_next_bathy - x_start_bathy;
        
        // displacements file
        if(err_val = nc_get_var1_float(ncid_displ, y_id_displ, &u0, &y_start_displ))
            cerr <<  nc_strerror(err_val) << endl;
        if(err_val = nc_get_var1_float(ncid_displ, y_id_displ, &u1, &y_next_displ))
            cerr <<  nc_strerror(err_val) << endl;
        
        if(err_val = nc_get_var1_float(ncid_displ, x_id_displ, &u0, &x_start_displ))
            cerr <<  nc_strerror(err_val) << endl;
        if(err_val = nc_get_var1_float(ncid_displ, x_id_displ, &u1, &x_next_displ))
            cerr <<  nc_strerror(err_val) << endl;
        
        y_delta_displ = y_next_displ - y_start_displ;
        x_delta_displ = x_next_displ - x_start_displ;
        
        
        return 0;
    }; 
    
private:
    // all variables for the bathymetric  data are marked with *_bathy
    //                       displacement data are marked with *_displ
    
    // file id
    int ncid_bathy, ncid_displ;
    
    // variable ids
    int x_id_bathy, y_id_bathy, z_id_bathy;
    int x_id_displ, y_id_displ, z_id_displ;
    
    // sizes of the grids in the nc-files
    size_t x_size_bathy, y_size_bathy;
    size_t x_size_displ, y_size_displ;
    
    // starting position and resolution of the bathymetry data
    float x_start_bathy, y_start_bathy;
    float x_delta_bathy, y_delta_bathy;
    
    // starting position and resolution of the displacement data
    float x_start_displ, y_start_displ;
    float x_delta_displ, y_delta_displ;
    
    
    /**
     * toGridCoordinates calculates the indices for accessing the bathymetry
     * and displacement data from the coordinates
     * 
     * @param source selects whether we want to get bathymetric or displacement data
     * @param x_coord x-coordinate
     * @param y_coord y-coordinate
     * @param x_index place to write the index for the x-axis
     * @param y_index place to write the index for the y-axis
     * @return true if the coordinate is inside the boundaries of the data, false otherwise
     */
    bool toGridCoordinates(DataSource source,
                           float x_coord, float y_coord,
                           size_t* x_index, size_t* y_index) {
        
        float y_start, x_start;
        float y_delta, x_delta;
        float y_size,  x_size;
        
        bool coordinate_valid = true;
        
        // select the correct data  
        if(source == BATHYMETRY) {
            y_start = y_start_bathy;
            x_start = x_start_bathy;
            y_delta = y_delta_bathy;
            x_delta = x_delta_bathy;
            y_size  = y_size_bathy;
            x_size  = x_size_bathy;
        }
        else {
            y_start = y_start_displ;
            x_start = x_start_displ;
            y_delta = y_delta_displ;
            x_delta = x_delta_displ;
            y_size  = y_size_displ;
            x_size  = x_size_displ;
        }
        
        *y_index = (size_t) ( ( (y_coord - y_start) / y_delta ) + 0.5f ); 
        *x_index = (size_t) ( ( (x_coord - x_start) / x_delta ) + 0.5f );
        
        // check if each of the indices is outside the boundaries
        if(*x_index < 0) {
            *x_index = 0;
            coordinate_valid = false;
        }
        else if(*x_index >= x_size) {
            *x_index = x_size - 1;
            coordinate_valid = false;
        }
        
        if(*y_index < 0) {
            *y_index = 0;
            coordinate_valid = false;
        }
        else if(*y_index >= y_size) {
            *y_index = y_size - 1;
            coordinate_valid = false;
        }
        
        return coordinate_valid;
    };
    
    
    /**
     * getOriginalBathymetry returns the bathymetry data without displacement
     * 
     * @param x x-coordinate
     * @param y y-coordinate
     * @return bathymetry at position (#x, #y) without displacement
     */
    float getOriginalBathymetry(float x, float y) {
        int err_val;
        size_t index[2];
        
        float result;
        
        if(!toGridCoordinates(BATHYMETRY, x, y, &index[0], &index[1])) {
            cerr << "Warning: bathymetry requested is outside the boundaries of the data!" << endl;
            cerr << "x: " << x << '\t' << (signed int) index[0] << endl;
            cerr << "y: " << y << '\t' << (signed int) index[1] << endl;
            
            // bathymetry outside of our data will be handled as landmass with minimum elevation
            return bath_min_zero_offset;
        }
        
        err_val = nc_get_var1_float(ncid_bathy, z_id_bathy, index, &result);
        if( err_val )
            cerr <<  nc_strerror(err_val) << endl;
        
        return result;
    };
};
