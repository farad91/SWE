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
        float height = (-1) * getOriginalBathymetry(x,y);
        
        // apply the minimum height
        if(height < bath_min_zero_offset)
            return bath_min_zero_offset;
        else
            return height;
    };
    
    
    float getDynamicBathymetry(float x, float y, float time) {
        int   err_val;
        float displacement = 0.f;
        float result;
        
        size_t index[3];        
        
        // try to get the indices for the displacement data
        if( !toGridCoordinates(DISPLACEMENT, x ,y, &index[2], &index[1]) || time == 0.f ) {
            // we are outside the area where the displacement is relevant
            result = getOriginalBathymetry(x,y);
        }
        else {
            if(DynamicDispl){
                toTimeCoordinates(time,&index[0]);
                // displacement data is available for this position
                err_val = nc_get_var1_float(ncid_displ, z_id_displ, index, &displacement);
                if( err_val )
                    cerr << "Error in getDynamicBathymetry: " << endl
                        << nc_strerror(err_val) << endl;
                    result = getOriginalBathymetry(x,y) + displacement;
                }
            else{
                 // displacement data is available for this position
                err_val = nc_get_var1_float(ncid_displ, z_id_displ, &index[1], &displacement);
                if( err_val )
                    cerr << "Error in getDynamicBathymetry: " << endl
                        << nc_strerror(err_val) << endl;
                if (time<180.f)
                    result = getOriginalBathymetry(x,y) + displacement*(time/180.f);
                else
                    result = getOriginalBathymetry(x,y) + displacement;
            }
               
        }
        
        // apply the minimum elevation or depth
        if( result > (-1) * bath_min_zero_offset && result < 0)
            // shallow water
            result = (-1) * bath_min_zero_offset;
        else if( result < bath_min_zero_offset && result >= 0)
            // very low elevation of the landmass
            result = bath_min_zero_offset;
        
        return result;
    };
    
    float getEruptionDuration() {
    if(DynamicDispl)
        return time_start_displ+(time_delta_displ*time_size_displ); 
    else
        return 180.f;
    };
    
    float getEruptionResolution() {
    if(DynamicDispl)
        return time_delta_displ; 
    else
        return 2.f;
    };
    
    float getBathymetry(float x, float y) {
        int   err_val;
        float displacement = 0.f;
        float result;
        
        size_t index[2];        
        
        // try to get the indices for the displacement data
        if( !toGridCoordinates(DISPLACEMENT, x ,y, &index[1], &index[0]) || DynamicDispl ) {
            // we are outside the area where the displacement is relevant
            result = getOriginalBathymetry(x,y);
        }
        else {
            // displacement data is available for this position
            err_val = nc_get_var1_float(ncid_displ, z_id_displ, index, &displacement);
            if( err_val )
                cerr <<"Error in getBathymetry: "<< nc_strerror(err_val) << endl;
            
            result = getOriginalBathymetry(x,y) + displacement;
        }
        
        // apply the minimum elevation or depth
        if( result > (-1) * bath_min_zero_offset && result < 0)
            // shallow water
            result = (-1) * bath_min_zero_offset;
        else if( result < bath_min_zero_offset && result >= 0)
            // very low elevation of the landmass
            result = bath_min_zero_offset;
        
        return result;
    };
    
    
    /**
     * getBoundaryPos will return the position of the boundary #i_edge
     * on the axis orthogonal to the boundary
     * 
     * @param i_edge the boundary we want to get the position of
     * @return the position of the boundary on the axis orthogonal to it
     */
    float getBoundaryPos(BoundaryEdge i_edge) {
        
        //! variable ID of the axis orthogonal to the boundary
        int var_id;
        //! position of the boundary on the axis
        size_t index;
        
        //! error value from the netcdf call
        int err_val;
        
        float result;
        
        size_t u0    = 0;
        size_t x_max = x_size_bathy-1;
        size_t y_max = y_size_bathy-1;
        
        // select the axis (x or y)
        if( i_edge == BND_RIGHT || i_edge == BND_LEFT )
            var_id = x_id_bathy;
        else
            var_id = y_id_bathy;
        
        // get the index of the boundary
        if( i_edge == BND_RIGHT )
            index = x_max;
        else if( i_edge == BND_TOP )
            index = y_max;
        else
            index = u0;
        
        // get the coordinate from the nc-file
        // or print an error message if that fails
        err_val = nc_get_var1_float(ncid_bathy, var_id, &index, &result);
        if( err_val )
            cerr << "Error reading boundary position:" << endl
                 << "\t" << nc_strerror(err_val) << endl;
        
        return result;
    };
    
    BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };

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
        int dim_x_bathy, dim_y_bathy, dim_x_displ, dim_y_displ, dim_time_displ;
        
        // open the bathymetry file
        if(err_val = nc_open(file_bathy, NC_NOWRITE, &ncid_bathy)) {
            cerr << "Error: cannot open bathymetry file!" << endl
                 << nc_strerror(err_val) << endl;
            return err_val;
        }
        
        // open the displacements file
        if(err_val = nc_open(file_displ, NC_NOWRITE, &ncid_displ)) {
            cerr << "Error: cannot open displacement file!" << endl
                 << nc_strerror(err_val) << endl;
            return err_val;
        }
        
        // get the variables IDs 
        // for the bathymetry
        err_val = nc_inq_varid(ncid_bathy, "x", &x_id_bathy );
        if( err_val )
            cerr << "Error getting variable #x for the bathymetry:" << endl
                 << nc_strerror(err_val) << endl;
        
        err_val = nc_inq_varid(ncid_bathy, "y", &y_id_bathy );
        if( err_val )
            cerr << "Error getting variable #y for the bathymetry:" << endl
                 << nc_strerror(err_val) << endl;
        
        err_val = nc_inq_varid(ncid_bathy, "z", &z_id_bathy );
        if( err_val )
            cerr << "Error getting variable #z for the bathymetry:" << endl
                 << nc_strerror(err_val) << endl;
        
        // for the displacement
        err_val = nc_inq_varid(ncid_displ, "x", &x_id_displ );
        if( err_val )
            cerr << "Error getting variable #x for the displacement:" << endl
                 << nc_strerror(err_val) << endl;
        
        err_val = nc_inq_varid(ncid_displ, "y", &y_id_displ );
        if( err_val )
            cerr << "Error getting variable #y for the displacement:" << endl
                 << nc_strerror(err_val) << endl;
        
        err_val = nc_inq_varid(ncid_displ, "z", &z_id_displ );
        if( err_val )
            cerr << "Error getting variable #z for the displacement:" << endl
                 << nc_strerror(err_val) << endl;
        
        err_val = nc_inq_varid(ncid_displ, "time", &time_id_displ );
        if( err_val ){
            DynamicDispl = false;
            cerr << "Error getting variable #time for the displacement probatly not existend:" << endl
                 << nc_strerror(err_val) << endl;
        }
        else{
            DynamicDispl = true;
        }
        // get the number of dimensions contained by 'x' and 'y'
        err_val = nc_inq_dimid(ncid_bathy, "x", &dim_x_bathy);
        if( err_val )
            cerr << "Error getting the dimension id of #x for the bathymetry:" << endl
                 << nc_strerror(err_val) << endl;
        
        err_val = nc_inq_dimid(ncid_bathy, "y", &dim_y_bathy);
        if( err_val )
            cerr << "Error getting the dimension id of #y for the bathymetry:" << endl
                 << nc_strerror(err_val) << endl;
        
        err_val = nc_inq_dimid(ncid_displ, "x", &dim_x_displ);
        if( err_val )
            cerr << "Error getting the dimension id of #y for the bathymetry:" << endl
                 << nc_strerror(err_val) << endl;
        
        err_val = nc_inq_dimid(ncid_displ, "y", &dim_y_displ);
        if( err_val )
            cerr << "Error getting the dimension id of #y for the bathymetry:" << endl
                 << nc_strerror(err_val) << endl;
        if(DynamicDispl){         
            err_val = nc_inq_dimid(ncid_displ, "time", &dim_time_displ);
            if( err_val )
                cerr << "Error getting the dimension id of #time for the bathymetry:" << endl
                     << nc_strerror(err_val) << endl;
        }
        
        // get length x, y;
        nc_inq_dimlen(ncid_bathy, dim_x_bathy, &x_size_bathy);
        nc_inq_dimlen(ncid_bathy, dim_y_bathy, &y_size_bathy);
        
        nc_inq_dimlen(ncid_displ, dim_x_displ, &x_size_displ);
        nc_inq_dimlen(ncid_displ, dim_y_displ, &y_size_displ);
        if(DynamicDispl)
            nc_inq_dimlen(ncid_displ, dim_time_displ, &time_size_displ);
        
        
        // get the distance of the first two values on the x- and y-axis and
        // calculate the resolution of the axes from them
        // (does only work for equidistant values in the grid)
        size_t u0 = 0;
        size_t u1 = 1;
        
        float x_next_bathy,  y_next_bathy;
        float x_next_displ,  y_next_displ, time_next_displ;
        
        // bathymetry file
        err_val = nc_get_var1_float(ncid_bathy, y_id_bathy, &u0, &y_start_bathy);
        if( err_val )
            cerr << nc_strerror(err_val) << endl;
        
        err_val = nc_get_var1_float(ncid_bathy, y_id_bathy, &u1, &y_next_bathy);
        if( err_val )
            cerr << nc_strerror(err_val) << endl;
        
        
        err_val = nc_get_var1_float(ncid_bathy, x_id_bathy, &u0, &x_start_bathy);
        if( err_val )
            cerr << nc_strerror(err_val) << endl;
        
        err_val = nc_get_var1_float(ncid_bathy, x_id_bathy, &u1, &x_next_bathy);
        if( err_val )
            cerr << nc_strerror(err_val) << endl;
        
        
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
        if(DynamicDispl){    
            if(err_val = nc_get_var1_float(ncid_displ, time_id_displ, &u0, &time_start_displ))
                cerr <<  nc_strerror(err_val) << endl;
            if(err_val = nc_get_var1_float(ncid_displ, time_id_displ, &u1, &time_next_displ))
                cerr <<  nc_strerror(err_val) << endl;
        
            y_delta_displ = y_next_displ - y_start_displ;
            x_delta_displ = x_next_displ - x_start_displ;
            time_delta_displ = time_next_displ - time_start_displ;
        }
        
#ifdef DEBUG
        printDebugInfo();
#endif
        
        return 0;
    }; 
    
private:
    // all variables for the bathymetric  data are marked with *_bathy
    //                       displacement data are marked with *_displ
    
    //type of bath file
    bool DynamicDispl;
    // file id
    int ncid_bathy, ncid_displ;
    
    // variable ids
    int x_id_bathy, y_id_bathy, z_id_bathy;
    int x_id_displ, y_id_displ, z_id_displ, time_id_displ;
    
    // sizes of the grids in the nc-files
    size_t x_size_bathy, y_size_bathy;
    size_t x_size_displ, y_size_displ, time_size_displ;
    
    // starting position and resolution of the bathymetry data
    float x_start_bathy, y_start_bathy;
    float x_delta_bathy, y_delta_bathy;
    
    // starting position and resolution of the displacement data
    float x_start_displ, y_start_displ, time_start_displ;
    float x_delta_displ, y_delta_displ, time_delta_displ;
    
    /**
    * toTimeCoordinates calculates the nearest time where a displacement exists
    * @param time
    * @param time_index place to write the index for the Timestep
    * @return true if Earthpuake is over
    */
    bool toTimeCoordinates(float time, size_t* time_index){
        *time_index = (size_t) ( ( (time - time_start_displ) / time_delta_displ ) + 0.5f );
        if(*time_index < 0)
            *time_index = 0;
        if(*time_index >= time_size_displ){
            *time_index = time_size_displ-1;
            return true;
        } 
        return false;
    };
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
        
        
        // check if any of the indices is outside of the boundaries
        if(*x_index < 0 || *x_index >= x_size)
            coordinate_valid = false;
        
        if(*y_index < 0 || *y_index >= y_size)
            coordinate_valid = false;
        
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
        
        if(!toGridCoordinates(BATHYMETRY, x, y, &index[1], &index[0])) {
            
            // move position to the next boundary
            if((int) index[0] < 0)
                index[0] = 0;
            else if(index[0] >= y_size_bathy)
                index[0] = y_size_bathy-1;
            
            if((int) index[1] < 0)
                index[1] = 0;
            else if(index[1] >= x_size_bathy)
                index[1] = x_size_bathy-1;
            
        }
        
        // get bathymetry from file
        err_val = nc_get_var1_float(ncid_bathy, z_id_bathy, index, &result);
        if( err_val )
            cerr << "getOrginalBathymetrie: " << nc_strerror(err_val) << endl;
        
        return result;
    };
    
    void printDebugInfo() {
        cerr << "x_start_bathy: " << x_start_bathy << endl;
        cerr << "y_start_bathy: " << x_start_bathy << endl;
        cerr << "x_delta_bathy: " << x_delta_bathy << endl;
        cerr << "y_delta_bathy: " << x_delta_bathy << endl;
        cerr << "x_size_bathy: " << x_size_bathy << endl;
        cerr << "y_size_bathy: " << y_size_bathy << endl;
        
        cerr << "x_start_displ: " << x_start_displ << endl;
        cerr << "y_start_displ: " << x_start_displ << endl;
        cerr << "x_delta_displ: " << x_delta_displ << endl;
        cerr << "y_delta_displ: " << x_delta_displ << endl;
        cerr << "x_size_displ: " << x_size_displ << endl;
        cerr << "y_size_displ: " << y_size_displ << endl;
        cerr << "time_delta_displ: " << time_delta_displ << endl;
        cerr << "time_size_displ: " << time_size_displ << endl;
        cerr << "time of earth movement: " << time_start_displ+(time_delta_displ*time_size_displ) << endl;
        
        
        cerr << "boundary TOP: \t" << getBoundaryPos(BND_TOP) << endl;
        cerr << "boundary BOTTOM: \t" << getBoundaryPos(BND_BOTTOM) << endl;
        cerr << "boundary LEFT: \t" << getBoundaryPos(BND_LEFT) << endl;
        cerr << "boundary RIGHT: \t" << getBoundaryPos(BND_RIGHT) << endl;
    }
};
