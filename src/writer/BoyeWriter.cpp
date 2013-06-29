/**
 * @file
 * This file is part of SWE.
 *
 * @author Thomas Blocher (blocher@in.tum.de)
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
 * A writer for the netCDF-format: http://www.unidata.ucar.edu/software/netcdf/
 */

#include "BoyeWriter.hh"
#include <string>
#include <vector>
#include <iostream>
#include <cassert>

/**
 * Create a netCdf-file
 * Any existing file will be replaced.
 *
 * @param i_baseName base name of the netCDF-file to which the data will be written to.
 * @param NumberOfBoyes contains the Number of Boyes Written in this File (if it is less 1 it's interpreted as Checkpoint to coninue)
 */
io::BoyeWriter::BoyeWriter( const std::string &i_baseName,
                                int NumberOfBoyes)
{
    if( NumberOfBoyes > 0){
        int status;
        timeStep = 0;
        initBoyes = 0;
        MaxBoyes = NumberOfBoyes;
        x_int = new int[NumberOfBoyes];
        y_int = new int[NumberOfBoyes];
        //create a netCDF-file, an existing file will be replaced
        status = nc_create((i_baseName+"_Boye.nc").c_str(), NC_NETCDF4, &dataFile);

        //check if the netCDF-file creation constructor succeeded.
        if (status != NC_NOERR) {
	        assert(false);
	        return;
        }

        //dimensions
        int l_timeDim, l_boyeDim;
        nc_def_dim(dataFile, "time", NC_UNLIMITED, &l_timeDim);
        nc_def_dim(dataFile, "Boye",  NumberOfBoyes , &l_boyeDim);


        nc_def_var(dataFile, "time", NC_FLOAT, 1, &l_timeDim, &timeVar);

        //variables, fastest changing index is on the right (C syntax), will be mirrored by the library
        int dims[] = {l_timeDim, l_boyeDim};
        nc_def_var(dataFile, "Data",  NC_FLOAT, 2, dims, &hVar);
        
        nc_def_var(dataFile, "x", NC_FLOAT, 1, &dims[1], &xVar);
        nc_def_var(dataFile, "y",  NC_FLOAT, 1, &dims[1], &yVar);
    }
    else{
        nc_open((i_baseName+"_Boye.nc").c_str(), NC_WRITE, &dataFile);
        
        //get dim id
        int l_timeDim,l_boyeDim;
        nc_inq_dimid(dataFile, "time", &l_timeDim);
        nc_inq_dimlen(dataFile, l_timeDim, &timeStep);
        nc_inq_dimid(dataFile, "Boye", &l_boyeDim);
        nc_inq_dimlen(dataFile, l_boyeDim, &initBoyes);
        // get the variables IDs 
        nc_inq_varid(dataFile, "Data",  &hVar);
        nc_inq_varid(dataFile, "x", &xVar);
        nc_inq_varid(dataFile, "y",  &yVar);
        nc_inq_varid(dataFile, "time", &timeVar);
        x_int = new int[initBoyes];
        y_int = new int[initBoyes];
        MaxBoyes = initBoyes;
    }
}
io::BoyeWriter::~BoyeWriter() {
	
    nc_close(dataFile);
}

/**
 * Initialize a Boye at given x and y position 
 *
 * @param X-Position of Boye
 * @param Y-position of boye
 */
void io::BoyeWriter::initBoye( float l_x, float l_y, SWE_DimensionalSplitting &block) {
	//Define Positon of Boye #number 
    if(initBoyes < MaxBoyes){
	    size_t pos = initBoyes;
	    nc_put_var1_float(dataFile, xVar, &pos, &l_x);
	    nc_put_var1_float(dataFile, yVar, &pos, &l_y);    
        x_int[initBoyes] = block.getXpos(l_x);
        y_int[initBoyes] = block.getYpos(l_y);
	    nc_sync(dataFile);
        initBoyes++;
    }
    else{
        float x;
        float y;
        for(int i = 0; i < MaxBoyes; i++){
            size_t pos = i;
	        nc_get_var1_float(dataFile, xVar, &pos, &x);
	        nc_get_var1_float(dataFile, yVar, &pos, &y);
            x_int[i] = block.getXpos(x);
            y_int[i] = block.getYpos(y);
        }
    }
}
/**
 * Write BoyeData write data of all initialized Boyes at given Time time 
 * @param Time of Data
 * @param Reference to Array containing Waterheight
 * @param Reference to Array containing Bathymetry
 */
void io::BoyeWriter::writeBoye( float time, const Float2D &h, const Float2D &b) {
	//Put waterhigh for Boye
    nc_put_var1_float(dataFile, timeVar, &timeStep, &time);
	size_t Pos[] = {0,0}; 
    for(int i = 0; i<initBoyes; i++){
        Pos[0] = timeStep;
        Pos[1] = i;
        int x = x_int[i];
        int y = y_int[i]; 
        float height = h[x][y]+b[x][y]; 
	    nc_put_var1_float(dataFile, hVar, Pos, &height);
    }
    nc_sync(dataFile);
    timeStep++;
}
