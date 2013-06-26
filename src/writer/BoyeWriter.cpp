/**
 * @file
 * This file is part of SWE.
 *
 * @author Thomas Blocher
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
 * @param NumberOfBoyes contains the Number of Boyes Written in this Filw
 */
io::BoyeWriter::BoyeWriter( const std::string &i_baseName,
                                int NumberOfBoyes)
{
    int status;
    timeStep = 0;
    initBoyes = 0;
    int x_int[NumberOfBoyes];
    int y_int[NumberOfBoyes];
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
io::BoyeWriter::~BoyeWriter() {
	
    nc_close(dataFile);
}

/**
 * Initialize A Boye give x and y position 
 *
 * @param X-Position of Boye
 * @param Y-position of boye
 * @param numberOfBoye
 */
void io::BoyeWriter::initBoye( float l_x, float l_y, SWE_DimensionalSplitting &block, int number) {
	//Define Positon of Boye #number 
	size_t pos = number;
	nc_put_var1_float(dataFile, xVar, &pos, &l_x);
	nc_put_var1_float(dataFile, yVar, &pos, &l_y);
    x_int[number] = block.getXpos(l_x);
    y_int[number] = block.getXpos(l_y);
    assert(x_int[number]>=0);
    assert(y_int[number]>=0);
	nc_sync(dataFile);
    initBoyes++;
}
/**
 * Write BoyeData 
 * @param Time of Data
 * @param Value of WaterHigh
 * @param Number of the Boye being writtend
 */
void io::BoyeWriter::writeBoye( float time, SWE_DimensionalSplitting &block) {
	//Put waterhigh for Boye
    nc_put_var1_float(dataFile, timeVar, &timeStep, &time);
	size_t Pos[] = {0,0}; 
    for(int i = 0; i<initBoyes; i++){
        Pos[0] = timeStep;
        Pos[1] = i;
        int x = x_int[i];
        int y = y_int[i]; 
        float height = block.getWaterHeight()[x][y]; 
	    nc_put_var1_float(dataFile, hVar, Pos, &height);
    }
    nc_sync(dataFile);
    timeStep++;
}
