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

#include "NetCdfWriter.hh"
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

    //create a netCDF-file, an existing file will be replaced
    status = nc_create((i_baseName+"_Boye").c_str(), NC_NETCDF4, &dataFile);

  //check if the netCDF-file creation constructor succeeded.
    if (status != NC_NOERR) {
	    assert(false);
	    return;
    }

    //dimensions
    int l_timeDim, l_xDim, l_yDim, l_boundDim;
    nc_def_dim(dataFile, "time", NC_UNLIMITED, &l_timeDim);
    nc_def_dim(dataFile, "Boye",  NumberOfBoyes , &l_boyeDim);

    //variables (TODO: add rest of CF-1.5)
    int l_xVar, l_yVar;

    nc_def_var(dataFile, "time", NC_FLOAT, 1, &l_timeDim, &timeVar);
    ncPutAttText(timeVar, "long_name", "Time");
    ncPutAttText(timeVar, "units", "seconds since simulation start"); // the word "since" is important for the paraview reader


    //variables, fastest changing index is on the right (C syntax), will be mirrored by the library
    int dims[] = {l_timeDim, l_BoyeDim};
    nc_def_var(dataFile, "Data",  NC_FLOAT, 2, dims, &hVar);
    
    nc_def_var(dataFile, "x", NC_FLOAT, 1, &dims[1], &xVar);
    nc_def_var(dataFile, "y",  NC_FLOAT, 1, &dims[1], &yVar);
    //set attributes to match CF-1.5 convention
    ncPutAttText(NC_GLOBAL, "Conventions", "CF-1.5");
    ncPutAttText(NC_GLOBAL, "title", "Computed tsunami solution");
    ncPutAttText(NC_GLOBAL, "history", "SWE");
    ncPutAttText(NC_GLOBAL, "institution", "Technische Universitaet Muenchen, Department of Informatics, Chair of Scientific Computing");
    ncPutAttText(NC_GLOBAL, "source", "Boye Data.");
    ncPutAttText(NC_GLOBAL, "references", "http://www5.in.tum.de/SWE");
    ncPutAttText(NC_GLOBAL, "comment", "SWE is free software and licensed under the GNU General Public License. Remark: In general this does not hold for the used input data.");
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
void io::BoyeWriter::initBoye( float l_x, float l_y, int number) {
	//Define Positon of Boye #number 
	size_t pos = number;
	nc_put_var_float(dataFile, xVar, pos, &l_x);
	nc_put_var_float(dataFile, xVar, pos, &l_x);
	nc_sync(dataFile);
}
/**
 * Write BoyeData 
 * @param Time of Data
 * @param Value of WaterHigh
 * @param Number of the Boye being writtend
 */
void io::BoyeWriter::writeBoye( float time, float waterhigh, int number) {
	//Put waterhigh for Boye
	size_t Pos[] = {time,number};
	nc_put_var_float(dataFile, BoyeVar, Pos, &l_x);
    nc_sync(dataFile);
}
