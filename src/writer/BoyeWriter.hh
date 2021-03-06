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

#ifndef BOYEWRITER_HH_
#define BOYEWRITER_HH_

#include <cstring>
#include <string>
#include <vector>
#include <netcdf.h>
#include "blocks/SWE_DimensionalSplitting.hh"
#include "tools/help.hh"
#ifdef MPI_INCLUDED_NETCDF
#undef MPI_INCLUDED
#undef MPI_INCLUDED_NETCDF
#endif

namespace io {
  class BoyeWriter;
}

class io::BoyeWriter{
private:
    /** netCDF file id*/
    int dataFile;

    /** Variable ids */
    int yVar, xVar, hVar, timeVar;
    /** Counter variables*/
    size_t initBoyes, MaxBoyes;
    size_t timeStep;
    /** Helper Arrays*/
    int *y_int, *x_int;

  public:
    BoyeWriter(const std::string &i_fileName, int NumberOfBoyes);
    virtual ~BoyeWriter();

    // Init a boye
    void initBoye( float x, float y, SWE_DimensionalSplitting &block);
    
    //write data for boye
    void writeBoye( float time,const Float2D &h, const Float2D &b);
};
#endif /* BOYEWriter_HH_ */
