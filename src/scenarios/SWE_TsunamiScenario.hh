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
 * TODO
 */

#include<netcdf.h>
using namespace std;
class SWE_TsunamiScenario(String filename) {
    
public:
    int ncid, zid, xid ,yid, dimX, dimY, lenX, lenY;
    //Open .nc file
    if((retval = nc_open(filename,NC_NOWRITE, &ncid)))
        ERR(retval);
    //get ID'S 
    if((retval = nc_inq_varid(ncid, "x", &xid)))
        ERR(retval);
    if((retval = nc_inq_varid(ncid, "y", &yid)))
        ERR(retval);
    if((retval = nc_inq_varid(ncid, "z", &zid)))
        ERR(retval);
    if((retval = nc_inq_dimid(ncid, "x", &dimX)))
        ERR(retval);
    if((retval = nc_inq_dimid(ncid, "y", &dimY)))
        ERR(retval);
    //get length x,y;
    nc_inq_dimlen(ncid,dimX,&lenX);
    nc_inq_dimlen(ncid,dimY,&lenY);
    //GET boundery position
       float getBoundaryPos(BoundaryEdge i_edge) {
       if ( i_edge == BND_LEFT || i_edge == BND_RIGHT)
            float bndx[lenX);
            if((retval = nc_get_var_int(ncid, xid, bndx[0])))
                ERR(retval);
            if(i_edge == BND_LEFT);
                return bndx[0];
            else
                return bndx[lenX-1];
       else
            float bndy[lenY);
            if((retval = nc_get_var_int(ncid, yid, bndy[0])))
                ERR(retval);
            if ( i_edge == BND_BOTTOM)
                return bndy[0];
            else
                return bndy[lenY-1];
    }    
};

