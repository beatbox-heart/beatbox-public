/**
 * Copyright (C) (2010-2016) Vadim Biktashev, Irina Biktasheva et al. 
 * (see ../AUTHORS for the full list of contributors)
 *
 * This file is part of Beatbox.
 *
 * Beatbox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beatbox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Beatbox.  If not, see <http://www.gnu.org/licenses/>.
 */


#define NDIMS 4
MPI_Datatype filetype;
MPI_Datatype sourceType;
int sizes[NDIMS],subsizes[NDIMS],starts[NDIMS];
/*---------------------------------------------------------------------------
 * This first datatype describes how the local space
 * contributes to the file (the global space)
 *---------------------------------------------------------------------------*/
/*  Global space dimensions. */
int xsize = (dev->s.global_x1 - dev->s.global_x0) + 1;
int ysize = (dev->s.global_y1 - dev->s.global_y0) + 1;
int zsize = (dev->s.global_z1 - dev->s.global_z0) + 1;
int vsize = (dev->s.v1 - dev->s.v0) + 1;

/*  Local space dimensions */
int local_xsize = (dev->s.x1 - dev->s.x0) + 1;
int local_ysize = (dev->s.y1 - dev->s.y0) + 1;
int local_zsize = (dev->s.z1 - dev->s.z0) + 1;
int local_vsize = vsize;

/*  Dimensions of the global space. */
sizes[0] = xsize;
sizes[1] = ysize;
sizes[2] = zsize;
sizes[3] = vsize;

/*  Dimensions of the local space. */
subsizes[0] = local_xsize;
subsizes[1] = local_ysize;
subsizes[2] = local_zsize;
subsizes[3] = local_vsize;

starts[0] = dev->s.x0 - dev->s.global_x0;
starts[1] = dev->s.y0 - dev->s.global_y0;
starts[2] = dev->s.z0 - dev->s.global_z0;
starts[3] = 0;

MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&filetype),"Couldn't define filetype.");
MPIDO(MPI_Type_commit(&filetype),"Couldn't commit filetype.");

/*---------------------------------------------------------------------------
 * The datatype defined below describes how which data will
 * be pulled from New to be written to file.
 *
 * With this defined, it's possible to write a single block
 * of data to file containing the whole space.
 *---------------------------------------------------------------------------*/

/*  Dimensions of New. */
sizes[0] = (local_xmax-local_xmin)+(ONE*2);
sizes[1] = (local_ymax-local_ymin)+(TWO*2);
sizes[2] = (local_zmax-local_zmin)+(TRI*2);
sizes[3] = vmax;

/*  Position of local space in New. */
starts[0] = (dev->s.x0 + ONE) - local_xmin;
starts[1] = (dev->s.y0 + TWO) - local_ymin;
starts[2] = (dev->s.z0 + TRI) - local_zmin;
starts[3] = dev->s.v0;

MPIDO(MPI_Type_create_subarray(NDIMS,sizes,subsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&sourceType),"Couldn't define sourceType.");
MPIDO(MPI_Type_commit(&sourceType),"Couldn't commit sourceType.");

#undef NDIMS
/*---------------------------------------------------------------------------*/
