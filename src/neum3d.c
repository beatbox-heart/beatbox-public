/**
 * Copyright (C) (2010-2023) Vadim Biktashev, Irina Biktasheva et al. 
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

/* 
 * Neumann boundary conditions for the 
 * variable diffusion devices diff2dv, diff2vh. 
 * Only for box geometry.
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system.h"
#include "beatbox.h"
#include "device.h"
#include "state.h"
#include "bikt.h"
#include "qpp.h"

typedef int STR;	/* dummy declaration; for compatibility only */

/****************/
RUN_HEAD(neum3d)
{
  int x, y, z, v;
  real *u;

  for(v=s.v0;v<=s.v1;v++) {
    for(z=s.z0;z<=s.z1;z++) for(x=s.x0;x<=s.x1;x++) {
      if (mpi_iy==0)        {y=s.y0; u=New+ind(x,y,z,v); u[-DY]=u[0];} /* fore */
      if (mpi_iy==mpi_ny-1) {y=s.y1; u=New+ind(x,y,z,v); u[+DY]=u[0];} /* hind */
    }
    /* limits extended by 1 to include corner cells */
    for(z=s.z0;z<=s.z1;z++) for(y=s.y0-1;y<=s.y1+1;y++) {
      if (mpi_ix==0)        {x=s.x0; u=New+ind(x,y,z,v); u[-DX]=u[0];} /* left */
      if (mpi_ix==mpi_nx-1) {x=s.x1; u=New+ind(x,y,z,v); u[+DX]=u[0];} /* right */
    }
    /* limits extended by 1 to include corner cells */
    for(y=s.y0-1;y<=s.y1+1;y++) for(x=s.x0-1;x<=s.x1+1;x++) {
      if (mpi_iz==0)        {z=s.z0; u=New+ind(x,y,z,v); u[-DZ]=u[0];} /* bottom */
      if (mpi_iz==mpi_nz-1) {z=s.z1; u=New+ind(x,y,z,v); u[+DZ]=u[0];} /* top */
    }
  }
}
RUN_TAIL(neum3d)

/****************/
DESTROY_HEAD(neum3d)
DESTROY_TAIL(neum3d)

/****************/
CREATE_HEAD(neum3d)
{
  /* No extra parameters to accept; just make some checks */
  if (dim!=3) MESSAGE("/* Warning: 3D device is used in %d dimensional simulation */",dim);
  if (GEOMETRY_ON) MESSAGE("/* Warning: this device is intended for boxes only */");
}
CREATE_TAIL(neum3d,1)

