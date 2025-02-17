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

/* 2D diffusion with variable diffusivity */

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

typedef struct {
  int vD;
  real hx;
} STR;

/****************/
RUN_HEAD(diff2dv)

  DEVICE_CONST(int,vD)
  DEVICE_CONST(real,hx)

  int x, y, z;
  real *u, *D, *out;
  real gam=0.5/(hx*hx);
  for(z=s.z0;z<=s.z1;z++) for(y=s.y0;y<=s.y1;y++) for(x=s.x0;x<=s.x1;x++) {
    u=New+ind(x,y,z,s.v0);
    D=New+ind(x,y,z,vD);
    out=New+ind(x,y,z,s.v1);
    *out = gam*(
      (D[+DX]+D[0])*u[+DX] +
      (D[-DX]+D[0])*u[-DX] +
      (D[+DY]+D[0])*u[+DY] +
      (D[-DY]+D[0])*u[-DY]
       - (D[+DX]+D[-DX]+D[+DY]+D[-DY]+4*D[0])*u[0]
    );
  }
RUN_TAIL(diff2dv)

DESTROY_HEAD(diff2dv)
DESTROY_TAIL(diff2dv)

CREATE_HEAD(diff2dv)
  ACCEPTI(vD,INONE,0,vmax-1);
  ACCEPTR(hx,RNONE,0.,RNONE);
  ASSERT( dev->s.v1 != dev->s.v0 );
  ASSERT( S->vD != dev->s.v0 );
  ASSERT( S->vD != dev->s.v1 );
CREATE_TAIL(diff2dv,1)

