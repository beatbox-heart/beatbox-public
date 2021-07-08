/**
 * Copyright (C) (2010-2021) Vadim Biktashev, Irina Biktasheva et al. 
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
 * Periodic boundary conditions along z axis 
 * Only for box geometry, sequential mode.
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
RUN_HEAD(torz)
  int x, y, v;
  real *u;

  for(v=s.v0;v<=s.v1;v++) for(y=s.y0;y<=s.y1;y++) for(x=s.x0;x<=s.x1;x++) {
    New[ind(x,y,s.z0-1,v)] = New[ind(x,y,s.z1,v)];
    New[ind(x,y,s.z1+1,v)] = New[ind(x,y,s.z0,v)];
  }
RUN_TAIL(torz)

DESTROY_HEAD(torz)
DESTROY_TAIL(torz)

CREATE_HEAD(torz) {
  Space s=dev->s;
  ASSERT(s.z0>1)
  ASSERT(s.z1<zmax-1)
} CREATE_TAIL(torz,1)

