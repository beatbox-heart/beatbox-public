
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
 * Experimental module, not for distribution. 
 * Periodic boundary conditions along y axis 
 *
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
RUN_HEAD(tory)
  int x, z, v;
  real *u;

  for(v=s.v0;v<=s.v1;v++) for(z=s.z0;z<=s.z1;z++) for(x=s.x0;x<=s.x1;x++) {
    New[ind(x,s.y0-1,z,v)] = New[ind(x,s.y1,z,v)];
    New[ind(x,s.y1+1,z,v)] = New[ind(x,s.y0,z,v)];
  }
RUN_TAIL(tory)

DESTROY_HEAD(tory)
DESTROY_TAIL(tory)

CREATE_HEAD(tory) {
  Space s=dev->s;
  ASSERT(s.y0>1)
  ASSERT(s.y1<ymax-1)
} CREATE_TAIL(tory,1)

