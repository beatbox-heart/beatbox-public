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

/* 2D advection term (no diffusion) */

/* d_x(C_x u)+d_y(C_yu) */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"

#include "device.h"
#include "state.h"
#include "bikt.h"
#include "rhs.h"
#include "qpp.h"

typedef struct {
  real Dxx, Dxy, Dyx, Dyy, Cx, Cy;
  real hx;
  int vi;
} STR;

RHS_HEAD(conv2d,1)
{
  DEVICE_CONST(real,Cx);
  DEVICE_CONST(real,Cy);
  DEVICE_CONST(real,hx);
  DEVICE_CONST(int,vi);
  real jx, jy;
  u[vi]=( Cx*(u[+DX]-u[-DX]) + Cy*(u[+DY]-u[-DY]) )/(2*hx);
}
RHS_TAIL(conv2d)

RHS_CREATE_HEAD(conv2d)
{
  ACCEPTP(Cx,0,0,RNONE);
  ACCEPTP(Cy,0,0,RNONE);
  ACCEPTR(hx,RNONE,RSUCC(0.),RNONE);;
  ACCEPTI(vi,INONE,0,INONE);
  ASSERT(vi != 0);
}
RHS_CREATE_TAIL(conv2d,1)


