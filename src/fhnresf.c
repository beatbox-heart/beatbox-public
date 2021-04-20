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

/* spiral wave in a rotating system of reference for 
 * response functions in the Cubic FitzHugh-Nagumo 
 * model in Winfree's notation [Chaos1(3):305,1991],
 * with some more parameters added
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
#include "rhs.h"
#include "qpp.h"

typedef struct {
  real
    a,  b, vc, 
    bet, gam,			/* Winfree's Kinetics A parameters */
    epsu, epsv;			/* may be required different */
  real 
    Iu,  Iv,  Iw; 			/* External influence */
} STR;

RHS_HEAD(fhnresf,3) {
  DEVICE_CONST(real,a)
  DEVICE_CONST(real,b)
  DEVICE_CONST(real,vc)
  DEVICE_CONST(real,bet)
  DEVICE_CONST(real,gam)
  DEVICE_CONST(real,epsu)
  DEVICE_CONST(real,epsv)
  DEVICE_CONST(real,Iu)
  DEVICE_CONST(real,Iv)
  DEVICE_CONST(real,Iw)
  real U, U3, V, dUdt;
  U = u[0];
  V = u[1];
  U3=U*U; U3*=U;
  du[0] = dUdt = (U-U3/3.0-V)/epsu + Iu;
  du[1] = (U+bet-gam*V)*epsv + Iv;
  du[2] = -b*dUdt*(V-vc) + a*Iw;  
} RHS_TAIL(fhnresf)

RHS_CREATE_HEAD(fhnresf)
  ACCEPTP(a,RNONE,RNONE,RNONE);
  ACCEPTP(b,RNONE,RNONE,RNONE);
  ACCEPTP(vc,0.,RNONE,RNONE);
  ACCEPTP(epsu,RNONE,RSUCC(0.),RNONE);
  ACCEPTP(epsv,RNONE,RSUCC(0.),RNONE);
  ACCEPTP(bet,RNONE,RSUCC(0.),RNONE);
  ACCEPTP(gam,RNONE,RSUCC(0.),RNONE);
  ACCEPTP(Iu,0,RNONE,RNONE);
  ACCEPTP(Iv,0,RNONE,RNONE);
  ACCEPTP(Iw,0,RNONE,RNONE);
RHS_CREATE_TAIL(fhnresf,3)
