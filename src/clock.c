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


/* SHOW/PRINT THE STATE OF A TIMER */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"
#include "device.h"
#include "state.h"
#include "qpp.h"
#include "bikt.h"

typedef int STR;

RUN_HEAD(clock)
#if MPI
	MESSAGE("\rt=%-6ld",t);
#else
  char l[80];
  sprintf (l,"t=%-6ld",t);
  crt_text(l,w.row0,w.col0,w.color);
#endif
RUN_TAIL(clock)

DESTROY_HEAD(clock)
DESTROY_TAIL(clock)

CREATE_HEAD(clock)
CREATE_TAIL(clock,0)


