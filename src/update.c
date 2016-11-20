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

#include <stdio.h>
#include <stdlib.h>
#include "system.h"
#include "beatbox.h"
#include "bgi.h"
#include "bikt.h"
#include "device.h"

int update_now=0;
typedef int STR;	/* dummy declaration; for compatibility only */

RUN_HEAD(update) {
  if (update_now) update_graph();
} RUN_TAIL(update)


DESTROY_HEAD(update)
DESTROY_TAIL(update)

CREATE_HEAD(update)
CREATE_TAIL(update,1)


