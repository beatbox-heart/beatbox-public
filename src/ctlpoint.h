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


#if MPI

CHAR_VAR("version",	&V_ ,	&VERSTRING,	VERLEN,		if (fatal()) return 0)
INT_VAR("xmax",   	&xm_,	&xmax,		1,			return 0)
INT_VAR("ymax",   	&ym_,	&ymax,		1,			return 0)
INT_VAR("zmax",   	&zm_,	&zmax,		1,			return 0)
INT_VAR("vmax",   	&vm_,	&vmax,		1,			if (fatal()) return 0)

#else

D("version",&V_ ,&VERSTRING,char,	   VERLEN,if (fatal()) return 0)
D("xmax",   &xm_,&xmax,    INT,			1,	       return 0)
D("ymax",   &ym_,&ymax,    INT,			1,	       return 0)
D("zmax",   &zm_,&zmax,    INT,			1,	       return 0)
D("vmax",   &vm_,&vmax,    INT,			1,if (fatal()) return 0)
D("field",  New ,New,     real,xmax*ymax*zmax*vm_,if (fatal()) return 0)

#endif
