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


_(Cm,	    200.e-6	,   uF)
_(IKx,	    1.0 	,   nA)
_(kmKi,	    10.		,   mM)
_(kmK,	    1. 		,   mM)
_(kmNa,	    40. 	,   mM)
_(kmCa,	    7.e-4 	,   mM)
_(Vsh,	    20. 	,   mV)
_(INaKx,    0.7 	,   nA)
_(GNa,	    2.5 	,   uS)
_(Gto,	    .005 	,   uS)
_(GbK,	    .0006 	,   uS)
_(GKi,	    1.	 	,   uS)
_(GbNa,	    .0006 	,   uS)
_(GbCa,	    .00025 	,   uS)
_(PCa,	    .25 	,   nA/mM)
_(PCaK,	    .002	,   1)
_(PCaNa,    .002	,   1)
_(Cao,	    2. 		,   mM)
_(Ko,	    4. 		,   mM)
_(Nao,	    140. 	,   mM)
_(kNaCa,    5.e-4	,   nA)
_(dNaCa,    0.		,   1)
_(gamma,    .5		,   1)
_(kcyca,    3.e-4	,   mM)
_(kxcs,	    .4		,   mM)
_(ksrca,    .5		,   mM)
_(F,	    96485.	,   Coul/mol)
_(R,	    8314.41	,   mJ/(mol*Kel))
_(T,	    310.	,   Kel)
_(Vecs,	    .4		,   1)
_(radius,   15.		,   um)
_(length,   80.		,   um)
_(Vcell,    1.e-9*M_PI*sqr(radius)*length,  1/*uL*/) /* NB: length and volume disagree */
_(Vi,	    (1.-Vecs-Vup-Vrel)*Vcell,	    1/*uL*/)
_(Vup,	    .01		,   1)
_(Vrel,	    .1		,   1)
_(Vsrup,    Vcell*Vup	,   1/*uL*/)
_(kmCad,    625.	,   nA/mM)
_(Mtrop,    .02		,   mM)
_(Ctrop,    .05		,   mM)
_(RT_F,	    R*T/F	,   1)
