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
 * IONIC description of the Courtemanche et al. 1998 model.
  Human atrial cell model as described in:
  Courtemanche M, Ramirez RJ, Nattel S. Ionic mechanisms underlying
  human atrial action potential properties: insights from a mathematical
  model. Am J Physiol (1998);275:H301-H321. 
  This code is based on the M.Courtemanche's `error-free' electronic 
  version of equations and parameters. 
  Changes made: the forward Euler updates replaced with ODE descriptions, 
  and Rush-Larsen updates replaced with alpha-beta descriptions
  when those depend only on V, and ODE descriptions otherwise. 
  NB all essential contents are stoved away into crn98_*.h files. 
*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"
#include "beatbox.h"

#include "device.h"
#include "state.h"
#include "bikt.h"
#include "ionic.h"
#include "qpp.h"

static real cub(real x) {return x*x*x;}

/* Define constants (unchangeable parameters) */
/* that can be used in tabulated functions. */
#include "crn98_const.h"

/* Enumerate all dynamic variables */
enum {
  #define _(n,i) var_##n,
  #include "crn98_other.h"
  #include "crn98_gate.h"
  #undef _
  NV		    		/* total number of variables */
};

/* Enumerate the other (non-gate) variables */
enum {
  #define _(n,i) other_##n,
  #include "crn98_other.h"
  #undef _
  NO		    		/* total number of other variables */
};

/* Enumerate the gates */
enum {
  #define _(n,i) gate_##n,
  #include "crn98_gate.h"
  #undef _
  NG		    		/* total number of gate variables */
};

/* Enumerate the transition rates and other tabulated functions */
enum {
#define _(n,i) _alp_##n,
  #include "crn98_gate.h"
  #undef _
#define _(n,i) _bet_##n,
  #include "crn98_gate.h"
  #undef _
  #define _(n) _##n,
  #include "crn98_fun.h"
  #undef _
  NTAB				/* total number of tabulated functions */
};

/* The structure containing the parameter values for this instance of the model */
typedef struct {                	/* First go the cell parameters as listed in M.C.'s code */
  #define _(name,default) real name;
  #include "crn98_par.h"
  #undef _
  real IV; 			/* The external current. */
  /* 
  ! NB: IV is in millivolt/millisecond rather than picoampere(/picofarad??),			!
  ! and is positive outword raher than positive inward, as the cellular current densities,	!
  ! so e.g. IV can be = D \nabla^2 V, where D is the voltage diffusion coefficient !!		!
  */
} STR;

/* Define all the functions to be tabulated, including gates' transition rates, */
/* for given function of the transmebrane voltage V. */
IONIC_FTAB_HEAD(crn98i) {
  /* /\* Declare the const pars and take their values from struct S==par (a formal parameter) *\/ */
  /* #define _(name,default) DEVICE_CONST(real,name); */
  /* #include "crn98_par.h" */
  /* #undef _ */
  #include "crn98_ftab.h"

  /* Copy the results into the output array values[]. */
  /* Care is taken that all, and only, tabulated functions are attended here. */
#define _(n,i) values[_alp_##n]=alp_##n;
  #include "crn98_gate.h"
  #undef _
  #define _(n,i) values[_bet_##n]=bet_##n;
  #include "crn98_gate.h"
  #undef _
  #define _(n) values[_##n]=n;
  #include "crn98_fun.h"
  #undef _
} IONIC_FTAB_TAIL(crn98i);	  

IONIC_FDDT_HEAD(crn98i,NO,NG) {
  /* Declare the const pars and take their values from struct S==par (a formal parameter) */
  #define _(name,default) DEVICE_CONST(real,name);
  #include "crn98_par.h"
  #undef _
  /* Declare and assign local variables for dynamic variables from state vector */
  /* ..., first for non-gate variables */
  #define _(name,initial) real name=u[var_##name];
  #include "crn98_other.h"
  #undef _
  #define _(name,i) real name=u[var_##name];
  /* ..., and then for gate variables */
  #include "crn98_gate.h"
  #undef _
  /* Get the values of the tabulable functions passed through the formal parameter */
  #define _(name) real name=values[_##name];
  #include "crn98_fun.h"
  #undef _
  /* Calculate the rates of non-gate variables */
  #include "crn98_fddt.h"
  /* Copy the calculated rates into the output array du[].  */
  /* Care is taken that all, and only, non-gating variables are attended here */
  #define _(name,initial) du[other_##name]=d_##name;
  #include "crn98_other.h"
  #undef _
} IONIC_FDDT_TAIL(brc);

/*
 **********************************************************************************
 * The function initialising an instance of the model.                            *
 * This includes: assigning values for all model parameters,                      *
 * i.e. reading from the script or keeping the default values otherwise;      *
 * and creating the vector of initial (steady-state) values of dynamic variables, *
 * which will be used by euler to initialise the whole medium.                    *
 **********************************************************************************
 */
IONIC_CREATE_HEAD(crn98i) {
  int ig; 	/* gates counter */
  double V0;	/* initial voltage */

  I->V_index=0; /* should this be specified here or elsewhere? */

  /* Here we assign the parameter values to the structure AND to namesake local variable */
  #define _(name,default) ACCEPTP(name,default,0,RNONE);
  // #include "crn98_par.h"
  //======================================
_(  fscale1 ,	1.0) //				( 1 )
_(  fscale2 ,	1.0) //				( 1 )
_(  vcell ,	 20100.0) //			( um3 )
_(  vi ,	 vcell*0.68) //			( 1 )
_(  vup ,	 fscale1*0.0552*vcell) //	( 1 )
_(  vrel ,	 fscale2*0.0048*vcell) //	( 1 )
/* _(  T ,	 	 310) //			( Celcius ) */
/* _(  Tfac ,	 3) //				( 1 )  */
_(  Csp ,	 1e+6) //			( pF/cm2 )
_(  Cm ,	 100.0) //			( pF )
/* _(  F ,	 	 96.4867) //			( coul/mmol ) */
/* _(  R ,	 	 8.3143) //			( J K-1 mol-1 ) */
_(  kb ,	 5.4) //			( mM )
_(  nab ,	 140) //			( mM )
_(  cab ,	 1.8) //			( mM )
_(  nac ,	 nab) //			( mM )
_(  cac ,	 cab) //			( mM )
_(  kc ,	 kb) //				( mM )
_(  gna ,	 7.8) //			( nS/pF ) 
_(  gto ,	 0.1652) //			( nS/pF ) 
_(  gkr ,	 0.029411765) //		( nS/pF )
_(  gks ,	 0.12941176) //			( nS/pF )
_(  gcaL ,	 0.12375) //			( nS/pF ) 
_(  ErL ,	 65.0) //			( mV )
_(  gk1 ,	 0.09) //			( nS/pF ) 
_(  gbna ,	 0.0006744375) //		( nS/pF ) 
_(  gbk ,	 0.0) //			( 1 )
_(  gbca ,	 0.001131) //			( nS/pF ) 
_(  inakbar ,	 0.59933874) //			( pA/pF )
_(  kmnai ,	 10.0) //			( mM )
_(  kmko ,	 1.5) //			( mM )
_(  icapbar ,	 0.275) //			( pA/pF )
_(  kmcap ,	 0.0005) //			( mM )
_(  knacalr ,	 1600.0) //			( pA/pF )
_(  kmnalr ,	 87.5) //			( mM )
_(  kmcalr ,	 1.38) //			( mM )
_(  ksatlr ,	 0.1) //			( 1 )
_(  gammalr ,	 0.35) //			( 1 )
_(  trpnbar ,	 0.070) //			( mM )
_(  cmdnbar ,	 0.050) //			( mM )
_(  csqnbar ,	 10) //				( mM )
_(  kmcsqn ,	 0.8) //			( mM )
_(  kmcmdn ,	 0.00238) //			( mM )
_(  kmtrpn ,	 0.0005) //			( mM )
_(  grelbar ,	 30.0) //			( ms-1 ) 
_(  kmup ,	 0.00092) //			( mM ) 
_(  iupbar ,	 0.005) //			( mM/ms )
_(  caupmax ,	 15.0) //			( mM )
_(  kupleak ,	 iupbar/caupmax) //		( ms-1 )
_(  tautr ,	 180.0) //			( ms )
_(  gkur_scale , 1.0) //			( 1 ) /* gkur_scale as in Marshall et al. Pflugers Arch. ??(2011):?? */
  //======================================
  #undef _

  /* Accept the non-cell parameter values by the standard macro */
  ACCEPTP(IV,0,RNONE,RNONE);		/* By default, no stimulation */

  /* Create the vector of initial (steady-state) values of dynamic (state) variables. */
  MALLOC(*u,(long int)NV*sizeof(real));

  /* Assign the initial values as given the following *.h files */
#define _(name,initial) (*u)[var_##name]=initial;
  //  #include "crn98_other.h"
  //======================================
_(V , -8.118e+01)
_(nai , 1.117e+01)
_(cai , 1.013e-04)
_(ki , 1.390e+02)
_(caup , 1.488e+00)
_(carel , 1.488e+00)
/* _(cmdn , 2.042e-03) */
/* _(trpn , 1.180e-02) */
/* _(csqn , 6.504e+00) */
_(fca , 7.755e-01)
_(uu, 2.350e-112)	/* u in M.C.'s code */
_(vv, 1.000e+00)	/* v in M.C.'s code */
_(ww, 9.992e-01)	/* w in M.C.'s code */
  //======================================
  #include "crn98_gate.h"
  #undef _
} IONIC_CREATE_TAIL(crn98i,NV)


/* /\* */
/*  ******************************************************************************* */
/*  * The function "computing the right-hand sides".                              * */
/*  * Actually, simply recalculates the value of the variables for the next step. * */
/*  * Made as an rhs module to take advantage of parameter substitution mechanism.* */
/*  ******************************************************************************* */
/*  *\/ */
/* RHS_HEAD(crn,N) { */
/*   /\* Declare the constant parameters and take their values from structure S==par (a formal parameter) *\/ */
/*   #define _(name,default,unit) DEVICE_CONST(real,name) */
/*   #include "crn98_par.h" */
/*   #undef _ */
/*   DEVICE_CONST(real,IV) */
/*   DEVICE_CONST(real,ht) */
/*   /\* Declare the dynamic variables and give them values from the vector u[] (a formal parameter) *\/ */
/*   #define _(name,initial) real name=u[crn98_##name]; */
/*   #include "crn98_var.h" */
/*   #undef _ */
/*   /\* The time step is called ht everywhere in Beatbox and dt in M.C.'s code *\/ */
/*   real dt=ht; */

/*   /\* Here go all actual computations of M.C.'s code *\/ */
/*   #include "crn98_step.h" */

/*   /\* Add the external current's effect *\/ */
/*   V += ht*IV; */

/*   /\* Copy the updated values of dynamic variables back to the vector u[].  *\/ */
/*   /\* Assign the derivatives to zero as the values are already fully updated *\/ */
/*   #define _(name,initial) u[crn98_##name]=name; du[crn98_##name]=0; */
/*   #include "crn98_var.h" */
/*   #undef _ */
/* } RHS_TAIL(crn) */
