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

/* Definitions related to rhs-type cell model definitions */

#ifndef _rhs
#define _rhs
typedef struct {		/* description of dependent parameters */
  int n;
  int *src;
  real **dst;
} Var;

#define RHSPROC(name) int name(real *u, real *du, Par par, Var var, int ln)
typedef RHSPROC(RhsProc);
#define RHS_HEAD(name,LN) \
RHSPROC(name) {           \
  STR *S = (STR *)par;	  \
  int ivar;               \
  assert(ln==LN);         \
  if (var.n) for(ivar=0;ivar<var.n;ivar++) *(var.dst[ivar])=u[var.src[ivar]];
#define RHS_TAIL(name)	\
  return 1;             \
}

typedef struct {
  real ht;		/* time step */
  RhsProc  *f;		/* function of right-hand sides */
  real *u;		/* vector of vars for steady state if any */
  real *du;		/* vector of derivatives + other working memory */
  Par p;		/* vector of ODE par-s */
  Name ode;		/* name of the ODE system */
  int rest;		/* how many steps to do to find resting values */
  Var var;		/* description of dependent parameters */
} rhs_str;

#define RHSCREATE(name) int name(Par *par, Var *var, char *w, real **u, int v0)
typedef RHSCREATE(RhsCreate);

#define RHS_CREATE_HEAD(name)				\
RHSCREATE(create_##name) {				\
  STR *S = (STR *)Calloc(1,sizeof(STR));        	\
  char *ptr=w;                                    	\
  int ivar=0;						\
  if (!S) ABORT("cannot create %s",#name);		\
  for(var->n=0;*ptr;var->n+=(*(ptr++)==AT));		\
  if(var->n){CALLOC(var->dst,var->n,sizeof(real *));	\
  CALLOC(var->src,var->n,sizeof(int));}			\
  else{var->src=NULL;(var->dst)=NULL;}

#define RHS_CREATE_TAIL(name,rc)    	\
  var->n=ivar;				\
  if(ivar){REALLOC(var->dst,1L*ivar*sizeof(real *));\
  REALLOC(var->src,1L*ivar*sizeof(int));}	\
  else{FREE(var->dst);FREE(var->src);}	\
  *par = S;				\
  return rc;				\
}

#endif

