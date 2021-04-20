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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
			    #ifdef gamma
			    #undef gamma
			    #endif

#include "system.h"
#include "beatbox.h"
#include "GPV.on"

#if GPV
#include "device.h"
#include "state.h"
#include "bikt.h"
#include "qpp.h"
/*#include "def.h"*/
/*#include "laplac.h"*/

static real sqr(real a) {return a*a;}
static real cub(real a) {return a*a*a;}
					/* UNITS */
#define _(n,m) static real n=1.;
#include "gpv_un.h"
#undef _				
					/* INDICES OF DYNAMIC VARIABLES */
enum {
  #define _(n,i) gpv_##n,
  #include "gpv_var.h"
  #undef _
  N		    /* total number of variables */
};
					/* PARAMETERS */
#define _(n,v,u) static real n;
#include "gpv_par.h"
#undef _
					/* GATING RATES */
#define alpm	(200.*(V+41.)/(1.-exp(-0.1*(V+41.))))
#define betm	(8000.*exp(-0.056*(V+66.)))
#define alph	(20.*exp(-0.125*(V+75.)))
#define beth	(2000./(1+320.*exp(-0.1*(V+75))))
#define alpd	(90.*(V+19.)/(1.-exp(-(V+19.)/4.)))
#define betd	(36.*(V+19.)/(exp((V+19.)/10.)-1.))
#define alpx	(0.5*exp(0.0826*(V+50.))/(1.+exp(0.057*(V+50.))))
#define betx	(1.3*exp(-0.06*(V+20.))/(1.+exp(-0.04*(V+20.))))
#define alpf	(3.125*(V+34.)/(exp((V+34.)/4.)-1.))
#define betf	(25./(1.+exp(-(V+34.)/4.)))
#define kq	(333.)
#define q0	(1./(1.+exp(-(V+4.)/5.)))
#define alpq	(kq*q0)
#define betq	(kq*(1.-q0))
#define alpr	(0.033*exp(-V/17.))
#define betr	(33./(1.+exp(-(V+10.)/8.)))

real gate_dt(real g0, real alp, real bet, real dt) {
  real lam=alp+bet;
  real ginf=alp/lam;
  return ginf+(g0-ginf)*exp(-lam*dt);    
}

int gpv_step(real u[], real dt) {
					/* DYNAMIC VARIABLES DEFINITIONS */
  #define _(n,i) real n=u[__LINE__-1];
  #include "gpv_var.h"
  #undef _
					/* CURRENTS DEFINITIONS */
  #define _(n,v) real n;
  #include "gpv_cur.h"
  #undef _
			    /*printf("gpv_step entered\n");*/
  
					/* GATES */
  #define gate(g) u[gpv_##g]=gate_dt(g,alp##g,bet##g,dt);
  gate(m); 
  gate(h); 
  gate(d);
  gate(x);
  gate(f);
  gate(q);
  gate(r);
					/* CURRENTS */
  #define _(n,v) n=v;
  #include "gpv_cur.h"
  #undef _
					/* VOLTAGE */
  #define D(n,a) u[gpv_##n]+=dt*(a)
  #define S(n) u[gpv_##n]=n
  D(V, (Ik + Iki + Ito + IsiK + IbK + INaK + INa + IbNa + IsiNa + INaCa + IsiCa + IbCa) / Cm );
  
					/* CONCENTRATIONS */
  D(Nai, - (INa + IbNa*Nao/140. + 3*INaK + 3*INaCa + IsiNa) / (Vi*F));
  D(Ki,    - (Ik + Iki + IsiK + IbK + Ito - 2*INaK) / (Vi*F));
  D(Caup,   Vi/Vsrup*Iup - Itr);
  D(Carel,   Vup/Vrel*Itr - Irel);
  
  Cai -= dt*( (IsiCa + IbCa - 2*INaCa)/(2*Vi*F) - Iup + Irel*(Vsrup*Vrel)/(Vi*Vup));
  do {
    real Catrop0 = Catrop;
    real Cacalmod0 = Cacalmod;
    real Catot=Cai+Catrop+Cacalmod;
    Catrop   = Ctrop * (Cacalmod - Catot) / (Catrop   - Catot - Ctrop + Cacalmod - 200. / 1.e5);
    Cacalmod = Mtrop * (Catrop   - Catot) / (Cacalmod - Catot - Mtrop + Catrop   -  50. / 1.e5);
    Cai	= Catot-Catrop-Cacalmod;
    if( fabs((Catrop0 - Catrop) / Catrop) + fabs((Cacalmod0 - Cacalmod) / Cacalmod)<1.e-3) break;
    /*if (Catrop0==Catrop && Cacalmod0==Cacalmod) break;*/
									/*printf(".");*/
									/*printf("%lf\n", Cai);*/
  } while(1);
									/*printf("\n\n\n\n\n\n");*/
  S(Cai);
  S(Cacalmod);
  S(Catrop);
  
  
  {
    real tmp=500.*sqr(Cai/(Cai+kmCa));
    D(fact,   (1. - fact - fprod) * tmp - fact*(tmp + 60.));
    D(fprod,   fact * (tmp + 60.) - fprod);
  }
  #undef gate
  #undef D
			    /*MESSAGE("gpv_step exiting");*/
  return 1;
}

static real D, ht, hx;

DEVICE_PROC(gpv){
  int x, y, z, v;
  real u[N];
			    /*printf("gpv entered\n");*/
  if (!New || !Aux) ERROR0("states were not allocated");
  #define vNew 0
  #define vAux 0
  Laplacian(X0,X1,Y0,Y1,Z0,Z1,vNew,vAux,D,hx,ht);
  for (z=Z0;z<=Z1;z++)
  for (y=Y0;y<=Y1;y++)
  for (x=X0;x<=X1;x++) {
			    /*printf("x=%d y=%d z=%d N=%d\n", x, y, z, N);*/
    for (v=0;v<N;v++) u[v]=New[ind(x,y,z,v)];
			    /*printf("calling gpv_step\n");*/
    if NOT(gpv_step(u,ht)) ERROR3("error calculating gpv_step at (%d,%d,%d)",x,y,z);
    for (v=0;v<=N;v++) New[ind(x,y,z,v)] = u[v];
    New[ind(x,y,z,vNew)] += Aux[ind(x,y,z,vAux)];
  }
			    /*printf("gpv exiting\n");*/
  return(1);
}

CREATE_DEVICE(create_gpv) {
  real gam;
  real u[N];
  int x, y, z, v;
  
  assert(vmax>=N); 
  assert(vaux>=1);
  if (step_already_created) ERROR0("step must be defined only once");
  step_already_created=1;
  ACCEPTR(ht ,RNONE,RSUCC(0.),RNONE);
  ACCEPTR(hx ,RNONE,RSUCC(0.),RNONE);
  ACCEPTR(D  , 0.04,0 /*RSUCC(0.)*/,RNONE);
  gam = ht * dim * D / (hx*hx);
  if (gam >= 1.0) {
    MESSAGE0("\nWarning: diffusion stability criterion violated: gam=");
    MESSAGE1(REAL,gam);   
  }
  #define _(n,v,u) ACCEPTR(n,v*u,0.,RNONE);
  #include "gpv_par.h"
  #undef _
  
  #define _(n,i) u[__LINE__-1]=i;
  #include "gpv_var.h"
  #undef _
  
  for(z=0;z<zmax;z++) for(y=0;y<ymax;y++) for(x=0;x<xmax;x++) for(v=0;v<N;v++) New[ind(x,y,z,v)]=u[v];
  
  dev->proc = gpv;
  dev->par  = NULL;
  strcpy(dev->name,"gpv");

  return(1);
}
#endif /* GPV */



