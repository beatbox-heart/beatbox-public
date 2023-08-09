/**
 * Copyright (C) (2010-2023) Vadim Biktashev, Irina Biktasheva et al. 
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

/* INITIALIZATION: READING PARAMETERS AND FORMING THE COMPUTATION LOOP */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "system.h"

#include "beatbox.h"
#include "device.h"
#include "screen.h"
#include "state.h"
#include "init.h"
#include "qpp.h"
#include "bikt.h"
#include "k_.h"

/* Declare list of all possible devices */
#define D(name) Create create_##name;
#define S(name) Create create_##name;
#include "devlist.h"
#undef D
#undef S

/* Imported variables */
extern FILE     *debug;
extern int      Verbose;
extern INT	Graph;
extern int      idev;
extern int      ndev;
extern Device   dev[];
extern char *device_name;
extern long int inf;

static int create (Name name, Create c, Device *d, char *rest)
{
  if (!accept_condition(&(d->c),rest)) return 0;
  if (!accept_space(&(d->s),rest)) return 0;
  if (!accepts("name=",&(d->n[0]),name,rest)) return 0;
  if (strcmp(d->n,name)) if(!def_dev(d)) return 0;
#if MPI
  d->sync      = 0;
  d->alwaysRun = 0;
#endif
  return(c(d,rest));
}

int init (void) {
  static char s[MAXSTRLEN];
  char *rest;
  char *w;
  int begun;
  int state_called;

  step_already_created=0;

  // this function defines tolerances, and also several functions.
  if (!init_const()) EXPECTED_ERROR("initializing constants"); 
  
  idev=0; begun=0; state_called=0;
  for(;;) {
    // this is what reads in the text line.
    if NOT(read_command(s,MAXSTRLEN)) EXPECTED_ERROR("reading command");
    if (s[0]=='\0') continue;               /* null command - OK */
    if (s[0]==(char)EOF) break;             /* top level EOF */
    rest=first_word(s,&w," \t\n\r;$"); /* extract command */
    device_name=w;
    DEBUG("\n#%d %s:",mpi_rank,w);
    #define CASE(a) else if (0==stricmp(w,a)) 
    if NOT(*w) continue;    /* empty cmd => eof - pass to next read_command */
    CASE("rem") continue;   /* comment command */
    CASE("if") {                    /* conditional statement */
      double condition;
      rest=first_word(rest,&w," \t\n\r;$");
      if (!calc(&condition,t_real,w)) goto ERR_CREATE;
      if (!condition) continue;
      /* this parses the rest of the command to headword + further rest */
      rest=first_word(rest,&w," \t\n\r;$"); 
    }
    CASE("def") {                   /* definintion */
      /* this reads definition of a k-variable in the string rest */
      if(!def(rest)) goto ERR_CREATE; 
      begun=1;
      continue;
    }
    CASE("state") {                 /* create the state arrays */
      if (begun) MESSAGE("\n");
      MESSAGE("state ");
      /* this creates the computational grid. */
      if(!state(rest)) goto ERR_CREATE; 
      MESSAGE("$");
      begun=0;
      state_called = 1;
      continue;
    }
    CASE("screen") {                /* create the BGI window */
#if MPI
      MESSAGE("\n/* The screen command is disabled when using MPI. Your simulation will continue without it. */");
#else
      if (Graph) {
	MESSAGE("\nscreen ");
	// this initialises the screen.
	if (makescreen(rest)) {
	  MESSAGE("$");
	} else {
	  MESSAGE("\n/* Will proceed without onscreen graphics */ $");
	}
      } else {
	MESSAGE("\n/* With nograph option, the 'screen' command is ignored */");
      }
      begun=0;
#endif
      continue;
    }    
    CASE("end") break;              /* end of input stream */
    
    if (idev>=MAXDEV) {
      MESSAGE("\n The number of devices in the input file exceeded the maximum of %d allowed in this compile.\n",MAXDEV);
      goto ERR_CREATE;
    }
    
#if MPI
    /*  Don't create devices on idle processes. */
    #define D(name)                                                     \
    CASE(#name) {                                                       \
      if(!I_AM_IDLE){                                                   \
        if(!state_called) {                                             \
          MESSAGE("\nDevices cannot be created before state is called."); \
          goto ERR_CREATE;                                              \
        }                                                               \
        MESSAGE("\n%s ",#name);                                        \
        if(!create(#name,create_##name,dev+idev,rest)) goto ERR_CREATE; \
        MESSAGE("$");                                                	\
        begun=0;                                                        \
      }                                                                 \
    }
    
    /*  Disable Sequential-only devices. */
    #define S(name)                                                     \
    CASE(#name) {                                                       \
      if(!I_AM_IDLE){                                                   \
        if(!state_called) {                                             \
          MESSAGE("\nDevices cannot be created before state is called."); \
          goto ERR_CREATE;                                              \
        }                                                               \
        MESSAGE("\n/* The %s device is disabled when using MPI. Your simulation will continue without it. */",#name); \
        begun=0;                                                        \
        continue;                                                       \
      }                                                                 \
    }
#else
    #define D(name)                                                     \
    CASE(#name) {                                                       \
      if(!state_called) {                                               \
        MESSAGE("\nDevices cannot be created before state is called."); \
        goto ERR_CREATE;                                                \
      }                                                                 \
      MESSAGE("\n%s ",#name);						\
      if(!create(#name,create_##name,dev+idev,rest)) goto ERR_CREATE;   \
      MESSAGE("$");                                                  	\
      begun=0;                                                          \
    }
      
    #define S(name)                                                     \
    CASE(#name) {                                                       \
      if(!state_called) {                                               \
        MESSAGE("\nDevices cannot be created before state is called."); \
        goto ERR_CREATE;                                                \
      }                                                                 \
      MESSAGE("\n%s ",#name);                                           \
      if(!create(#name,create_##name,dev+idev,rest)) goto ERR_CREATE;   \
      MESSAGE("$");                                                  	\
      begun=0;                                                          \
    }
#endif

    #include "devlist.h"
    #undef D
    #undef S
    #undef eq
    else { MESSAGE("\nUnknown device name '%s'\n",w); goto ERR_CREATE; }
    idev++;
  } /*  for(;;) */
  DEBUG("\n#%d end of input file\n",mpi_rank);

  MESSAGE("\nend of input file $");
  MESSAGE("\nLoop of %d devices created:",ndev=idev);
  for (idev=0;idev<ndev;idev++) {
    if (idev%5==0) MESSAGE("\x01""\n  ");
        MESSAGE(" (%d)%s",idev,dev[idev].n);
  }
  MESSAGE(" $");

  return(1);

ERR_CREATE:
  MESSAGE(
    "\nError occurred in input file \"%s\" command \"%s\" ending at line %ld pos %ld",
    inname[depth-1], w, inpline[depth-1], inppos[depth-1]
  );
  return(0);
}

void term(void) {               /* free all dynamical objects */
  Device d;
  if (ndev) for (idev=ndev-1;idev>=0;idev--) {
    d=dev[idev];
#if MPI
    d.d(d.s,d.par,d.sync,d.alwaysRun);
#else
    d.d(d.s,d.par);
#endif
  }
  state_free();                 /* free the state array(s) */
  term_const();                 /* clear the symbol table and free variables */
}


