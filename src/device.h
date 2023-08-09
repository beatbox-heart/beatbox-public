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

#ifndef _DEVICE_H_
#define _DEVICE_H_
#ifdef _OWN
#undef _OWN
#define OWN
#endif
#include "extern.h"
#include "error.h"

#if MPI
#include <mpi.h>
#endif

extern FILE *res;
extern FILE *debug;

#define MAXDEV 1024    		/* Maximum number of Devices that can be used */
#define MAXSTRLEN 8192 		/* length of various buffers, including that for reading input file */
EXTERN char buf[MAXSTRLEN];	/* Common general-purpose buffer */

typedef double *Condition; 	/* Condition when the Device will operate */

typedef struct {
  int x0, x1,   /* Set the limits on which the Device will operate */
      y0, y1,
      z0, z1,
      v0, v1,		
      global_x0, global_x1, /* Local == global in seq mode */
      global_y0, global_y1,
      global_z0, global_z1,
#if MPI
      runHere,
#endif
      nowhere;
} Space;

/* Graphics display information for the Device */
typedef struct {int row0,col0,row1,col1,color,area;} BGIWindow;

/* Internal persistent storage used by the Device */
typedef void *Par;

#if MPI
/**********************/
/* PARALLEL */
#define PROC(name) int name(Space s,Par par,int sync,int alwaysRun)
#else
#define PROC(name) int name(Space s,Par par)
#endif
typedef PROC(Proc);

/* Device name */
typedef char Name[32];

/* Device */
typedef struct {	
    Condition c; /* When the device is to be run */
    Space s;     /* Coordinate bounds on which the device operates */
    Par par;     /* Store for persistent data within the device */
    Proc *p;     /* Pointer to the device's run function */
    Proc *d;     /* Pointer to the device's destroy function */
    Name n;      /* Device name */
#if MPI
    int sync;
    int alwaysRun;
#endif
} Device;

typedef int Create (Device *dev, char *w);

/* Shortcuts for std pieces of device codes */

#define DEVICE_CONST(type,name) type name=S->name;
#define DEVICE_VAR(type,name)   type *name=&(S->name);
#define DEVICE_ARRAY(type,name) type *name=&(S->name[0]);
#define DEVICE_PAR(type,name)   type name=(S->name##ptr)?(S->name=*(type *)(S->name##ptr)):(memcpy(&(S->name),execute(S->name##code),sizeof(type)),S->name);
#define PAR(type,name) type name; type *name##ptr; pp_fn name##code;

#if MPI
#define RUN_HEAD(name)		 \
PROC (run_##name) {		 \
   if (sync) {haloSwap();}	 \
   if (s.runHere || alwaysRun) { \
      STR *S = (STR *) par;

#define RUN_TAIL(name)		 \
  } /* if s.runHere */		 \
  return 1;			 \
}

#define DEVICE_REQUIRES_SYNC dev->sync = 1;

#define DEVICE_ALWAYS_RUNS dev->alwaysRun = 1;

/* Sync is set to 0 to prevent haloSwaps in delegated devices
 * This isn't just for efficiency, but also avoids deadlocks if
 * any device's runHere is false.
 */
#define DELEGATE_TO_DEVICE(name) run_##name(s,w,par,0,alwaysRun);

#else
/**********************/
/* SEQUENTIAL */

#define RUN_HEAD(name)	\
PROC (run_##name) {	\
  STR *S = (STR *) par;

#define RUN_TAIL(name)	\
  return 1;		\
}

#define DEVICE_REQUIRES_SYNC
#define DEVICE_ALWAYS_RUNS

#define DELEGATE_TO_DEVICE(name) run_##name(s,w,par);

#endif
/**********************/

#define DESTROY_HEAD(name)	\
PROC (destroy_##name) {		\
  STR *S = (STR *) par;

#define DESTROY_TAIL(name)	\
  FREE(par);			\
  return 1;			\
}

#define SAFE_CLOSE(file)        \
  if ((file) != NULL && 	\
      (file) != res && 		\
      (file) != debug && 	\
      (file) != stdout && 	\
      (file) != stderr && 	\
      ftell(file) !=-1 		\
      ) {			\
    fseek((file),0,SEEK_END);	\
    fflush(file);		\
    fclose(file);		\
  }				\
  (file)=NULL;


#define CREATE_HEAD(name)			\
int create_##name (Device *dev, char *w)  {	\
  STR *S = (STR *) Calloc(1,sizeof(STR));	\
  if (!S) ABORT("cannot create %s",#name);

#define CREATE_TAIL(name,areaval)	\
  dev->p = (Proc *) run_##name;		\
  dev->d = (Proc *) destroy_##name;	\
  dev->par  = S;			\
  return(1);				\
}

#define DEVICE_IS_SPACELESS	    \
if (spaceParametersExist(w)) {	    \
  MESSAGE("WARNING: %s does not use space parameters (x0,x1,y0,y1,z0,z1).\nThe parameters provided will be ignored.",dev->n); \
}

#define DEVICE_HAS_DEFAULT_SPACE	  \
if (spaceParametersExist(w)) {		  \
  EXPECTED_ERROR("%s insists on using the default space.\nPlease try again without space parameters (x0,x1,y0,y1,z0,z1).\n",dev->n); \
}

#define DEVICE_MUST_BE_NOWHERE		\
if (dev->s.nowhere != 1) {		\
  EXPECTED_ERROR("%s requires that nowhere is equal to 1. Try again, adding 'nowhere=1' to the device parameters.\n",dev->n); \
}

#define DEVICE_OPERATES_ON_A_SINGLE_POINT	\
if (						\
  (find_key("x1=",w)) ||			\
  (find_key("y1=",w)) ||			\
  (find_key("z1=",w)) ||			\
  (find_key("v1=",w))				\
) {					       	\
  MESSAGE("WARNING: %S does not use the x1,y1,z1 or v1 space parameters",dev->n); \
  MESSAGE(" -- it only operates on a single point.");			\
  MESSAGE(" The parameter(s) provided will be ignored.\n");		\
}

/*  Some handy utilities for devices running on MPI. */
#if MPI

/*  Error handling variables. */
EXTERN MPI_Status 	status;
EXTERN int 		mpi_errno;
EXTERN char 		error_string[MPI_MAX_ERROR_STRING];
EXTERN int 		error_length;

#define CHECK_MPI_SUCCESS(error_message)                      \
  if(mpi_errno != MPI_SUCCESS) {                              \
    MPI_Error_string(mpi_errno, error_string, &error_length); \
    ABORT("Process %d: %s\n\tMPI Message: %s\n", mpi_rank, error_message, error_string); \
  }

/* this defines its own buf in case there is another one outside */
#define MPIDO(cmd,...)							\
  mpi_errno = cmd;							\
  if (mpi_errno != MPI_SUCCESS) {					\
    MPI_Error_string(mpi_errno, error_string, &error_length);		\
    char buf[MAXSTRLEN];						\
    sprintf(buf,__VA_ARGS__);						\
    ABORT("%s\nMPI Message: %s\n",buf,error_string);			\
  }

/* 	Define new group and communicator for I/O, 
 *	including only instances with runHere = true.	*/
int deviceCommunicator(int runHere, MPI_Comm *new_comm);
int deviceCommunicatorWithFirstRank(int runHere, MPI_Comm *new_comm, int *first);

/*  	Exchange internal boundaries. */
int haloSwap(void);

#endif /*  MPI */
#undef EXTERN
#endif /* end of include guard: _DEVICE_H_ */

