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

/* Record values in the given 4D cube to a text file */

#if MPI
#include <mpi.h>
#endif
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "system.h"
#include <limits.h>

#include "beatbox.h"
#include "state.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"
#include "mpi_io_choice.h"

#define RECORD_SEPARATOR_LENGTH 2
typedef struct {
  char filename[MAXPATH];
  int append;
  char filehead[1024];
  char format[1024];
  int formatwidth;
  char vsep[RECORD_SEPARATOR_LENGTH +1];
  char zsep[RECORD_SEPARATOR_LENGTH +1];
  char ysep[RECORD_SEPARATOR_LENGTH +1];
  char xsep[RECORD_SEPARATOR_LENGTH +1];
  char recordsep[RECORD_SEPARATOR_LENGTH +1];
  int timestamp;
#if MPI
  MPI_File file;
  MPI_Datatype V_Type;
  int value_length;
#else
  FILE *file;
#endif
} STR;

RUN_HEAD(record)
  DEVICE_ARRAY(char,vsep)
  DEVICE_ARRAY(char,zsep)
  DEVICE_ARRAY(char,ysep)
  DEVICE_ARRAY(char,xsep)
  DEVICE_ARRAY(char,recordsep)
  DEVICE_ARRAY(char,format)
  DEVICE_CONST(int,formatwidth)
  DEVICE_CONST(char,timestamp)
  int v, z, y, x;
  char *nextSep;

#if MPI
  DEVICE_CONST(MPI_File, file)
  DEVICE_CONST(MPI_Datatype, V_Type)
  int buf_length = S->value_length + 1;
  char buf[buf_length];
#else
  DEVICE_CONST(FILE *,file)
  if (!file) return 1;
#endif

/*  Write data, selecting hierarchical separator */
/*  Record->X->Y->Z->V (most->least significant). */

  if (timestamp) fprintf(file, "%-6ld\t", (long)t);

  for (x=s.x0; x<=s.x1; x++) {
    for (y=s.y0; y<=s.y1; y++) {
      for (z=s.z0; z<=s.z1; z++) {
	for (v=s.v0; v<=s.v1; v++) {
	  if (v!=s.v1) {
	    nextSep = vsep;
	  } else {
	    if (z!=s.z1) {
	      nextSep = zsep;
	    } else {
	      if (y!=s.y1) {
		nextSep = ysep;
	      } else {
		/* if (x!=(xmax-1)) {} // xmax is global, so thats how many entries in the 1D you can make. This is a special case, you cannot choose to record a "patch" of the 1D or any other. VNB 2012/12/27: I found this cryptic and rewrote as before */
		if (x!=s.x1) {
		  nextSep = xsep;
		} else {
		  nextSep = recordsep;
		}
	      }
	    }
	  }
	  sprintf(buf, format, New[ind(x,y,z,v)]);
	  strcat(buf, nextSep);
#if MPI
	  mpi_errno = MPI_File_write(file, buf, 1, V_Type, MPI_STATUS_IGNORE);
	  CHECK_MPI_SUCCESS("Couldn't write to file.");
#else
	  fprintf(file, "%s", buf);
#endif
	}
      }
    }
  }
#if MPI==0
  FFLUSH(file);
#endif
RUN_TAIL(record)

DESTROY_HEAD(record)
#if MPI
  if (s.runHere) MPI_File_close(&(S->file));
#else
  if (S->file) fclose(S->file); S->file=NULL;
#endif
DESTROY_TAIL(record)


/* ----------------------------------------------------------------- */
CREATE_HEAD(record)
  ACCEPTI(append,1,0,1);
  ACCEPTI(timestamp,0,0,1);
  ACCEPTS(filehead,"");
  /* Add newline to filehead if it is nonempty and does not have one already */
  if (S->filehead[0]!='\0' && S->filehead[strlen(S->filehead)-1]!='\n') {
    strcat(S->filehead,"\n");
    MESSAGE("/* amended filehead=\"%s\" */\n",S->filehead);
  }
	
  #define CHECK_SEPARATOR(s)						\
    if (strlen(S->s) > RECORD_SEPARATOR_LENGTH) {			\
      MESSAGE("WARNING from %s: %s is longer than the required separator size of %d and will be truncated.\n", dev->n, #s, RECORD_SEPARATOR_LENGTH); \
    } else if(strlen(S->s) < RECORD_SEPARATOR_LENGTH) {			\
      int i;								\
      for (i=strlen(S->s); i<RECORD_SEPARATOR_LENGTH; i++) 		\
	S->s[i]=' ';							\
    }									\
    S->s[RECORD_SEPARATOR_LENGTH] = '\0';

  ACCEPTS(vsep,"  "); CHECK_SEPARATOR(vsep)
  ACCEPTS(zsep,", "); CHECK_SEPARATOR(zsep)
  ACCEPTS(ysep,"; "); CHECK_SEPARATOR(ysep)
  ACCEPTS(xsep," \n"); CHECK_SEPARATOR(xsep)
  ACCEPTS(recordsep," \n"); CHECK_SEPARATOR(recordsep)
  
  #undef CHECK_SEPARATOR

  ACCEPTS(format,"%38f");
  ACCEPTI(formatwidth,38,2,INONE);
  /* check the consistency of the two on some arbitrary chosen numbers */
  {
    int i, l;
    char s[1024];
    real testvalue[6]={0, 1, 3.14159265358, -1, 0.1, -1.e6};
    for (i=0;i<6;i++) {
      sprintf(s,S->format,testvalue[i]);
      l=strlen(s);
      if (S->formatwidth!=l) {
	MESSAGE("/* Warning: test value "REALF" printed by format '%s' produces string of length %d not %d\n",
		S->format,l,S->formatwidth);
	break;
      }
    }
  }

#if MPI
  if (!accepts("file=",&(S->filename[0]),"",w)) return(0);

  /*  Temporary Datatypes for file I/O */
  MPI_Datatype V_Type, Point_Type, Z_Type, Y_Type, X_Type;

  /*  For defining MPI datatypes */
  int *blocklengths = calloc(3, sizeof(int));
  MPI_Aint *disp = calloc(3, sizeof(MPI_Aint));
  MPI_Datatype *old_types = calloc(3, sizeof(MPI_Datatype));

  /*  Extents */
  MPI_Aint char_extent,v_extent,point_extent,z_extent,y_extent,x_extent,lb;
  int num_vars = (dev->s.v1 - dev->s.v0) + 1;
  int count;

  /*  Lengths */
  int filehead_length = strlen(S->filehead);
  int value_length = S->formatwidth + RECORD_SEPARATOR_LENGTH;
  MPI_Type_get_extent(MPI_CHAR,&lb,&char_extent);

  /*  Represents a single value with separator. */
  count = value_length;
  mpi_errno = MPI_Type_contiguous(count,MPI_CHAR,&V_Type);
  CHECK_MPI_SUCCESS("Couldn't create V_Type.")
  mpi_errno = MPI_Type_commit(&V_Type);
  CHECK_MPI_SUCCESS("Couldn't commit V_Type.")
  MPI_Type_get_extent(V_Type, &lb, &v_extent);

  /*  Represents one point, i.e. num_vars values with separators. */
  count = num_vars;
  mpi_errno = MPI_Type_contiguous(count,V_Type,&Point_Type);
  CHECK_MPI_SUCCESS("Couldn't create Point_Type.");
  mpi_errno = MPI_Type_commit(&Point_Type);
  CHECK_MPI_SUCCESS("Couldn't commit Point_Type.");
  MPI_Type_get_extent(Point_Type, &lb, &point_extent);

  /*  Struct representing local contribution to ZV plane */
  blocklengths[0] = 1;
  blocklengths[1] = (dev->s.z1 - dev->s.z0) + 1;
  blocklengths[2] = 1;

  disp[0] = 0;
  disp[1] = (dev->s.z0 - dev->s.global_z0) * point_extent;
  disp[2] = ((dev->s.global_z1 - dev->s.global_z0) + 1) * point_extent;

  old_types[0] = MPI_LB;
  old_types[1] = Point_Type;
  old_types[2] = MPI_UB;

  mpi_errno = MPI_Type_create_struct(3, blocklengths, disp, old_types, &Z_Type);
  CHECK_MPI_SUCCESS("Couldn't create Z_Type.")
  mpi_errno = MPI_Type_commit(&Z_Type);
  CHECK_MPI_SUCCESS("Couldn't commit Z_Type.")
  MPI_Type_get_extent(Z_Type, &lb, &z_extent);

  /*  Struct representing local contribution to YZV hyperplane */
  blocklengths[0] = 1;
  blocklengths[1] = (dev->s.y1 - dev->s.y0) + 1;
  blocklengths[2] = 1;

  disp[0] = 0;
  disp[1] = (dev->s.y0 - dev->s.global_y0) * z_extent;
  disp[2] = ((dev->s.global_y1 - dev->s.global_y0) + 1) * z_extent;

  old_types[0] = MPI_LB;
  old_types[1] = Z_Type;
  old_types[2] = MPI_UB;

  mpi_errno = MPI_Type_create_struct(3, blocklengths, disp, old_types, &Y_Type);
  CHECK_MPI_SUCCESS("Couldn't create Y_Type.");
  mpi_errno = MPI_Type_commit(&Y_Type);
  CHECK_MPI_SUCCESS("Couldn't commit Y_Type.");
  MPI_Type_get_extent(Y_Type, &lb, &y_extent);

  /*  Struct representing local contribution to XYZV hypercube */
  blocklengths[0] = 1;
  blocklengths[1] = (dev->s.x1 - dev->s.x0) + 1;
  blocklengths[2] = 1;

  disp[0] = 0;
  disp[1] = (dev->s.x0 - dev->s.global_x0) * y_extent;
  disp[2] = ((dev->s.global_x1 - dev->s.global_x0) + 1) * y_extent;
  
  old_types[0] = MPI_LB;
  old_types[1] = Y_Type;
  old_types[2] = MPI_UB;

  mpi_errno = MPI_Type_create_struct(3, blocklengths, disp, old_types, &X_Type);
  CHECK_MPI_SUCCESS("Couldn't create X_Type.");
  mpi_errno = MPI_Type_commit(&X_Type);
  CHECK_MPI_SUCCESS("Couldn't commit X_Type.");
  MPI_Type_get_extent(X_Type, &lb, &x_extent);
  /*---------------------------------------------------------------------*/

  /*  Communicator for active dump instances. */
  MPI_Comm comm_record;
  int firstRank; /*  Rank of first processor in communicator, for writing headers. */
  
  /*  File */
  MPI_File file;
  int mode;
  int displacement=0;
  MPI_Offset filesize; /*  existing file size, for appending. */
  
  if (!deviceCommunicatorWithFirstRank(dev->s.runHere, &comm_record, &firstRank))
    ABORT("Could not create communicator.");
	
  if (dev->s.runHere) { /*  Active devices only!	 */
    /*  Set access mode */
    if (S->append) {
      mode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND;
    } else {
      mode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
    }
	
    /*  Open file and set view */
    mpi_errno = MPI_File_open(comm_record, S->filename, mode, MPI_INFO_NULL, &file);		
    
    /*  If requested by user, append to EOF */
    if (S->append) {
      MPI_File_get_size(file, &filesize);
      displacement = (int)filesize;
    }
	
    /*  Write file header, if we have one. */
    if (*S->filehead) {
      mpi_errno = MPI_File_set_view(file, displacement, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
      
      if (mpi_rank == firstRank) {
	MPI_File_write(file, S->filehead, filehead_length, MPI_CHAR, MPI_STATUS_IGNORE);
      }
      displacement += filehead_length;
    }
		
    /*  Have to reset view for X_Type. */
    mpi_errno = MPI_File_set_view(file, displacement, MPI_CHAR, X_Type, "native", MPI_INFO_NULL);
    
    /*  Assign the file to the parameter structure */
    S->file = file;
    S->value_length = value_length;
    S->V_Type = V_Type;
  }

#else
  /*  Sequential version */
  ACCEPTF(file,S->append?"at":"wt",NULL);
  if (*S->filehead) {
    fprintf(S->file,"%s",S->filehead);
    FFLUSH(S->file);	
  }
#endif
CREATE_TAIL(record,0)

#undef RECORD_SEPARATOR_LENGTH
