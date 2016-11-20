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


#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "beatbox.h"
#include "error.h"
#include "state.h"
#include "device.h"
#define OWN
#include "geometry.h"
#undef OWN
#include "mpi_io_choice.h"

extern int Verbose;

/*  Shifts required to put first tissue point at ONE,TWO,TRI. */
static int x_offset,y_offset,z_offset;

/* Read one line from the geometry file (defined later) */
static int readGeomLine(FILE *geomFile,int *x,int *y,int *z, int *status,double *k1,double *k2,double *k3);

/* Validates vector components. Returns 1 if the components describe a unit vector, 0 otherwise. */
static int isUnitVector(double v1, double v2, double v3, double tolerance){
  double error = fabs(sqrt((v1*v1) + (v2*v2) + (v3*v3)) -1 );
  if (error > tolerance || isnan(error)) return 0;
  return 1;
}

/* Normalises vector components to produce a unit vector. */
static int normaliseVector(double *v1, double *v2, double *v3){
  double length = sqrt( (*v1)*(*v1) + (*v2)*(*v2) + (*v3)*(*v3) );
  if(length <= UNIT_VECTOR_TOLERANCE) return 0;
  *v1 = *v1 / length;
  *v2 = *v2 / length;
  *v3 = *v3 / length;
  return 1;
}

/* Loads geometry data from file into Geom.                                  */
/* This is to be called when this processes is already allocated a subdomain */
int populateMesh(FILE *geomFile, const char *geomFileName, int normaliseVectors) {
  long line = 0;	/* input line counter */
  long numPoints = 0;	/* total number of points */
  int status;		/* read-in point status */
  int x_in,y_in,z_in; 	/* read-in point coords */
  int x,y,z;		/* corrected coords, based on offsets */
  int x0,x1,y0,y1,z0,z1;/* this subdomain's bounds */
  double k1, k2, k3;	/* read-in direction cosines */
  int rank;		/* counter of procs, for doing stats */
  int invalidVectors = 0; /* flag signalling there were invalid vectors in the file */
#if MPI
  long fillNums[num_active_procs];		/* used for stats of subdomains */
  long fillNum, maxFill, sumFill, numEmpty;	/*   filling numbers            */
#endif
  
  /* Set bounds for seq/MPI versions */
#if MPI
  /*  Include halos if any, since we'll need to know if they're tissue too. */
  x0 = local_xmin - ONE;
  x1 = local_xmax + ONE;
  y0 = local_ymin - TWO;
  y1 = local_ymax + TWO;
  z0 = local_zmin - TRI;
  z1 = local_zmax + TRI;
#else
  x0 = 0;
  x1 = xmax;
  y0 = 0;
  y1 = ymax;
  z0 = 0;
  z1 = zmax;
#endif
  
  /* Zero all points, since files are sparse.                             */
  /* This is an unnecessary precaution as Geom=geom was created by CALLOC */
  for (x=x0; x<x1; x++) {
    for (y=y0; y<y1; y++) {
      for (z=z0; z<z1; z++) {
	Geom[ geom_ind(x, y, z, GEOM_STATUS  ) ] = GEOM_VOID;
	Geom[ geom_ind(x, y, z, GEOM_FIBRE_1 ) ] = 0.0;
	Geom[ geom_ind(x, y, z, GEOM_FIBRE_2 ) ] = 0.0;
	Geom[ geom_ind(x, y, z, GEOM_FIBRE_3 ) ] = 0.0;
      }
    }
  }
  
  if (Verbose){
    MESSAGE("/* Loading geometry data"); 
    if (normaliseVectors) {MESSAGE(" and normalising vectors");}
    MESSAGE("...");
  }
  rewind(geomFile);
  while (EOF!=readGeomLine(geomFile,&x_in,&y_in,&z_in,&status,&k1,&k2,&k3)) {
    
    /* Validate status.                                                    */
    /* This is subject to review. Currently only 0 and 1 are used, but ... */
    if (!(status==GEOM_VOID || status==GEOM_TISSUE || status==GEOM_TISSUE2))
      EXPECTED_ERROR(
        "Problem reading geometry file. Invalid status value %d given on line %d of %s\n.", 
	status, line, geomFileName
      );
    
    /* Shift coords */
    x = x_in + x_offset;
    y = y_in + y_offset;
    z = z_in + z_offset;
    
    /* Ensure all adjusted coords are inside containing box.                */
    /* These assertions must be true if getMeshDimensions worked correctly. */
    assert(x>=0 && x<xmax);
    assert(y>=0 && y<ymax);
    assert(z>=0 && z<zmax);
    
    /* Can a sscanf really produce a NaN ?? */
    assert(!(isnan(k1)));
    assert(!(isnan(k2)));
    assert(!(isnan(k3)));
    
    /* Don't need to do anything if the point is void */
    if (status != GEOM_VOID) {
      
      /* Ensure no tissue points in global halos.                             */
      /* These assertions must be true if getMeshDimensions worked correctly. */
      if (dim>=1) assert(x!=0 && x!=xmax-1);
      if (dim>=2) assert(y!=0 && y!=ymax-1);
      if (dim==3) assert(z!=0 && z!=zmax-1);

      /* Check and optionally normalise fibre vectors if anisotropy is on. */
      if (normaliseVectors && ANISOTROPY_ON) {
	if (!normaliseVector(&k1,&k2,&k3)) {
	  EXPECTED_ERROR(
            "Error normalising vector at (%d,%d,%d): [%lf,%lf,%lf]. Vector is too short.",
	    x_in,y_in,z_in,k1,k2,k3
	  );
	}
      }
      
      /* Check vectors, even if we're normalising. */
      /* 2D meshes may cause errors here, if the only non-zero vector component is along the z axis. */
      if (ANISOTROPY_ON) {
	if (!isUnitVector(k1,k2,k3,UNIT_VECTOR_TOLERANCE)) {
	  MESSAGE("\nFibre vector at point (%d,%d,%d): [%lf,%lf,%lf] is not a valid unit vector.\n",
		x_in,y_in,z_in,k1,k2,k3);
	  invalidVectors = 1;
	}
      }	
      
      /* If it belongs to the local subdomain assign the point */
      if (x>=x0 && x<x1 &&
	  y>=y0 && y<y1 &&
	  z>=z0 && z<z1) {
	/* NB here integer status converts to real */
	Geom[ geom_ind(x, y, z, GEOM_STATUS ) ] = (real) status;
	Geom[ geom_ind(x, y, z, GEOM_FIBRE_1) ] = (real) k1;
	Geom[ geom_ind(x, y, z, GEOM_FIBRE_2) ] = (real) k2;
	Geom[ geom_ind(x, y, z, GEOM_FIBRE_3) ] = (real) k3;
	numPoints++;
      }
    } /* status != GEOM_VOID */
    
    line++;
  } /* while (readGeomLine) */

  if (invalidVectors)
    EXPECTED_ERROR("Invalid fibre vectors encountered. "
		   "Try using 'normaliseVectors=1' in state to correct them.\n");

  if (Verbose) MESSAGE("done. */\n");

#if MPI
  if (Verbose) {
    FFLUSH(stdout);
    MPI_Barrier(ALL_ACTIVE_PROCS);
    URGENT_MESSAGE("/* Process %d: %d tissue points */\n", mpi_rank, numPoints);
    MPI_Barrier(ALL_ACTIVE_PROCS);
  }
  /* We only need the filling numbers in root to do statistics, */
  /* so no need for MPI_Allgather.                              */
  /* NB we set root at rank 0; this is debatable. 		*/
  mpi_errno = MPI_Gather(&numPoints, 1, MPI_LONG, &fillNums, 1, MPI_LONG, 0, ALL_ACTIVE_PROCS);
  CHECK_MPI_SUCCESS("Couldn't gather filling numbers.");
  if (mpi_rank==0) {
    maxFill=0;
    sumFill=0;
    numEmpty=0;
    for (rank=0;rank<num_active_procs;rank++) {
      /* printf("i=%d\t",i); */
      fillNum=fillNums[rank];
      if (maxFill<fillNum) maxFill=fillNum;
      sumFill += fillNum;
      if (fillNum==0) numEmpty++;
      /* printf("/\* fillNum=%d maxFill=%d sumFill=%d numEmpty=%d *\/\n", */
      /* 	     fillNum, maxFill, sumFill, numEmpty); */
      /* fflush(stdout); */
    }
    if (Verbose) {
      MESSAGE("/*----------------------------------------------------*/\n");
      MESSAGE("/* Filling quality %.1f%%, %d=%.1f%% processes with empty subdomains */\n",
	      (sumFill*100.0)/(float)(num_active_procs*maxFill),
	      numEmpty, (numEmpty*100.0)/(float)num_active_procs
	      );
    }
  }
#else 
  if (Verbose)
    MESSAGE("/* %d tissue points out of total volume of %d point, %.1f%% full */\n",
	  numPoints,xmax*ymax*zmax,numPoints*100.0/(xmax*ymax*zmax));
#endif
  return 1;
}

/* Gets primary information about the given geometry:               */
/* bulk dimensions and the non-void points.                         */
/* Explicit return: bulk dimensions (k-variables).                  */
/* Side effects: gpoints array, to be used for subdomain allocation */
int getMeshDimensions (FILE *geomFile,INT *mesh_xmax,INT *mesh_ymax,INT *mesh_zmax) {
  int first,
    x,y,z,status,
    x_in, y_in, z_in,
    min_x,max_x,min_y,max_y,min_z,max_z,
    geom_xlen, geom_ylen, geom_zlen;
  double k1, k2, k3;
  
  first = 1;
  rewind(geomFile);
  while (EOF!=readGeomLine(geomFile,&x,&y,&z,&status,&k1,&k2,&k3)) {
    if (first) {
      min_x = max_x = x;
      min_y = max_y = y;
      min_z = max_z = z;
      first = 0;
    } else if (status != GEOM_VOID) {
      if (x<min_x) min_x=x;
      if (x>max_x) max_x=x;
      if (y<min_y) min_y=y;
      if (y>max_y) max_y=y;
      if (z<min_z) min_z=z;
      if (z>max_z) max_z=z;
    }
  }
  if (Verbose) MESSAGE("/* Geometry is contained in a box [%d:%d]x[%d:%d]x[%d:%d] */\n",min_x,max_x,min_y,max_y,min_z,max_z);
  
  geom_xlen = ((max_x - min_x) + 1);
  geom_ylen = ((max_y - min_y) + 1);
  geom_zlen = ((max_z - min_z) + 1);
  
  /* How should we shift the incoming coords to 
   * put the first tissue point at 1,1,1?	*/
  x_offset = ((geom_xlen>1?1:0)-min_x);
  y_offset = ((geom_ylen>1?1:0)-min_y);
  z_offset = ((geom_zlen>1?1:0)-min_z);
  if (Verbose) MESSAGE("/* Geometry coordinates will be offset by (%d,%d,%d) */\n",x_offset,y_offset,z_offset);

  *mesh_xmax = geom_xlen + (geom_xlen>1?2:0); /* Num points + bounds. */
  *mesh_ymax = geom_ylen + (geom_ylen>1?2:0); /* Num points + bounds. */
  *mesh_zmax = geom_zlen + (geom_zlen>1?2:0); /* Num points + bounds. */
  
  FFLUSH(stdout);

/* To do: consider using this in sequential mode too? */
#if MPI
  /* Now we need to create the gpoints array for subdomain allocation. */

  /* !!! We rely that xmax,... are synonymous to *mesh_xmax,... */
  /* as formal vs factual arguments of this function.           */
  /* This is potentially dangerous, subject to review.          */
  CALLOC(gpoints,xmax*ymax*zmax,sizeof(unsigned char));
  rewind(geomFile);
  while (EOF!=readGeomLine(geomFile,&x_in,&y_in,&z_in,&status,&k1,&k2,&k3)) {
    x = x_in + x_offset;
    y = y_in + y_offset;
    z = z_in + z_offset;
    gpoints[gind(x,y,z)] = (status!=GEOM_VOID) ? 1 : 0;
  }
#endif

  return 1;
}



/* List of characters that can separate fields in a geometry line, other than simple blank */
/* Perhaps this should be extended to include all characters not used in numbers           */
char GEOMSEPS[]=",:;\t\r\b!$=";

/* Read buffer: should be distinct from qpp buffer as they are processed "in parallel" */
#define MAXSTRLEN 8192
static char buf[MAXSTRLEN];
int readGeomLine(FILE *geomFile,int *x,int *y,int *z, int *status,double *k1,double *k2,double *k3) {
  char *p, *c;
  if (NULL==(fgets(buf,MAXSTRLEN,geomFile))) return EOF;
  /* Convert all separator characters to blanks */
  for (c=&(GEOMSEPS[0]);*c;c++) {
    while (NULL!=(p=strchr(buf,*c))) *p=' ';
  }
  return sscanf(buf,"%d %d %d %d %lf %lf %lf\n",x,y,z,status,k1,k2,k3);
}
