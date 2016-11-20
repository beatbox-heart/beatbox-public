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
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "system.h"
#include "beatbox.h"
#include "k_.h"
#include "device.h"
#include "qpp.h"
#include "bikt.h"
#include "decomp.h"
#include "error.h"
#include "geometry.h"
#define OWN
#include "state.h"
#undef OWN
#include "mpi_io_choice.h"

/* macro for accepting global values, unlike local one in qpp.h */
#define ACCEPTLG(b,c,d,e) if (!acceptl(#b"=",&(S->b),c,d,e,w)) return(0); b=S->b

extern int Verbose;

static real  *state_u1;
static real  *geom;
static int already_created=0;
typedef struct {
  /* Geometry file */
  char geometryname[MAXPATH];
  FILE *geometry;
  /* Geometry options */
  int normaliseVectors;
  int correctBoundaries;
  int anisotropy;
  
  /*  Global bounds: "allow sign" to suppress compiler warnings */
  ssize_t xmax, ymax, zmax, vmax;
} STR;

/* state is performed once just when reading from parameter file */
int state (char *w) {
  int i; /*  Loop iterator */
  int base; /* Index from which shifts are measured. */
  size_t num_values; /* Number of values in the state (or in this process's subgrid) */

  STR *S= (STR *) Calloc(1,sizeof(STR));
  
  if (!S) ABORT("not enough memory");
	
  if (already_created) EXPECTED_ERROR("state must be defined only once");
  
  ACCEPTF(geometry,"r","");
  GEOMETRY_ON = (S->geometry!=NULL);
  ACCEPTI(anisotropy,0,0,1);
  ANISOTROPY_ON = (anisotropy && GEOMETRY_ON); /* Anisotropy is a geometry-only feature */
	
  if (GEOMETRY_ON) {

    ACCEPTI(normaliseVectors,0,0,1);
    
    if (stateDimensionsExist(w)) {
      MESSAGE("\n/**************************************************** WARNING ****************************************************\n");
      MESSAGE("Dimensions of the simulation medium (xmax,ymax,zmax) cannot be specified when geometry is in use.\n");
      MESSAGE("The dimensions provided will be ignored.\n");
      MESSAGE("*****************************************************************************************************************/\n\n");
    }

    if (Verbose && ANISOTROPY_ON) MESSAGE("\n/* Anisotropy is on. */");

    if (Verbose) MESSAGE("\n/* Finding minimal enclosing box for geometry from %s ... */\n",
			 S->geometryname);
    getMeshDimensions(S->geometry,&xmax,&ymax,&zmax);
    S->xmax = xmax;
    S->ymax = ymax;
    S->zmax = zmax;
    geom_vmax = 4; /* STATUS + 3 VECTOR COMPONENTS */
    
  } else {


    ACCEPTLG(xmax,LNONE,1L,LNONE);
    ACCEPTLG(ymax,LNONE,1L,LNONE);
    ACCEPTLG(zmax,   1L,1L,LNONE);
    geom_vmax = 0; /* extra insurance against unintended use of Geom */
    if (anisotropy) {
      MESSAGE("\n/**************************************************** WARNING ****************************************************\n");
      MESSAGE("Anisotropy can only be activated when geometry is in use.\n");
      MESSAGE("We'll continue without for now.\n");
      MESSAGE("*****************************************************************************************************************/\n\n");
    }

  }

  ACCEPTLG(vmax,     2L,2L,LNONE);
  if (Verbose) MESSAGE("\n/* Grid size (%d x %d x %d x %d) */",xmax,ymax,zmax,vmax);
	
  /* If any dimension is exactly 2, complain.
   * Since it's not a single point, we'll need diffusion,
   * but we don't have enough points for boundaries.
   */
  if (xmax==2||ymax==2||zmax==2) {
    EXPECTED_ERROR("One or more dimensions (xmax,ymax,zmax) is equal to 2. Dimensions must be either 1 or >=3.");
  }
  
  /* If the dimension is in use, flag it so we can 
     treat the first and last points as boundaries.*/	
  TRI = zmax >=3 ? 1 : 0;
  TWO = ymax >=3 ? 1 : 0;
  ONE = xmax >=3 ? 1 : 0;
  dim = ONE + TWO + TRI;
  
  if(dim==1 && !ONE)		EXPECTED_ERROR("1-Dimensional simulations must be defined on the x axis.");
  if(dim==1 && GEOMETRY_ON)	EXPECTED_ERROR("Geometry requires a 2- or 3-dimensional simulation medium.");
  if(dim==2 && !(ONE && TWO))	EXPECTED_ERROR("2-Dimensional simulations must be defined on the x and y axes.");

/**************************************************************************
 ************************** MPI ONLY **************************************
 **************************************************************************/
#if MPI
  /*  Get total number of processes. */
  MPIDO(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size),"Could not get number of MPI processes.");
	
  /* Compute domain decomposition.                 */
  /* Also, define bounds of subdomains along axes. */
  num_subdoms = decompose(mpi_size, &mpi_nx, &mpi_ny, &mpi_nz);
  if (num_subdoms<=0) EXPECTED_ERROR("Could not decompose the domain\n");
  if (Verbose) MESSAGE("\n/* Domain decomposed as (%03d,%03d,%03d). */",mpi_nx,mpi_ny,mpi_nz);

  /* Allocate subdomains to processes */
  num_active_procs = decomp_allocateSubdomains();
  if (num_active_procs<=0) 
    EXPECTED_ERROR("could not allocate subdomains per processes: num_active_procs=%d\n",num_active_procs);
  if (Verbose) MESSAGE("\n/* %d processes will be active, %d will be idle. */",num_active_procs,mpi_size-num_active_procs);
  num_empty_subdoms = num_subdoms-num_active_procs;
  if (Verbose) MESSAGE("\n/* %d subdomains are empty. */",num_empty_subdoms);


  /*  Define communicator over active processes. If we're idle, we're done. */
  I_AM_IDLE = decomp_globalCommunicator (mpi_rank, &ALL_ACTIVE_PROCS, &mpi_rank);

  if (I_AM_IDLE) {
    if (Verbose) URGENT_MESSAGE("\n/* Process %d will be idle */",mpi_rank);
    FREE(S);
    return (1);
  }
	
  decomp_getSuperindices(mpi_rank, &mpi_ix, &mpi_iy, &mpi_iz);
  decomp_getSubdomainDimensions();
  decomp_getSubdomainSize();
  decomp_defineHaloTypes();

  /* Locate processes of neighbouring subdomains. */
  /* These will return NORANK for empty neighbours. */
  xn_neighbour=getProcessRank(mpi_ix-1,mpi_iy  ,mpi_iz  ); /* neighbour to LEFT	*/
  xp_neighbour=getProcessRank(mpi_ix+1,mpi_iy  ,mpi_iz  ); /* neighbour to RIGHT	*/
  yn_neighbour=getProcessRank(mpi_ix  ,mpi_iy-1,mpi_iz  ); /* neighbour ABOVE	*/
  yp_neighbour=getProcessRank(mpi_ix  ,mpi_iy+1,mpi_iz  ); /* neighbour BELOW	*/
  zn_neighbour=getProcessRank(mpi_ix  ,mpi_iy  ,mpi_iz-1); /* neighbour in FRONT	*/
  zp_neighbour=getProcessRank(mpi_ix  ,mpi_iy  ,mpi_iz+1); /* neighbour BEHIND	*/

  /* Report subdomain details of this processes. */
  if (Verbose) {
    FFLUSH(stdout);
    MPI_Barrier(ALL_ACTIVE_PROCS);
    URGENT_MESSAGE(
      "\n/* Process %d=(%d,%d,%d) domain size (%dx%dx%d) limits [%d:%d]x[%d:%d]x[%d:%d] neighbours %2d:%2d, %2d:%2d, %2d:%2d */",
      mpi_rank,mpi_ix,mpi_iy,mpi_iz,xlen,ylen,zlen,
      local_xmin, local_xmax-1, local_ymin, local_ymax-1, local_zmin, local_zmax-1,
      xn_neighbour,xp_neighbour,yn_neighbour,yp_neighbour,zn_neighbour,zp_neighbour
    );
    MPI_Barrier(ALL_ACTIVE_PROCS);
  }

#else
/**************************************************************************
 ************************** SEQUENTIAL ONLY *******************************
 **************************************************************************/
  /*  Size of the medium, including boundaries. */
  xlen = xmax;
  ylen = ymax;
  zlen = zmax;	
#endif
/**************************************************************************
 ************************** COMMON ****************************************
 **************************************************************************/

  /* This is the actual declaration of memory (used by this subdomain) */
  num_values =  xlen*ylen*zlen*vmax;
  CALLOC(state_u1,num_values,sizeof(real));
	
  /*  This is unnecessary if CALLOC works correctly ? */
  for (i=0;i<num_values;i++) state_u1[i] = 0.0;
	
  /*  Get constants for index computation. */
  vmax_zmax = vmax * zlen;
  vmax_zmax_ymax = vmax * zlen * ylen;
  vmax_ymax = vmax * ylen;
	
  /*  Get constants for geometry index computation. */
  if (GEOMETRY_ON) {
    geom_vmax_zmax= geom_vmax * zlen;
    geom_vmax_zmax_ymax= geom_vmax * zlen * ylen;
    geom_vmax_ymax = geom_vmax * ylen;
  }

/*  Shift increments */
#if MPI
  /*  Computed from local points to avoid upsetting ind() */
  base = ind(local_xmin,   local_ymin,   local_zmin,   0);
  DX   = ind(local_xmin+1, local_ymin,   local_zmin,   0) - base;
  DY   = ind(local_xmin,   local_ymin+1, local_zmin,   0) - base;
  DZ   = ind(local_xmin,   local_ymin,   local_zmin+1, 0) - base;
  DV   = ind(local_xmin,   local_ymin,   local_zmin,   1) - base;
#else /*  Sequential */
  base = ind(0,0,0,0);
  DX   = ind(1,0,0,0) - base;
  DY   = ind(0,1,0,0) - base;
  DZ   = ind(0,0,1,0) - base;
  DV   = ind(0,0,0,1) - base;
#endif

  New = state_u1;
	
  if (GEOMETRY_ON) {
    int num_geom_values = xlen*ylen*zlen*geom_vmax;
    CALLOC(geom,num_geom_values,sizeof(real));
    Geom = geom;
    if (!populateMesh(S->geometry,S->geometryname,S->normaliseVectors)) return 0;
    fclose(S->geometry);
  }
	
  FREE(S);
  return (1);
	
} /*  state() */

void state_free(void) {
  if(GEOMETRY_ON) FREE(Geom);
  FREE(New);
}
/*------------*/


/**************************************************************************
 ************************** MPI FUNCTIONS *********************************
 **************************************************************************/
#if MPI

/* This is a slow and careful version intended for use */
/* at start time rather than run time                  */
int getProcessRank (int ix,int iy,int iz) {
  /* These checks here are unnecessary if done in the calling function */
  if (ix<0) return NORANK; if (ix>=mpi_nx) return NORANK;
  if (iy<0) return NORANK; if (iy>=mpi_ny) return NORANK;
  if (iz<0) return NORANK; if (iz>=mpi_nz) return NORANK;
  /* The line below is enough if checks done in the calling function. */
  /* If the subdomain is empty it still will return NORANK. */
  return ProcessRank[sind(ix,iy,iz)];
}
/*------------*/

/* This version exploits axes partitions rather than going through all subdomains. */
/* To do: optimize this using details of the partitioning algorithm */
int getRankContainingPoint(int x,int y,int z) {
  int i, ix, iy, iz;
  ix=iy=iz=-1;
  for (i=0;i<mpi_nx;i++) {if (all_xmin[i]<=x && x<all_xmax[i]) {ix=i; break;};} if (ix<0) return NORANK;
  for (i=0;i<mpi_ny;i++) {if (all_ymin[i]<=y && y<all_ymax[i]) {iy=i; break;};} if (iy<0) return NORANK;
  for (i=0;i<mpi_nz;i++) {if (all_zmin[i]<=z && z<all_zmax[i]) {iz=i; break;};} if (iz<0) return NORANK;
  /* No need for extra checks now. */
  /* If the subdomain is empty it still will return NORANK. */
  return ProcessRank[sind(ix,iy,iz)];
}
/*------------*/

#endif /* MPI */
