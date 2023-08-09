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

#ifndef _STATE_H_
#define _STATE_H_ 

#ifdef OWN
#undef OWN
#define _OWN
#endif
#include "extern.h"
#include "k_.h"
#ifdef MPI
#include <mpi.h>
#endif /*	MPI	*/
#ifdef _OWN
#undef _OWN
#define OWN
#endif
#include "extern.h"
/* as some state vars are to be k-vars */

/*
 * All items to be saved for continuing calculations later are here,
 * except k-vars
 */

/* Pointer to current state Ordered V,Z,Y,X, so that steps between
 * neighbouring cells in the x-axis are largest.  This means that a
 * y/z plane is represented by a contiguous block of memory.  This
 * declaration means that New is a global EXTERN variable. Need to
 * "gatherv this variable on master node.  my_rank = 0 to be treated
 * as master node.
 */

EXTERN real *New;		

/**
 *	Pointer to geometry data.
 *	Separated from New to prevent confusions over vmax.
 */
EXTERN real *Geom;

EXTERN INT t;				/* device loops counter */
EXTERN INT advance;			/* nonzero "before the time" */
EXTERN INT vmax,xmax,ymax,zmax;		/* the sizes of a state */
EXTERN size_t xlen,ylen,zlen;		/* Lengths of axes in memory (including halos) */
EXTERN INT GEOMETRY_ON;			/* is geometry on?		*/
EXTERN INT geom_vmax;			/* number of geom variables */

#if MPI
EXTERN unsigned char *gpoints;		/* "skeleton" of the full geometry   */
EXTERN int *ProcessRank;		/* table subdomain -> process */
EXTERN int *rank_ix;			/* table process -> subdomain x-superindex */
EXTERN int *rank_iy;			/* table process -> subdomain y-superindex */
EXTERN int *rank_iz;			/* table process -> subdomain z-superindex */
#define NORANK (-1)			/* a no-process to which an empty subdomain is allocated */
#endif

EXTERN INT ANISOTROPY_ON;		/* is anisotropy on?		*/

EXTERN size_t vmax_zmax, vmax_zmax_ymax;		/* for index computation */
EXTERN size_t geom_vmax_zmax, geom_vmax_zmax_ymax;	/* for geometry index computation */

EXTERN size_t DX, DY, DZ, DV;		/* index shifts */
EXTERN int dim,ONE,TWO,TRI;		/* used to control dimension */

/*
 *	Device Space defaults.
 */
#define SPACE_DEFAULT_X0 ONE
#define SPACE_DEFAULT_X1 ((int)xmax-(1+ONE))
#define SPACE_DEFAULT_Y0 TWO
#define SPACE_DEFAULT_Y1 ((int)ymax-(1+TWO))
#define SPACE_DEFAULT_Z0 TRI
#define SPACE_DEFAULT_Z1 ((int)zmax-(1+TRI))

/**************************************************************************
 ************************** GEOMETRY **************************************
 **************************************************************************/

/*	Layers for geometry data.	*/
#define GEOM_STATUS  0
#define GEOM_FIBRE_1 1
#define GEOM_FIBRE_2 2
#define GEOM_FIBRE_3 3

/*	Point statuses.	*/
#define GEOM_VOID   0
#define GEOM_TISSUE 1
#define GEOM_TISSUE2 2

/**************************************************************************
 ************************** MPI ONLY **************************************
 **************************************************************************/
#ifdef MPI
EXTERN int mpi_size;			/* Total number of processes in MPI_COMM_WORLD */
EXTERN int mpi_rank; 			/* Rank of this process among MPI_COMM_WORLD and later among ALL_ACTIVE_PROCS */
EXTERN int num_subdoms;			/* Number of subdomains */
EXTERN int num_active_procs;		/* Number of active (non-idle) processes. Could be < num_subdoms. */
EXTERN MPI_Comm ALL_ACTIVE_PROCS;	/* Communicator between active processes */
EXTERN int I_AM_IDLE;			/* This process is idle or empty */
EXTERN long mpi_nx,mpi_ny,mpi_nz;	/* Number of partitions for each axis */
EXTERN long mpi_ix, mpi_iy, mpi_iz;	/* Superindices for this process */
typedef struct {			/* Limits of a subdomain, "netto"; */
  int local_xmin, local_xmax;		/*   allowed upper values are local?_max-1 */
  int local_ymin, local_ymax;
  int local_zmin, local_zmax;
} Subdomain;
EXTERN Subdomain *subdomain_limits;	/* Local minima and maxima for all subdomains */
EXTERN int *all_xmin, *all_xmax;	/* Local minima and maxima */
EXTERN int *all_ymin, *all_ymax;	/*   of all coordinate     */
EXTERN int *all_zmin, *all_zmax;	/*   axes' partitions      */
EXTERN int num_empty_subdoms;		/* Number of subdomains not allocated to processes */

/***********************************/
/* This subdomain's details        */

/* Local bounds for this parallel process */
EXTERN int local_xmin, local_xmax, local_ymin, local_ymax, local_zmin, local_zmax;
/* Datatypes for sending hyperplanes */
EXTERN MPI_Datatype XN_Type, XP_Type, YN_Type, YP_Type, ZN_Type, ZP_Type;
/* Datatypes for receiving hyperplanes */
EXTERN MPI_Datatype XN_Halo_Type, XP_Halo_Type, YN_Halo_Type, YP_Halo_Type, ZN_Halo_Type, ZP_Halo_Type;

/* Macros for calculating neighbouring processes
 * e.g. XN_NEIGHBOUR is neighbour in negative direction along x-axis,
 * 		YP_NEIGHBOUR is neighbour in positive direction along y-axis. */
#define XN_NEIGHBOUR xn_neighbour /* neighbour to LEFT	*/
#define XP_NEIGHBOUR xp_neighbour /* neighbour to RIGHT	*/
#define YN_NEIGHBOUR yn_neighbour /* neighbour ABOVE	*/
#define YP_NEIGHBOUR yp_neighbour /* neighbour BELOW	*/
#define ZN_NEIGHBOUR zn_neighbour /* neighbour in FRONT	*/
#define ZP_NEIGHBOUR zp_neighbour /* neighbour BEHIND	*/

/* Ranks of neighbours are pre-calculated at start time */
EXTERN int xn_neighbour;
EXTERN int xp_neighbour;
EXTERN int yn_neighbour;
EXTERN int yp_neighbour;
EXTERN int zn_neighbour;
EXTERN int zp_neighbour;

/* .. using this function, mapping supercoords -> process rank */
int getProcessRank(int ix,int iy,int iz); 

/* end of this subdomain's details */
/***********************************/

/* Get the rank of the process containing a      */
/*    given internal or absolute boundary point. */
int getRankContainingPoint(int x,int y,int z);

/* Macro to compute the index of a variable in New by its coordinates.
 * The MPI version has to be adjusted to match the continuous coordinate 
 * system to the local subdomain stored in New.
 *
 * Takes account of local bounds and whether or not boundary halos exist 
 * on each axis.
 */
/* To do: optimize this */
#define ind(x,y,z,v) ( \
  ((x + ONE) - local_xmin)*vmax_zmax_ymax + \
  ((y + TWO) - local_ymin)*vmax_zmax + \
  ((z + TRI) - local_zmin)*vmax + \
  (v) \
)

static size_t indfun(size_t x,size_t y,size_t z,int v) {
  return ((x + ONE) - local_xmin)*vmax_zmax_ymax +	
  ((y + TWO) - local_ymin)*vmax_zmax + 
  ((z + TRI) - local_zmin)*vmax + 
  (v);
}

/* TODO: check this is equivalent to, and which is faster */
/* #define ind(x,y,z,v) ((v)+vmax*((z)+TRI-local_zmin+zmax*((y)+TWO-local_ymin+ymax*((x)+ONE-local_xmin)))) */

/* Macro to compute the index in the global geometry array(s) */
/* used by domain decomposition routines, by the global coordinates */
#define gind(x,y,z) (((x)*ymax+(y))*zmax+(z))

/* Macro to compute the index in the subdomains' supergrid, by supercoordinates */
#define sind(ix,iy,iz) (((ix)*mpi_ny+(iy))*mpi_nz+(iz))

/* Macro to compute the index of a variable in Geom by its coordinates.
 * The MPI version has to be adjusted to match the continuous coordinate system to
 * the local subdomain stored in New.
 *
 * Takes account of local bounds and whether or not boundary halos exist on each axis.
 */
/* To do: optimize this by hand */
#define geom_ind(x,y,z,v) ( \
  ((x + ONE) - local_xmin)*geom_vmax_zmax_ymax + \
  ((y + TWO) - local_ymin)*geom_vmax_zmax + \
  ((z + TRI) - local_zmin)*geom_vmax + \
  (v) \
)

#else
/**************************************************************************
 ************************** SEQUENTIAL ONLY *******************************
 **************************************************************************/
#define ind(x,y,z,v) ((x)*vmax_zmax_ymax+(y)*vmax_zmax+(z)*vmax+(v))
static size_t indfun(size_t x,size_t y,size_t z,int v) {
  return ((x)*vmax_zmax_ymax+(y)*vmax_zmax+(z)*vmax+(v));
}

#define geom_ind(x,y,z,v) ((x)*geom_vmax_zmax_ymax+(y)*geom_vmax_zmax+(z)*geom_vmax+(v))
/* and these definition can help avoiding unnecessary #if(MPI) clauses: */
#define mpi_size 1
#define mpi_rank 0
EXTERN long mpi_nx,mpi_ny,mpi_nz;	/* ==1 */
#define mpi_ix 0L
#define mpi_iy 0L
#define mpi_iz 0L
#define local_xmin 0
#define local_xmax xmax
#define local_ymin 0
#define local_ymax ymax
#define local_zmin 0
#define local_zmax zmax
#endif
/**************************************************************************
 ************************** COMMON ****************************************
 **************************************************************************/

/* Tests the geometry status (TISSUE/VOID) of a point.
 * Takes (x,y,z) coordinates and evaluates to true or false.
 */
#if MPI
/* newer version: uses gpoints */
#define isTissue(x,y,z)	( \
  ( \
    x>=SPACE_DEFAULT_X0 && x<=SPACE_DEFAULT_X1 && \
    y>=SPACE_DEFAULT_Y0 && y<=SPACE_DEFAULT_Y1 && \
    z>=SPACE_DEFAULT_Z0 && z<=SPACE_DEFAULT_Z1 \
  ) && ( \
    !GEOMETRY_ON || \
    gpoints[gind(x,y,z)] \
  ) \
)
#define istissue(x,y,z)	(!GEOMETRY_ON || gpoints[gind(x,y,z)])
#else
/* older version: uses Geom */
#define isTissue(x,y,z)	( \
  ( \
    x>=SPACE_DEFAULT_X0 && x<=SPACE_DEFAULT_X1 && \
    y>=SPACE_DEFAULT_Y0 && y<=SPACE_DEFAULT_Y1 && \
    z>=SPACE_DEFAULT_Z0 && z<=SPACE_DEFAULT_Z1 \
  ) && ( \
    !GEOMETRY_ON || \
    Geom[geom_ind(x,y,z,GEOM_STATUS)] == (real)GEOM_TISSUE || \
    Geom[geom_ind(x,y,z,GEOM_STATUS)] == (real)GEOM_TISSUE2 \
  ) \
)
#define istissue(x,y,z)	(!GEOMETRY_ON || Geom[geom_ind(x,y,z,GEOM_STATUS)] != GEOM_VOID)
#endif

/* Macro for calculating index of screen arrays by 2D pixel coordinates */
#define indxy(x,y) ((x)+xmax*(y))

/* Flag showing that some step procedure was already defined */
EXTERN int step_already_created;

/* Function defining the state */
int state(char *w);

/* Function destroying the state */
void state_free(void);

#undef EXTERN

#endif /* _STATE_H_ */
