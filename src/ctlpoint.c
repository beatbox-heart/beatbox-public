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

/* Control point device */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <stdarg.h>
#include "system.h"
#include "beatbox.h"
#include "k_.h"
#include "bikt.h"
#include "device.h"
#include "qpp.h"
#include "state.h"
#include "mpi_io_choice.h"

extern int  Verbose;
extern char VERSTRING[];

typedef struct {
  char file[MAXPATH];
  int enrich;
  long last;
  char debugname[MAXPATH];
  FILE *debug;
#if MPI
  MPI_Datatype k_var_Type;
  MPI_Datatype filetype;
  MPI_Datatype sourceType;
  MPI_Aint     filetype_extent;
#endif
} STR;

typedef struct {
  char          nm[maxname];
  unsigned char tp;
  unsigned char val[256];
} k_var;

static k_var empty_var={"", 0, ""};

static int fatal(void) {
  return 1;
}

#if MPI
/**********************/
/* PARALLEL VARIATION */
RUN_HEAD(ctlpoint)
  DEVICE_ARRAY(char,file)
  DEVICE_CONST(int,enrich)
  DEVICE_VAR(long,last)
  DEVICE_CONST(FILE *,debug)
  #define DB if(debug) {fprintf(debug, 
  #define BD ); FFLUSH(debug);}

  DEVICE_CONST(MPI_Datatype, k_var_Type)
  DEVICE_CONST(MPI_Datatype, filetype)
  DEVICE_CONST(MPI_Datatype, sourceType)
  DEVICE_CONST(MPI_Aint, filetype_extent)
  MPI_File f;
  MPI_Offset offset;
  int x,y,z,v;
  int mode;
  int file_exists;
  FILE *temp_file;
  
  char V_[32], *p;
  INT xm_,ym_,zm_,vm_;
  int i;
  k_var var;
  
  DB "\nctlpoint called at %ld\n",t BD
    
  if (*last==LNONE) {                      /* never saved, try to restore */
    DB "\nreading\n" BD
	
    #define D(name,new,old,type,size,erraction)				\
      if(mpi_errno != MPI_SUCCESS){					\
	MPI_Error_string(mpi_errno, error_string, &error_length);	\
	MESSAGE("Process %d: error reading from %s at t=%d --- %s", mpi_rank, file, t, error_string); \
	FFLUSH(stdout);							\
	/*erraction;*/							\
      } else if (new!=old && memcmp(new,old,sizeof(type)*size)) {	\
	MESSAGE("P%02d: wrong %s read from %s: %d(x%0x) instead of %d(x%0x) \n", \
		mpi_rank,name,file,*(long *)new,*(long *)new,*(long *)old,*(long *)old); \
	/*erraction;*/							\
      } else if (debug)							\
	fprintf(debug,"\nread %s[0]=%f\n",name,(float)(*(type *)(new)));
      
    #define CHAR_VAR(name, new, old, size, erraction)			\
      mpi_errno = MPI_File_read(f, new, size, MPI_CHAR, MPI_STATUS_IGNORE); \
      D(name, new, old, char, size, erraction)
      
    #define INT_VAR(name, new, old, size, erraction)			\
      mpi_errno = MPI_File_read(f, new, size, MPI_LONG, MPI_STATUS_IGNORE); \
      D(name, new, old, INT, size, erraction)
	
    #define K_VAR(name, new, old, size, erraction)                          \
      mpi_errno = MPI_File_read(f, new, size, k_var_Type, MPI_STATUS_IGNORE); \
      D(name, new, old, k_var, size, erraction)
	
    /*      Use sequential method to check file existance.
     *      Prevents MPI errors all over the shop.  */
    if NOT(temp_file=fopen(file,"rb")) {
	file_exists = 0;
	MESSAGE("Cannot read control point %s at t=%d.\n", file, t);
      } else {
      file_exists = 1;
      fclose(temp_file);
    }
      
    if (file_exists==1) {
      /*  READ FROM DUMP FILE */
      mpi_errno = MPI_File_open(ALL_ACTIVE_PROCS, file, MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
      mpi_errno = MPI_File_set_view(f, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
      CHECK_MPI_SUCCESS("Couldn't set view.")
	
      /*  Read first block of global data. */
      #include "ctlpoint.h"
	  
      #undef CHAR_VAR
      #undef INT_VAR
      #undef K_VAR
                         
      MPI_File_get_position(f, &offset);
      mpi_errno = MPI_File_set_view(f, offset, sourceType, filetype, "native", MPI_INFO_NULL);
      mpi_errno = MPI_File_read(f, New, 1, sourceType, MPI_STATUS_IGNORE);
      CHECK_MPI_SUCCESS("Couldn't read from file.")
      haloSwap();

      offset += filetype_extent;
      mpi_errno = MPI_File_set_view(f, offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
      
      for (;;) {
	mpi_errno = MPI_File_read(f, &var, 1, k_var_Type, MPI_STATUS_IGNORE);
	if(mpi_errno != MPI_SUCCESS){
	  MPI_Error_string(mpi_errno, error_string, &error_length);
	  printf("Process %d: error reading k_var from %s at t=%ld--- %s\n", mpi_rank, file, t, error_string);
	  FFLUSH(stdout);
	  break;
	}
	
	if NOT(var.nm[0]) break;                /* eof mark */
	if NOT(i=tb_find(deftb,var.nm)) {
	    if (!enrich) continue;
	    if NOT(i=tb_insert_abstract(
					deftb,var.nm,var.tp,(p_vd)Calloc(1,sizetable[var.tp]),0,f_vb|f_rs
                                        )) EXPECTED_ERROR("cannot define variable %s",var.nm);
	  } /*  if not defined. */
	if (tb_type(deftb,i)!=var.tp) EXPECTED_ERROR("variable %s of wrong type",var.nm);
	if ((tb_flag(deftb,i)&f_vb)==0) EXPECTED_ERROR("%s is not a variable",var.nm);
	memcpy(tb_addr(deftb,i),&(var.val),sizetable[var.tp]);
      } /*  for */
      mpi_errno = MPI_File_close(&f);
      CHECK_MPI_SUCCESS("Couldn't close file after reading.")
      f = MPI_FILE_NULL;                      
      DB "\nread ok\n" BD
    } /*  if(file_exists==1) */
    #undef D
  } /*  if(*last==LNONE) */
  else {                          /* saved, try to save */
    DB "\nwriting\n" BD
    vm_=vmax;
      
    #define D(name,new,old,type,size,erraction)				\
      if(mpi_errno != MPI_SUCCESS){					\
	MPI_Error_string(mpi_errno, error_string, &error_length);	\
	MESSAGE("Process %d: error writing %s to %s at t=%ld --- %s\n",mpi_rank,name,file,t, error_string); \
	FFLUSH(stdout);							\
	ABORT("\n");							\
      }
    
    #define CHAR_VAR(name, new, old, size, erraction)			\
      mpi_errno = MPI_File_write_all(f, old, size, MPI_CHAR, MPI_STATUS_IGNORE); \
      D(name,new,old,char,size,erraction)        
    
    #define INT_VAR(name, new, old, size, erraction)			\
      mpi_errno = MPI_File_write_all(f, old, size, MPI_LONG, MPI_STATUS_IGNORE); \
      D(name,new,old,INT,size,erraction)
    
    #define K_VAR(name, new, old, size, erraction)				\
      mpi_errno = MPI_File_write_all(f, &var, size, k_var_Type, MPI_STATUS_IGNORE); \
      D(name,new,old,k_var,size,erraction)
    
    strcpy(buf,file);
    for (p=buf+strlen(buf)-1;p>=buf&&*p=='~';p--);
    if (p<buf) {MESSAGE("cannot generate tempname for %s",file);return(0);}
    *p='~';
    
    mode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
    MPIDO(MPI_File_open(ALL_ACTIVE_PROCS, buf, mode, MPI_INFO_NULL, &f),
	  "cannot open dump file for writing");
    
    /*  Set view with no displacement (top of file). */
    MPIDO( MPI_File_set_view(f, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL),
	   "Couldn't set view.");
      
    if (mpi_rank == 0) {
      /*  Write first block of global data. */
      #include "ctlpoint.h"   
    }
    
    MPIDO(MPI_File_get_position(f, &offset),
	  "Could not get file position");
    MPIDO(MPI_Bcast(&offset, 1, MPI_INT, 0, ALL_ACTIVE_PROCS),
	  "Could not broadcast offset");
    MPIDO(MPI_File_set_view(f,offset,sourceType,filetype,"native",MPI_INFO_NULL),
	  "Could not set view");
    MPIDO(MPI_File_write_all(f, New, 1, sourceType, MPI_STATUS_IGNORE),
	  "Could not write to file.");
      
      
    /*  Adjust file pointer to after New. */
    offset += filetype_extent;
    MPIDO(MPI_File_set_view(f,offset,MPI_CHAR,MPI_CHAR,"native",MPI_INFO_NULL),
	  "Could not set view.");
    
    if (mpi_rank == 0) {
      /*  Write more global data... */
      for (i=1;i<=maxtab;i++) {
	if(!*tb_name(deftb,i)) continue;
	if((tb_flag(deftb,i)&f_vb)==0) continue;
	memcpy(var.nm,tb_name(deftb,i),maxname);
	var.tp=tb_type(deftb,i);
	memcpy(&(var.val),tb_addr(deftb,i),sizetable[var.tp]);
	K_VAR(var.nm,&var,&var,1,return 0);
      } /*  for i-->maxtab */
      K_VAR("empty_var",&empty_var,&empty_var,1,return 0); /* eof mark */
    } /*  if(mpi_rank == 0) */
    
    /*
      TODO Supress error message when closing the file.
      */
    MPIDO(MPI_File_close(&f),"Couldn't close ctlpoint file.");
    f = MPI_FILE_NULL;
      
    frename(buf,file);
                
    sprintf(buf,"dumped at %ld to %s\n",t,file);
    /* if(!MESSAGE(buf)) crt_text(buf,w.row0,w.col0,15); */
    #undef D
    
    DB "\nwritten ok\n" BD
      
    #undef CHAR_VAR
    #undef INT_VAR
    #undef K_VAR
  }
        
  *last = t;
RUN_TAIL(ctlpoint)
#else
/*************************/
/*  SEQUENTIAL VARIATION */
RUN_HEAD(ctlpoint)
  DEVICE_ARRAY(char,file)
  DEVICE_CONST(int,enrich)
  DEVICE_VAR(long,last)
  DEVICE_CONST(FILE *,debug)
  #define DB if(debug) {fprintf(debug, 
  #define BD ); FFLUSH(debug);}
 
  FILE *f;

  char V_[32], *p;
  INT xm_,ym_,zm_,vm_;
  int i;
  k_var var;

  DB "\nctlpoint called at %ld\n",t BD

  if (*last==LNONE) {                      /* never saved, try to restore */
    DB "\nreading\n" BD

    #define D(name,new,old,type,size,erraction)				\
      if(size!=fread(new,sizeof(type),size,f)) {			\
	MESSAGE("error reading %s from %s at t=%ld\n",name,file,t);	\
	erraction;							\
      } else if (new!=old && memcmp(new,old,sizeof(type)*size)) {	\
	MESSAGE("wrong %s read from %s: %d(x%0x) instead of %d(x%0x) \n", \
		name,file,*(long *)new,*(long *)new,*(long *)old,*(long *)old);	\
	erraction;							\
      } else if (debug)							\
	fprintf(debug,"\nread %s[0]=%f\n",name,(float)(*(type *)(new)));
      if NOT(f=fopen(file,"rb")) {
	  MESSAGE("\nWARNING: cannot read control point %s\n", file);
	  *last = t;
	  return 1;
	}
    #include "ctlpoint.h"
    for (;;) {
      if(1!=fread(&var,sizeof(k_var),1,f)) {				\
	MESSAGE("error reading k_var from %s at t=%ld\n",file,t);	\
	break;								\
      } 
      if NOT(var.nm[0]) break;                /* eof mark */
      DB "\n\t%s",var.nm BD
	if NOT(i=tb_find(deftb,var.nm)) {
	    if (!enrich) continue;
	    if NOT(i=tb_insert_abstract(
					deftb,var.nm,var.tp,(p_vd)Calloc(1,sizetable[var.tp]),0,f_vb|f_rs
					)) EXPECTED_ERROR("cannot define variable %s",var.nm);
	  }
      
      if (tb_type(deftb,i)!=var.tp) EXPECTED_ERROR("variable %s of wrong type",var.nm);
      DB "\t(%s) ", typenametable[var.tp] BD
	if ((tb_flag(deftb,i)&f_vb)==0) EXPECTED_ERROR("%s is not a variable",var.nm);
      memcpy(tb_addr(deftb,i),&(var.val),sizetable[var.tp]);
      DB "\t= %s",prt(tb_addr(deftb,i),var.tp) BD
	}
    
    MESSAGE(buf,"loaded at %ld from %s\n",t,file);
#undef D
    DB "\nread ok\n" BD
  } else {                                /* saved, try to save */
    DB "\nwriting\n" BD
    vm_=vmax;
    #define D(name,new,old,type,size,erraction)                             \
      if(size!=fwrite(old,sizeof(type),size,f)) {			\
	MESSAGE("error writing %s to %s at t=%ld\n",name,file,t);	\
	return 0;							\
      }
    
    strcpy(buf,file);
    for(p=buf+strlen(buf)-1;p>=buf&&*p=='~';p--);
    if(p<buf) {MESSAGE("cannot generate tempname for %s",file);return(0);}
    *p='~';
    
    if NOT(f=fopen(buf,"wb")) {
	MESSAGE("cannot open dump file %s",buf);
	perror(0);
	return(0);
    }

    #include "ctlpoint.h"
    
    for(i=1;i<=maxtab;i++) {
      if(!*tb_name(deftb,i)) continue;
      if((tb_flag(deftb,i)&f_vb)==0) continue;
      memcpy(var.nm,tb_name(deftb,i),maxname);
      var.tp=tb_type(deftb,i);
      memcpy(&(var.val),tb_addr(deftb,i),sizetable[var.tp]);
      
      D(var.nm,&var,&var,k_var,1,return 0);
    } /*  for i-->maxtab */
    D("empty_var",&empty_var,&empty_var,k_var,1,return 0); /* eof mark */
    
    if(0!=fclose(f)) {
      MESSAGE("error closing %s at t=%ld\n",buf,t);
      return(0);
    }
    
    f=NULL;
    
    frename(buf,file);
    
    MESSAGE(buf,"dumped at %ld to %s\n",t,file);
    #undef D
    DB "\nwritten ok\n" BD
  }

  *last = t;
RUN_TAIL(ctlpoint)
#endif

/********************/
DESTROY_HEAD(ctlpoint)
        if (S->debug) fclose(S->debug); S->debug=NULL;
DESTROY_TAIL(ctlpoint)

/*******************/
CREATE_HEAD(ctlpoint)
  DEVICE_IS_SPACELESS
        
  ACCEPTS(file,NULL);
  ACCEPTI(enrich,0,0,1);
  ACCEPTF(debug,"wt","");
  S->last=LNONE;
  
#if MPI
  
  MPI_Datatype k_var_Type;
  k_var var;
  
  int i;
  int lengths[3] = {maxname, 1, 256};
  MPI_Aint *displacements = calloc(3, sizeof(MPI_Aint));
  MPI_Datatype old_datatypes[3] = {MPI_CHAR, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR};

  MPI_Get_address(var.nm,         &displacements[0]);
  MPI_Get_address(&var.tp,        &displacements[1]);
  MPI_Get_address(var.val,        &displacements[2]);

  for (i = 2; i>=0; i--) {
    displacements[i] -= displacements[0];
  }

  MPI_Type_create_struct(3, lengths, displacements, old_datatypes, &k_var_Type);
  MPI_Type_commit(&k_var_Type);

  /*  Global space must be whole medium (inc absolute boundaries.) */
  dev->s.runHere = 1; /*  Must run on all processes. */
  dev->s.global_x0 = 0;
  dev->s.global_x1 = xmax-1;
  dev->s.global_y0 = 0;
  dev->s.global_y1 = ymax-1;
  dev->s.global_z0 = 0;
  dev->s.global_z1 = zmax-1;
  dev->s.v0 = 0;
  dev->s.v1 = vmax-1;

  if (mpi_ix == 0) {
    dev->s.x0 = dev->s.global_x0;
  } else {
    dev->s.x0 = local_xmin;
  }
  if (mpi_ix == mpi_nx-1) {
    dev->s.x1 = dev->s.global_x1;
  } else {
    dev->s.x1 = local_xmax-1;
  }
  if (mpi_iy == 0) {
    dev->s.y0 = dev->s.global_y0;
  } else {
    dev->s.y0 = local_ymin;
  }
  if (mpi_iy == mpi_nx-1) {
    dev->s.y1 = dev->s.global_y1;
  } else {
    dev->s.y1 = local_ymax-1;
  }
  if (mpi_iz == 0) {
    dev->s.z0 = dev->s.global_z0;
  } else {
    dev->s.z0 = local_zmin;
  }
  if (mpi_iz == mpi_nz-1) {
    dev->s.z1 = dev->s.global_z1;
  } else {
    dev->s.z1 = local_zmax-1;
  }

  /*  Define types for New. */
  #include "dump_types.h"
  MPI_Aint filetype_extent;
  MPI_Type_extent(filetype, &filetype_extent);
  S->k_var_Type = k_var_Type;
  S->filetype = filetype;
  S->sourceType = sourceType;
  S->filetype_extent = filetype_extent;
#endif
CREATE_TAIL(ctlpoint,0)

