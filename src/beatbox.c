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


/* MAIN PROGRAM */

/* 
 * The next two files will be autogenerated in the parent (trunk)
 * directory after you compile for the first time.
 */
#include <config.h>
#include <subversion.h>

#if MPI
#include <mpi.h>
#endif

#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAIN
#include "system.h"
#include "beatbox.h"
#include "init.h"
#include "device.h"
#include "state.h"
#include "bikt.h"
#include "windraw.h"
#include "qpp.h"
#include "mpi_io_choice.h"

/* Exported variables */

/* options */
#if MPI
/*  Adjust for lack of user interaction when running with MPI. */
  INT Graph=0;                /* don't use graphical output with MPI */
#else
  INT Graph=1;                /* use graphical output */
#endif
int Profile=0;              /* Output timing for each device. */
int Append=0;               /* append to output file */
int Verbose=0;              /* verbose messages in output file */
int Mute=0;                 /* no default output to stdout */
int Optoff=0;               /* switch off options parsing, */
                            /* to allow for args that begin with `-' etc */
FILE *debug=0;              /* debug output file  */
#if MPI
int Decomp_version=1;	    /* default: automatic decomposition by RMF method */
char *Decomp_string=NULL;   /* explicit decomposition formula */
#endif

/* others */
int narg=0;                 /* the additional arguments */
char **arg=NULL;            /* exported for def only !! */
int ReturnCode=EXIT_SUCCESS;/* return code: success by default */
int ExitAlone=0;	    /* if nonzero, only one process is aborting */
char *WMODE[2]={"wt","at"}; /* write/append (text) modes */
int ndev=0;                 /* num of devices in current run */
Device dev[MAXDEV];         /* the array of devices - exported for init only!!*/
int idev=-1;           	    /* device counter */
char *device_name;	    /* name of the device currently parsed or executed */

/* Static variables */
static char HEADLINE[256];	/* First line of output with package version variation and date */
static FILE *res;               /* Log file, separate from stdout in sequential mode */
static char logname[MAXPATH], debugname[MAXPATH];
static int iarg, argc_new;      /* argument counters */

/*  Timing variables */
static time_t start, finish;
static double run_start_time, run_end_time, device_start_time, device_end_time;

/* The Wallclock time, measured differently in different modes */
static double Wtime(void) {
#if MPI
    return MPI_Wtime();
#else
    return (double)clock()/(double)CLOCKS_PER_SEC;
#endif
}

/* Generic message procedure, writes both to stdout and log file. */
/* Works from any thread only if expected that only one thread may generate it. */
/* Outputs to stdout and log file separately if in sequential mode, and all to stdout in MPI mode. */
void ANY_MESSAGE(int urgent, char *fmt, ...) {
  char s[4096], *p;
  va_list argptr;
  if (mpi_rank==0 || urgent) { /*  Avoid repeated messages on MPI unless urgent */
    va_start(argptr, fmt);
    vsprintf(s, fmt, argptr);                 /* form the message */
    va_end(argptr);
    p=s;
    if (s[0]=='\x01') p++;                    /* \x01 - stdout mutes */
#if MPI
    fputs(p,stdout);
    FFLUSH(stdout);
#else
    if (res) {				    /* to log file - always */
      fputs(p,res);
      FFLUSH(res);
    } else {
      fputs(p,stdout);
      FFLUSH(stdout);
    }
    if (p==s && Mute==0) {		    /* to stdout */
      fputs(p,stdout); 
      FFLUSH(stdout);
    }
#endif
  }
}



void DEBUG(char *fmt, ...) {
  va_list argptr;
  if NOT(debug) return;
  va_start(argptr, fmt);
  vfprintf(debug, fmt, argptr);             /* form the message */
  va_end(argptr);
  /* FFLUSH(debug); */ /* debug is about thoroughness, not performance! */
  fflush(debug);
}

int nofflush(void *f) {return 0;}

/* Abnormal termination by one process only */
jmp_buf errjmp;
void ABORT(char *fmt, ...) {
  char s[4092], *p;
  va_list argptr;
  va_start(argptr, fmt);
  vsprintf(s, fmt, argptr);                     /* form the message */
  va_end(argptr);
  p=s;
  if (s[0]=='\x01') p++;                        /* \x01 - stdout mutes */
  if (res) {     /* to log file */
#if MPI
    fprintf(res,"\nProcess %d: %s",mpi_rank,p);
#else
    fprintf(res,"\n%s",p);
#endif
    FFLUSH(res);
    FFLUSH(res);
  }
  if (p<=s) 
#if MPI
    printf("\nProcess %d: %s",mpi_rank,p); /* to stdout */
#else
    printf("\n%s",p); /* to stdout */
#endif
  ExitAlone=1;
  longjmp(errjmp,1);
}

/* Print a usage message. */
#define HELPMSG \
"\nUsage: \n"\
"   Beatbox [options] [--] <parameter_file> [<log_file> [additional_args..]]"\
"\nOptions: \n" \
"   -- : switches the parsing of further options off.\n"\
"   -append : appends to the output file.\n"\
"   -debug <one of stdout, stderr or a filename> : to output debugging info to.\n"\
"   -help : print this message to stdout and exit.\n"\
"   -log <filename> : to log to, else if no filename defaults to \"parameter_file.log\".\n"\
"   -mute : mutes output to stdout.\n"\
"   -nograph : do not produce graphical output (off for MPI-based jobs).\n"\
"   -profile : outputs timings for each device.\n"\
"   -verbose : produces verbose messages in the output file.\n"

#define HELP() {MESSAGE("%s%s",HEADLINE,HELPMSG); goto Exit;}

#ifdef UNIX
#include <unistd.h>
int    _argc; 
char **_argv;
extern char WindowName[], IconName[];
#endif

/**********************
 * Start main program *
 **********************/
int main (int argc, char **argv) {

  int id;		/* secondary device loops counter */
  int rc;		/* return code of a device */

  /* Profiling stuff */
  int *timesCalled;	/* number of times each device was executed */
  double spent;		/* time spent for a device call */
  double *totalSpent;	/* total time spent by each device */
  double grandtotal;	/* total time spent by all devices (less than wall clock time by the overheads */
  int tc;		/* times this device was executed */
  double ts;		/* total time spent by this device */
  #if MPI
  int mpi_on=0;		/* if 1 then MPI has been initiated so MPI_Finalize/MPI_Abort is required */
  int localrc;		/* return code of a device in this thread */
  double totalspent;	/* time spent by a device call in all threads */
  #endif
  #undef MARK_ERROR
  #define MARK_ERROR   goto Exit;

  #if MPI
    /*  Start the MPI processes. */
    mpi_errno = MPI_Init(NULL, NULL);
    if(mpi_errno != MPI_SUCCESS) {
      MPI_Error_string(mpi_errno, error_string, &error_length);
      fprintf(stderr,"Could not initialise MPI.\n\tMPI Message: %s\n",error_string);
      exit(1);
    }
    MPI_Comm_set_errhandler(MPI_COMM_WORLD,  MPI_ERRORS_RETURN);
    mpi_on = 1;
    mpi_errno = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if(mpi_errno != MPI_SUCCESS) {
      MPI_Error_string(mpi_errno, error_string, &error_length);
      fprintf(stderr,"Could not get my rank.\n\tMPI Message: %s\n",error_string);
      exit(1);
    }
  #endif

  /*  Start timing */
  run_start_time = Wtime();

  sprintf(HEADLINE,
	  /* "#! %s %s version compiled %s %s\n" */
	  /* "#------------------------------------------------------------------------------------\n" */
	  /* ,VERSTRING,VARIATION,__DATE__,__TIME__); */
	  "#! %s revision %s, %s compile %s %s\n"
	  "#------------------------------------------------------------------------------------"
	  ,VERSTRING,SUBVERSION,VARIATION,__DATE__,__TIME__);

  /* Prepare for possible abort during start-up */
  t=-1;
  device_name="Beatbox start-up";
  if (setjmp(errjmp)) { ReturnCode=EXIT_FAILURE; goto Exit;} /* nonzero return code */

  
       
  #ifdef UNIX /* Generate a title for the X window */
    _argc=argc; _argv=argv;
    gethostname(IconName,1024); 
    strcat(IconName,":"); 
    strcat(IconName,argv[0]);
    gethostname(WindowName,1024); 
    strcat(WindowName,":");
    for(iarg=0;iarg<argc;iarg++) {
       strcat(WindowName,"  ");
       strcat(WindowName,argv[iarg]);
    }
  #endif

  /*********************/
  /* Parse arguments   */
  /*********************/

  /* First, extract options */
  argc_new=1;
  for (iarg=1;iarg<argc;iarg++) {
    if (argv[iarg][0] && strchr("-",argv[iarg][0]) && !Optoff) {
      #define CASE(a) else if (0==stricmp(argv[iarg]+1,a))
      if (argv[iarg][1]=='?')   HELP()
      CASE("-")         Optoff   = 1;
      CASE("append")    Append   = 1;
      CASE("debug")     {
	iarg++;
	if (iarg>=argc)
	  EXPECTED_ERROR("debug file name is not specified after -debug option")
	if (0==strcmp(strcpy(debugname,argv[iarg]),"stdout"))
	  debug=stdout;
	else if (0==strcmp(debugname,"stderr"))
	  debug=stderr;
	else if (0!=strcmp(debugname,"log")) {
	  if NOT(debug=fopen(debugname,"wt"))
		  ABORT("could not open '%s' for writing\n",debugname);
	} else
	  debug=NULL; /* for now; will be identified with log file */
      }
#if MPI
      CASE("decomp")	{
	iarg++;
	if (iarg>=argc)
	  EXPECTED_ERROR("decomposition method is not specified after -decomp option");
	Decomp_string=strdup(argv[iarg]);
	if (strchr(Decomp_string,'x')) {
	  Decomp_version=0;
	} else {
	  if (1!=sscanf(Decomp_string,"%d",&Decomp_version))
	    EXPECTED_ERROR("cannot read '%s' as decomposition version",Decomp_string);
	}
      }
#endif
      CASE("help")      HELP() 
      CASE("log")       {
	iarg++;
	if (iarg>=argc)
	  EXPECTED_ERROR("log file name is not specified after -log option");
	strcpy(logname, argv[iarg]);
      }
      CASE("mute")      Mute    = 1;
      CASE("nograph")   Graph   = 0;
      CASE("profile")   Profile = 1;
      CASE("verbose")   Verbose = 1;
      else              {MESSAGE("Unknown option %s",argv[iarg]); HELP();}
    } else {
      argv[argc_new++] = argv[iarg];
    }
  }
  argc = argc_new;

  /* Check number of arguments */
  if (argc<2) HELP();
  if (argv[1][0]=='?') HELP();

  /* First argument - input file */
  if (!qopen(argv[1])) EXPECTED_ERROR("Cannot open input file %s", argv[1]);

#if MPI
  /* In MPI, everything goes to stdout only */
  if (debug) {
    if (debug!=stdout && debug!=stderr && debug!=res)
      fclose(debug);
    debug=stdout;
    strcpy(debugname,"stdout");
  }
  res=NULL;
  if (*logname)
    MESSAGE("Log file '%s' will not be used in this MPI run; all messages will be to standard output\n",logname);
  *logname='\0';
#else
  /* In sequential mode, log file is separate from stdout. */
  /* Log file name is made from script name by default  */
  if NOT(*logname) fputext(argv[1],logname,"log");  /* default */
  if (0==strcmp(logname,null)) {
    res=NULL; 
    *logname='\0';
  } else if NOT(res=fopen(logname,WMODE[Append])) {
    EXPECTED_ERROR("Cannot open log file %s", logname);
  }
  if (0==strcmp(debugname,"log"))
    debug=res;
#endif

  /* Any further arguments are additional parameters */
  if (argc>2) {narg = argc-2; arg = argv+2;}

  /**********************/
  /* Header of log file */
  /**********************/
  MESSAGE(HEADLINE);
  time(&start);
  MESSAGE("\n# Execution begin at %24.24s ",ctime(&start));
  MESSAGE("\n# Input file %s ",argv[1]);
  if (narg) {
    MESSAGE("with additional argument");
    if(narg>1) MESSAGE("s");
    for (iarg=0;iarg<narg;iarg++)
      MESSAGE(" \"%s\"",arg[iarg]);
  } else {
    MESSAGE("without additional arguments");
  }

  MESSAGE("\n# with options:");
  if(Append) MESSAGE(" append"); else MESSAGE(" noappend");
  if(debug) MESSAGE(" debug=%s",debugname); else MESSAGE(" nodebug");
  if(Graph) MESSAGE(" graph");else MESSAGE(" nograph");
  if(Profile) MESSAGE(" profile");else MESSAGE(" noprofile");
  if(Verbose) MESSAGE(" verbose");else MESSAGE(" noverbose");
  if(*logname) MESSAGE(" logname=%s",logname);else MESSAGE(" logname=default");
  MESSAGE("\n");
  
  /**************************************/
  /* Read input file and build the task */
  /**************************************/
  device_name="init";
  if (!init()) EXPECTED_ERROR("Failure in reading input file");
  if (!ndev) EXPECTED_ERROR("No devices were created");
  MESSAGE("\n#------------------------------------------------------------------------------------\n\n");

  /********************/
  /* Perform the task */
  /********************/
  #if MPI
    /*  Inactive processes can wait for the others at the exit barrier. */
    if (I_AM_IDLE) goto Exit;
  #endif

  t = 0; /* Set the main loop counter to 0 */

  /*  Initialize profiling */
  if (Profile) {
    CALLOC(totalSpent,ndev, sizeof(double));
    CALLOC(timesCalled,ndev, sizeof(int));
    DEBUG("============================ PROFILING DATA ============================\n");
    DEBUG("t       ,");
    for (id=0; id<ndev; id++) {
      DEBUG("%8s",dev[id].n);
      if (id<(ndev-1)) {
        DEBUG(",");
      } else {
        DEBUG("\n");
      }
    }
  }
  
  /*********************************************************/
  /* Time loop. Exit for loop when any Device returns a 0. */
  /*********************************************************/
  for (;;) { 
    DEBUG("t=%ld:",t);
    for (idev=0; idev<ndev; idev++) {  /* Loop over the active devices */
      #define d (dev[idev])
      device_name=d.n;
      if (Profile) {
        /*  Add timestep to profile. */
        if(idev==0){
          DEBUG("%08d,",t);
        }
      }
      if (*(d.c)) { /* is the device condition met? */
	DEBUG(" %d(%s",idev,d.n);
        if (Profile) device_start_time = Wtime();
	#if MPI
	  rc=d.p(d.s, d.w, d.par, d.sync, d.alwaysRun);
          if (Profile) device_end_time = Wtime();
	#else
          if (Graph && (d.w.row1>d.w.row0 || d.w.col1>d.w.col0)) {
            opengraph(); /* if X unavailable, proceed without */
            if (!d.w.drawn) {
               SetWindow(d.w);
               Frame(d.s); 
               d.w.drawn=1;
            }
          }
          /*
            The same structure is used in the serial and parallel codes.
            So the structure dev[idev].XXXX has some dangling components.
          */
          rc=d.p(d.s, d.w, d.par);
	  DEBUG(")");
          if (Profile) device_end_time = Wtime();
  	#endif

        if (Profile) {
          spent = (device_end_time - device_start_time);
	  #if MPI
	  mpi_errno = MPI_Reduce(&spent, &totalspent, 1, MPI_DOUBLE, MPI_SUM, 0, ALL_ACTIVE_PROCS);
	  spent=totalspent/num_active_procs;
	  #endif
          DEBUG("%8f",spent);
          totalSpent[idev] += spent;
          timesCalled[idev] ++;
        }

	if (rc==0) /* This only happens in "organized termination" say by stop device, */
	  break;   /* so will be collective. Individual crashes are handles via ABORT. */

      } else { /*  this device is not run now */
        if (Profile) DEBUG("--------");
      }
      
      /*  Separator */
      if (Profile) DEBUG((idev<(ndev-1))?",":"\n");

    } /* End the device loop */

    if (rc==0) {
      break; 
    }
    t++; /* increment the loop counter */
    DEBUG("\n");
 
   }; /* End main loop */

/*************************/
/* Terminate the program */
/*************************/
Exit: /* this the target for the long jump */

#if MPI
  /* The barrier has to be for all threads not just active ones. */
  /* No point in waiting here if only the process where the error occurred */
  /* has jumped to here. Will try to close files, so too early to abort now. */
  if (mpi_on && ReturnCode==EXIT_SUCCESS) MPI_Barrier(MPI_COMM_WORLD);
  /* MESSAGE("\nall threads are through the exit barrier\n"); */

  if (!I_AM_IDLE) {
#endif
    /* more than 1 proc is hopefully active: run_end_time is for each proc. */
    run_end_time = Wtime();
    /* but it is reported only for rank = 0 */
    MESSAGE("\n\nTIMING INFO:\nThis run took %f seconds.\n\n", 
            run_end_time-run_start_time);
    
    if (Profile) {
      MESSAGE("Profiling summary:\n");
      MESSAGE("%-4s%-20s %8s %4s %16s %16s\n","No","Device name","# calls","%","ms per step","ms per call");
      grandtotal=0;
      for(id=0;id<ndev;id++)
	grandtotal+=totalSpent[id];
      for(id=0;id<ndev;id++) {
	tc=timesCalled[id];
	ts=totalSpent[id];
	MESSAGE("%3d %-20s %8d %4.1f %16.5f %16.5f\n",
		id,
		dev[id].n,
		tc,
		ts*100.0/grandtotal,
		1.e3*ts/(double)t,
		1.e3*ts/(tc?tc:1)
		);
      }	
      free(timesCalled);
      free(totalSpent);
    }

    if (graphon) {
      usleep(ReturnCode==EXIT_SUCCESS?200000:0);
      closegraph();
      graphon=0;
    }
    MESSAGE("\n%s finished at t=%ld by device %d \"%s\"\n",
            VERSTRING,t,idev,device_name);
    time(&finish);
    MESSAGE("%24.24s\n",ctime(&finish));
    MESSAGE("======================================================\n");
    if(res) fclose(res);
    if(debug) fclose(debug);

    /* Terminate the other devices cleanly */
    term();

#if MPI
  } /* If(!I_AM_IDLE) */

  /* Stop the MPI processes */
  if (mpi_on) {
    /* printf("%d: about to exit, rc=%d\n",mpi_rank,ReturnCode); */
    if (ExitAlone) {
      printf("\nProcess %d will abort\n",mpi_rank);
      fflush(stdout);
      /* MPI_Abort(ALL_ACTIVE_PROCS,ReturnCode); */
      MPI_Abort(MPI_COMM_WORLD, ReturnCode);
      /* printf("%d: aborted\n",mpi_rank); */
    } else {
      printf("\nProcess %d will finalize with return code %d\n",mpi_rank,ReturnCode);
      fflush(stdout);
      MPI_Finalize();
      /* printf("%d: finalized\n",mpi_rank); */
    }
  }
#endif

  /* printf("%d: will return with code %d\n",mpi_rank,ReturnCode); */
  fflush(stdout);
  return(ReturnCode);
} /* End of main */

