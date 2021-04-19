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

/*! Expands a string macro. 
 * \param name The name of the macro to expand.
 * \return Flag indicating success. 1 if successful, 0 otherwise.
 */
int qopen(char *name);

/*
	TODO Comment for Doxygen 
*/
int qclose(void);

/*! Points w to first word in s, and return pointer to next after
 *	\param s String to search.
 *	\param w Pointer to the first word.
 *      \param delim 
 *	\return Pointer to the rest of the search string.
 */
char *first_word (char *s, char **w,char *delim);

/*!
 *	Read to next terminator or end of file.
 *      \param s
 *      \param strlen
 *	\return Number of bytes read, or 0 if end of file (error) 
 *              is the first char.
*/
int read_command (char *s, int strlen);

/*!
 *	Evaluates a script expression.
 *	\param resaddr The address into which the result will be written.
 *	\param restype The datatype of the result.
 *	\param expr The expression to be evaluated.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int calc (void *resaddr, int restype, char *expr);

/*!
 *	Defines a k_ variable in the default user table (deftb).
 *	Parses a script expression such as "int cheese 5;", allocates memory for
 *	the result and inserts it into the symbol table.
 *	\param s The string to be parsed.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int def (char *s);

/*!
 *	Defines a k_ variable in the given symbol table.
 *	Parses a script expression such as "int cheese 5;", allocates memory for
 *	the result and inserts it into the symbol table.
 *	\param s The string to be parsed.
 *	\param table The symbol table into which the variable is to be inserted.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int def_local (char *s, p_tb table);

/*!
 *	Define a device name.
 *	\param d The device whose name is to be defined.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int def_dev (Device *d);

/*!
 *	Obtain a reference to a device from its name.
 *	\param n The name of the device to be found.
 *	\param d Address of the pointer that should be pointed to the named 
 *               device, if found.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int get_dev (Name n, Device **d);

/*!
 *  Find key in buffer, return to the next char after key or NULL if not 
 *  \param key
 *  \param s
 *  \return Pointer to the next character after the key or NULL if none.
 */
char *find_key (const char *key, char *s);

/*!
 *	Accepts an integer value, identified by name, from a given 
 *      parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param var Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" 
 *                   is found. If INONE, a script parameter MUST be found.
 *	\param minv Minimum allowed value. If INONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If INONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int accepti(const char *name,int *var,int deflt, int minv, int maxv, char *w);

/* Same, for (long) integer parameters linked to k-variables or k-expressions */
int acceptik(const char *name, INT **varptr, INT *varval, pp_fn *varcode, INT deflt, INT minv, INT maxv, char *w);

/*
 * Accepts an entry of an integer parameter array.
 */
int acceptie(const char *mask,int *arr,int i,int deflt,int minv,int maxv,char *w);

/*!
 *	Accepts a long value, identified by name, from a given parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param var Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" 
 *                   is found. If LNONE, a script parameter MUST be found.
 *	\param minv Minimum allowed value. If LNONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If LNONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int acceptl(const char *name,long *var,long deflt,long minv,long maxv,char *w);

/*!
 *	Accepts a real (double) value, identified by name, from a given 
 *      parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param var Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" is 
 *                   found. If RNONE, a script parameter MUST be found.
 *	\param minv Minimum allowed value. If RNONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If RNONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int acceptr(const char *name,real *var,real deflt,real minv,real maxv,char *w);

/* Same, for real parameters linked to k-variables or k-expressions */
int acceptrk(const char *name, REAL **varptr, REAL *varval, pp_fn *varcode, real deflt, real minv, real maxv, char *w);

/*!
 *	Accepts a real (double) value, identified by name with a number, from a given 
 *      parameter string.
 *	\param name Mask for the name of the parameter whose value is to be read.
 *	\param arr Pointer to beginning of the array where value is to be assigned.
 *	\param int Index of the array entry = number in the parameter name. 
 *	\param deflt Default value, in case no parameter called "name#number" is 
 *                   found. If RNONE, a script parameter MUST be found.
 *	\param minv Minimum allowed value. If RNONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If RNONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int acceptre(const char *mask,real *arr,int i,real deflt,real minv,real maxv,char *w);

#if defined(_rhs) || defined(_ionic)
/*!
 *	RHS ONLY
 *	Accepts a parameter or dependent parameter from the script.
 *	\param name  Name of the parameter to be found in the script.
 *	\param var   Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" 
 *                  is found. If RNONE, a script parameter MUST be found.
 *	\param minv Minimum allowed value. If RNONE, no minimum is prescribed.
 *	\param maxv Maximum allowed value. If RNONE, no maximum is prescribed.
 *	\param w Parameter string to search.
 *	\param v Address of the dependent-parameter-describing Var structure.
 *	\param iv Address of the variable indicating the index in v of the 
 *                current dependent parameter.
 *	\param v0 Lower variable bound for the calling device's Space 
 *                structure. Used to ensure that the v's src fields point to 
 *                variable 0.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int acceptp(const char *name,real *var,real deflt,real minv,real maxv,char *w,Var *v,int *iv,int v0);
#endif

/*!
 *	Accepts a string, identified by name, from a given parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param var Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" 
 *             is found. If NULL, a script parameter MUST be found.
 *	\param w Parameter string to search.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int accepts(const char *name,char *var,const char *deflt, char *w);
int acceptsn(const char *name,char *var,int len,const char *deflt, char *w);

/*!
 *	Accepts a file handle, identified by parameter name, from a given 
 *      parameter string.
 *	\param name Name of the parameter whose value is to be read.
 *	\param mode Access mode in which the file is to be opened. Syntax is 
 *                  identical to a fopen() call:
 *		- r or rb Open existing ﬁle for reading. 
 *		- w or wb Create ﬁle or wipe existing ﬁle before writing. 
 *		- a or ab Append to end of existing ﬁle, creating if necessary. 
 *		- rt or rbt or rtb Open existing ﬁle for update—reading and 
 *                writing. 
 *		- wt or wbt or wtb Create ﬁle or wipe existing ﬁle before 
 *                updating. 
 *		- at or abt or atb Append—Open or create ﬁle for update, 
 *                writing at end of ﬁle. 
 *	\param deflt Default filename, in case no parameter called "name" is 
 *             found. If NULL, a filename MUST be found in the script.
 *	\param fname Address at which to store the filename.
 *	\param f Pointer to the open file.
 *	\param w Parameter string to search.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int acceptf (const char *name,const char *mode,const char *deflt,char *fname,FILE **f, char *w);

/*!
 *	Parses a code-string; an expression in Beatbox script syntax
 *	that may be evaluated using calc(). A code-string is fetched
 *	as a complete expression rather than as a reduced value so
 *	that it may be computed at runtime using k_ variables.
 *	\param name Name of the parameter whose value is to be read.
 *	\param type Data type returned by the code-string when executed. 
 *	\param type Address of variable where value is to be assigned.
 *	\param deflt Default value, in case no parameter called "name" is 
 *             found. If NULL, a script parameter MUST be found.
 *	\param codestring Address at which the executable code should be stored.
 *	\param code Pointer to the compiled code block.
 *	\param w Parameter string to search.
 *	\return Flag indicating success. 1 if successful, 0 otherwise.
 */
int acceptc (const char *name,const k_type type,const char *deflt,char *codestring,void **code, char *w);

/*!
 *	Accepts a block of Beatbox script code from a given parameter string.
 *	A block is any number of lines of code contained in braces.
 *	\param name Name of the script parameter to which the block is 
 *             assigned, e.g. name={BLOCK}.
 *	\param var Address of the variable to which the code block is to 
 *             be assigned.
 *	\param w Parameter string to search.
 */
int acceptb(const char *name,char *var, char *w);

/*
 *	UNSURE This just prints a message about separators. What's
 *	going on? Vestigial?
 */
void close_block(void);

/*!
 *	Accepts a device's condition from a given parameter string.
 *	The condition is evaluated in beatbox.c to decide whether or not a 
 *      device should be run at the current time step.
 *	\param c Address at which condition (real/double value) should be 
 *             stored.
 *	\param w Parameter string to search for "when" key.
 */
int accept_condition(Condition *c, char *w);

/*!
 *	Accepts a real k_variable from a given parameter string.
 *	\param v Address at which condition (real/double value) should be 
 *               stored.
 *      \param vname Name of the variable found.
 *	\param name Name of the variable to be found.
 *	\param w Parameter string to search for name key.
 */
int accept_real_variable(REAL **v, char **vname, char *name, char *w);

/*!
 *	Tests for state dimensions (xmax,ymax,zmax) in the parameter string 
 *      provided.
 *	\param parameterString The parameter string to search.
 *	\return True if state dimensions are found, false otherwise.
 */
int stateDimensionsExist(char *parameterString);

/*
	TODO Improve documentation.
*/
int spaceParametersExist(char *parameterString);

/*!
 *	Accepts a device's Space constraints from a given parameter string.
 *	\param s Address of Space structure to be initialised.
 *	\param w Parameter string to search.
 */
int accept_space(Space *s, char *w);

/*#ifdef IBMPC*/
int accept_window(BGIWindow *s, char *w);
/*#endif*/

double _u(double _x, double _y, double _z, double _v);

/*!
 *	Shortcut for accepting an integer parameter via accepti()
 *	If accepti() returns false, the function in which ACCEPTI() is called will also return false.
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 *	\sa accepti()
 */
#define ACCEPTI(b,c,d,e) if (!accepti(#b"=",&(S->b),c,d,e,w)) return(0); int b=S->b

/*!
 *      Same, for a parameter (potentially) linked to a k-variable or k-expression
 */
#define ACCEPTIK(b,c,d,e) if (!acceptik(#b"=",&(S->b##ptr),&(S->b),&(S->b##code),c,d,e,w)) return(0); 

/*!
 *	Shortcut for accepting an integer array element via acceptie()
 *	If acceptie() returns false, the function in which ACCEPTIE() is called 
 *      will also return false.
 *	\param m Mask for the enumerated parameter name to be found (with %d or something like that).
 *	\param a The pointer to the accepted array. 
 *	\param i The index of the accepted array element, also used for enumeration of parameter name. 
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 *	\sa acceptre()
 */
#define ACCEPTIE(m,a,i,c,d,e) if (!acceptie(m"=",a,i,c,d,e,w)) return(0)

/*!
 *	Shortcut for accepting a long parameter via acceptl()
 *	If acceptl() returns false, the function in which ACCEPTL() is called 
 *      will also return false.
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 *	\sa acceptl()
 */
#define ACCEPTL(b,c,d,e) if (!acceptl(#b"=",&(S->b),c,d,e,w)) return(0); long b=S->b

/*!
 *	Shortcut for accepting a real (double) parameter via acceptr()
 *	If acceptr() returns false, the function in which ACCEPTR() is called 
 *      will also return false.
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 *	\sa acceptr()
 */
#define ACCEPTR(b,c,d,e) if (!acceptr(#b"=",&(S->b),c,d,e,w)) return(0); real b=S->b


/*!
 *      Same, for a parameter (potentially) linked to a k-variable or k-expression
 */
#define ACCEPTRK(b,c,d,e) if (!acceptrk(#b"=",&(S->b##ptr),&(S->b),&(S->b##code),c,d,e,w)) return(0); 

/*!
 *	Shortcut for accepting a real array element via acceptre()
 *	If acceptre() returns false, the function in which ACCEPTRE() is called 
 *      will also return false.
 *	\param m Mask for the enumerated parameter name to be found (with %d or something like that).
 *	\param a The pointer to the accepted array. 
 *	\param i The index of the accepted array element, also used for enumeration of parameter name. 
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 *	\sa acceptre()
 */
#define ACCEPTRE(m,a,i,c,d,e) if (!acceptre(m"=",a,i,c,d,e,w)) return(0)

#if defined _rhs || defined _ionic
/*!
 *	RHS ONLY:
 *	Shortcut for accepting a parameter or dependent parameter via acceptp()
 *	If acceptp() returns false, the function in which ACCEPTP() is called 
 *      will also return false.
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\param d Minimum value.
 *	\param e Maximum value.
 *	\sa acceptp()
 */
#define ACCEPTP(b,c,d,e) if (!acceptp(#b"=",&(S->b),c,d,e,w,var,&ivar,v0)) return(0); real b=S->b
#endif

/*!
 *	Shortcut for accepting a string parameter via accepts()
 *	If accepts() returns false, the function in which ACCEPTS() is called 
 *      will also return false.
 *	\param b Parameter name/key to be found.
 *	\param c Default value.
 *	\sa accepts()
 */
#define ACCEPTS(b,c)     if (!accepts(#b"=",&(S->b[0]),c,w)) return(0); char *b=&(S->b[0])
#define ACCEPTSN(b,n,c)     if (!acceptsn(#b"=",&(S->b[0]),n,c,w)) return(0); char *b=&(S->b[0])

/*!
 *	Shortcut for accepting a file via acceptf()
 *	If acceptf() returns false, the function in which ACCEPTF() is called 
 *      will also return false.
 *	\param b Parameter name/key to be found.
 *	\param c Access mode.
 *	\param d Default value.
 *	\sa acceptf()
 */
#define ACCEPTF(b,c,d)   if (!acceptf(#b"=",c,d,&(S->b##name[0]),&(S->b),w)) return(0); FILE *b=S->b; char *b##name=&(S->b##name[0])

/*!
 *	Shortcut for accepting a code string via acceptc()
 *	If acceptc() returns false, the function in which ACCEPTC() is called will also return false.
 *	\param b Parameter name/key to be found.
 *	\param c k_type to be returned by the code string.
 *	\param d Default code string.
 *	\sa acceptc()
 */
#define ACCEPTC(b,c,d)   if (!acceptc(#b"=",c,d,&(S->b[0]),&(S->b##code),w)) return(0)

/*!
 *	Shortcut for accepting a real k_variable via accept_real_variable()
 *	If accept_real_variable() returns false, the function in which 
 *      ACCEPTV() is called will also return false.
 *	\param a Variable name/key to be found
 *	\sa accept_real_variable()
 */
/* 8< mpi_acceptv_macro */
#define ACCEPTV(a)   if (!accept_real_variable(&(S->a),&(S->a##name),#a,w)) return(0)
/* >8 mpi_acceptv_macro */

/*!
 *	Shortcut for accepting a code block via acceptb()
 *	If acceptb() returns false, the function in which BEGINBLOCK() is 
 *      called will also return false.
 *	\param a Parameter name/key to be found.
 *	\param b Address to which code block should be assigned.
 *
 *	\sa acceptb()
 */
#define BEGINBLOCK(a,b)    if (!acceptb (a,b,w)) return(0)

/*
	UNSURE Why would you want to use close_block()?
*/

/*!
 *	Shortcut for close_block()
 *	\sa close_block()
 */
#define ENDBLOCK	   close_block()

/* initiation of constants used in checking limits etc */
/*
 *	UNSURE Dunno what's happening here. Looks like "minimum, but not quite".
 */
real RSUCC(real val);
real RPRED(real val);


/*! Initiation of constants used in checking limits etc */
int init_const(void);
/* Clear the symbol table and deallocate variables */
int term_const(void);

#include "extern.h"
			
			/***** Exported variables *****/
/*! Maximum depth of nested includes in script files.*/
#define MAXDEPTH 10

/*! Current depth of include. */
EXTERN int depth;		

/*! Names of includes. */
EXTERN char inname[MAXDEPTH][MAXPATH];

/*! Coordinates in file of includes. */
EXTERN long inpline[MAXDEPTH], inppos[MAXDEPTH]; 

/*! Integer null value (-MAXINT+1). */
EXTERN int  INONE;

/*! Long null value (-MAXLONG+1). */
EXTERN long LNONE;

/*! Real (double) null value (-MAXREAL). */
EXTERN real RNONE;

/*! Real (double) "infinite" value (+MAXREAL). */
EXTERN REAL real_inf;

EXTERN real macheps;

/*! Null String */
#define null "null"

/*! Similar to assert() but cannot be screened */
/* #define ASSERT(p) { if(0==(p)) EXPECTED_ERROR("Assertion failed:\n %s",#p); } - defined in error.h */
#undef EXTERN

/*! Command terminators */
#define TERMINATOR ";$"

/*! Escape character */
#define ESCAPE "#"

/*! Token separators
decimal point ".", and colon ":" are used in numbers and file ids, "," in functions */
#define SEPARATORS " ;\t\r\n!$"

/* reference to a file or stdin */
/*! Marks the beginning of a reference to a file or stdin. */
#define INCLUDEBEGIN '<'
/*! Marks the end of a reference to a file or stdin. */
#define INCLUDEEND '>'

/* reference to a string macro */
/*! marks the beginning of a reference to a string macro. */
#define PASTEBEGIN '['
/*! marks the end of a reference to a string macro. */
#define PASTEEND ']'

/* reference to a system command output */
/*! marks the beginning of a reference to a system command output. */
#define CATCHBEGIN '`'
/*! marks the end of a reference to a system command output. */
#define CATCHEND '`'

/* Block of parameters */
/*! marks the beginning of a block of parameters. */
#define BLOCKBEGIN '{'
/*! marks the end of a block of parameters. */
#define BLOCKEND '}'

/* String with blanks */
/*! marks the end of a string with blanks. */
#define STRBEGIN '\"'
/*! marks the end of a string with blanks. */
#define STREND '\"'

/*! reference to an additional argument passed to the program
 * #define AMPERS '&'
 * Linking parameter to a level of dynamic variables */
#define AT '@'

/* links device parameter to a k-variable */
#define LINK '~'

/*! how default values are shown in res file */
#define DFLT "!"

