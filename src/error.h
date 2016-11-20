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


/* Generic error reporting routines */

#ifndef __ERROR_H
#define __ERROR_H

/* Most urgent fail: */
/* Print msg on error and do a long jump to main */
/* This is required when the failure may occur not in all processes */
void ABORT(char *fmt, ...);


/* Asks Beatbox to fail gracefully:                                */
/* print msg on error and perform MARK_ERROR defined in user pgm.  */
/* This looks like a half-measure between ABORT and EXPECTED_ERROR */
/* for which there is no legitmate place.                          */
/* #define ERROR(...) {					\ */
/*   printf("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   printf(__VA_ARGS__);MARK_ERROR} */
/* #define ERROR0(msg) {						\ */
/*   printf("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   printf(msg);MARK_ERROR} */
/* #define ERROR1(msg,a1) {					\ */
/*   printf("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   printf(msg,a1);MARK_ERROR} */
/* #define ERROR2(msg,a1,a2) {					\ */
/*   printf("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   printf(msg,a1,a2);MARK_ERROR} */
/* #define ERROR3(msg,a1,a2,a3) {					\ */
/*   printf("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   printf(msg,a1,a2,a3);MARK_ERROR} */
/* #define ERROR4(msg,a1,a2,a3,a4) {				\ */
/*   printf("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   printf(msg,a1,a2,a3,a4);MARK_ERROR} */
/* #define ERROR5(msg,a1,a2,a3,a4,a5) {				\ */
/*   printf("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   printf(msg,a1,a2,a3,a4,a5);MARK_ERROR} */
/* #define ERROR6(msg,a1,a2,a3,a4,a5,a6) {				\ */
/*   printf("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   printf(msg,a1,a2,a3,a4,a5,a6);MARK_ERROR} */
/* #define ERROR8(msg,a1,a2,a3,a4,a5,a6,a7,a8) {			\ */
/*   printf("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   printf(msg,a1,a2,a3,a4,a5,a6,a7,a8);MARK_ERROR} */

/* Lower class of error. 
 * Asks Beatbox to fail gracefully, but can get by
 * with only process 0 producing a message.
 *
 * Useful for user errors in the script, etc.
 */
#define EXPECTED_ERROR(...) {					\
  MESSAGE("\nError in %s at line %u: ",__FILE__,__LINE__);       \
  MESSAGE(__VA_ARGS__);MARK_ERROR}
/* #define EXPECTED_ERROR0(msg) {						\ */
/*   MESSAGE("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   MESSAGE(msg);MARK_ERROR} */
/* #define EXPECTED_ERROR1(msg,a1) {					\ */
/*   MESSAGE("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   MESSAGE(msg,a1);MARK_ERROR} */
/* #define EXPECTED_ERROR2(msg,a1,a2) {					\ */
/*   MESSAGE("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   MESSAGE(msg,a1,a2);MARK_ERROR} */
/* #define EXPECTED_ERROR3(msg,a1,a2,a3) {					\ */
/*   MESSAGE("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   MESSAGE(msg,a1,a2,a3);MARK_ERROR} */
/* #define EXPECTED_ERROR4(msg,a1,a2,a3,a4) {				\ */
/*   MESSAGE("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   MESSAGE(msg,a1,a2,a3,a4);MARK_ERROR} */
/* #define EXPECTED_ERROR6(msg,a1,a2,a3,a4,a5,a6) {				\ */
/*   MESSAGE("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   MESSAGE(msg,a1,a2,a3,a4,a5,a6);MARK_ERROR} */
/* #define EXPECTED_ERROR8(msg,a1,a2,a3,a4,a5,a6,a7,a8) {			\ */
/*   MESSAGE("\nError in %s at line %u: ",__FILE__,__LINE__);       \ */
/*   MESSAGE(msg,a1,a2,a3,a4,a5,a6,a7,a8);MARK_ERROR} */

#define ASSERT(p) { if(0==(p)) EXPECTED_ERROR("Assertion failed:\n %s\n",#p); }

#endif
