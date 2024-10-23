/*
*		Useful Extensions to C language:
*
*	Copyright 1999 by The University at Stony Brook, All rights reserved.
*/

#if !defined(_CDECS_H)
#define _CDECS_H

#if defined(float)
#define __double_precision__
#undef float
#endif /* defined(float) */

#if !defined(USE_OVERTURE) 
#if defined(linux)
#  if !defined(_GNU_SOURCE)
#    define _GNU_SOURCE
#  endif /*!defined(_GNU_SOURCE)*/
#endif /* defined(linux) */
#endif /* if !defined(USE_OVERTURE) */  

#include <stdlib.h>
#include <stdio.h>
#if defined(__cplusplus)
/* Change to standard C++ header
#  include <strings.h>
*/
#  include <cstdlib> // 09/14/2017 added
#  include <string>
#  include <iostream>
#  include <fstream>
#  include <cmath>
#  include <iomanip> // 09/14/2017 added
#  include <list>    // 09/14/2017 added
/* added for geodesic */
#  include <vector>
#  include <limits>
/*END: added for geodesic */

#  if defined(float)
#    undef float
#    include <algorithm>
#    define float double
#  else
#    include <algorithm>
#  endif
#endif /* defined(__cplusplus) */
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#if defined(__INTEL_COMPILER)
#include <stdint.h>
#endif /* defined(__INTEL_COMPILER) */
#include <limits.h>
#include <float.h>
#include <errno.h>

#if defined(__MPI__)
#   include <mpi.h>
#endif /* defined(__MPI__) */

#if defined(__GNUC__)
typedef u_int64_t uint64_t;
#elif defined(__alpha)
typedef long long int int64_t;
typedef unsigned long long uint64_t;
#endif /*defined(linux) || defined(__alpha)*/

#define ptr2ull(p) u_ptr2ull((void*)(p))

#if defined(__double_precision__)
#define float double
#undef __double_precision__
#endif /* defined(__double_precision__) */

#if defined(mips) || defined(__GNUC__) || defined(linux) || __SUNPRO_CC>=0x500 || __SUNPRO_C>=0x500

#define HasGen 1

#else /* #define HasGen defined(mips) || defined(__GNUC__) || defined(linux) || __SUNPRO_CC>=0x500 || __SUNPRO_C>=0x500 */

#define HasGen 0

#endif /* #define HasGen defined(mips) || defined(__GNUC__) || defined(linux) || __SUNPRO_CC>=0x500 || __SUNPRO_C>=0x500 */

#if HasGen
#   include <libgen.h>
#endif /* HasGen */
#include <limits.h>

#if defined(USE_OVERTURE)
/* Including header files of overture into FronTier */
#include "FT_overture.h"  
#endif /* defined(USE_OVERTURE) */

#include "my_vector.h"

#undef HUGE
#define HUGE 1.e+18

		/* Machine Dependent Quantities: */

#if defined(cray)
#   define isnan(x)     (NO)
#   define isnanf(x)    (NO)
#   define isnand(x)    (NO)
#elif defined(sun) && !(defined(__SUNPRO_C) || defined(__SUNPRO_CC))
#   define isnanf(x)    isnan(x)
#   define isnand(x)    isnan(x)
#elif defined(__GNUC__) || defined(__PGI__) || defined(__INTEL_COMPILER)
#   if !defined(isnan)
#     if !defined(float)
#       define isnan(x) isnanf(x)
#     endif /* !defined(float) */
#   endif /*!defined(isnan)*/
#elif !defined(__alpha) && !defined(__hpux) && !defined(linux) && !defined(_AIX)
#   if defined(float)
#       undef float
#       include <ieeefp.h>
#       define float double
#       define isnan(x)    isnand(x)
#   else /* defined(float) */
#       include <ieeefp.h>
#       define isnan(x)    isnanf(x)
#   endif /* defined(float) */
#endif /* defined(cray) */

#if defined(float)
#   define MACH_EPS DBL_EPSILON
#else /* defined(float) */
#   define MACH_EPS FLT_EPSILON
#endif /* defined(float) */

#if defined(__hpux)
#if defined(isfinite)
#define finite(x) isfinite(x)
#endif /* defined(isfinite) */
#endif /* defined(__hpux) */

/*
*		Control of single vs. Double Precision: 
*
*	Single precision, double precision or mixed precision storage
*	of variables can be controlled by setting the macros `float'
*	and `REAL'.  `float' may be undefined (in which case
*	it retains its  system declaration of
*		`typedef float single'
*	or defined to be `double'.  REAL should be set to either `float' or
*	`double'. The default is single precision storage.  `float' and `REAL'
*	should be defined in the appropriate `makefile' in order to
*	achieve non-single-precision code as indicated below.
*
*	To get a fully double precision code:
*		add -Dfloat=double to the makefile and use libdutil.a
*
*	To get a hybrid version:
*		add -DREAL=double to the makefile and use libutil.a
*
*	To get a fully single precision code:
*		leave float and REAL undefined and libutil.a
*
*	In the hybrid case all variables declared to be of type REAL will
*	be double, while all others will have their declared types.
*/

#if defined(float)
#   undef float
    typedef float   truefloat;
#   define float   double
#else  /* defined(float) */
    typedef float   truefloat;
#endif /* defined(float) */

#if defined(float)

#   if defined(REAL)
#      undef REAL
#   endif /* defined(REAL) */
#   define REAL double

#elif !defined(REAL)

#   define REAL float
#   define RealEqualsFloat

#endif /* defined(float) */

#if defined(float) || defined(cray)
#   define SFMT "22s"
#   define FFMT "22.16g"
#else /* defined(float) || defined(cray) */
#   define SFMT "14s"
#   define FFMT "14.8g"
#endif /* defined(float) || defined(cray) */

typedef	double 	ALIGN;	/* Strictest alignment type of all data types */
#define	num_aligns(size)	(((size) + sizeof(ALIGN)-1)/sizeof(ALIGN))
enum {
	Gets_BUF_SIZE = 513
};

		/* Macros for REAL variables and Vectors: */

typedef unsigned char byte;	/* Generic object of size 1, sizeof(byte) = 1 */

#if defined(__cplusplus)
#define NO                 false
#define FUNCTION_FAILED    false
#define YES		   true
#define FUNCTION_SUCCEEDED true
#else /* defined(__cplusplus) */
enum _bool { false              = 0,
		NO                 = false,
		FUNCTION_FAILED    = false,
		true               = 1,
		YES                = true,
		FUNCTION_SUCCEEDED = true};
typedef enum _bool bool;
#endif /* defined(__cplusplus) */

#if defined(ERROR)
#  undef ERROR
#endif /* defined(ERROR) */
enum {ERROR=-1};

typedef  void	    *POINTER;	/* Pointer to Unknown Data Types. */
typedef  const void *CPOINTER;	/* Constant Pointer to Unknown Data Types. */

#define  ERROR_FLOAT  -HUGE_VAL	/* Returned by Float Functions on Errors */


		/* stuff for parallel computers */

enum {
	IO_NODE_ID = 0,    /* processor node used as IO node*/
		/* identifier tags for timed global op  */
	ALL_GATHER_TAG,
	GLOBAL_IOR_TAG,
	GLOBAL_ISUM_TAG,
	GLOBAL_SUM_TAG,
	GLOBAL_IMAX_TAG,
	GLOBAL_MAX_TAG,
	GLOBAL_IMIN_TAG,
	GLOBAL_MIN_TAG,
	SYNC_TAG,
	GLOBAL_STATUS_TAG,
	USER_MIN_MESG_ID = GLOBAL_STATUS_TAG + 10
};

#define is_io_node(id)  ((id) == IO_NODE_ID)

/*      Re-direct standard I/O to files:
*       This structure is used on both host and node: normally two
*       entirely different machies.
*/
 
typedef struct {
	char    stderr_file[80];
	char    stdin_data[80];
	char    stdin_argv[80];
	char    stdout_file[80];
} REDIRECT_FILES;

#define	pp_isend(tag,buf,len,node,request)				\
	u_pp_isend(tag,buf,len,node,request,__FILE__,__LINE__)

#define	pp_send(tag,buf,len,node)					\
	u_pp_send(tag,buf,len,node,__FILE__,__LINE__)

#define	pp_send_all(tag,buf,len)					\
	u_pp_send_all(tag,buf,len,__FILE__,__LINE__)

#define	pp_irecv(tag,source,buf,len,request)				\
	u_pp_irecv(tag,source,buf,len,request,__FILE__,__LINE__)

#define	pp_recv(tag,source,buf,len)					\
	u_pp_recv(tag,source,buf,len,__FILE__,__LINE__)

#define	pp_recv_any(tag,buf,len)					\
	u_pp_recv_any(tag,buf,len,__FILE__,__LINE__)

#define	pp_all_gather(sendbuf,sendcount,recvbuf,recvcount)		\
	u_pp_all_gather(sendbuf,sendcount,recvbuf,recvcount,__FILE__,__LINE__)

#define pp_bcast(root,buf,len)						\
	u_pp_bcast(root,buf,len,__FILE__,__LINE__)

	/* Enclosed Code Compiled iff Nonzero */
#define  DONT_COMPILE  0

		/*  Inline Functions - Be Careful: */

#if defined(__cplusplus) && !defined(__hpux)

using std::min;
using std::max;

#else /* defined(__cplusplus) && !defined(__hpux) */

#if !defined(max)
#  define     max(a,b)     (((a) > (b)) ? (a) : (b))
#endif /* !defined(max) */

#if !defined(min)
#define	    min(a,b)     (((a) < (b)) ? (a) : (b))
#endif /* !defined(min) */

#endif /* defined(__cplusplus) && !defined(__hpux) */

#define     sqr(x)       ((x)*(x))

#define     cub(x)       ((x)*(x)*(x))

#define     quad(x)      ((x)*(x)*(x)*(x))


/* Structure copying; num is the number of chars to be copied */

//#define assign(l,r,num)	(void) memcpy((void*)(l),(const void*)(r),num)

/* Structure zeroing; num is the number of chars to be zeroed */

#define zero_scalar(s,num)	(void) memset((POINTER)(s),0,num)

enum {
	READ	  = 0,
	WRITE,
	READ_WRITE
};

#define   IMPORT        extern
#define   LOCAL         static
#define   EXPORT      
#define   LIB_LOCAL

/****
*******************Proper use of IMPORT,  LOCAL,  and EXPORT ******************

The explanation of the LOCAL vs.  static is due to the split nature of
the construction in C.  The static declaration really does two things,
first it declares a variable to be static,  that is its value is retained
from one function call to another,  the second effect is that the variable
becomes local,  that is its access is restricted to the enclosing block.
Actually in C all variables declared within a block (brackets) are local,
so the second case only applies to external variables in a file.  The
default is that externals are global (ie known to all files),  so that to
restrict the scope of an external variable the static is declaration is
used.  Note that all external variable are by default static and global.
The defines EXPORT,  IMPORT,  and LOCAL were designed to explicitly
declare the scope of a given procedure or variable,  thus leaving the
declaration static to perform its first (primary) function.  Thus the
declarations should be used as follows

EXPORT - any variable or procedure with global scope,  this includes
all functions to be called outside of their file of definition,  as well
as any transfile global variables (these are strongly discouraged).

IMPORT - any variable or procedure which is defined elsewhere,  IMPORT
is really just another name for extern and serves exactly the same function.
IMPORTED variable are just that,  variables or procedures that are imported
from somewhere else.

LOCAL - any variable or procedure whose scope is to be limited to its
file of definition.

static - any variable (presumably within a block) whose value is static and
will not change between function calls to that block.

Even though in C LOCAL and static are the same,  we see that the division
of this single variable into two types provides a useful documentation
function by separating the static declaration from the access restriction.

*****/


#if defined(cray)
#   include <fortran.h>
#   if defined(__cplusplus)
#       define FORTRAN IMPORT "C"
#   else /* defined(__cplusplus) */
#       define	FORTRAN	fortran
#   endif /* defined(__cplusplus) */
#else /* defined(cray) */
#   if defined(__cplusplus)
#       if defined(FORTRAN)
#           undef FORTRAN
#       endif /* defined(FORTRAN) */
#       define FORTRAN IMPORT "C"
#   else /* defined(__cplusplus) */
#       define	FORTRAN	IMPORT
#   endif /* defined(__cplusplus) */
#endif /* defined(cray) */

#if defined(cray) || defined(_AIX)
#   define   FORTRAN_NAME(a)    a
#else /* defined(cray) || defined(_AIX) */
#   define   FORTRAN_NAME(a)    a ## _
#endif /* defined(cray) || defined(_AIX) */
#if defined(__GNUC__) && !defined(LAHEY_FORTRAN)
#   define   SFORTRAN_NAME(a)   a ## __  /* GCC has a different naming
					  * convention for FORTRAN routines
					  * with character strings as
					  * arguments.
					  */
#else /* defined(__GNUC__) */
#   define   SFORTRAN_NAME(a)   FORTRAN_NAME(a)
#endif /* defined(__GNUC__) */

#if defined(LAHEY_FORTRAN)
#   define C_MAIN_PROGRAM MAIN__
#   if defined(__cplusplus)
        extern "C" int MAIN__(int,char**);
#   endif /* defined(__cplusplus) */
#else /* defined(LAHEY_FORTRAN) */
#   define C_MAIN_PROGRAM main
#endif /* defined(LAHEY_FORTRAN) */

#if defined(__cplusplus)
extern "C" {
#endif /* defined(__cplusplus) */
    typedef void (LSODE_FUNC)(int*,float*,float*,float*);
    typedef void (LSODE_JAC)(int*,float*,float*,int*,int*,float*,int*);
#if defined(__cplusplus)
}
#endif /* defined(__cplusplus) */

#define	Null(x)		(!(x))
#if !defined(PI)
#   define	PI		3.14159265358979
#endif /* !defined(PI) */
#define	degrees(ang)	((ang)*180.0/PI)
#define radians(ang)    ((ang)*PI/180.0)

struct _Error {
	const char *filename;
	int        line_number;
	int        number;
	const char *message;
	struct _Error *next;
};

#define Error(__num,__mess)  log_error(__FILE__,__LINE__,__num,__mess)

/* Opaque holder for storing location in file for output */

typedef void OUTPUT;

#define	fclose(file)	Fclose(file)

	/* Simpleio Macros */

#define	print_double_matrix(title, rows, cols, mtrx, fmt) \
	fprint_double_matrix(stdout,title, rows, cols, mtrx, fmt)

#define	print_double_vector( title, n, vec, fmt) \
	fprint_double_vector(stdout, title, n, vec, fmt)

#define	print_double_vector_as_matrix( title, length, cols, vec, fmt) \
	fprint_double_vector_as_matrix(stdout, title, length, cols, vec, fmt)

#define	print_float_matrix( title, rows, cols, mtrx, fmt) \
	fprint_float_matrix(stdout, title, rows, cols, mtrx, fmt)

#define	print_float_vector( title, n, vec, fmt) \
	fprint_float_vector(stdout, title, n, vec, fmt)

#define	print_float_vector_as_matrix( title, length, cols, vec, fmt) \
	fprint_float_vector_as_matrix(stdout, title, length, cols, vec, fmt)

#define	print_matrix( title, rows, cols, mtrx, fmt) \
	fprint_matrix(stdout, title, rows, cols, mtrx, fmt)

#define	print_int_matrix( title, rows, cols, mtrx, fmt) \
	fprint_int_matrix(stdout, title, rows, cols, mtrx, fmt)

#define	print_int_vector_as_matrix( title, length, cols, vec, fmt) \
	fprint_int_vector_as_matrix(stdout, title, length, cols, vec, fmt)

#define	print_vector( title, n, vec, fmt) \
	fprint_vector(stdout, title, n, vec, fmt)

#define	print_vector_as_matrix( title,length, cols, vec, fmt) \
	fprint_vector_as_matrix(stdout, title,length, cols, vec, fmt)

#define	print_vector_of_floats( num, f) \
	fprint_vector_of_floats(stdout, num, f)


 /* True iff  x  lies between  y  and  z  */
#define  Between(x,y,z)   ( ((x)>=(y) && (x)<=(z)) || ((x)<=(y) && (x)>=(z)) )

#if DONT_COMPILE
 /* Alternate definition of Between, may be faster for floating comparisons but
  * is subject to inaccurate answers due to floating point degeneracies. May
  * also give inaccurate results for integer comparisons. */
#define  Between(x,y,z)   ( ((x)-(y)) * ((z)-(x)) >= 0. )
#endif /* DONT_COMPILE */

struct _Prompt_type {
	const char *prompt;	/* Full prompt name */
	const char *select;  	/* Abbreviated name for input */
	int	   ncmp;	/* # of chars to uniquely identify select */
	union {int itype; const char *ctype;} type;
};
typedef struct _Prompt_type Prompt_type;


	/* Possible values of variable  debug_mode: */

enum _DEBUG_MODE {
	PROMPT_FOR_DEBUG=-1,	/*Initiate debug prompting*/
	NONE=             0,	/* Indicates no debugging requested */
	SOME=             1,	/* Indicates debugging requested not with all */
	ALL=              2,	/* Indicates all was the first request */
	TRACE_ONLY=	  3	/* Trace debug calls but don't print messages*/
};
typedef enum _DEBUG_MODE DEBUG_MODE;

struct _DEBUG_PARAMS {
	DEBUG_MODE	_debug_mode;
};
typedef struct _DEBUG_PARAMS DEBUG_PARAMS;
#define	debug_params(prms)	((DEBUG_PARAMS*)(prms))
#define debug_mode(prms)	(debug_params(prms))->_debug_mode

/*
* Basic initialization structure.
*/

/*
* Machine Endian type
*/

enum _FT_ENDIAN {
	FT_LITTLE_ENDIAN  = -1, /* Little endian byte ordering */
	FT_UNKNOWN_ENDIAN =  0, /* Undetermined byte ordering */
	FT_BIG_ENDIAN     =  1 /* Big endian byte ordering */
};
typedef enum _FT_ENDIAN FT_ENDIAN;

/*
* Return status from quadrature programs
*/

enum _QUADRATURE_STATUS {
        ACCURATE_INTEGRAL    = 0,
	INACCURATE_INTEGRAL  = 1,
	INVALID_EPSILON      = 6
};
typedef enum _QUADRATURE_STATUS QUADRATURE_STATUS;

/* Reverse bytes of argument to convert endian sex */
#define byte_reverse(_x_)	reverse_string((char*)&_x_,sizeof(_x_))

/*
* IO conversion specifications
*/

struct _IO_TYPE {
	FILE        *file;
	size_t      read_float_size;
	size_t      cpu_float_size;
	FT_ENDIAN   read_endian;
	FT_ENDIAN   cpu_endian;
	bool        reverse_endian;
};
typedef struct _IO_TYPE IO_TYPE;

struct _INIT_DATA {
	const char   *_title;
	bool	     _interactive_prompting;
	DEBUG_PARAMS *_dbparams;	/*PROMPTED*/
        char         bname[300];
};
typedef struct _INIT_DATA INIT_DATA;
#define init_data(init) ((INIT_DATA *)(init))
#define	interactive_prompting(init) init_data(init)->_interactive_prompting
#define	dbparams(init) init_data(init)->_dbparams
#define title(init) init_data(init)->_title


#include "fnamedebug.h"

#include "uprotos.h"

#endif /* !defined(_CDECS_H) */
