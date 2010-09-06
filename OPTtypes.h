/*--------------------------------------------------------------------------*/
/*--------------------------- File OPTtypes.h ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Standard types and constants definitions: by changing the
 * definitions here, the dimension/precision of all the numbers
 * used within the programs can be customized.
 *
 * Also, small classes are provided for:
 * - reading the time of a code;
 * - generating random numbers.
 *
 * The classes can be adapted to different environments setting a
 * compile-time switch in this file.
 *
 * Additionally, a function is provided for safely reading parameters
 * out of a stream.
 *
 * \version 2.61
 *
 * \date 15 - 10 - 2004
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright 1994 - 2004
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __OPTtypes
 #define __OPTtypes  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*----------------------------- MACROS -------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTTYPES_MACROS Compile-time switches in OPTtypes.h
    These macros control how the classes OPTTimers and OPTrand are
    implemented; choose the appropriate value for your environment,
    or program a new version if no value suits you.
    Also, namespaces can be eliminated if they create problems.
    @{ */


/*----------------------- OPT_USE_NAMESPACES -------------------------------*/

#define OPT_USE_NAMESPACES 1

/**< Setting OPT_USE_NAMESPACES == 0 should instruct all codes that use
   OPTtypes stuff to avoid using namespaces; to start with, the common
   namespace OPTtypes_di_unipi_it, that contains all the types defined
   herein, is *not* defined. */

/*---------------------------- OPT_TIMERS ----------------------------------*/

#define OPT_TIMERS 1

/**< The class OPTtimers is defined below to give an abstract interface to the
   different timing routines that are used in different platforms. This is
   needed since time-related functions are one of the less standard parts of
   the C[++] library. The value of the OPT_TIMERS constant selects among the
   different timing routines:

   - 1 = Unix, using sys/times.h

   - 2 = as 1 but uses sys/timeb.h (typical for Microsoft(TM) compilers)

   - 3 = timing routines still use sys/times.h, but return wallclock time
         rather than CPU time

   - 4 = as 3 but uses sys/timeb.h (typical for Microsoft(TM) compilers)

   - 5 = return the user time obtained with ANSI C clock() function; this
         may result in more accurate running times w.r.t. "1", "2", "3" and
	 "4", but it is limited to ~ 72 hours on systems where ints are
	 32bits.

   Any unsupported value would simply make the class to report constant
   zero as the time.

   The values 1 .. 4 rely on the constant CLK_TCK for converting between
   clock ticks and seconds; for the case where the constant is not defined by
   the compiler -- should not happen, but it does -- or it is defined in a
   wrong way, the constant is re-defined below. */

/*---------------------------- OPT_RANDOM ---------------------------------*/

#define OPT_RANDOM 0

/**< The class OPTrand is defined below to give an abstract interface to the
   different random generators that are used in different platforms. This is
   needed since random generators are one of the less standard parts of the
   C[++] library. The value of the OPT_RANDOM constant selects among the
   different timing routines:

   - 0 = an hand-made implementation of a rather good random number generator
         is used; note that this assumes that long ints >= 32 bits

   - 1 = standard rand() / srand() pair, common to all C libreries but not
         very sophisticated

   - 2 = drand48() / srand48(), common on Unix architectures and pretty good.

   Any unsupported value would simply make the functions to report constant
   zero, which is not nice but useful to quickly fix problems if you don't
   use random numbers at all. */

/*@} -----------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <float.h>

/* For
   DBL_MAX , FLT_MAX           (largest double , float )
   DSIGNIF , FSIGNIF           (number of significant bits in double , float)
   */

#include <limits.h>

/* For
   UCHAR_MAX , CHAR_MAX             (largest [unsigned] char)
   USHRT_MAX , UINT_MAX , ULONG_MAX (largest unsigned short , int , long)
   */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( OPT_RANDOM )
 #include <stdlib.h>

 /* For random routines, see OPTrand() below. */
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#include <math.h>

/* For the function double pow( double b , double e ) */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( ( OPT_TIMERS > 0 ) && ( OPT_TIMERS <= 5 ) )
 #include <time.h>

 #if( OPT_TIMERS <= 4 )
  #include <sys/types.h>

  #if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 3 ) )
   #include <sys/times.h>
  #else
   #include <sys/timeb.h>
  #endif
 #endif
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#include <iostream>

/* For istream and the >> operator, used in DfltdSfInpt(). */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#include <assert.h>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace OPTtypes_di_unipi_it
{
 /** @namespace OPTtypes_di_unipi_it
     The namespace OPTtypes_di_unipi_it is defined to hold all the data
     types, constants, classes and functions defined here. It also
     comprises the namespace std. */
#endif

 using namespace std;  // I know it's not elegant, but ...

/*--------------------------------------------------------------------------*/
/*----------------------------- NULL ---------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NULL
 #define NULL 0
#endif

/* Curiously enough, some compilers do not always define NULL; for instance,
   sometimes it is defined in stdlib.h. */

/*--------------------------------------------------------------------------*/
/*---------------------------- CLK_TCK -------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef CLK_TCK
 #define CLK_TCK 100
#endif

/* CLK_TCK is the constant used to transform (integer) time values in seconds
   for times()-based routines. Its normal value is 100. */

/*--------------------------------------------------------------------------*/
/*------------------------ DBL_EPS and FLT_EPS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup EPS DBL_EPS and FLT_EPS
    These two constants represent the machine precision of doubles and floats,
    i.e., the relative error that one should expect from performing arithmetic
   operations with doubles and floats.
   They can be automatically calculated if the number of significant bits for
   doubles and floats is known: however, not all versions of values.h return
   that information. In the unlucky case that your one is among them, you
   have to provide the numbers by yourself: the ones above are for IEEE
   standard 64- and 32-bits floating point numbers.
    @{ */

#ifdef DSIGNIF
 const double DBL_EPS = pow( 2.0 , -DSIGNIF );
#else
 const double DBL_EPS = 1.1102230246251565e-16;
#endif

#ifdef FSIGNIF
 const double FLT_EPS = pow( 2.0 , -FSIGNIF );
#else
 const double FLT_EPS = 5.96046448e-08;
#endif

/*@} -----------------------------------------------------------------------*/
/*---------------------------- TYPES ---------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTTYPES_TYPES Types in OPTtypes.h
    @{ */

/*--------------------------------------------------------------------------*/
/** @defgroup OPTTYPES_INTREAL INT_TYPE and REAL_TYPE
    Some type will have an associated macro, that must be set to INT_TYPE if
    it is short, int, long ... and to REAL_TYPE if it is float, double ...
    so that codes using that type can distinguish the two cases.
    @{ */

#define INT_TYPE  1
///< This is an "int" type

#define REAL_TYPE 0
///< This is a "float" type

/*@} -----------------------------------------------------------------------*/
/*------------------- General-purpose type definitions ---------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTTYPES_GENPUR General-purpose type definitions
    @{ */

typedef unsigned char   BOOL;            ///< booleans
typedef BOOL           *Bool_Vec;        ///< vector of booleans
typedef Bool_Vec       *Bool_Mat;        ///< matrix of booleans

typedef const BOOL     cBOOL;            ///< a read-only boolean
typedef cBOOL         *cBool_Vec;        ///< read-only array


static cBOOL TRUE = 1;                          ///< TRUE
static cBOOL FALSE = 0;                         ///< FALSE


/*--------------------------------------------------------------------------*/

typedef unsigned int    Index;           ///< Index in a vector ( >= 0 )
typedef Index          *Index_Set;       ///< set (array) of indices
typedef Index_Set      *Index_Mat;       ///< set of set of indices

typedef const Index    cIndex;           ///< a read-only Index
typedef cIndex        *cIndex_Set;       ///< read-only array

cIndex InINF = UINT_MAX;                 ///< the largest Index

/*--------------------------------------------------------------------------*/

typedef int             SIndex;         ///< A Signed Index, for when an index
                                        ///< may have a "direction"
typedef SIndex         *SIndex_Set;     ///< set (array) of s. indices
typedef SIndex_Set     *SIndex_Mat;     ///< set of set (matrix) of s. indices

typedef const SIndex   cSIndex;         ///< a read-only Signed Index
typedef cSIndex       *cSIndex_Set;     ///< read-only array

cSIndex SInINF = INT_MAX;               ///< the largest SIndex

/*--------------------------------------------------------------------------*/

typedef long int        INum;           ///< integer numbers
typedef INum           *IRow;           ///< vectors of integers
typedef IRow           *IMat;           ///< matrices (vectors of vectors)
                                        ///< of integers

typedef const INum     cINum;           ///< a read-only integer
typedef cINum         *cIRow;           ///< read-only array

cINum IntINF = LONG_MAX;                ///< the largest INum

/*--------------------------------------------------------------------------*/

typedef double          Number;         ///< "normal" floating point numbers
typedef Number         *Row;            ///< "normal" (fp) array
typedef Row            *Mat;            ///< "normal" (fp) matrix

typedef const Number   cNumber;         ///< a read-only Number
typedef cNumber       *cRow;            ///< read-only array

cNumber eM = DBL_EPS;                   ///< machine precision of Number
cNumber INF = DBL_MAX;                  ///< the largest Number

/*--------------------------------------------------------------------------*/

typedef double          HpNum;          ///< "finer" floating point numbers
typedef HpNum          *HpRow;          ///< "finer" (fp) array
typedef HpRow          *HpMat;          ///< "finer" (fp) matrix

typedef const HpNum    cHpNum;          ///< a read-only HpNum
typedef cHpNum        *cHpRow;          ///< read-only array

cHpNum HpeM = DBL_EPS;                  ///< machine precision of HpNum
cHpNum HpINF = DBL_MAX;                 ///< the largest HpNum

/*@} -----------------------------------------------------------------------*/
/*----------- Type definitions for subgradient-based algorithms ------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTTYPES_SUBGT Type definitions for subgradient-based algorithms
    @{ */

typedef double          SgNum;          ///< subgradient entries
typedef SgNum          *SgRow;          ///< a subgradient
typedef SgRow          *SgMat;          ///< a bundle (set of subgradients)

typedef const SgNum    cSgNum;          ///< a read-only subgradient entry
typedef cSgNum        *cSgRow;          ///< a read-only subgradient

const SgNum SgINF = DBL_MAX;            ///< the largest subgradient entry

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

typedef double          QuNum;          ///< numbers in Q ( G{i}{T} * G{j} )
typedef QuNum          *QuRow;          ///< row of Q
typedef QuRow          *QuMat;          ///< Q (itself)

typedef const QuNum    cQuNum;          ///< a read-only number in Q
typedef cQuNum        *cQuRow;          ///< a read-only row of Q

cQuNum QuINF = DBL_MAX;                 ///< the largest number in Q

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

typedef double          LMNum;          ///< a Lagrangean Multiplier
typedef LMNum          *LMRow;          ///< a vector of Lagrangean Multipliers
typedef LMRow          *LMMat;          ///< a matrix of Lagrangean Multipliers

typedef const LMNum    cLMNum;          ///< a read-only Lagrangean Multiplier
typedef cLMNum        *cLMRow;          ///< a read-only vector of LMs

cLMNum LMeM = DBL_EPS;                  ///< machine precision of LMNum
cLMNum LMINF = DBL_MAX;                 ///< the largest Lagrangean Multiplier

/*@} -----------------------------------------------------------------------*/

/* @} end( group( OPTTYPES_TYPES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------- OPT_TIMERS -----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTTYPES_CLASSES Classes in OPTtypes.h
    @{ */

#if( OPT_TIMERS )

/** Provides a common interface to the different timing routines that are
    available in different platforms. */

class OPTtimers {

 public:  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /// constructor of the class
  OPTtimers( void )
  {
   ReSet();
   }

  /// start the timer
  void Start( void )
  {
   if( ! ticking )
   {
    #if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 2 ) )
     times( &buff );
     t_u = buff.tms_utime;
     t_s = buff.tms_stime;
    #elif( ( OPT_TIMERS == 3 ) || ( OPT_TIMERS == 4 ) )
     t_u = times( &buff );
    #elif( OPT_TIMERS == 5 )
     t_u = clock();
    #endif

    ticking = 1;
    }
   }

  /// stop the timer
  void Stop( void )
  {
   if( ticking )
   {
    Read( u , s );
    ticking = 0;
    }
   }

  /** Return the elapsed time. If the clock is ticking, return the *total*
    time since the last Start() without stopping the clock; otherwise,
    return the total elapsed time of all the past runs of the clock since
    the last ReSet() [see below]. */

  double Read( void )
  {
   double tu = 0;
   double ts = 0;
   Read( tu , ts );
   return( tu + ts );
   }

  /// As Read( void ) but *adds* user and system time to tu and ts.
  void Read( double &tu , double &ts )
  {
   if( ticking )
   {
    #if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 2 ) )
     times( &buff );
     tu += ( double( buff.tms_utime - t_u ) ) / double( CLK_TCK );
     ts += ( double( buff.tms_stime - t_s ) ) / double( CLK_TCK );
    #elif( ( OPT_TIMERS == 3 ) || ( OPT_TIMERS == 4 ) )
     tu += ( double( times( &buff ) - t_u ) ) / double( CLK_TCK );
    #elif( OPT_TIMERS == 5 )
     tu += ( double( clock() - t_u ) ) / double( CLOCKS_PER_SEC );
    #endif
    }
   else { tu += u; ts += s; }
   }

  /// reset the timer
  void ReSet( void )
  {
   u = s = 0; ticking = 0;
   }

 private:  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 double u;      // elapsed *user* time, in seconds
 double s;      // elapsed *system* time, in seconds
 char ticking;  // 1 if the clock is ticking

 #if( ( OPT_TIMERS > 0 ) && ( OPT_TIMERS <= 5 ) )
  clock_t t_u;

  #if( OPT_TIMERS <= 4 )
   struct tms buff;

   #if( ( OPT_TIMERS == 1 ) || ( OPT_TIMERS == 3 ) )
    clock_t t_s;
   #endif
  #endif
 #endif

 };  // end( class OPTtimers );

#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ OPTrand() ---------------------------------*/
/*--------------------------------------------------------------------------*/

/** Provide a common interface to the different random generators that are
    available in different platforms. */

class OPTrand {

 public:  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  /// constructor of the class
  OPTrand( void )
  {
   #if( OPT_RANDOM == 0 )
    A[ 0 ] = -1;
    gb_fptr = A;
   #else
    OPTrand::srand( long( 1 ) );
   #endif
   }

  /** Returns a random number uniformly distributed in [0, 1).
      \note each object of class OPTrand has its own sequence, so that
      multiple OPTrand objects being used within the same program do not
      interfere with each other (as opposed to what C random routines
      would do). */

  double rand( void )
  {
   #if( OPT_RANDOM == 0 )
    long nmbr = ( *gb_fptr >= 0 ? *gb_fptr-- : gb_flip_cycle() );
    return( double( nmbr ) / double( (unsigned long)0x80000000 ) );
   #elif( OPT_RANDOM == 1 )
    ::srand( myseed );
    myseed = ::rand();
    return( double( myseed ) / double( RAND_MAX ) );
   #elif( OPT_RANDOM == 2 )
    return( erand48( myseed ) );
   #else
    return( 0 );  // random, eh?
   #endif
   }

  /// Seeds the random generator for this instance of OPTrand.
  void srand( long seed )
  {
   #if( OPT_RANDOM == 0 )
    register long prev = seed;
    register long next = 1;
    seed = prev = mod_diff( prev , 0 );
    A[ 55 ] = prev;

    for( register long i = 21 ; i ; i = ( i + 21 ) % 55 )
    {
     A[ i ] = next;
     next = mod_diff( prev , next );
     if( seed & 1 )
      seed = 0x40000000 + ( seed >> 1 );
     else
      seed >>= 1;

     next = mod_diff( next , seed );
     prev = A[ i ];
     }

    gb_flip_cycle();
    gb_flip_cycle();
    gb_flip_cycle();
    gb_flip_cycle();
    gb_flip_cycle();
   #elif( OPT_RANDOM == 1 )
    myseed = int( seed );
   #elif( OPT_RANDOM == 2 )
    long *sp = (long*)( &myseed );
    *sp = seed;                // copy higher 32 bits
    myseed[ 2 ] = 0x330E;      // initialize lower 16 bits
   #endif
   }

 private:  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( OPT_RANDOM == 0 )
  long A[ 56 ];
  long *gb_fptr;

  long mod_diff( long x , long y )
  {
   return( ( ( x ) - ( y ) ) & 0x7fffffff );
   }

  long gb_flip_cycle( void )
  {
   register long *ii, *jj;
   for( ii = &A[ 1 ] , jj = &A[ 32 ] ; jj <= &A[ 55 ] ; ii++ , jj++ )
    *ii = mod_diff( *ii , *jj );

   for( jj = &A[ 1 ] ; ii <= &A[ 55 ] ; ii++ , jj++ )
    *ii = mod_diff( *ii , *jj );

   gb_fptr = &A[ 54 ];

   return A[ 55 ];
   }
 #elif( OPT_RANDOM == 1 )
  int myseed;
 #elif( OPT_RANDOM == 2 )
  unsigned short int myseed[ 3 ];
 #endif

 };  // end( class( OPTrand ) )

/* @} end( group( OPTTYPES_CLASSES ) ) */
/*--------------------------------------------------------------------------*/
/*--------------------------- DfltdSfInpt() --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OPTTYPES_FUNCTIONS Functions in OPTtypes.h
    @{ */

/** Template function for reading parameters from a istream. The function is
   "safe" because it works also if the istream is not given, is not be long
   enough or contains erroneous things.

   Given a &istream (possibly NULL), DfltdSfInpt() attempts to read Param out
   of it, skipping any line that begins with the comment carachter (defaulted
   to '#'), any blank line and any line starting with anything that can not
   be interpreted as a `T'. If, for any reason, the read operation fails,
   then the parameter is given the default value `Dflt'. Otherwise, all the
   rest of the line up to the nearest newline ('\n') carachter is flushed.

   \note lines should not be longer than 255 carachters. */

template<class T>
inline void DfltdSfInpt( istream *iStrm , T &Param , const T Dflt ,
                         const char cmntc = '#' )
{
 static char buf[ 255 ];

 if( iStrm && ( ! ( ! (*iStrm) ) ) )
  // the "! ! stream" trick is there to force the compiler to apply the
  // stream -> bool conversion, in case it is too dumb to do it by itself
  for( char c ;; (*iStrm).getline( buf , 255 ) )
  {
   if( ! ( (*iStrm) >> ws ) )         // skip whitespace
    break;
   
   if( ! ( (*iStrm).get( c ) ) )      // read first non-whitespace
    break;

   if( c != cmntc )                   // that's not a comment
   {
    (*iStrm).seekg( -1 , ios::cur );  // backtrack
    if( ( (*iStrm) >> Param ) )       // try reading it
     (*iStrm).getline( buf , 255 );   // upon success, skip the rest of line

    return;                           // done
    }
   }

 Param = Dflt; 

 }  // end( DfltdSfInpt )

/* @} end( group( OPTTYPES_FUNCTIONS ) ) */
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace OPTtypes_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* OPTtypes.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File OPTtypes.h -------------------------------*/
/*--------------------------------------------------------------------------*/
