/*--------------------------------------------------------------------------*/
/*---------------------------- File QPBundle.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--   Implementation of the (QP) direction finding subproblem for the    --*/
/*-- Bundle class using the specialized (QP) solvers implemented in the   --*/
/*-- [B]MinQuad class (refere to [B]MinQuad.h and references therein).    --*/
/*--                                                                      --*/
/*--                            VERSION 1.00                              --*/
/*--                           14 - 10 - 2004                             --*/
/*--                                                                      --*/
/*--                   Original Idea and Implementation by:               --*/
/*--                                                                      --*/
/*--                           Antonio Frangioni                          --*/
/*--                                                                      --*/
/*--                        Operations Research Group                     --*/
/*--                       Dipartimento di Informatica                    --*/
/*--                           Universita' di Pisa                        --*/
/*--                                                                      --*/
/*--                         Copyright 1994 - 2004                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef _QPBundle
 #define _QPBundle  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*------------------------------ HV_NNVAR ----------------------------------*/

#define HV_NNVAR 1

/* If no Multipliers are constrained to be NonNegative (NN), the QPBundle
   class can inherit directly from the base MinQuad class rather than from
   the derived class BMinQuad, saving time and space; in fact, if ! HV_NNVAR
   the BMinQuad.* files need not to be compiled and linked with the code.
   Obviously, tring to solve a problem with NN Multipliers and HV_NNVAR == 0
   will result in erratic behaviour. */

/*------------------------------ G_IMPLM -----------------------------------*/

#define G_IMPLM 4

/* The two computationally expensive tasks of the (QP) algoritmh are the
   calculation of the entries of the Hessian matrix Q[][] (each time a new
   item is inserted) and the calculation of the entries of the tentative
   direction d[] (at each iteration, and possibly more than once). How these
   two tasks are accomplished is in large part controlled by the switches
   LAZY_Q in CMinQuad.h and LAZY_D in BMinQuad.h: however, another important
   issue is how the set of items (Bundle) is implemented.

   There are essentially two ways: either by columns or by rows. A third
   possibility is clearly that of having both representations. Furthermore,
   in the latter case another degree of freedom is available, that is which
   of the two forms have to be used for each task.

   Each choice has its pros and cons, depending on the instance (max dim. of
   the Bundle vs number of Multipliers), on other algorithmic options (wheter
   or not a LVG is used) and on the target architecture (caching etc.). 
   All these issues are controlled by this switch, that has the following
   four possible values:

   0  =>  the Bundle is implemented by columns only, i.e. as a set of
          NumVar-vectors each one being an item;

   1  =>  the Bundle is implemented by rows *and* columns, that is both forms
          are kept: however, the *column-wise* implementation is used for
          calculating the entries of Q[][] while the *row-wise* implementation
          is used for calculating the entries of d[];

   2  =>  as for 1, the Bundle is implemented by rows *and* columns; however,
          in this case the *row-wise* implementation is used for calculating
          the entries of Q[][] while the *column-wise* implementation is used
          for calculating the entries of d[];

   3  =>  the Bundle is implemented by rows only, i.e. as a set of
          (max Bundle dim.)-vectors each one containing one entry of all the
          items: in this case, LAZY_Q *must* be 0 [see CMinQuad.h].

   4  =>  as 3, plus memory is allocated and deallocated on the fly to keep
          the number of allocated rows to the bare minimum.

   Clearly, if only the column-wise implementation is available (G_IMPLM == 0)
   or only the row-wise implementation is available (G_IMPLM >= 3), they are
   used for both the tasks. */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Bundle.h"

#if( HV_NNVAR )
 #include "BMinQuad.h"
#else
 #include "MinQuad.h"
#endif

/*--------------------------------------------------------------------------*/
/*---------------------- REQUIRED MACRO SETTINGS ---------------------------*/
/*--------------------------------------------------------------------------*/

/* Some relevant switches are defined in CMinQuad.h and BMinQuad.h.
   *mandotory switches* (the code won't compile/run otherwise) are

   CONSTR        1 if ADD_CNST > 0                (CMinQuad.h)
   LAZY_Q        0 if G_IMPLM > 1                 (CMinQuad.h)
   SIGNAL_MBCHG  0                                (BMinQuad.h)
   SIGNAL_B2CHG  0                                (BMinQuad.h)
   SPRBL_FI      0                                (Bundle.h)

   Obviously, switches in BMinQuad.h have no influence if HV_NNVAR == 0. */

/*--------------------------------------------------------------------------*/
/*-------------------- STRONGLY ADVISED MACRO SETTINGS ---------------------*/
/*--------------------------------------------------------------------------*/

/* Although the compile-time switches are generally thought to be hortogonal
   to each other, certain combinations make no sense and should be avoided.
   Some *strongly advised* settings are

   CONSTR        0 if ADD_CNST == 0              (CMinQuad.h)
   TWOSIDED      0                               (BMinQuad.h) */

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace Bundle_di_unipi_it
{
 using namespace MinQuad_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

const Index SltSize = 50;

/* In the column-wise implementation of the matrix, not all the columns are
   allocated from the beginning: rather, they are allocated in "slots" of the
   above size, i.e. SltSize subgradients are allocated/deallocated in one
   blow when needed. The same is done (even though the column-wise form of
   the matrix is not available) for the memory of MinQuad. */

#if( G_IMPLM > 0 )
 #if( G_IMPLM < 4 )
  const Index GTSltSize = 50;  // same as SltSize for the rows of the Bundle
 #else
  const Index GTSltSize = 1;   // keep memory to the bare minimum
 #endif
#endif

#if( ! HV_NNVAR )
 const HpNum EpsD = 1e-6;      // value for EpsilonD()
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/

class QPBundle : public Bundle ,
                 #if( HV_NNVAR )
                  protected BMinQuad
                 #else
                  protected MinQuad
                 #endif
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods and data are the actual interface of the      --*/
/*--  class: the standard user should use these methods and data only.    --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

   QPBundle( istream *iStrm , cIndex NV , cIndex ENV = InINF );

/* Constructor of the class: see Bundle.h for the meaning of most parameters.
   However, the parameter

    MPar2    (unused by the base Bundle class) here is used to initialize
             the "break" value for the pricing in the base class MinQuad
             [see SetPricing() in CMinQuad.h]: positive values are passed
             untouched, while any negative value is turned into HpINF. */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   void SetUC( cIndex Strt , cIndex Stp , cIndex_Set Which = NULL );

/*--------------------------------------------------------------------------*/

   void SetRHS( cSgRow tRHS );

   void GetRHS( SgRow tRHS );

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   cHpRow ReadMult( void )
   {
    return( MinQuad::ReadMult() );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   cIndex_Set ReadBase( void )
   {
    return( MinQuad::ReadBase() );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   Index ReadBDim( void )
   {
    return( MinQuad::ReadBDim() );
    }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   virtual ~QPBundle();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- INTERFACE WITH THE BUNDLE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

   SgRow GetItem( cIndex i );

   BOOL SetItem( cHpNum DFi , cHpNum Tau , HpNum &Ai , HpNum &ScPri ,
		 cIndex_Set SGBse );

/*--------------------------------------------------------------------------*/

   void RmvItem( cIndex i )
   {
    ResetBundle( i );
    }

   void RmvItems( void )
   {
    ResetBundle();
    }

/*--------------------------------------------------------------------------*/

   Index BSize( void )
   {
    return( ActBDim );
    }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( ADD_CNST )

   Index BCSize( void )
   {
    return( ActCNum );
    }

#endif

/*--------------------------------------------------------------------------*/

   BOOL IsItem( cIndex i )
   {
    return( IsThere( i ) );
    }

/*--------------------------------------------------------------------------*/

   void SetVars( void );

/*--------------------------------------------------------------------------*/

   char *GetNN( void )
   {
    #if( HV_NNVAR )
     return( GetVars() );
    #else
     return( NULL );
    #endif
    }

/*--------------------------------------------------------------------------*/

   void AddVar( cIndex i , SgRow NwVar , cLMNum lb );

/*--------------------------------------------------------------------------*/

   void RemoveVar( cIndex i );

/*--------------------------------------------------------------------------*/

   void SubstVar( cIndex i );

/*--------------------------------------------------------------------------*/

   BStatus SolveSubP( cHpNum tt , HpNum &dNrm , HpNum &Sgm );

/*--------------------------------------------------------------------------*/

   cLMRow Price( void );

/*--------------------------------------------------------------------------*/

   HpNum CheckAlfa( cBOOL All = FALSE );

/*--------------------------------------------------------------------------*/

   void ChgAlfa( cIndex i , cHpNum DeltaAlfai );

   void ChgAlfa( cHpRow DeltaAlfa , BOOL Mde );

/*--------------------------------------------------------------------------*/

   void ChangeCurrPoint( cLMRow DLambda , cHpNum DFi );

   void ChangeCurrPoint( cHpNum Tau , cHpNum DFi );

/*--------------------------------------------------------------------------*/

   cHpRow ReadAlfa( void );

/*--------------------------------------------------------------------------*/

   void ChgRHSi( cIndex i , cSgNum DeltaRHSi );

/*--------------------------------------------------------------------------*/

   void MakeLambda1( cHpNum Tau );

/*--------------------------------------------------------------------------*/

   HpNum EpsilonD( void );

/*--------------------------------------------------------------------------*/

   void SensitAnals( HpNum &lp , HpNum &cp );

/*--------------------------------------------------------------------------*/

#if( DO_AGGR )

   void AggregateZ( cIndex where );

#endif

/*--------------------------------------------------------------------------*/
/*------------------ "REAL" PROTECTED PART OF THE CLASS --------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The standard user should not care about the following part: users   --*/
/*--  who need to extend the code by deriving a new class may use these   --*/
/*--  methods and data structures. It is *dangerous* to *modify* the      --*/
/*--  data structures, while it safe to read them.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

#if( LAZY_Q )

   QuNum GiTGj( cIndex i , cIndex j );

#else

   void GiTG( cIndex i , register QuRow Qi , cIndex iMax );

#endif

/* The "standard" interface with the (QP) solver, common to both the MinQuad
   and BMinQuad classes. */

/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )

   cSgRow GiTilde( cIndex i );

/* This method of the BMinQuad-specific part of the interface with the (QP)
   solver returns a pointer to the i-th row of the matrix of subgradients:
   if such a row is not available for free, GiTilde() constructs it in some
   temporary memory (the vector TSGk) and returns a pointer to it. */

#endif

/*--------------------------------------------------------------------------*/

   LMNum CalculateZ( cIndex h );

   void CalculateZ( register cIndex_Set Wh , register LMRow z );

   void CalculateZ( register LMRow z );

/* These methods are used to construct the tentative direction. If HV_NNVAR,
   one of the CalculateZ()'s (depending on LAZY_D) is called by the BMinQuad
   class at each iteration; otherwise they are called by the methods of the
   QPBundle class when the direction needs to be computed.

   If PPar2 > 0, the above methods only calculate the entries of z
   corresponding to *existing* variables. */

/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )

   LMNum GiTLB( cIndex h , register cLMRow l );

/* This method of the BMinQuad-specific part of the interface with the (QP)
   solver needs to be defined only if HV_NNVAR > 0. */

#endif

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

   //!! private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

#if( ( G_IMPLM > 0 ) && ( G_IMPLM < 3 ) )

   inline void GetNewGTRow( cIndex i );

   inline void SetGTRow( cIndex i );

   inline void PutGTRow( cIndex i );

/* GetNewGTRow() and PutGTRow() can *only* be called if PPar2 != 0, i.e.,
   Buffer[] is non-NULL. */

#endif

/*--------------------------------------------------------------------------*/

#if( G_IMPLM < 1 )

  inline SgRow TSGi( cIndex i );

#endif

/*--------------------------------------------------------------------------*/

   inline void FullZ( register LMRow z , cIndex_Set tB , cHpRow tM );

/* Construct the full Z into z, using tB and tM as the Base and Mult. */

/*--------------------------------------------------------------------------*/

   inline void MakeZ( void );

/* If necessary, construct the *existing* entries of d[]. */

/*--------------------------------------------------------------------------*/

   cLMRow CompleteZ( cIndex_Set tB , cHpRow tM );

/* If necessary, complete a partly calculated d[] to a "full" d[]; returns a
   pointer to the full d[]. */

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

  Index DimMinQuad;    // dimension of the Bundle "seen" by [B]MinQuad
  Index Insrtd;        // the item to be inserted with SetItem

  HpNum MinNewAlfa;    // the min. among the Alfa's of "new" items

  char ZCmptd;         // 0 = the (QP) has not been (correctly) solved
                       // 1 = no entries of d[] have been calculated
                       // 2 = only the "defined" entries for NN variables
                       // 3 = only the "defined" entries for *all* variables
                       // 4 = *all* the entries

  SgRow RHSp;          // the Right Hand Side of the variables

  #if( G_IMPLM < 3 )
   SgMat SubG;         // the Bundle (itself)

   #if( G_IMPLM < 1 )
    SgRow TSGk;        // a temporary
   #endif
  #endif

  #if( G_IMPLM > 0 )
   SgMat TSubG;        // Transpose of G (GT)
   Index NrAllRws;     // how many rows of GT have been allocated so far

   #if( ( G_IMPLM > 0 ) && ( G_IMPLM < 3 ) )
    SgMat Buffer;      // Temporary buffer of rows of G
    Index BFCntr;      // current size of Buffer
   #endif
  #endif

  #if( G_IMPLM >= 3 )
   SgRow tmpG1;        // the current subgradient
  #endif

  #if( ! HV_NNVAR )
   LMRow dir;          // the direction
  #endif

  #if( DO_AGGR )
   HpNum ZSigma;       // the aggregated linearization error
   HpNum ZxD;          // the scalar product between z[] and d[]

   #if( ADD_CNST )
    Index CBZDim;      // number of constraints in the last optimal Base
   #endif

   #if( ADD_CNST || ( DO_AGGR == 2 ) )
    Index_Set ZBase;   // the "base" for doing aggregation
    HpRow ZMult;       // its multipliers
    Index ZBDim;       // its size
   #endif
  #endif

/*--------------------------------------------------------------------------*/

 };  // end( class QPBundle )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace Bundle_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* QPBundle.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File QPBundle.h ----------------------------*/
/*--------------------------------------------------------------------------*/
