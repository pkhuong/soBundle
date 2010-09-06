/*--------------------------------------------------------------------------*/
/*-------------------------- File MinQuad.h --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Implementation of the TT (Tall and Thin) algorithm for solving the   --*/
/*-- Quadratic Problems arising as tentative descent direction finding    --*/
/*-- subproblem within Bundle algorithms for the (linearly constrained)   --*/
/*-- minimization of NonDifferentiable convex functions.                  --*/
/*--                                                                      --*/
/*-- The user is assumed to be familiar with the kind of problems that    --*/
/*-- are solved by this code: for a description of the algorithm and the  --*/
/*-- basic notations, refere to                                           --*/
/*--                                                                      --*/
/*--   A. Frangioni "Solving Semidefinite Quadratic Problems Within       --*/
/*--   Nonsmooth Optimization Algorithms" Computers & O.R. 23(1),         --*/
/*--   p. 1099-1118, 1996                                                 --*/
/*--                                                                      --*/
/*-- available at ftp://ftp.di.unipi.it/pub/techreports/TR-95-10.ps.Z.    --*/
/*--                                                                      --*/
/*--                            VERSION 4.10      			  --*/
/*--                	       18 - 10 - 2004			       	  --*/
/*--                                                                      --*/
/*--                   Original Idea and Implementation by:               --*/
/*--                                                                      --*/
/*--                           Antonio Frangioni                          --*/
/*--                                                                      --*/
/*--                        Operations Research Group                     --*/
/*--                       Dipartimento di Informatica                    --*/
/*--                           Universita' di Pisa                        --*/
/*--                                                                      --*/
/*--                         Copyright 1992 - 2004                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MinQuad
 #define __MinQuad  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

/*------------------------------- LOG_MQ -----------------------------------*/

#define LOG_MQ 0

/* This macro controls how (if any) MinQuad produces a log of its activities
   on a ostream object set with the method SetMQLog() [see below]:

   0  =>  no log at all (SetMQLog() is even removed from the interface);

   1  =>  "basic" log: only the ERRORs are reported;

   2  =>  as 1, plus the following indices are kept updated and *succintly*
          (in a row, tab-separated) reported when the object is destroyed:
          - the number of calls
          - the average dimension of the base (in each step)
          - the average dimension of the Bundle (in each call)
          - the total time spent by SolveQP() (if TIMERS_MQ > 0);

   3  =>  as 2, but performances reports are more verbose and the "Faults" are
          also reported: Faults can happen in the normal run of the algorithm,
          but might also depend on erroneous settings of the parameters;

   4  =>  a detailed step-by-step log the algorithm is displayed;

   5  =>  as 4, plus some further debug informations are printed.

   6  =>  as 5, plus the data structures are checked for correctness at each
          iteration (slows down the code considerably). */

/*----------------------------- TIMERS_MQ ----------------------------------*/

#define TIMERS_MQ 0

/* If TIMERS_MQ > 0, then timing of the code is done: the specific type of
   timing routines used is decided in OPTtypes.h */

/*------------------------------ VARCOEFF ----------------------------------*/

#define VARCOEFF 0

/* The typical (QP) that has to be solved has the linear equality constraint

     Sum{ i } Mult[ i ] == 1.

   However, by setting VARCOEFF > 0 it is turned into the more general

     Sum{ i } cf[ i ] * Mult[ i ] == 1

   where cf[ i ] are user-provided (nonnegative) coefficients, that can be
   changed between two subsequent iterations. */

/*--------------------------- VASTE_MEM ------------------------------------*/

#define VASTE_MEM 2

/* This macro controls the handling of memory within an object of class
   MinQuad. There are three possible values:

   0 = the minimum amount of memory is allocated, and no calls to new() and
       delete() are done within the normal execution of the algorithm: some
       tasks require considerably more computational effort;

   1 = almost the minimum amount of memory is allocated, but some calls to
       new() and delete() are done (i.e. the amount of memory used can grow
       a little bit during a call to a method): this speeds up calculations;

   2 = more memory is used, but calculations are even faster: this option can
       be convenient for small problems and/or if the C++ RTS is slow in
       allocating and deallocating memory. */

/*------------------------------ LAZY_Q ------------------------------------*/

#define LAZY_Q 0

/* If LAZY_Q > 0, Q[ i , j ] is calculated, with a call to GiTGj( i , j )
   [see below], only the first time - if any - that it is required. This may
   save the calculation of some scalar products, but slows down the access to
   Q[]. If LAZY_Q == 0 instead, all the scalar products corresponding to item
   i are required in one blow with a call to GiTG( i ... ) [see below] when i
   is inserted in the Bundle [see Add[SubGrad/Constr]() below]. */

/*------------------------------- EXACT ------------------------------------*/

#define EXACT 1

/* If EXACT > 0, the scalar products x * y for x and y in {z1, z2, Lh, Lj} are
   calculated from scratch at each iteration rather than kept updated in O(1):
   this is costly, but it avoids accumulation of rounding errors that may lead
   to faults (and therefore to costly refactorizations) or even to errors. */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OPTtypes.h"

#if( LOG_MQ )
 #include <iostream>
#endif

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace MinQuad_di_unipi_it
{
 using namespace OPTtypes_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS MinQuad -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- GENERAL NOTES -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- The code works with "virtual" items in the bundle (subgradients or   --*/
/*-- constraints): it never deals directly with any implementation of     --*/
/*-- items (such as the "naive" one, a vector of SgNum's), but it allows  --*/
/*-- the user to choose it. In fact, the code does not even need to know  --*/
/*-- the dimension of space, that can be dynamically changed [see         --*/
/*-- [Add/Cut]SGSpaceDim() below]. All the operations on the bundle are   --*/
/*-- performed by referring to "names" [see Add[SubGrad/Constr]() and     --*/
/*-- ResteBundle() below].                                                --*/
/*--                                                                      --*/
/*-- MinQuad is an *abstract class*, since it has a *pure virtual method* --*/
/*-- [see the "pure virtual methods" part of the protected interface      --*/
/*-- below]: instances of MinQuad cannot be built, and a derived class    --*/
/*-- has to be defined where that method is actually implemented.         --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

class MinQuad {

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
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   MinQuad( void );

/* Constructor: initialize the object, but does not allocate memory. This is
   done is SetMaxDim() [see below], because only there inside the size of the
   various data structures is known. */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

   virtual void SetMaxDim( Index m , Index n , Index SDim = 0 );

/* This method must be called prior to any other method in the class (but the
   constructor) to tell the MinQuad object the size of the problem to be
   solved. The size is specified by the three numbers

     m    = max. dimension of the Bundle;
     n    = max. ACTIVE dimension of the Bundle (n <= m);
     SDim = dim. (number of entries) of the items [optional].

   This method can be called more than once, but each call causes a complete
   reset of the internal state of the object, and requires quite a lot of
   memory allocation/deallocation. Calling the method with m == 0 makes all
   the memory to be deallocated, and the object to sit down quietly in its
   corner, waiting for either a new call to SetMaxDim() or its destruction.

   Note on memory requirements: m, n and SDim controls the way in which memory
   is allocated by the object.

   - if SDim is != 0, it is taken as the dimension of the items: no more than
     SDim + 3 elements of matrix L are allocated, and therefore their max.
     length is SDim + 3; SDim has no impact if SDim + 3 >= m, in fact we
     consider SDim = min( SDim + 3 , m );

   - m cannot be changed: this is the actual max. dimension of the Bundle,
     and 4 vectors of dimension m (+ 5 of dim. SDim) are allocated until the
     object is destroyed;

   - n can be changed (but must always be <= m), see ChangeCrrBDim() below:
     at any time, no more than (the current value of) n items can be present,
     since only the memory for those items is allocated.

   If VASTE_MEM > 1, then O( l * n + n * n / 2 ) memory is allocated, 
   otherwise O( n * n ) memory is allocated, where l = min( SDim + 3 , m ). */

/*--------------------------------------------------------------------------*/

   void SetCrrBDim( cIndex Newn );

/* Sets the max. *active* dimension of the Bundle, that must always be <= the
   max. *absolute* dimension. The max. active dimension can be both increased
   or decreased: if it is decreased and there are items with "name" >= Newn
   in the Bundle, they are simply deleted.
   Since this procedure has no other cost than the allocation/deallocation of
   the memory, the caller is responsible for a proper memory management in
   case of memory shortage. */

/*--------------------------------------------------------------------------*/

   inline void SetLB( HpNum LB = - HpINF );

/* Sets a lower bound on the value of the variable v in the QP (the predicted
   value of the Cutting Plane model). This lower bound is automatically
   updated when the current point is changed [see MoveAlongD() below]. Giving
   - HpINF as the lower bound means that no such bound is available.

   A finite LB is equivalent to inserting in the bundle an all-0 subgradient
   with linearization error equals to - LB (LB must be negative), but it is
   directly handled by the class without requiring to waste space. When a LB
   is "in", the multipliers corresponding to subgradients returned by
   ReadMult() [see below] may sum to a quantity *strictly smaller* than 1:
   1 - Sum{ multipliers } is the optimal multiplier that the all-0 subgradient
   would get if it would be actually inserted in the bundle, and it is
   returned by ReadLBMult() [see below]. */

/*--------------------------------------------------------------------------*/

   void SetEpsilonR( HpNum NeweR = 0 );

/* In the code, a number eR is used as a measure of the relative precision
   required from the algorithm: a number x is considered to be zero <=>

     ABS( x ) <= eR * Lmu * size( x ).

   In some cases, slightly different formulae give better results, and they
   are used instead.

   eR is automatically changed (either increased or decreased) by SolveQP()
   to try to face badly conditioned problems: this method allows to "manually"
   set its current value. Passing 0 (clearly, an impossible value) resets the
   value of eR to some "default" value. */

/*--------------------------------------------------------------------------*/

   inline void SetPricing( cHpNum Brk = HpINF );

/* At each "outer" iteration, i.e. when a restricted optima is reached, the
   non-basic items are checked for insertion in the base. The MinQuad class
   implements a "Steepest - Reverse Entering Order with Early Break" pricing,
   that is:

   - in principle, all non-basic items are checked and the item with the *most
     negative* directional derivative (the mostly violated dual constraint),
     if any, is put in the base: (unlikely) ties are break since the items are
     (essentially) checked in Reverse Entering Order, i.e. items that have
     entered first are checked last;

   - however, pricing is interrupted as soon as an item with derivative <
     - Brk * BDim * | f( < current point > ) | is found (the BDim factor is
     because it is assumed that no multiplier will be larger than 1 / BDim,
     hence Brk * BDim is an upper bound on the expected decrease in f()).

   Setting Brk == 0 (the default if SetPricing() is not called) means that
   the pricing is in fact a "pure" Reverse Entering Order: the most recently
   entered item with negative directional derivative is selected. Note that
   this avoids to calculate the directional derivative for all the items
   entered before that one, hence it is less expensive. On the other hand,
   setting Brk == HpINF means using a "pure" Steepest pricing, where just the
   item with the most negative directional derivative is selected: this *may*
   help in reducing the overall iterations count, but is more expensive, and
   therefore an intermediate value of Brk may help.

   The value of Brk can be changed at any time. */

/*--------------------------------------------------------------------------*/

   inline void SetMinFVal( cHpNum NewMinFVal );

/* In typical conditions, reaching f() = 0 means terminating: hence, the extra
   stopping condition "f() <= MinFVal" can help in avoiding useless operations
   and the numerical problems that can occur when f() is too near to the
   machine precision. This method allow to set MinFVal to an user-provided
   value: if it is never called, MinFVal is set to - HpINF. */

/*--------------------------------------------------------------------------*/

#if( LOG_MQ )

   inline void SetMQLog( ostream *log );

/* The output of the code is directed onto the ostream object pointed by log:
   if this method is never called, 'clog' (the standard error bufferized) is
   used as default. Under Unix-like environments, it can be redirected to a
   file by using ">&" from the command shell. */

#endif

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR ADDING / REMOVING / CHANGING DATA -------------*/
/*--------------------------------------------------------------------------*/

   virtual void AddSubGrad( cIndex n , cHpNum alfan );

   virtual void AddConstr( cIndex n , cHpNum alfan );

   inline BOOL IsThere( cIndex n );

   inline BOOL IsASubG( cIndex n );

   inline BOOL IsAConst( cIndex n );

   inline Index MaxItemN( void );

/* Add to the current Bundle one more subgradient / constraint (in general, a
   new "item") whose "name" in 'n'. Note that nothing is known about the item,
   except that somewhere "outside" there is a subgradient [constraint] whose
   name is 'n' and whose linearization error [Right Hand Side] is alfan.

   IsThere( n ) tells if there is an item with name 'n'. IsASubG( n ) and
   IsAConst( n ) tell if the item with "name" 'n' is respectively a
   subgradient or a constraint.

   There are no restrictions on the use of names, except that they must be
   0 <= n < CrrBDim and unique. MaxItemN() returns 1 + the maximum name of
   an item currently in the Bundle if the Bundle is nonempty, and 0 if the
   Bundle is empty; that is, MaxItemN() returns the first name for which
   calling IsThere() is useless because it can only return FALSE.

   IMPORTANT ADVICE: in constrained problem, a wrong scaling of constraints
   may have an adverse influence on the behaviour of the algorithm. In fact,
   even though the constraints C * d <= e and [ k * C ] * d <= k * e (for k a
   positive scalar) are equivalent in theory, using a vector C whose norm is
   much different from the "average" norm of the subgradients in the Bundle
   may result in a badly conditioned problem. Since error-recovering proedures
   are costly, and may cause significant losses of time/precious information,
   the caller should scale the constraints appropriatedly. */

/*--------------------------------------------------------------------------*/

   void ResetBundle( cIndex which );

/* Eliminates the item named 'which' from the Bundle: if which is in Base,
   then all the data structures (L, z1, z2...) are properly updated. */

/*--------------------------------------------------------------------------*/

   void ResetBundle( void );

/* Complete reset of the Bundle: all items are eliminated. */

/*--------------------------------------------------------------------------*/

   virtual inline void ChangeAlfa( cIndex i , cHpNum NewAlfai );

/* Changes the Alfa (linearization error of a subgradient, or RHS of a linear
   constraints) of the item with name 'i', setting it to NewAlfai. */

/*--------------------------------------------------------------------------*/

   virtual void ChangeAlfa( cHpNum DeltaAlfa );

/* Performs a "relative" change of all the vector Alfa[], i.e.

                  / DeltaAlfa   'i' is a subgradient
     Alfa[ i ] += |
                  \     0       'i' is a constraint

   This is the typical correction that ought to be done in Bundle algorithms
   when "nonexact" functions have to be handled, and can be managed faster
   than the "general" case above. */

/*--------------------------------------------------------------------------*/

   virtual void MoveAlongD( cHpNum Tau , cHpNum DeltaFi );

/* Has the same effect of ChangeAlfa() when the next current point is

     next current point = current point - Tau * last direction d

   where d is intended without the "- ti" factor (Tau is the absolute step
   taken along the convex combination of the items). The method checks for the
   case Tau == ti (of the latest call to SolveQP()) and in case exploits this
   to save some computations.

   Note that a step Tau > ti leads to a point violating the "active"
   constraints, if any, i.e. the constraints in the optimal "base" [see
   ReadBase() below]; hence, the new point is feasible if and only if the
   optimal base does not contain constraints.

   The OBTAINED increase Fi( next curr. point ) - Fi( curr. point ) has to be
   provided in DeltaFi: MoveAlongD[oft]() efficiently computes the new Alfa's
   (O(1) for each of them), but it is faster in recalculating some of the
   internal data structures.

   MoveAlongD[oft]() uses the well-known Alfa updating formula

     Alfa[ i ] = Alfa[ i ] + Tau * ( d * G_i ) + DeltaFi

   where in this case d is the latest direction (without the (-ti) factor),
   hence it needs the scalar product between d and the item 'i' to perform
   the task. Such a scalar product is available for free as a byproduct of
   the latest optimization when i was already in the Bundle at that time, but
   it is clearly *not* available if i has been inserted after the latest call
   to SolveQP().

   Therefore, in order to make MoveAlongD[oft]() work properly, one needs

   - either that no new subgradients [constraints] be added between the end of
     SolveQP() and the call to MoveAlongD[oft]();

   - or that SetGTd( i , .. ) [see below] is invoked after
     Add[Subgrad/Constr]( i , .. ) to provide the required information;

   - that Alfa has not been changed in the meantime, i.e. that no calls to
     *any version of* ChangeAlfa() have occurred;

   - that a previous direction d exists, i.e. SolveQP() have already been
     called and was successfull. */

/*--------------------------------------------------------------------------*/

   virtual void ReadAlfa( HpRow NewAlfa );

   virtual inline cHpRow ReadAlfa( void );

/* Reads the current vector of Alfas, with the same "naming convenction" of
   ChangeAlfa(), either by writing it in a vector or by returning a read-only
   pointer to it. */

/*--------------------------------------------------------------------------*/

   inline void SetGTd( cIndex i , cHpNum GTdi );

/* Used to provide the scalar product between the latest direction d (without
   the (-t) factor) and the newely inserted item 'i' prior to a call to
   MoveAlongD() [see above]. */

/*--------------------------------------------------------------------------*/

   void AddSGSpaceDim( cSgRow NewDim );

   void AddSGSpaceDim( cSgRow NewDim , cHpNum lbh );

/* Increases the dimension of the space of items, implicitly "adding" one
   entry to all the items (i.e. a dual variable of the problem).

   - NewDim[ i ] is supposed to be G[ i ][ h ], i.e. the "missing" component
     of the item whose "name" is 'i': the matrix Q[][] is updated for all the
     already calculated entries (everywhere if ! LAZY_Q);

   - The current base is properly updated: this can at most decrease (by one)
     Dependent, hence it cannot create problems.

   The version with the lbh extra parameter is provided for use of MinQuad as
   "core solver" for the box constrained case: it handles the case when a
   box-constrained active variable have to be *deleted* from the Base, where
   lbh is its lower bound (or - its upper bound). */

/*--------------------------------------------------------------------------*/

   void CutSGSpaceDim( cSgRow OldDim );

   void CutSGSpaceDim( cSgRow OldDim , cHpNum lbh );

/* Decreases the dimension of the space of items, "cutting away" one of the
   entries of the items (i.e. a dual variable of the problem).

   - OldDim[ i ] is supposed to be G[ i ][ h ], i.e. the "undesidered"
     component of the item whose "name" is 'i': the matrix Q[][] is updated
     for all the already calculated entries (everywhere if ! LAZY_Q);

   - The current base is properly updated: this can increase (by one)
     Dependent, so that if Dependent is == 2 it may become == 3; in this
     case, the last item in base is discarded.

   The version with the lbh extra parameter is provided for use of MinQuad as
   "core solver" for the box constrained case: it handles the case when a
   box-constrained active variable have to be *inserted* in the Base, where
   lbh is its lower bound (or - its upper bound). */

/*--------------------------------------------------------------------------*/

   void ChangeQ( cBOOL KeepBase = TRUE );

/* The current matrix Q is substituted with an entirely new Q, whose elements
   are requested with calls to GiTG[j]() [see below]. This can be a cheaper
   alternative to the use of [Add/Cut]SGSpaceDim() [see above] in the case of
   many changes to the Q matrix, or be simply used to insert an entirely
   uncorrelated problem.
   If KeepBase = TRUE (the default), the latest optimal base is retained and
   the data structures are rebuilt in order to allow a "warm start" of the
   method. This can clearly be wasteful in certain circumnstances, and it is
   avoided if KeepBase = FALSE.

   This method can also be used to recover from failures of SolveQP() due to
   accumulation of numerical errors in *the matrix Q*. These errors may
   happen if QuNum is a floating point type and [Add/Cut]SGSpaceDim() are
   used to update Q many times. MinQuad implements a sophisticated machinery
   for handling this problem in the lower trapezoidal factor L of Q, but Q is
   assumed to be "exact": one possibility for reacting to a failure of
   SolveQP() is therefore to call ChangeQ(), for having Q "refreshed". */

/*--------------------------------------------------------------------------*/

#if( VARCOEFF )

   void ChangeCoefficient( cIndex i , cHpNum ci );

/* If VARCOEFF == 1, the primal linear equality constraint is 

     Sum{ i } cf[ i ] * Mult[ i ] == 1

   where the cf[ i ] are user-provided (nonnegative) coefficients: they are
   defaulted to 1, but they can be changed with this method.

   The coefficient of the primal variable relative to the item whose "name"
   is 'i' is set to ci: note that ci == 0 ==> 'i' is a constraint, and it is
   actually possible to change a subgradient into a constraint and vice-versa
   by setting ci == 0 (> 0) with this method. */

#endif

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

   enum MQError { kOK =  0 ,

                  kNegativeDelta ,
                  kIncreasingStep ,
                  kNegVarInBase ,
                  kInvalidV ,
                  kLoop ,

                  kQPPrimUnbndd ,

                  kFatal
                  };

   virtual MQError SolveQP( HpNum ti );

/* Solves the problem with the current data, i.e. the current Bundle, Alfa
   and ti. The second level of dynamic tolerances adjustment is implemented
   here inside: the relative precision parameter eR [see EpsilonR() above]
   can be changed to accomplish with increase/decrease "requests".

   The method returns:

   kOK    = Normal termination: a finite optimal solution has been found;

   kQPPrimUnbndd = the problem is primal unbounded: this can happen only if
                   there are constraints with *negative* right hand side in
		   the Bundle - in this case, the solution returned by
		   ReadMult() [see below] is an *infinite direction* for the
		   original Lagrangean problem;

   kFatal = even adjusting eR was not sufficient to handle the problems,
            i.e. either eR is increased/decreased out of the bounds or two
            contrasting requests (increase eR and then decrease it or
            vice-versa) have been detected. */

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR READING RESULTS ---------------------*/
/*--------------------------------------------------------------------------*/

   virtual inline HpNum ReaddNorm( void );

   virtual inline HpNum ReadSigma( cBOOL IncldCnst = TRUE );

   inline HpNum Readv( void );

/* After (correct) termination of SolveQP(), the optimal value of the (dual
   and primal) objective function can be read with these methods, as

     ReaddNorm() = Mult{T} * Q * Mult = || d ||^2

     ReadSigma() = Alfa * Mult = Sigma

   so that f( Mult ) = (1/2) * || d ||^2 + Sigma / ti, ti being the one of the
   latest call to SolveQP(). ReadSigma() returns HpINF if SolveQP() returns
   kQPPrimUnbndd.

   If IncldCnst == FALSE, the value of Sigma "restricted" to only those items
   that are actually subgradients - i.e. excluding the constraints - is
   returned.

   Readv() returns the predicted decrease of the function in the tentative
   point, i.e. the decrease in the Cutting Plane model. */

/*--------------------------------------------------------------------------*/

   inline cHpRow ReadMult( void );

   inline cIndex_Set ReadBase( void );

   inline Index ReadBDim( void );

   inline BOOL IsInBase( cIndex n );

   inline Index ReadCBDim( void );

/* After (correct) termination of SolveQP(), pointers to read-only vectors
   describing the the optimal primal solution can be read with these methods.

   Only the non-zero entries are given: the number of such entries is returned
   by ReadBDim(), and ReadMult()[ i ] contains the multiplier relative to the
   item whose "name" is ReadBase()[ i ] ( i = 0 .. ReadBDim() - 1 ). Hence,
   the dual optimal solution is

    d = - ti * Sum{ i = 0 .. ReadBDim() - 1 }
                    G_( ReadBase()[ i ] ) * ReadMult()[ i ]
   
   where G_h is the item whose "name" is 'h'. Note that the vector pointed by
   ReadBase() is InINF-terminated, i.e. ReadBase()[ ReadBDim() ] == InINF.

   It is required that no ResetBundle( i ) (where i is in OptBase[]) calls are
   done between the latest call to SolveQP() and this call to ReadMult().

   If SolveQP() returns kQPPrimUnbndd, then what is returned by the above
   methods is rather an *infinite direction* for the problem, i.e. a vector
   v such that x + beta * v is feasible for any feasible x and any beta > 0.
   Clearly, v must be zero on all the entries relative to subgradients: for
   technical reasons, however, there can subgradients in ReadBase()[], hence
   in this case the returned multipliers can be == 0.

   IsInBase( n ) returns TRUE if the "name" 'n' is somewhere in OptBase[].

   ReadCBDim() returns the number of basic items that are constraints: note
   that ReadCBDim() <= ReadBDim(), and the inequality holds strictly (i.e. at
   least one subgradient must be in Base) whenever there are subgradients in
   the Bundle and the solution is finite. */


/*--------------------------------------------------------------------------*/

   inline HpNum ReadLBMult( void );

/* Returns the optimal multiplier that the all-0 subgradient corresponding to
   the Lower Bound [see SetLB() above] would get if it would be inserted in
   the bundle. */

/*--------------------------------------------------------------------------*/

   void ReadGTd( HpRow g );

   cHpRow ReadGTd( void );

/* After (correct) termination of SolveQP(), writes in the caller-supplied
   vector g (or returns a read-only pointer vector to) the vector GTd, i.e.

     g[ i ] = d * G_i       where as usual G_i = item whose "name" is i

   and d[] is the same as in the comments to ReadMult() [see above], but
   without the "- t" factor. */

/*--------------------------------------------------------------------------*/

   virtual void SensitAnals1( HpNum &v1 , HpNum &v2 , HpNum &v3 );

/* Let x( t ) be the optimal primal solution of the problem as a function of
   the parameter t, and let v( t ), d( t ) be the corresponding dual
   solutions: there exists an interval [tMin, tMax] [see SensitAnals3()] in
   which the currently optimal base stays optimal.
   This function reports three numbers v1, v2 and v3 such that within that
   interval
 
     v( t )               =    v1 + ( 1 / t ) * v2

     || d( t ) ||^2       =  - v1 + ( 1 / t )^2 * v3

     Alfa * x( t )        =  - v2 - ( 1 / t ) * v3

   Outside the "stability interval", the above expressions give *lower
   approximations* of the real figures. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void SensitAnals2( HpNum v1 , HpNum v2 , HpRow x1 , HpRow x2 );

/* Given v1 and v2 reported by SensitAnals1(), calculates the two vectors
   x1 and x2 such that x( t ) = x1 + ( 1 / t ) * x2 is the primal optimal
   solution of the problem as a function of t, within the "stability interval"
   [tMin, tMax] [see SensitAnals3()]. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   void SensitAnals3( cHpNum v1 , cHpNum v2 , HpRow x1 , HpRow x2 ,
                      HpNum &tMin , HpNum &tMax );

/* Given v1 and v2 reported by SensitAnals1() and the two vectors x1 and x2
   calculated by SensitAnals2(), calculates the "stability interval"
   [tMin, tMax] in which the currently optimal base stays optimal, so that
   v( t ), x( t ) etc. as reported by SensitAnals[1/2]() are "exact". */

/*--------------------------------------------------------------------------*/

#if( TIMERS_MQ )

   inline void TimeMQ( double &tu , double &ts );

/* Returns total user and sistem time (in seconds) spent so far by SolveQP().
   */

#endif

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   inline HpNum EpsilonR( void );

/* Returns the current value of eR (actually, the returns eR * Lmu * BDim),
   see SetEpsilonR() above. */

/*--------------------------------------------------------------------------*/

   inline QuNum LowQ( cIndex i , cIndex j );

   inline QuNum LowQ( cIndex i );

/* Q[ i ][ j ]: with the "lazy" calculation method, the first time (if any)
   that Q[ i ][ j ] is requested the appropriate scalar product is calculated,
   by invoking GiTG[j]() [see below]. Note that invoking LowQ() with either
   'i' or 'j' being names of items currently not in the Bundle is possible,
   but clearly wrong. LowQ( i ) = LowQ( i , i ). */

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   virtual ~MinQuad();

/* Memory deallocation. Statistics (if any) are printed. */

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PURE VIRTUAL METHODS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods are a part of the protected interface of the  --*/
/*--     class, but they *must* be implemented in the derived classes.    --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

  #if( LAZY_Q )

   virtual QuNum GiTGj( cIndex i , cIndex j ) = 0;

  #else

   virtual void GiTG( cIndex i , register QuRow Qi , cIndex iMax ) = 0;

  #endif

/* If Q is calculated "lazyly", calculation of its elements is performed
   "asincronously" when (and only if) it is needed by a call to

     GiTGj( i , j ) = Sum{k} G_i[ k ] * G_j[ k ], where

   G_h = item whose "name" is h.

   In the simplest case (there exist a "naive" form of items as n-vectors of
   SgNum's, and the item whose "name" is 'i' is just the i-th row of the
   matrix G) the above scalar product can be calculated with [see OPTvect.h]

     ScalarProduct( G[ i ] , G[ j ] , n )

   However, the objects of class MinQuad doesn't need to know anything about
   the "real" implementation of the items: therefore, the implementation of
   GiTG[j]() is *not* provided, and it's due to the derived class. The
   rationale is that the items may have some structure (e.g. sparsity) that
   can be exploited to fasten the calcultation of the scalar product.

   To ease the implementation, it is guaranteed that GiTGj( i , j ) will
   always be called with i <= j.

   Conversely, if Q is not calculated "lazyly", all the elements of row Q[ i ]
   are computed together when the item of "name" 'i' is inserted in the Bundle
   (with Add[SubGrad/Constr]()): hence, in this case the derived class is
   rather required to implement a method that gives them alltogether - this
   may fasten the implementation, e.g. in the case that the matrix G is stored
   by rows (variables-wise) instead of columns (items-wise).

   The required effect of GiTG( i , Qi , iMax ) is described by the following
   pseudo-code:

   for( j = 0 ; j < iMax ; j++ )
    if( < an item of name "j" is in the Bundle > )
     Qi[ j ] = GiTGj( i , j );

   where GiTGj( i , j ) here means "whatever GiTGj() whould have returned
   if called with i , j". This method is called each time a new item is
   inserted in the Bundle, and in this case iMax will be just the largest
   "name" among all items currently in the Bundle (comprised i). However,
   the method can also be called during a "from-scratch" recomputation of
   the matrix Q [see ChangeQ() above], and in this case iMax == i. */

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
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

   inline void AlfaChanged( void );

/* For efficiency, the Alfa vector is protected and it can be modified by
   derived classes: however, when the Alfa[] of "basic" items are changed,
   AlfaChanged() must be called to warn the base class. */

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED DATA STRUCTURES --------------------------*/
/*--------------------------------------------------------------------------*/

  Index MaxBDim;        // Max. dimension of the Bundle
  Index ActBDim;        // Number of items currently in the Bundle
  Index ActCNum;        // <= ActBDim: how many of the items currently in the
                        // Bundle are constraints  
  Index NxtBIdx;        // Next available name for a new item == max. index
                        // of an item in the Bundle + 1
  MQError QPStatus;     // termination code of SolveQP()

  HpRow Mult;           // Primal variables
  Index_Set Base;       // Base (set of strictly positive variables)
  Index BDim;           // n. of items in Base
  HpNum Quad;           // Mult{T} * Q * Mult
  HpNum Lin;            // Alfa{T} * Mult

  HpRow GTd;            // Scalar product between the (optimal) direction and
                        // the items in the Bundle
  HpRow Alfa;           // Linearization errors

  HpNum PrvsTi;         // Value of ti in the latest call of CalcOptDir()
  HpNum eR;             // Relative error for "== 0" tests
  HpNum MinFVal;        // Stop if f() <= MinFVal

  HpRow tmpa;           // Temporary vectors, at least as long as the max
  HpRow tmpv;           // Bundle dimension: can be freely used by derived
                        // classes, since their content need not to be
                        // conserved between two calls to any method

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Nobody should ever look at this part: everything that is under this  --*/
/*-- advice may be changed without notice in any new release of the code. --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

  inline void GjTGBk( cIndex j , register cIndex_Set Bk , register Index Dim ,
                      register HpRow v );

/*--------------------------------------------------------------------------*/

  void RebuildBase( void );

/*--------------------------------------------------------------------------*/

  inline void CalcOptDir( cHpNum ti );

/*--------------------------------------------------------------------------*/

  inline BOOL Feasible( void );

/*--------------------------------------------------------------------------*/

  inline HpNum vPerAlfa( register cHpRow v );

/*--------------------------------------------------------------------------*/

  inline HpNum CalculateDelta( cIndex i );

/*--------------------------------------------------------------------------*/

  inline void AddSubGradToBase( Index BD , cHpNum Delta );

/*--------------------------------------------------------------------------*/

  inline void ChkNwLinInd( void );

/*--------------------------------------------------------------------------*/

  void CutSubGradFromBase( cIndex s );

/*--------------------------------------------------------------------------*/

  inline void MoveSubgradToLastPos( cIndex s );

/*--------------------------------------------------------------------------*/

  inline void UpdateLmu( cIndex s );

/*--------------------------------------------------------------------------*/

  inline void MemDealloc( void );

/*--------------------------------------------------------------------------*/

  #if( LOG_MQ > 5 )

   inline void CheckBaseValidity( void );

   inline void CheckQValidity( void );

  #endif

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

  Index_Set Order;      // Order[ 0 .. ActBDim - 1 ] contains the "names" of
                        // all the items currently in the Bundle
  char *G;              // Which position contains what in the Bundle

  QuMat Q;              // Q = G{T}G (the lower triangular part)

  HpMat L;              // L{B}*L{B}T = Q{B} (L may be trapezoidal)

  HpRow z1;             // L{-1} * e(S) ( L * z1 = e(S) )
  HpRow z2;             // L{-1} * Alfa{B}
  #if( VARCOEFF )
   HpRow cf;            // Coefficients of the linear equality constraint
  #endif

  HpNum f;              // f( x ) = (1/2) * Quad + Lin / ti
  HpNum Minf;           // Min. value of f() got in this run
  HpNum LastRo;         // Latest value of ro = - v
  HpNum z1Tz2;          // z1 * z2
  HpNum z1Norm;         // || z1 ||^2
  HpNum LhTz1;          // L[ BDim - 1 ] * z1  \ used by SensytAnalsX()
  HpNum LhTz2;          // L[ BDim - 1 ] * z2  / if Dependent == 1
  #if( ! EXACT )
   HpNum LhNorm;        // || Lh ||^2  \ set by CalculateDelta, on the base
   HpNum LjNorm;        // || Lj ||^2  / of the current value of Dependent
  #endif
  HpNum LwrBnd;         // the Lower Bound on v
  HpNum LBMult;         // the multiplier corresponding to LwrBnd

  Index Dependent;      // How many items in Base are lin. dep. (0, 1 or 2)
  Index LBIsInB;        // 1 <=> if the Lower Bound is "in Base"
  Index SgSpDim;        // Max. dimension of a base
  Index CrrBDim;        // Max. *active* dimension of the Bundle

  Index Ename;          // The name of the item that has entered last
  Index Ewher;          // The position in Base of the item Ename
  Index Eordr;          // The position in Order of the item Ename
  Index TLDim;          // The number of "taboo" items
  Index LiMax;          // Index of the max ...
  Index LiMin;          // ... and of the min. diagonal element of L
  HpNum Lmu;            // Lmu = L[ iMax ][ iMax ] / L[ iMin ][ iMin ], i.e.
                        // the condition number of L
  HpNum EpsRo;          // The "steepest" pricing is interrupted if
                        // df < - EpsRo * BDim * | f |
  Index CBDim;          // <= BDim - 1: how many of the items in Base are NOT
                        // subgradients (i.e. are constraints)

  BOOL AddedOne;        // TRUE if some new item has been added to the Bundle
                        // since the last call
  BOOL LocalOptima;     // TRUE if CalcOptDir() ended in a local optima

  #if( ! LAZY_Q )
   QuRow tempQ;         // Temporary of QuNum's
  #endif

  #if( ! VASTE_MEM )
   HpRow tmpr;          // Extra temporary for MoveSubgradToLastPos()
  #endif

  #if( LOG_MQ )
   ostream *MQLog;              // the output stream object

   #if( LOG_MQ > 1 )
    unsigned long int Calls;    // Calls counter
    unsigned long int Step;     // Step counter
    unsigned long int Success;  // Succesfull calls counter
    float SumBDim;              // Sum{ all the steps i } BDim{i}
    float SumActBDim;           // Sum{ all the calls i } ActBDim{i}
    float SumGSBase;            // Sum( all the steps i } CBDim{i}/BDim{i}
    float SumGSBundle;          // Sum( all the calls i } ActCNum{i}/ActBDim{i}
    unsigned long int Insert;   // Total n. of insertions
    unsigned long int Delete;   // Total n. of deletions
    unsigned long int SGSIncr;  // Total n. of calls to AddSGSpaceDim()
    unsigned long int SGSDecr;  // Total n. of calls to CutSGSpaceDim()
   #endif
  #endif

  #if( TIMERS_MQ )
   OPTtimers *MQt;
  #endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class MinQuad )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void MinQuad::SetLB( HpNum LB )
{
 LwrBnd = - LB;
 if( LBIsInB )
 {
  Lin = HpINF;
  if( LwrBnd == HpINF )
   LBIsInB = 0;
  }
 }

/*--------------------------------------------------------------------------*/

inline void MinQuad::SetPricing( cHpNum Brk )
{
 EpsRo = Brk;
 }

/*--------------------------------------------------------------------------*/

inline void MinQuad::SetMinFVal( cHpNum NewMinFVal )
{
 MinFVal = NewMinFVal;
 }

/*--------------------------------------------------------------------------*/

#if( LOG_MQ )

 inline void MinQuad::SetMQLog( ostream *log )
 {
  MQLog = log;

  #if( LOG_MQ > 2 )
   *MQLog << endl << "MinQuad: #Bundle = " << MaxBDim << " ~ active = "
	  << CrrBDim << "." << endl << endl; 
  #endif
  }

#endif

/*--------------------------------------------------------------------------*/

inline void MinQuad::ChangeAlfa( cIndex i , cHpNum DeltaAlfai )
{
 if( G[ i ] && DeltaAlfai )
 {
  Alfa[ i ] += DeltaAlfai;

  if( G[ i ] & 4 )       // IsInBse == 4
   z1Tz2 = Lin = HpINF;  // Alfa (=> z2) is changed
  }
 }  // end( MinQuad::ChangeAlfa( Index , HpNum ) )

/*--------------------------------------------------------------------------*/

inline BOOL MinQuad::IsThere( cIndex n )
{
 return( BOOL( G[ n ] ) );  // Void == 0
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline BOOL MinQuad::IsASubG( cIndex n )
{
 return( G[ n ] & 1 );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline BOOL MinQuad::IsAConst( cIndex n )
{
 return( G[ n ] & 2 );  // IsACnst == 2
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline Index MinQuad::MaxItemN( void )
{
 return( NxtBIdx );
 }

/*--------------------------------------------------------------------------*/

inline cHpRow MinQuad::ReadAlfa( void )
{
 return( Alfa );
 }

/*--------------------------------------------------------------------------*/

inline void MinQuad::SetGTd( cIndex i , cHpNum GTdi )
{
 GTd[ i ] = GTdi;
 }

/*--------------------------------------------------------------------------*/

inline cHpRow MinQuad::ReadGTd( void )
{
 return( GTd );
 }

/*--------------------------------------------------------------------------*/

inline HpNum MinQuad::ReaddNorm( void )
{
 return( Quad );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline HpNum MinQuad::ReadSigma( cBOOL IncldCnst )
{
 if( IncldCnst || ( ! CBDim ) )
  return( Lin );
 else
 {
  register HpNum tS = 0;
  register cHpRow tM = Mult;
  register cIndex_Set tB = Base;
  for( register Index h ; ( h = *(tB++) ) < InINF ; tM++ )
   if( ! IsAConst( h ) )
    tS += (*tM) * Alfa[ h ];

  return( tS );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline HpNum MinQuad::Readv( void )
{
 HpNum tv = - LastRo;

 if( tv < HpINF )
  tv *= PrvsTi;

 return( tv );
 }

/*--------------------------------------------------------------------------*/

inline cHpRow MinQuad::ReadMult( void )
{
 return( Mult );
 }

/*--------------------------------------------------------------------------*/

inline cIndex_Set MinQuad::ReadBase( void )
{
 return( Base );
 }

/*--------------------------------------------------------------------------*/

inline Index MinQuad::ReadBDim( void )
{
 return( BDim );
 }

/*--------------------------------------------------------------------------*/

inline Index MinQuad::ReadCBDim( void )
{
 return( CBDim );
 }

/*--------------------------------------------------------------------------*/

inline HpNum MinQuad::ReadLBMult( void )
{
 return( LBIsInB ? LBMult : 0 );
 }

/*--------------------------------------------------------------------------*/

inline BOOL MinQuad::IsInBase( cIndex n )
{
 return( G[ n ] & 4 );  // IsInBse == 4
 }

/*--------------------------------------------------------------------------*/

#if( TIMERS_MQ )

 inline void MinQuad::TimeMQ( double &tu , double &ts )
 {
  MQt->Read( tu , ts );
  }

#endif

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

inline HpNum MinQuad::EpsilonR( void )
{
 return( eR * Lmu * ( BDim ? BDim : 1 ) );
 }

/*--------------------------------------------------------------------------*/

inline QuNum MinQuad::LowQ( cIndex i , cIndex j )
{
 #if( LAZY_Q )
  if( i >= j )
  {
   if( Q[ i ][ j ] == -QuINF )
    Q[ i ][ j ] = GiTGj( j , i );

   return( Q[ i ][ j ] );
   }
  else
  {
   if( Q[ j ][ i ] == -QuINF )
    Q[ j ][ i ] = GiTGj( i , j );

   return( Q[ j ][ i ] );
   }
 #else
  return( i >= j ? Q[ i ][ j ] : Q[ j ][ i ] );
 #endif

 }  // end( LowQ( i , j ) )

/*--------------------------------------------------------------------------*/

inline QuNum MinQuad::LowQ( cIndex i )
{
 #if( LAZY_Q )
  if( Q[ i ][ i ] == -QuINF )
   Q[ i ][ i ] = GiTGj( i , i );
 #endif

 return( Q[ i ][ i ] );

 }  // end( LowQ( i ) )

/*--------------------------------------------------------------------------*/

inline void MinQuad::AlfaChanged( void )
{
 z1Tz2 = Lin = HpINF;  // Alfa (=> z2) is changed
 }

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace MinQuad_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MinQuad.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File MinQuad.h -------------------------------*/
/*--------------------------------------------------------------------------*/
