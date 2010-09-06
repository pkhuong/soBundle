/*--------------------------------------------------------------------------*/
/*---------------------------- File Bundle.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--    Minimization of a convex NonDifferentiable (possibly nonexact)    --*/
/*--       function over a polyhedral set via a "Bundle" algorithm.       --*/
/*--                                                                      --*/
/*-- The user is assumed to be familiar with the kind of problems that    --*/
/*-- are solved by this code: for a description of the algorithm and the  --*/
/*-- basic notations, refere to                                           --*/
/*--                                                                      --*/
/*-- A. Frangioni "Dual-Ascent Methods and Multicommodity Flow Problems"  --*/
/*-- TD 5/97 (Ph.D. Thesis), Dipartimento di Informatica, Universita' di  --*/
/*-- Pisa, 1997                                                           --*/
/*--                                                                      --*/
/*-- available at http://www.di.unipi.it/~frangio/thesis.html.            --*/
/*--                                                                      --*/
/*-- The code is capable of minimizing any NonDifferentiable convex       --*/
/*-- function, subject to polyhedral constraints, on a finite variables   --*/
/*-- space: however, the special case where the function                  --*/
/*--                                                                      --*/
/*--        Fi : LMNum[ NumVar ] -> HpNum                                 --*/
/*--                                                                      --*/
/*-- is a Lagrangean function is often referred to and discussed.         --*/
/*-- In this case, a "primal problem"                                     --*/
/*--                                                                      --*/
/*--  (P)   sup{ c*x : A*x [<]= b , x .in. X }                            --*/
/*--                                                                      --*/
/*-- is assumed to exist, where X is some set, and the function is        --*/
/*-- thought to be                                                        --*/
/*--                                                                      --*/
/*--   Fi( Lambda ) = sup{ ( c - Lambda*A )*x : x .in. X } - Lambda*b     --*/
/*--                                                                      --*/
/*-- (a subgradient of the function is Gi( Lambda ) = b - A*x( Lambda ),  --*/
/*--  where x( Lambda ) is any optimal solution of the subproblem).       --*/
/*--                                                                      --*/
/*-- The function inputs a vector of LMnum's (intended to be the type of  --*/
/*-- Lagrangean Multipliers, see OPTtypes.h) and outputs a value in the   --*/
/*-- user-defined type HpNum.                                             --*/
/*--                                                                      --*/
/*-- This module is parametric over the type of direction finding         --*/
/*-- subproblem employed: in the protected interface of the class there   --*/
/*-- is a collection of pure virtual methods representing the interface   --*/
/*-- with the subproblem solver, that must be properly implemented.       --*/
/*--                                                                      --*/
/*--                            VERSION 2.12                              --*/
/*--                           18 - 10 - 2004                             --*/
/*--                                                                      --*/
/*--                   Original Idea and Implementation by:               --*/
/*--                                                                      --*/
/*--                           Antonio Frangioni                          --*/
/*--                                                                      --*/
/*--                        Operations Research Group                     --*/
/*--                       Dipartimento di Informatica                    --*/
/*--                           Universita' di Pisa                        --*/
/*--                                                                      --*/
/*--                         Copyright 1993 - 2004                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef _Bundle
 #define _Bundle  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*------------------------------- LOG_BND ----------------------------------*/

#define LOG_BND 1

/* This macro controls if Bundle produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetBLog() [see below]. */

/*------------------------------- TIMERS_B ---------------------------------*/

#define TIMERS_B 0

/* If TIMERS_B > 0, then timing of the code is done: the specific type of
   timing routines used is decided in OPTtypes.h */

/*------------------------------ ADD_CNST ----------------------------------*/

#define ADD_CNST 0

/* If ADD_CNST == 1, methods for handling general linear constraints on the
   Lambda space [see AddConstraint() below] are provided, and automatic
   generation of constraints is done for cases where the feasible set is not
   know a priori [see Fi() and SetGi() / GetGi() below].

   From a "primal" viewpoint, constraints are usually needed when the set X is
   *not* compact, and therefore the Lagrangean subproblem can be *unbounded*:
   in this case, constraints correspond to *extreme rays* of X. */

/*------------------------------ DO_AGGR ----------------------------------*/

#define DO_AGGR 1

/* The normal B-strategy is to count how many consecutve iterations a given
   item has been "out of base" (has had a 0 multiplier in the solution of the
   subproblem), and to eliminate it after BPar1 such iterations. If the Bundle
   is "small" this may not be sufficient to create the space for the new item
   that has to be inserted at each step. If the Bundle is full, the code finds
   one of the items having not been "in base" for more iterations and
   eliminates it: however, if the item is "in base" at the current iteration,
   this may lead (also in practice) to cycling.

   If DO_AGGR == 0, nothing is done to prevent this from happening.

   If DO_AGGR == 1, each time a "basic" item must be eliminated another item
   (obviously, a basic one) is substituted with the "aggregated subgradient"
   Z. It is known from the theory that convergence is retained if a (QP)
   subproblem is used.

   If DO_AGGR == 2, aggregation is also used to face *failures* in the
   subproblem solver. In fact, the standard method for handling such failures
   is to discard the "removable items" in base one by one until the solver
   gets OK: clearly, there is a danger of cycling. When DO_AGGR == 2, the
   optimal Z of the previous iteration (where the solved did not failed) is
   inserted in the Bundle: in the worst case, only Z will remain (and the
   problem should be solvable with only one subgradient), but convergence is
   still ensured if a (QP) subproblem is used.

   If DO_AGGR > 0, an extra "pure virtual" protected method AggrPrimalSol() 
   [see below] is defined to give the derived class an "hook" for doing the
   corresponding "primal" aggregation. Also, a method AggregateZ() [see below]
   is added to the (protected) interface between the Bundle class and the
   subproblem solver.

   Aggregating usually slows down convergence, hence it may not be convenient
   if it is possible to keep the Bundle large enough. Also, it may not be
   possible: one example is the case when SgNum is an *integer* data type,
   and a corresponding problem is posed by the aggregation of the *primal
   solutions* X[ i ] [see AggrPrimalSol() below]. */

/*----------------------------- SPRBL_FI ----------------------------------*/

#define SPRBL_FI 0

/* In many cases, the function Fi( Lambda ) to be minimized is in fact the
   sum of k independent (nondifferentiable) functions, i.e.

       Fi( Lambda ) = Fi[ 0 ]( Lambda ) + ... + Fi[ k - 1 ]( Lambda )

   for each of which a "subgradient function" Gi[ j ]( Lambda ) is available.
   In this case, there are two alternatives:

   - either the Bundle contains the usual "aggregated" information, i.e.
       Gi( Lambda ) = Gi[ 0 ]( Lambda ) + ... + Gi[ k - 1 ]( Lambda );

   - or the Bundle contains each Gi[ j ]( Lambda ) as a separate element.

   Obviously, the subproblems to be solved in the two cases are (slightly)
   different: typically, in the second case the primal multipliers are
   partitioned into k sets, each of which sums to 1. With SPRBL_FI == 1, the
   "disaggregated" subgradients gives a more "accurate description" Fi(): the
   cost is of about a factor of k in memory, and a potentially much larger
   (=> expensive) subproblem to be solved at each step. Hence, SPRBL_FI == 0
   (to be clearly used if Fi() is "non separable"), may be more convenient in
   several applications even if Fi() is in fact separable.

   Setting SPRBL_FI == 1 changes some details in the public interface of the
   Bundle class and in the (protected) interface between Bundle and the
   subproblem solver.

   SPRBL_FI == 1 is currently *incompatible* with DO_AGGR > 0. */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OPTtypes.h"

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
namespace Bundle_di_unipi_it
{
 using namespace OPTtypes_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/

class Bundle
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

   Bundle( istream *iStrm , cIndex NV , cIndex ENV = InINF );

/* Constructor of the class.  The parameter `iStrm', if provided, istaken as
   a pointer to a istream from which the algorithmic parameters for the Bundle
   algorithm are sequentially read in the following order.  Each parameter
   must be placed at the beginning of a separate line, max 128 carachters
   long, with all the rest of the line up to the first newline carachter
   '\n' (apart from a separating whitespace) being available for comments. 
   If NULL is passed, the file ends before reaching a given parameter, or
   some parameter is in the wrong format, each non-specified parameter is
   given a default value.
   The meaning of each parameter is:

    MaxIter  Maximum number of iterations.

    BPar1    If an item has not been in the optimal base for the last BPar1
             steps, it is eliminated. If BPar1 is too small you risk to lose
             precious informations, but keeping the "bundle" small obviously
             makes the subproblem faster to solve.

    BPar2    Maximum dimension of the Bundle. Have more or less the same
             "problems" than BPar1, but if the latter is "well chosen" then
             the max bundle dimension can be kept big while the "reset
             strategy" keeps the actual number of subgradients low. A "small"
             BPar2 can affect the convergence of the algorithm, in theory as
             well as in practice, if aggregation is not performed.

    Incr     Even then t is enlarged, it is always kept <= current_t * Incr.
    Decr     As above, for reductions of t ( >= current_t * Decr ).

    m1       SS condition: if DeltaFi >= m1 * v, then a SS is done.
             The usual value for this parameter is .1: sometimes, 0 happens
             to work slightly better.
    m2       If MPar1 == 1 (the "soft" long-term t-strategy is used), then
             t decreases are inhibited whenever v < m2 * EpsU * | Fi |.
    m3       A nevely obtained subgradient is useless if Alfa >= m3 * Sigma:
             in this case, if a NS has to be done, t is decreased.
             This parameter is mostly critical: if no "long-term" t-strategy
             (see MPar1 below) is used, values < 2/3 usually make t to
             decrease quite fast to tMinor (see below), possibly making the
             algorithm to perform very short steps and therefore to converge
             very slowly. nversely, when a "long-term" t-strategy is used .9
             is a good value.

    tStar    The optimality test is dNorm * tStar + Sigma <=
    EpsLin       EpsLin * < best value of Fi() found so far >

    tMaior   Maximum, ...
    tMinor   ... minimum and ...
    tInit    ... initial value of t. 
             These parameters are potentially critical, but they are not very
             difficult to set. Typically, there is one "right" order of
             magnitude for t, that is usually the one that is guessed by the
             heuristics during most of the run, even though the starting
             value is very different. Hence, a good setting for tInit is in
             that order of magnitude, while tMaior and tMinor should be set
             respectively large and small enough to never enter into play.
             Note that a good value for tStar (i.e. one that usually ensures
             that the stopping point be actually EpsLin-optimal) is usually
             one or two orders of magnitude larger than such a tInit.

    MPar1    Select the t-strategy used:
               0 : only heuristic t-strategies;
	       1 : heuristics + "soft" long-term t-strategy;
	     > 1 : heuristics + "hard" long-term t-strategy; each time
	           "convergence" is attained, EpsU is divided by MPar1
                   (see the referred literature for details).
    MPar2    Currently unused by the base Bundle class.
    MPar3    If MPar1 != 0 (some long-term t-strategy is used), then the
             initial value of EpsU is set to EpsLin * MPar3.

    PPar1    Parameters controlling the Lagrangean Multipliers Generator:
             "price in" (discover if new Multipliers have to be added) is
	     done all the first PPar1 iterations ...
    PPar2    ... and then every PPar2 iterations; if PPar2 == 0, all the
             Lagrangean Multipliers are present from the beginning and
	     throughout all the execution (PPar1 is ignored if PPar2 == 0).
	     Note that the "price in" is forced anyway each time convergence
	     is detected.
    PPar3    A Multiplier that has been inactive for the last PPar3 pricings
             (this one included) is eliminated: note that the "price out"
	     operation is done every PPar2 iterations, so that a Multiplier
	     that is eliminated is likely to have been inactive for (about)
	     PPar2 * PPar3 iterations. For PPar3 == 1, a Multiplier is
	     eliminated in the very pricing in which it is discovered to
	     be zero (and the direction saying that it would stay zero).
	     If PPar3 == 0, Multipliers are *never* removed. PPar3 is
	     ignored if PPar2 == 0.

   NV is the *maximum* number of Lagrangean Multipliers, while ENV is the
   number of "existing" Lagrangean Multipliers: only multipliers 0 .. ENV - 1
   are significative at the beginning [but see [Add/Remove]Variable() below
   for how to change this]. ENV == InINF (the default) means that all the
   Multipliers exist.

   Each Lagrangean Multiplier can be either constrained to be NonNegative (NN)
   or UnConstrained in sign (UC): by default all the Multipliers are NN, but
   SetUC() [see below] can be used to change this. */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

#if( SPRBL_FI )

   void SetK( cIndex k );

/* If Fi() [see below] is "separable", i.e.

     Fi( Lambda ) = Fi[ 0 ]( Lambda ) + ... + Fi[ k - 1 ]( Lambda )

   this method must be used to specify the number k of components to be dealt
   with explicitly. SetK() *must* be called *exactly once*, *immediately
   after* the constructor. */

#endif

/*--------------------------------------------------------------------------*/

   virtual void SetUC( cIndex Strt , cIndex Stp , cIndex_Set Which = NULL );

/* Each Lagrangean Multiplier can be either constrained to be NonNegative (NN)
   or UnConstrained in sign (UC): this method tells that all the Multipliers

   - whose name is Strt <= i < Stp, *or*

   - whose name is in the vector Which (InINF-terminated)

   are UC. Which == NULL means "the empty set". By default, all Multiplierss
   are NN: this is equivalent to an automatic call to SetUC( 0 , 0 ).

   This method *must* be called once, *immediately after* the constructor. It
   can be called whenever the Bundle is *empty* (hence Solve() [see below] is
   not running), provided that evenctually SetLambda() [see below] is called
   to ensure that the current point is feasible w.r.t. the current setting of
   the signs of variables.

   It is "pure virtual" because derived classes implementing the subproblem
   solver may need this information also: if the method is redefined in a
   derived class, however, it is required to "throw" the method of its parent
   class as well. This is of no concern for the end user. */

/*--------------------------------------------------------------------------*/

   void SetLambda( cLMRow nLambda = NULL );

/* Sets the starting point of the Bundle algorithm; if NULL is given, or if
   the method is *not* called, the vector of all zeroes is taken. Note that,
   by default, Solve() repotimizes from the results of the previous call (if
   any) using the last curent point as the starting point; so, SetLambda()
   can be used to reset the starting point.

   No calls to SetLambda() are allowed while Solve() is running; the starting
   point can only be changed between two calls to Solve(), and this possibly
   has a nontrivial computational cost.

   Note that tLambda *must* be feasible, i.e., tLambda[ i ] must *not* be < 0
   if i is declared to be a Non-Negative (NN) variable. */

/*--------------------------------------------------------------------------*/

   virtual void SetRHS( cSgRow tRHS ) = 0;

   virtual void GetRHS( SgRow tRHS ) = 0;

/* In some cases, the subgradients/constraints *all* have the form g = u - x,
   where u is a fixed vector and x varies. In the case of Lagrangean
   optimization, u is the Right Hand Side of the relaxed constraints; more in
   general, this happens every time that the function Fi( Lambda ) to be
   maximized is the sum of a linear function u * Lambda and a nonlinear
   function (i.e. this is a special case of SPRBL_FI, see above).

   SetRHS() allows to set the Right Hand Side to the given vector: only the
   entries of tRHS corresponding to "existing" variables are significative.
   tRHS can be == NULL, meaning that the Right Hand Side is the all-0 vector:
   this is also assumed if SetRHS() is never called.
   
   GetRHS() allows to read it back.

   These methods are "pure virtual" because they are directly implemented by
   the subproblem solver: this is of no concern for the "end user", i.e. they
   must *not* be re-defined by the derived class implementing Fi() and
   SetGi() / GetGi() [see below]. */

/*--------------------------------------------------------------------------*/

   Index SetLowerBound( cHpNum LwBnd );

   void SetLBMode( const char LBMd );

/* SetLowerBound() sets a Lower Bound on the min. value of the function to be
   minimized [see Fi() below]: however, it depends on the value passed to
   SetLBMode() how this value is interpreted and used.

   0  [the default] means that such a value may not really be a LB, because
      Fi() can be unbounded below. In this case, LwBnd is interpreted as "the
      - Infinity", i.e. if Lambda is found s.t. Fi( Lambda ) < LwBnd then Fi()
      is declared unbounded below and Solve() immediately exits. This is the
      only way in which the Bundle code can detect unboundedness, unless Fi()
      returns - HpINF. If Fi() is a Lagrangean function for a problem that is
      not known to have a feasible solution, a value suitable for the task can
      be a LB (in case of a maximization problem) on the o.f. value of any
      feasible solution.

   1  means that any value > - HpINF passed to SetLowerBound() is intended to
      be a "true" LB, i.e. the function is not expected to be unbounded below.
      In this case, the code is terminated with success if
        Fi( Lambda ) - LwBnd <= EpsLin * | Fi( Lambda ) |;
      yet, unboundedness is reported if Fi( Lambda ) < LwBnd (meaning that the
      LB was not good). In case of Lagrangean optimization, a Lower Bound is
      equivalent to a *feasible primal solution*: hence, if the algorithm is
      terminated by the above control on LwBnd, then *that* feasible primal
      solution is EpsLin * | Fi( Lambda ) |-optimal. However, this is not
      correctly reflected by the vectors reported by ReadPSol(), because there
      is no subgradient "representing" that solution in the direction finding
      subproblem. It is therefore due to the calling code to recognize that
      solution as the "optimal" one instead of the "usual" convex combination.

   2  as 1, but a "floor" for the Cutting Plane model of the function (i.e. an
      "horizontal" all-0 subgradient) is explicitly constructed in the
      direction finding subproblem: this makes the CP model always lower
      bounded. This has a cost, i.e. the Bundle contains one more element, but
      it can be useful to drive the search for optimality, especially if the
      bound is tight. In the case of Lagrangean optimization, if the "floor
      subgradient" (whose "name" is returned by the method) happen to be in
      the optimal base, or in a base corresponding to an aggregation [see
      AggrPrimalSol()], then the corresponding (feasible) primal solution has
      to be used. Hence, there is no more the problem with ReadPSol() in the
      case 1. If there is no space available in the Bundle, no "floor" is
      constructed and InINF is returned. The "floor" is updated each time
      SetLowerBound() is called, and it is removed if LwBnd == - HpINF or if
      LBMd < 2 in some subsequent call to SetLBMode(). If the mode is switched
      again to 2, SetLowerBound() must be called again for the "floor" to be
      reconstructed. This can be repeated any number of times.

   Be careful with the setting of the LB if the function Fi() is not "exact":
   the Bundle has no way of knowing if and how much the value is wrong, and
   may quitely terminate in nonoptimal points if the information it has been
   given tell so. However, with LBMd == 2 the termination will go through the
   standard stopping criteria based on the outcome of the direction finding
   subproblem, hence EveryIteration() [see below] is a good place for checking
   if such a thing is happening and take your countermeasures. With LBMd < 2
   it is wiser to keep LwBnd == - HpINF and to check reliable stopping
   conditions directly inside EveryIteration(). */

/*--------------------------------------------------------------------------*/

   inline void SetStrctNZ( cBOOL SNZ = TRUE );

/* This method decides whether small negative Lagrangean Multipliers are
   allowed: if SNZ == TRUE, all the negative Multipliers among those that are
   "declared" nonnegative are turned into zeroes. If SetStrctNZ() is never
   called, or after a SetStrctNZ( FALSE ) call, small negative Multipliers
   can be found due to numerical errors. */

/*--------------------------------------------------------------------------*/

#if( LOG_BND )

   void SetBLog( ostream *log = NULL , const char lvl = 0 );

/* The output of the code is directed onto the ostream object pointed by log:
   lvl controls the "level of verbosity" of the code, as

   0  =>  no log at all (also assumed if log = NULL);

   1  =>  "basic" log: only the errors are reported;

   2  =>  a detailed step-by-step log of the algorithm is displayed;

   3, 4   unused, available to derived classes

   5  =>  as 2, plus the activity of the Lagrangean Multipliers Generator
          [see PPar* in the comments of Bundle() below] is also logged. */
#endif

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

   typedef enum { kOK = 0 ,

		  kUnbounded ,
		  kUnfeasible ,

		  kAbort ,
		  kMaxIter ,

		  kError

                  } BStatus;

   virtual BStatus Solve( void );

   inline BStatus GetBStatus( void );

/* Tries to minimize the function Fi() [see below]. If ADD_CNST > 0, the
   minimization is divided into two phases: in Phase 0 a feasible point is
   searched for, and in Phase 1 the function is minimized moving on feasible
   points only.

   Returns       if

   kOK           optimization has been carried out succesfully: a solution
                 that is "optimal" (w.r.t. the current parameters settings)
                 has been found;

   kUnbounded    there has been an error in the Fi() calculation (Fi() has
                 returned - HpINF) or the function is unbounded below: the
		 latter case can be detected only if a lower bound on the
		 min. value of the function has been provided [see above];

   kUnfeasible   the polyhedral set defined by the constraints is empty (this
                 can clearly happen only if ADD_CNST > 0): in this case, the
		 primal optimal solution [see ReadMult() ... below] is an
		 unbounded *extreme ray* for the primal problem;

   kMaxIter      the max. number of iterations has been exhausted;

   kAbort        Solve() has been stopped by EveryIteration() [see below];

   kError        There was an error in the subproblem solver during the
                 optimization, and this condition has not been corrected by
		 the elimination of items [see comments to DO_AGGR above].

   Note that, whatever the exit condition be, the current point is always
   available by calling ReadSol(), and its Fi() value by calling ReadFiVal()
   [see below]. If ADD_CNST > 0 and kMaxIter has been returned, Phase 0 may
   *not* have finished yet: hence, the current point may *not be feasible*,
   so that ReadFiVal() may return + HpINF.

   GetBStatus() returns the status of the latest call to Solve().

   Solve() is "virtual" in order to allow derived classes to implement
   different "main" strategies: this can be easily done by exploiting the
   methods FormD(), FormLambda1(), FiAndGi(), GotoLambda1() and UpdtCntrs()
   in the "non pure-virtual" protected interface of the class [see below]. */

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

   inline cLMRow ReadSol( cIndex_Set &I , Index &D );

   inline void ReadSol( register LMRow L );

   inline HpNum ReadFiVal( void );

   inline HpNum ReadBestFiVal( void );

/* The first method returns a read-only pointer to the current point: the
   format is the same of Fi() [see below], i.e. if I == NULL (=> D == NumVar)
   then the pointer points to an array whose i-th entry is Lambda[ i ],
   otherwise the only (possibly) nonzero elements are the first D entries of
   the vector, each being Lambda[ I[ i ] ].

   The second method writes the current point in L, as a full NumVar-vector.

   If Solve() has returned a kOK and the stopping parameter tStar has been
   properly set, this is an EpsLin-optimal point: in the case of Lagrangean
   optimization, an "almost" feasible primal solution also exists [see below].

   ReadFiVal() returns the Fi() value of the point returned by ReadSol(),
   while ReadBestFiVal() returns the best value of Fi() ever found by the
   algorithm: the two may not be the same due to the fact that a Bundle
   algorithm may decide not to perform a "slightly" decreasing step. Note
   that it is always possible to force ReadFiVal() == ReadBestFiVal() by just
   setting m1 == 0 in the parameters. */

/*--------------------------------------------------------------------------*/

   inline BOOL CurrentIsLast( void );

/* The point returned by ReadSol() [see above], called the `current point',
   may not be the point corresponding to the last call to Fi() [see below],
   because a number of `Null Steps' may have performed done after the last
   `Serious Step'.

   This method returns TRUE if the current point is actually the last point
   for which Fi() was called, and FALSE otherwise. This may be useful in some
   cases. */

/*--------------------------------------------------------------------------*/

   virtual cHpRow ReadMult( void ) = 0;

   virtual cIndex_Set ReadBase( void ) = 0;

   virtual Index ReadBDim( void ) = 0;

/* If Solve() [see above] returns kOK, the "aggregated subgradient"

   Sum{ i = 0 .. ReadBDim() - 1 } SubG[ ReadBase()[ i ] ] * ReadMult()[ i ]

   is ~= 0, where SubG[ i ] is the (last) subgradient with name 'i' that has
   been obtained by a call to SetGi(). If NN Multipliers are present, the
   above relation holds for the *projected* aggregated subgradient over the
   nonnegativity constraints active in the current point [see ReadSol()].

   In the case of Lagrangean optimization, a "primal" solution x[ i ] is
   associated with SubG[ i ]: the aggregated primal solution

   x* = Sum{ i = 0 .. ReadBDim() - 1 } x[ ReadBase()[ i ] ] * ReadMult()[ i ]

   is (almost) feasible w.r.t. the relaxed constraints A*x [<]= b and (almost)
   optimal for the Lagrangean subproblem w.r.t. the point returned by
   ReadSol(), i.e. (almost) optimal for the primal problem (P).

   If SPRBL_FI > 0, typically the feasible polyhedron X is the Cartesian
   product of k polihedra X[ 0 ] .. X[ k - 1 ], each one taking a subset of
   variables: in this case, the Lagrangean subproblem is separable into k
   smaller subproblems and the "disaggregated" subgradient of the j-th
   component of Fi() is Gi[ j ]( Lambda ) = b - A*x[ j ]( Lambda ), where
   x[ j ]( Lambda ) is any optimal solution of j-th the subproblem, i.e.
   belongs to the j-th sub-polyhedron X[ j ].

   If ADD_CNST > 0, some of the names in ReadBase() may be of *constraints*:
   they typically correspond to *extreme rays* of the feasible set of the
   Lagrangean subproblem. For obtaining the primal optimal solution x*, the
   above formula still works if x[ i ], for 'i' a name of a constraint,
   contains the corresponding extreme ray. Note that the multipliers relative
   to *subgradients* alone sum to 1, all the others being just nonnegative.

   Note that, if SPRBL_FI > 0 and ADD_CNST > 0, while subgradients are
   associated to this or that component, constraints are "global" to all the
   components: in fact, they are on the dual variables space and the dual
   variables cannot be partitioned.

   If Solve() returns kUnfeasible, the problem is *dual unfeasible*: this
   means that it is either *primal unbounded* or *primal empty*. In fact, in
   this case the x* obtained as above with the multipliers returned by
   ReadMult()[] is a *feasible ascent extreme ray*, that is c*(x*) > 0 ,
   A*(x*) [<]= 0 and x + beta*(x*) .in. X for each x .in. X and beta > 0.
   The existence of such an x* guarantees that if X is nonempty then (P) is
   unbounded, but it does not tell which of the two cases.

   These methods are "pure virtual" because they are directly implemented by
   the subproblem solver: this is of no concern for the "end user", i.e. these
   methods must *not* be re-defined by the derived class implementing Fi() and
   SetGi() / GetGi() [see below]. */

/*--------------------------------------------------------------------------*/

   inline HpNum ReaddNorm( void );

   inline HpNum ReadSigma( void );

/* These methods return respectively the norm of the (projected) direction,
   (in a Lagrangean problem, this is the norm of the constraint violation
   of the current continuous solution) and the value of the aggregate
   linearization error (in a Lagrangean problem, the current continuous
   solution is ReadSigma()-optimal for the Lagrangean subproblem for the
   current value of the multipliers).

   If the algorithm terminates with success, both these values should be
   "small". */

/*--------------------------------------------------------------------------*/

   inline Index FiEval( void );

   inline Index GiEval( void );

/* Return respectively the number of times that Fi() and [Set/Get]Gi() have
   been called: if called from within Fi() or [Set/Get]Gi(), the present call
   is excluded. */

/*--------------------------------------------------------------------------*/
   
   inline Index NrCalls( void );

   inline Index NrIters( void );

/* Return respectively the number of times that Solve() has been called (if
   called from within Solve(), the  present call is included) and the number
   of iterations whithin the present call. */

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

   inline Index GetNumVar( void );

/* Returns the current number of variables (Lagrangean Multipliers) of the
   function. */

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

#if( ADD_CNST )

   Index AddConstraint( cSgRow C , HpNum V , cSIndex Cntr = SInINF ,
			cIndex_Set CBse = NULL );

/* This method adds to the Bundle the linear constraint on the Lambda space

      C * Lambda <= V

   which is treated exactly like a subgradient, i.e. it occupies a "slot" in
   the Bundle and receives a "name" 'n' that is returned by AddConstraint()
   (a return value of InINF means that there was no space in the Bundle to
   acommodate for the constraint).
   If CBse == NULL, the vector C is intended to be in "dense" format, i.e.
   C[ i ] is the entry of the constraint relative to variable 'i' for i in
   0 .. GetNumVar() - 1. If CBse != NULL, it must point to the vector of
   indices of nonzero entries of the constraint (ordered in increasing sense
   and InINF-terminated), so that the "real" C is

               / C[ j ]    if exists j s.t. CBse[ j ] = i (before the inINF)
    C[ i ]  =  |
               \   0       otherwise

   Constraint names may appear in the optimal base [see ReadPSol()]: they are
   in fact usually associated to "primal extreme rays" of the Lagrangean
   subproblem, i.e. directions along which it is unbounded. Conversely, they
   *never* appear in "aggregation" bases [see AggrPrimalSol() below], since
   only subgradients are aggregated, not constraints.
   Note that constraints do not "belong" to any component [see SPRBL_FI].

   The constraint is guaranteed not to be removed for at least Cntr iterations
   (if Cntr == SInINF this is forever), unless RemoveItem( n ) [see below] is
   called to explicitly eliminate it: afterwards, the constraint is treated
   like a subgradient, i.e. it is discarded after having not been in base for
   BPar1 consecutive iterations.

   This method is not the only form for adding constraints to the Bundle: it
   can also be done by means of SetGi() / GetGi() [see below]. In fact, this
   method must *not* be called while Solve() [see above] is running, and the
   other way must be used instead.

   It *is* allowed that the current point (as returned by ReadSol()) be *not*
   feasible w.r.t. the newly inserted constraint. */

#endif

/*--------------------------------------------------------------------------*/

   void RemoveItem( cIndex Name );

   void RemoveItems( void );

/* These methods are provided for allowing the external code to remove
   subgradients and/or constraints (items) from the Bundle. The form first
   removes just the item with name 'Name' [see SetGi() below], while the
   second removes all them in one blow, comprised the "perpetual" constraints
   [see AddConstraint() above] and the "floor" of the Cutting Plane model [see
   SetLBMode() above]. If the "floor" of the Cutting Plane model is eliminated
   by a call to one of these methods, it won't be inserted again until the
   next call to SetLowerBound().

   Be careful when using RemoveItem(): deleting informations without a need
   may slow down or impair convergence. Removing *all* the subgradients from
   the Bundle while Solve() is running is forbidden: hence, RemoveItems()
   can only be called between two calls of Solve(). Also, if DO_AGGR > 0 it
   is required that no items in the currently optimal "base" [see ReadMult()
   below] be removed when the aggregation has to be performed.

   Removing *all the subgradients* from the Bundle, either with one call to
   RemoveItems() or by calling RemoveItem() for each one, essentially resets
   the whole internal state of the Bundle object, so that Solve() can be
   called again. This is a "cold restart", where evenctually the function to
   be minimized can be completely different, albeit on the same variables.
   Note that it is not necessary to explicitly remove the "floor" (the class
   already takes care of it), and that it is allowed to leave constraints in
   the Bundle - even though it might be more efficient to call RemoveItems()
   and then insert them again, depending on their number. */

/*--------------------------------------------------------------------------*/

   Index AddVariable( SgRow NwVar , cLMNum IV = 0 , cBOOL NNV = TRUE );

/* Adds a new variable to the problem: the name of the new variable will be
   the total number of variables (after the call) - 1, and it is returned by
   the method (note that the set of variable "names" is always contiguous,
   i.e. 0 .. GetNumVar() - 1). IV is taken as the initial value for the newly
   created variable (defaulted to 0). NNV == TRUE means that the new variable
   is NN, and UC otherwise.

   For each item in the Bundle with name 'i', its entry corresponding to the
   new variable must be provided in NwVar[ i ]. If a non-NULL (not all-0)
   Right Hand Side has been set [see SetRHS() above], the entry of the RHS
   corresponding to the new variable must be in NwVar[ BPar2 ] (the maximum
   Bundle dimension, a protected field inherited from BundlePar): otherwise,
   0 is obviously assumed.

   If a "floor" to the Cutting Plane model has been set [see SetLBMode()
   above], upon return NwVar[ name of the floor ] will be set == to the Right
   Hand Side of the new variable: all other entries of NwVar are untouched.
   That is, it is *not* necessary for the caller to calculate that entry.
   The name of the floor is the latest (< InINF) returned by SetLowerBound()
   [see above]. */

/*--------------------------------------------------------------------------*/

   Index RemoveVariable( cIndex i );

/* Removes the variable with name 'i'. To keep the set of variables names
   contiguous, the variable with largest "name" is "renamed" 'i', all the
   others being left untouched. The old name of the variable having become
   'i' after the call to RemoveVariable() is returned by the method: note that
   this is exactly the current number of variables returned by GetNumVar().
   If 'i' was the name of the last variable, InINF is returned (and variable
   'i' is just deleted). */

/*--------------------------------------------------------------------------*/

   void ChangeRHSi( cIndex i , cSgNum DeltaRHSi );

/* Upon termination of Solve(), this method can be used to update the Right
   Hand side of the variable with name 'i' [see SetRHS() above]. This can
   be done also if an all-0 RHS has been set: in this case, the all-0
   structure of the RHS is mantained, and DeltaRHSi is simply added to the
   i-th entry of each item (subgradient or constraint) in the Bundle.

   This operation is typically useful in the case of Lagrangean optimization
   to change the RHS of the i-th relaxed constraint: by using this method, all
   the subgradients in the Bundle are turned into subgradients of the modified
   problem, so that they can be used for a "warm" restart.
   It is easy to see that the linearization errors need not to be updated,
   while the value of Fi() in the current point changes. */

/*--------------------------------------------------------------------------*/

   void TranslateSubgradients( cHpRow Deltas );

/* Upon termination of Solve(), this method can be used to update the items
   in the Bundle by translating the subgradient with "name" i of Deltas[ i ]
   in the Fi() axe. This operation is useful, for instance, in the case of
   Lagrangean optimization when the costs in the original problem (hence, in
   the Lagrangean subproblem) are changed and the problem has to be solved
   again, exploiting the information gathered in the previous run(s). In this
   case, the subgradient with "name" i still defines a face of the function,
   but the face has to be translated of the difference between the objective
   function value of the corresponding primal solution with the old costs and
   the value with the new costs (it is easy to see that the Lagrangean
   Multipliers play no role here). In this way, the Bundle can be warm-started
   with a nonempty Cutting Plane model, hopefully leading to faster
   convergence. Note that the "floor" for the Cutting Plane model (if it has
   been constructed, see SetLBMode() above) is one possible subgradient and
   must therefore be considered among the others - after all, in the case of
   Lagrangean optimization it corresponds to one particular (feasible) primal
   solution whose cost changes as well. This is equivalent to changing the
   value of the Lower Bound with SetLowerBound() (remember that the "name" of
   the floor is reported by that method). 

   The method can be used as well for changing the RHS of constraints (if
   ADD_CNST > 0, see AddConstraint() above): Deltas[ i ] is added to the
   RHS of the constraint with "name" i. The two operations can (and should,
   if both are required) be performed with an unique call to the method:
   note that it is required to explicitly put a 0 in Deltas[ i ] if the
   RHS of the i-th constraint has not to be changed or the i-th subgradient
   has not to be translated. */

/*--------------------------------------------------------------------------*/

#if( TIMERS_B )

   inline void TimeBAll( double &t_us , double &t_ss );

   inline void TimeBFi( double &t_us , double &t_ss );

/* These methods return user and sistem time (in seconds) spent by the
   Bundle algorithm: the first returns the overall time spent by Solve(),
   while the second the time spent by Fi() and [Set/Get]Gi(). */

#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

  virtual ~Bundle();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PURE VIRTUAL METHODS --------------------------*/
/*--------------------------------------------------------------------------*/
/*--									  --*/
/*--  The following methods are a part of the protected interface of the  --*/
/*--     class, but they *must* be implemented in the derived classes.    --*/
/*--									  --*/
/*--------------------------------------------------------------------------*/

   virtual HpNum Fi( cLMRow Lam1 , cIndex_Set LamB , Index LamBd ) = 0;

   virtual void SetGi( SgRow SubG , cIndex Name ) = 0;

   virtual SIndex GetGi( HpNum &Eps1 , cIndex_Set &SGBse ) = 0;

/* These methods represent the (only) interface between the Bundle algorithm
   and the function to be minimized, and they have obviusly to be implemented
   by the derived classes - in fact, they are "pure virtual" methods.

   At each iteration, Fi() is called once: Lam1 (and LamB, LamBd - see below)
   contains the vector of Lagrangean Multipliers where the function has to be
   calculated, and the value Fi( Lam1 ) (an HpNum) of the (NonDifferentiable)
   function in Lam1 has to be returned.

   The "format" of Lam1 depends on LamB: if LamB == NULL then Lam1 is just the
   NumVar-vector of the Lagrangean Multipliers (NumVar is a protected field of
   type Index of class Bundle), while if LamB != NULL then Lam1 contains only
   the *nonzero entries* of the vector. There are LamBd of such entries
   (<= NumVar, and == NumVar if LamB == NULL), and their indices are in LamB:
   that is, the "real" Lam1 is

                  / Lam1[ j ]    if exists j < LamBd s.t. LamB[ j ] = i
    Lam1[ i ]  =  |
                  \   0          otherwise.

   Note that even some of the Lam1[ j ] may be == 0.
   To help in the implementation, LamB is ordered in increasing sense, i.e.
   LamB[ i ] < LamB[ i + 1 ] for each i < LamBd, and it is InINF-terminated,
   i.e. LamB[ LamBd ] = InINF. LamB is always == NULL if PPar2 == 0.

   In case an error occurs in the calculation of Fi(), the algorithm can be
   stopped by just returning - HpINF

   If ADD_CNST > 0, it is also possible to work even if the set Y where Fi()
   is defined is not known a priori: if Lam1 is outside Y, Fi() must return
   + HpINF and the subsequent call(s) to GetGi() must return (at least one)
   linear constraint that is valid for Y (i.e. does not cut away any point of
   it) but that cuts away Lam1 - in this way, an outer approximation of Y is
   automatically constructed (see below for more details). Remarkably, the
   algorithm can work even if no feasible point has been generated so far:
   by tightening the outer approximation of Y, a feasible point will be
   evenctually found during the run. Hence, the initial point need not to
   be feasible [see SetLambda() above].

   If SPRBL_FI > 0, however, the value of Fi( Lam1 ) is not enough: for each
   "component" Fi[ j ](), the value of Fi[ j ]( Lam1 ) need to be known. To
   avoid changes in the interface, these values must be written in the
   protected field FiLambda1k (of type HpRow), that is after termination of
   Fi() it must be FiLambda1k[ j ] = Fi[ j ]( Lam1 ).

   The two functions SetGi() and GetGi() are used to gather subgradient
   information. SetGi() provides to the Fi-solver (the user-defined derived
   class) a pointer to a NumVar-vector of SgNum where the subgradient have
   to be written, together with the "name" that the subgradient will have.
   That "name" has several possible uses: for instance, it is the one to be
   passed to RemoveItem() [see above] to eliminate that particular item.
   Also, for Lagrangean optimization, it is in general possible (at the end)
   to obtain a "primal" optimal solution as a convex combination of the primal
   solutions of the Lagrangean subproblems that are solved within Fi() [see
   ReadPSol() and AggrPrimalSol() above]: there is a primal solution for each
   of the subgradients inserted in the Bundle, and Name is the "name" that
   is given to the subgradient (hence to the corresponding primal solution).
   Note that subgradients can be deleted, and hence names can be used more
   than once. The max number of subgradients is BPar2 (protected field of the
   class, of type Index, inherited from BundlePar), and "names" are the
   integers from 0 to BPar2 - 1.

   SetGi() is called for the first time *before* Fi() is called, so that
   the Fi-solver can exploit the memory: producing the subgradient during the
   calculation of Fi() can be more efficient than constructing it afterwards.
   However, the subgradient is assumed to be ready in the vector only after
   that GetGi(), called immediately after Fi(), returns. In fact, since the
   returned vector can more in general be an epsilon-subgradient of Fi() in
   Lam1, GetGi() must return in Eps1 the epsilon (>= 0). Also, GetGi() must
   return in SGBse the information relative to the "format" of SubG. If SGBse
   == NULL, SubG is intended to be in "dense" format, i.e. SubG[ i ] is the 
   entry of the subgradient relative to variable 'i' for i in 0 .. NumVar - 1.
   If SGBse != NULL, it must point to the vector of indices of nonzero entries
   of the subgradient (ordered in increasing sense and InINF-terminated), so
   that the "real" SubG is

                  / SubG[ j ]    if exists j s.t. SGBse[ j ] = i
                  |              (before the inINF)
    SubG[ i ]  =  |
                  \   0       otherwise

   Note that *all the entries of the subgradient must be computed*, even
   those not corresponding to "active" variables (if LamB != NULL): in fact,
   those values will be used if the variable later becomes active.

   If Fi() returns + HpINF, the pair(s) (SubG, Eps1) returned by the call(s)
   of SetGi() / GetGi() - until the next call to Fi() - are rather taken as
   constraints on the space of variables [see AddConstraint() above]

      SubG * Lambda <= Eps1

   for all the rest being treated as any other subgradient: for the algorithm
   to work it is clearly required that SubG * Lam1 > Eps1, i.e. that Lam1 be
   "cut away" by the new constraint (otherwise the Bundle will cycle). Note
   that SetGi() / GetGi() are the only way of inserting constraints while
   Solve() is running, since AddConstraint() must *not* be called. However,
   it is not possible to insert subgradients and constraints together after
   the same call to Fi(). Yet, if Fi( Lam1 ) < HpINF but "inactive"
   constraints have to be added, it is possible to "fake" a Fi() == HpINF and
   insert them: at the subsequent iteration the same Lam1 will be produced.

   The returned value 'v' of GetGi() has several meanings:

   - v < 0 means that other (epsilon-)subgradients of Fi() in Lam1 are
     available and should be put in the Bundle: in this case, SetGi() and
     GetGi() are called again in sequence until a v >= 0 is returned.
     Actually, SetGi() will be called again only if the space for storing the
     other subgradients is available, i.e. it is not guaranteed that SetGi()
     will always be called whenever v < 0: however, it is guaranteed that
     GetGi() will always be called right after that SetGi() is;

   - if Fi() has returned something < + HpINF (=> SubG is a subgradient) and
     SPRBL_FI > 0, SubG belongs to Fi[ ABS( v ) - 1 ](), i.e. v tells to which
     "component" SubG referes to: note that (for obvious reasons) v = 0 is a
     possible return value only if SPRBL_FI == 0;

   - if Fi() has returned HpINF (=> SubG is a constraint), ABS( v ) tells for
     how many steps (possibly SInINF == forever) SubG must surely be kept in
     the Bundle prior to being treated like the subgradients, as the 'Cntr'
     parameter in AddConstraint() : note that constraints do not "belong" to
     any "component" of Fi().

   Lam1 and LamB do not change for all the calls of SetGi() / GetGi() until
   the next call to Fi(), hence the derived class is allowed to retain a
   copy of them and keep using the pointed vectors.

   The very first time that Fi() is evaluated, some initializations might be
   needed: this can be checked from FiEval() [see above]. It may also be
   useful to know whether or not LamB has changed since the latest call: the
   protected variable BHasChgd is set to TRUE when LamB is *not* the same as
   in the latest call to Fi(). Note that it is *never* set to FALSE: this
   must be done by the derived class, once it has done whatever it needs to
   do. If PPar2 == 0, BHasChgd is *always* FALSE; otherwise, it is TRUE at
   the very first iteration only if the initial set of Multipliers is *not*
   the all-zero vector [see SetLambda() above].

   The *last* subgradient inserted is taken as the "representative", and it
   is used for some tasks within the algorithm, such as the heuristics to
   determine the new t. */

/*--------------------------------------------------------------------------*/

#if( DO_AGGR )

   virtual void AggrPrimalSol( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm ,
			       cIndex NwNm ) = 0;

/* If the Bundle is "small", it might be necessary to "make space" for new
   items by [see DO_AGGR] substituting some of them with the "aggregate
   subgradient" Z. This method is an "hook" for derived classes: it is called
   each time the item with name NwNm is replaced by Z in the Bundle.

   From the primal point of view, this means that the corresponding element
   of the set of the primal solutions { x[ i ] } (where i is the "name" given
   by SetGi()) has to be substituted with a convex combination of the Dm
   primal solutions whose names are in NmSt, using Mlt as multipliers, i.e.

     x[ NwNm ] = Sum{ i = 0 .. Dm - 1 } Mlt[ NmSt[ i ] ] * x[ i ]

   If SPRBL_FI == 0, Sum{ i = 0 .. Dm - 1 } Mlt[ i ] = 1; otherwise, Mlt[]
   is partitioned into k subsets each of which sums to 1 [see DO_AGGR].
   Note that NwNm *can* be an element in NmSt[].
   NmSt[] is InINF-terminated, i.e. NmSt[ Dm ] == InINF.

   Obviosly, the actual method has to be implemented in the derived class
   (in fact, this is a "pure virtual" method), but this has to be done only
   if a primal solution is desired at the end. In case DO_AGGR > 0 but the
   primal aggregation need not to be done, it is necessary to define an
   AggrPrimalSol() in the derived class doing nothing. */

#endif

/*--------------------------------------------------------------------------*/
/*---------------------- HOOKS FOR DERIVED CLASSES -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods are called by the solver to give the user a   --*/
/*--  better control over the optimization process.                       --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

   typedef enum { kEINorm = 0,

		  kEIAbort ,
		  kEILoopNow ,
		  kEIContAnyway

                  } EIStatus;

   virtual EIStatus EveryIteration( void );

/* This method is an "hook" for derived classes: it is called at Every
   Iteration, between the computation of the tentative direction and the calls
   to Fi() and [Set/Get]Gi(). It can serve to various purposes, primarly
   checking extra stoping conditions or interfering with the usual stopping
   conditions of the Bundle code: however, any kind of operation can be
   performed here inside, comprised adding/removing constraints and/or
   subgradients [see RemoveItem() and AddConstraint() above], setting/removing
   the LB on the optimal value of Fi() [see SetLowerBound() above] and adding
   or removing variables from the problem [see [Add/Remove]Variable() above].
   More in general, this method can be used to merge the main cycle of the
   Bundle method within any other however complex code: the Bundle gives out
   the control at this time, and resumes its operations when EveryIteration()
   returns. The returned value influences the behaviour of the Bundle for the
   current iteration:

   kEINorm        the current iteration is continued normally;

   kEIAbort       the whole algorithm is aborted, and Solve() is immediately
                  terminated returning kAbort: this is useful for instance
                  to enforce new termination criteria;

   kEILoopNow     the current iteration is aborted, i.e. the stopping
                  conditions are *not* checked, and Fi() is *not* called: the
		  next iteration is immediately started, but the iterations
		  count is *not* increased. This is useful e.g. if something
		  has been changed in the data of the problem that advices to
		  try a new direction, like a new "active" constraint [see
		  AddConstraint() above] or a new variable [see AddVariable()
		  above] has been inserted;

   kEIContAnyway  the current iteration is continued normally but for the
                  fact that the usual stopping condition
                              vStar <= EpsLin * ABS( FiLambda )
                  is *not* checked (all the terms in the previous expression
                  are protected fields of the class and can be accessed by
                  derived classes). This is useful for instance if Fi() is
                  not calculated exactly, and termination would occur with
                  the present "inexact" information, while calculating Fi()
                  again with a better degree of accuracy would generate new
                  information that might make the code to continue.

  As an example of its possible use, Bundle::EveryIteration() (the version
  implemented by the base class) is the place where the "Hard" Long-Term
  t-strategy is implemented: if the maximum expected improvement is "small",
  t is increased and kEILoopNow is returned. Hence, classes re-implementing
  EveryIteration() that want to use the "Hard" Long-Term t-strategy must
  properly call Bundle::EveryIteration(). */

/*--------------------------------------------------------------------------*/

   virtual BOOL LooseSubgradients( void )
   {
    return( FALSE );
    }

/* In some cases, the function Fi() cannot be calculated exactly, and/or
   0-subgradients to Fi() are too difficult to obtain. There are essentially
   two possible ways of facing the problem.

   Either a "wrong" value of the function is given, together with at least one
   0-subgradient per iteration. In the case of Lagrangean optimization, this
   is what happens if the Lagrangean subproblem is not solved to optimality,
   but the suboptimal solution is treated as if it were optimal. In this case,
   a lower approximation of the "real" Fi() is being minimized instead of the
   original function, and some "incongruences" in the data may arise in the
   course of the algorithm. Typically, negative linearization errors may be
   detected, that show how the current value of the function is "wrong" and
   must be increased of (at least) some easily computable quantity - this is
   automatically done by the code. The problem with this approach occurs with
   termination: the algorithm terminates when the approximation to Fi(), not
   the real function, is (approximately) minimized. However, this can be
   detected inside EveryIteration() [see above]: if some control is possible
   on the "level of accuracy" with which Fi() is computed, termination can be
   avoided for the current iteration, increasing in the same time the accuracy
   in such a way that at the following call to Fi() new information is
   generated that allow the Bundle to continue. Note that, when using this
   strategy in the case of Lagrangean optimization, the resulting "optimal"
   Fi() value coming out of the minimimzation is potentially *not* a valid
   Lower Bound to the original problem, unless the error in the Fi()
   calculation is driven to 0 (i.e. at least the Lagrangean subproblem
   corresponding to the "optimal" point has been solved to optimality).

   Or a "correct" value of the function is computed, but no 0-subgradient is
   available. In the case of Lagrangean optimization, this corresponds to the
   case where a conveniently tight Lower Bound on the optimal value of the
   Lagrangean subproblem (evenctually, the optimal value itself) is known, but
   the best solution (hence subgradient) that can be retrieved is not optimal.
   Hence, at some iterations no 0-subgradient is available, but only epsilon-
   subgradients for some epsilon > 0. In this case, an upper approximation of
   Fi() (evenctually, Fi() itself) is minimized, so that there are no more
   problems with negative linearization errors and stopping to nonoptimal
   points; also, in the case of Lagrangean optimization, the best value of
   Fi() found during the minimization *is* a valid Lower Bound for the
   original problem. However, in this case a different problem may arise, that
   is the algorithm may cycle. This can happen if the none of the epsilon-
   subgradients found at a given iteration has epsilon small enough so that,
   when added to the Bundle, it makes the previous solution of the direction
   finding subproblem suboptimal. In this case, the same tentative point may
   be generated at the next iteration, and the algorithm cycles (assuming that
   the calculation of Fi() is a deterministic process).

   LooseSubgradients() is exactly provided for allerting derived classes about
   the possibility that cycling occurs, and giving them the chance to react.
   It is called, immediately after the last call to GetGi() returns, if

   i)  none of the nevly obtained epsilon-subgradients will make the solution
       of the direction finding subproblem change, and

   ii) the current point is not going to be changed, i.e. a Null Step is
       going to be performed.

   In this case, the algorithm *may* cycle: it won't necessarily do if t
   happens to change in the current iteration, due to the t-strategies used,
   but this is not very likely to happen. However, the process is
   deterministic, i.e. if the same t is used for two consecutive iterations
   where LooseSubgradients() is called, then the algorithm is cycling.
   This method will never be called if at least one 0-subgradient is always
   available, but possibly for the very last iteration (i.e. immediately
   before that the Bundle terminates reporting optimality).

   On a return value of TRUE, Fi() and [Set/Get]Gi() will be immediately
   called again on the same point Lam1, giving the possibility to provide a
   more accurate estimate of Fi() and/or epsilon-subgradients with smaller
   epsilons. This can be repeated any number of (but hopefully finitely many)
   times within the same iteration. However, note that *all* the subgradients
   accumulated are in principle kept, unless removed by calls to RemoveItem()
   from the derived class: it is responsibility of the derived class to
   ensure that there is space enough in the Bundle to accommodate for all
   them (or at least for the "important ones").

   Note that the Bundle class does not automatically react to a potential
   danger of cycling, i.e. if the derived class does not do nothing about it,
   the algorithm will simply cycle until the max number of iterations is
   reached. This is done for allowing the maximum flexibility in the reaction
   from the derived class. Of course, one (often the only) possible reaction
   is simply to terminate the run by returning the appropriate value at the
   next call to EveryIteration().

   Finally, note that this mechanism cannot be used for *constraints*. That
   is, when Fi() returns + HpINF then at least one of the constraints
   returned by GetGi() *must* cut away the point Lam1 on which Fi() has been
   called. This already ensures that the next trial point will be different,
   preventing the cycles. */

/*--------------------------------------------------------------------------*/

   virtual void DeletedI( cIndex i ) {}

/* To obtain a primal solution at the end of the optimization process, the
   user must keep all the primal solutions corresponding to items in the
   Bundle: this can obviously be memory consuming.
   When an item is deleted from the Bundle, its "name" is possibly re-used
   afterwards: seeing the "name" of a primal solution be used again is what
   tells the user that a particular primal solution can be discarded (but
   another one will immediately take its place).
   However, there may be a gap between the moment in which an item is
   actually discarded and the moment in which the name is re-used: it may be
   useful for the user to be timely informed that one particular primal
   solution is no longer needed.

   DeletedI( i ) is called each time the item with "name" i is deleted,
   unless this is done to make space for another item: in this case, in fact,
   the name 'i' will be immediately re-used however. */

/*--------------------------------------------------------------------------*/
/*------------------ "REAL" PROTECTED PART OF THE CLASS --------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The standard user should not care about the following part: users   --*/
/*--  who need to really extend the functionality of the Bundle class     --*/
/*--  (e.g. by implementing a different subproblem solver) may use these  --*/
/*--  methods and data structures. It is *dangerous* to *modify* the      --*/
/*--  data structures, while it safe to read them.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------- "NORMAL" PROTECTED METHODS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods can be used to implement variants of the      --*/
/*--  Solve() method, i.e. different "main" Bundle strategies.            --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

   void FormD( void );

/* When no Lagrangean Multipliers Generation is done (PPar2 == 0), FormD()
   just calls SolveSubP() once and calculates the direction d: however, it
   also implements some strategies to survive to "fatal" failures in the
   subproblem solver, typically eliminating some of the Bundle items.

   Set the protected field Result to kOK if (evenctually after some "fatal"
   failure) a tentative descent direction could be found, to kUnfeasible if
   the subproblem is dual unfeasible and to kError if this was returned by
   SolveSubP(): in the latter cases, the whole algorithm must abort.

   If Lagrangean Multipliers Generation is done (PPar2 > 0), this is where
   the corresponding strategies are implemented: in this case, SolveSubP()
   can be called more than once within the same call to FormD(), since the
   resulting direction has to be optimal w.r.t. all the current "active set"
   of Lagrangean Multipliers. */

/*--------------------------------------------------------------------------*/

   void FormLambda1( cHpNum Tau );

/* After a (succesfull) call to FormD(), sets the new tentative point Lambda1
   (a protected field of type LMRow) as Lambda1 = Lambda + ( Tau / t ) * d. */

/*--------------------------------------------------------------------------*/

   void FiAndGi( void );

/* Computes Fi( Lambda1 ), inserting the obtained items (subgradients or
   constraints) in the Bundle. */

/*--------------------------------------------------------------------------*/

   void GotoLambda1( void );

/* Move the current point to Lambda1. */

/*--------------------------------------------------------------------------*/

   void UpdtCntrs( void );

/* Updates the out-of-base counters for all items in the Bundle, evenctually
   eliminating outdated items. */

/*--------------------------------------------------------------------------*/
/*------------------ INTERFACE WITH SUBPROBLEM SOLVER ----------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  The following methods are the interface between the Bundle class    --*/
/*--  and the solver of the direction finding subproblems. These methods  --*/
/*--  are "pure virtual", since Bundle is independent from the choice of  --*/
/*--  the subproblem, whose implementation is let to a proper derived     --*/
/*--  class. However, the "final" user (who have to define Fi() and       --*/
/*--  [Set/Get]Gi()) does not need to bother with these details.          --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

/* General notes on the interaciotn between the Bundle class and the
   direction finding subproblem solver.

   The solver must alllocate the memory for the subgradients. There are at
   most BPar2 (protected field, inherited from ParBundle) subgradients or
   constraints (items), whose "names" go from 0 to BPar2 - 1.

   The Multipliers are at most MaxNumVar (protected field): in each moment,
   there are NumVar (protected field) "declared" Multipliers, with "names" in
   0 .. NumVar - 1. Of these, NNVars (protected field) are NN and the
   remaining NumVar - NNVars are UC (NNVars is calculated inside SetUC(), so
   be sure to invoke Bundle::SetUC() before counting on that value).

   If PPar2 > 0 (protected field, inherited from ParBundle), the protected
   field LamBase (Index_Set) points to the vector of the names of "active"
   Multipliers: the vector is ordered in increasing sense, InINF-terminated
   and has LamDim (protected field) elements. PPar2 == 0 <=> LamBase == NULL.
   The methods SetVars(), GetNN(), RemoveVar() and Price() [see below] are
   *only* called if PPar2 > 0: furthermore, RemoveVar() is only called if
   PPar3 (protected field, inherited from ParBundle) is *also* > 0 and
   GetNN() is only called if also ( NNVars > 0 ) && ( NNVars < NumVar ), that
   is if there is actually a doubt on the sign of a variable.

   The derived class may count on the above invariants to simplify or make
   more efficient its implementation. */

/*--------------------------------------------------------------------------*/

   virtual SgRow GetItem( cIndex i ) = 0;

   virtual BOOL SetItem( cHpNum DFi , cHpNum Tau , HpNum &Ai ,
			 HpNum &ScPri , cIndex_Set SGBse ) = 0;

/* These two methods are used to insert a new item (subgradient or constraint)
   in the Bundle. GetItem( i ) returns a pointer to a NumVar-long vector of
   SgNum's where the new item with name 'i' has to be written; i must be the
   name of an "empty" item.

   Then, SetItem() is called. The new item has been obtained in the point
   Lambda1 = Lambda + ( Tau / t ) * d, where d is the (scaled) optimal
   solution of the direction finding subproblem.
   If Tau != 0, the scalar product between the (scaled) direction (1 / t) * d
   and the item must be calculated and returned in ScPri.

   DFi is assumed to be Fi( Lambda1 ) - Fi( Lambda ): if DFi == HpINF then
   Fi( Lambda1 ) == HpINF ==> the item is a constraint, and Ai contains its
   Right Hand Side scaled w.r.t. the current point Lambda.
   Otherwise, the item is a subgradient and Ai contains its linearization
   error w.r.t. *Lambda1* (i.e., it is an Ai-subgradient in Lambda1).
   If Tau != 0, Ai must be updated to the linearization error of the
   subgradient w.r.t. *Lambda* (the current point, different from the
   tentative point Lambda1), and this is the lin. err. that must be stored
   in the subproblem; otherwise, the scalar product need not be calculated,
   and Ai is directly taken as the lin. err. of the item in the subproblem.

   The method returns TRUE if the new item changes the solution of the
   direction finding subproblem, FALSE otherwise.

   If SPRBL_FI > 0, the name of the component to which the new item belongs
   can be found in the protected field whtsG1 (of type Index). In this case,
   DFi holds the change in the value of *that component* of Fi(), that is
     Fi[ whtsG1 ]( Lambda1 ) - Fi[ whtsG1 ]( Lambda ).

   SetItem() is always called exactly once after any call to GetItem(), and
   before any subseguent call to GetItem().

   The "format" of what is written in SubG depends on SGBse: if SGBse == NULL
   then SubG is in "dense" format, while if SGBse != NULL then SubG is in
   "sparse" format [see comments to GetGi() for details: in fact, SGBse is
   exactly the pointer that has been returned by GetGi()].

   Note that *all the first NumVar entries of the vector are meaningful*, and
   the subproblem solver is assumed to retain that information in case it
   turns out to be useful later. That is, even if some variable i < NumVar
   (which is a protected field of type Index) is not defined in the direction
   finding subproblem, the corresponding entries of the items must be kept in
   case the variable is generated afterwards [see AddVar() below]. */

/*--------------------------------------------------------------------------*/

   virtual void RmvItem( cIndex i ) = 0;

   virtual void RmvItems( void ) = 0;

/* The first form removes the item with name 'i' from the Bundle; the second
   removes *all* items. Items removed with any of the two methods are "empty"
   and their names can be used again in GetItem(). Initially, all items are
   empty. */

/*--------------------------------------------------------------------------*/

   virtual Index BSize( void ) = 0;

#if( ADD_CNST )

   virtual Index BCSize( void ) = 0;

#endif

/* These two methods must return respectively the current number of items in
   the Bundle, and the number of these items that are constraints. */

/*--------------------------------------------------------------------------*/

   virtual BOOL IsItem( cIndex i ) = 0;

/* Returns TRUE if item `i' is in the Bundle. */

/*--------------------------------------------------------------------------*/

   virtual void SetVars( void ) = 0;

/* If PPar2 > 0, this method is used to change the current set of "active"
   Multipliers in the direction finding subproblem: the set of names must be
   taken from LamBase [see above]. If this method is *not* called, the set of
   variables is just *empty*.

   The "sign" of each Multiplier (whether it is NN or UC) is assumed to be
   the one set by SetUC() [see above]. For NN Multipliers, the lower bound
   on the Multiplier `LamBase[ i ]' is `- Lambda[ i ]' (Lambda is a protected
   field of type LMRow).

   This method can be called more than once: if there were a previous set of
   active Multipliers, it is deleted and LamBase becomes the set of *all and
   only* the active Multipliers. */

/*--------------------------------------------------------------------------*/

    virtual char *GetNN( void ) = 0;

 /* If PPar2 > 0, the Bundle code needs to check the "sign" of all Multipliers
    as set by SetUC() [see above]. Since the subproblem solver must keep this
    information internally anyway, it is better to avoid duplicating it.
    This method shall return a pointer v to a vector of 'char' such that the
    *first bit* of v[ i ] is 1 <=> 'i' is a NN Multiplier. */

/*--------------------------------------------------------------------------*/

   virtual void AddVar( cIndex i , SgRow NwVar , cLMNum lb ) = 0;

/* Adds the new variable with name 'i' to the direction finding subproblem: i
   must be the name of an unused variable. If NwVar == NULL, the entries of
   all the items in the Bundle corresponding to the variable i are taken from
   the i-th positions of the vectors passed with [Get/Set]Item(); otherwise,
   the entry corresponding to the item with name 'i' is taken from NwVar[ i ].

   lb is the lower bound for the new variable: lb == - LMINF means that 'i'
   is *not* constrained. However, this is only checked if NwVar != NULL (this
   is a "new" variable): otherwise, the sign of the variable stay the one set
   by SetUC() [see above].

   If a non-NULL (not all-0) Right Hand Side has been inserted, and NwVar !=
   NULL, then variable 'i' is "new" and the corresponding entry of the R.H.S.
   must be provided: it is written in NwVar[ BPar2 ] (the first "always free"
   position in NwVar). This can only happen if PPar2 > 0.

   If a "floor" to the Cutting Plane model has been set, i.e. FlrNme != InINF
   (FlrNme is a protected field), NwVar[ FlrNme ] *must* be set == to the
   Right Hand Side of the new variable: this is responsibility of the
   *derived class*, that however can write into NwVar[] (provided, obviously,
   that it is != NULL). */

/*--------------------------------------------------------------------------*/

   virtual void RemoveVar( cIndex i ) = 0;

/* Removes the Multiplier 'i' from the direction finding subproblem: 'i'
   is guaranteed to be a defined Multiplier. Note that this removal may be
   only temporary, i.e., the Multiplier may be re-inserted afterwards, hence
   the information about this Multiplier (the corresponding entries of all
   the items in the Bundle) must be kept. */

/*--------------------------------------------------------------------------*/

   virtual void SubstVar( cIndex i ) = 0;

/* Removes the Multiplier with name 'i', substituting it with the Multiplier
   with name 'NumVar' (NumVar is decremented *just before* the call to this
   method, hence it is the name of the last variable). The "original copy" of
   the latter is removed; hence, when called with i == NumVar this method
   simply removes the last Multiplier.
   The (previous) Multiplier 'i' is not assumed to be re-inserted again, hence
   the information about it (the entries ...) need *not* to be kept.

   Note that, when this method is called, both 'i' and 'NumVar' may *not*
   have been defined in the subproblem (if PPar2 > 0). */

/*--------------------------------------------------------------------------*/

/* Note that the methods for handling the Right Hand Side

   virtual void SetRHS( cSgRow tRHS ) = 0;

   virtual void GetRHS( SgRow tRHS ) = 0;

   defined in the public interface of Bundle must in fact be implemented
   by the derived class implementing the subproblem solver. */

/*--------------------------------------------------------------------------*/

   virtual BStatus SolveSubP( cHpNum tt , HpNum &dNrm , HpNum &Sgm ) = 0;

/* Attempts to solve the direction finding subproblem with t = tt.
   Note that the solver is required to conform to the following standards in
   some "extreme" cases:

   - if the Bundle is *empty*, the optimal (both primal and dual) solution is 
     all-0;

   - if the Bundle contains *no subgradients* (only constraints) the dual
     variable v (the primal linear equality constraint) must be ignored, and
     the "pure feasibility" problem solved.

   It must return:

   - kOK is the solver found an optimal solution;

   - kUnfeasible if the subproblem is dual unfeasible (primal unbounded or
     empty, see ReadMult() above);

   - kError in case of a failure in the solver.

   If PPar2 > 0, the subproblem to be solved is a *restriction* of the *full*
   subproblem where (possibly many, or *all*) dual variables are forced to 0:
   hence, the restricted subproblem may be *empty* (so that the solver might
   have to return kUnfeasible) while the full problem has a nonempty feasible
   set. Hence, after a kUnfeasible return the "price in" is done to find if
   new variables have to be added to the subproblem in order to try to make
   it nonempty: for this purpose, ReadMult() ... [see above] must return the
   multipliers that make a *feasible ascent extreme ray*. In fact, the
   direction d[] returned by Price() with those multipliers is nonzero
   (strictly positive for constrained variables) if and only if that variable
   have to be added. Clearly, this can happen only if ADD_CNST > 0.

   After (succesfull) termination, dNrm must hold the norm (1,2 or Infinity
   according to the type of subproblem solved) of the optimal solution d and
   Sgm the "aggregate linearization error". */

/*--------------------------------------------------------------------------*/

/* Note that the methods for reading the optimal primal multipliers

   virtual cHpRow ReadMult( void ) = 0;

   virtual cIndex_Set ReadBase( void ) = 0;

   virtual Index ReadBDim( void ) = 0;

   defined in the public interface of Bundle must in fact be implemented
   by the derived class implementing the subproblem solver. */

/*--------------------------------------------------------------------------*/

   virtual cLMRow Price( void ) = 0;

/* At each iteration, an optimal solution d for the direction finding
   subproblem must be found. If PPar2 > 0 not all the Multipliers are defined
   in the subproblem, and only the corresponding entries of d need to be
   calculated. Every PPar2 iterations, however, the Bundle algorithm needs to
   check the remaining Multipliers to see whether some of them should be
   "created" (a "price in").

   If PPar2 > 0, this method is (every PPar2 iterations or each time
   "convergence" is detected) called by Bundle *after* SolveSubP() to ask the
   subproblem solver for the entries of d corresponding to *all* Multipliers,
   even those that are not defined. The returned (read-only) pointer must
   point to a vector d such that d[ i ] contains the entry of the direction
   relative to the variable 'i' for i in [ 0 .. NumVar - 1 ].

   If MakeLambda1() is called in the same iteration, it is called *after*
   Price() - that is, the "new" entries of d can be assumed to have been
   already calculated. */

/*--------------------------------------------------------------------------*/

   virtual HpNum CheckAlfa( cBOOL All = FALSE ) = 0;

/* Checks the existence of negative linearization errors in the subproblem.
   If they are found, it means that the value of Fi( Lambda ) is incorrect:
   in order to regain nonnegativity, a (negative) number must be subtracted
   from Fi( Lambda ).

   CheckAlfa() must look for this condition and evenctually change the
   linearization errors: the returned value must be the (negative) number to
   be subtracted from Fi( Lambda ).

   If All == FALSE, only the "newly inserted" subgradients need to be checked
   (that is, all those inserted *after* the previous call to CheckAlfa());
   otherwise, all subgradients in the Bundle need to be checked. */

/*--------------------------------------------------------------------------*/

   virtual void ChgAlfa( cIndex i , cHpNum DeltaAlfai ) = 0;

/* Changes the Alfa (linearization error of a subgradient, or RHS of a linear
   constraints) of the i-th item, setting adding it DeltaAlfai. */

/*--------------------------------------------------------------------------*/

   virtual void ChgAlfa( cHpRow DeltaAlfa , BOOL Mde ) = 0;

/* Performs a "relative" change of the vector Alfa[], depending on Mde:

   - if Mde == TRUE, Alfa[ i ] += DeltaAlfa[ i ] for each i = name of an item
     in the Bundle;

   - if Mde == FALSE changes *only* the Alfa[] for the *subgradients*, i.e.

                  / DeltaAlfa[ j ]   i is a subgradient belonging to the
		  |                  j-th "component" of Fi()
     Alfa[ i ] += |
                  \     0            i is a constraint

     if SPRBL_FI == 0, DeltaAlfa[ 0 ] contains the change for the unique
     component of Fi(). */

/*--------------------------------------------------------------------------*/

   virtual void ChangeCurrPoint( cLMRow DLambda , cHpNum DFi ) = 0;

   virtual void ChangeCurrPoint( cHpNum Tau , cHpNum DFi ) = 0;

/* These two methods updates the data of the direction finding subproblem
   in response to a change of the "current point". In the first case,

     < new current point > = < old current point > + DLambda

   where DLambda is a NumVar-long vector containing the displacement of the
   variable with name 'i' in DeltaLambda[ i ]. The displacement will be 0
   (zero) if variable 'i' is not defined. DFi must be the change in the
   value of the function between the two points, i.e.

     Fi( < new current point > ) - Fi( < old current point > )

   The second case is equivalent for DLambda = ( Tau / t ) * d, where d is
   the optimal solution of the direction finding subproblem (Tau is a
   *relative step* along d, relative w.r.t. t). */

/*--------------------------------------------------------------------------*/

   virtual cHpRow ReadAlfa( void ) = 0;

/* Must return a read-only pointer to a vector containing in position i the
   current value of the linearization error/right hand side of item 'i'. */

/*--------------------------------------------------------------------------*/

   virtual void ChgRHSi( cIndex i , cSgNum DeltaRHSi ) = 0;

/* Changes the i-th entry of the R.H.S. to RHSi. If PPar2 > 0 and the
   Multiplier 'i' is *not* defined in the subproblem prior to the call, it
   *is defined* within the call like with AddVar( i , 0 , NULL ): note that
   'i' is inserted into LamBase only *after* the call to ChgRHSi().
   i is guaranteed to be < NumVar. */

/*--------------------------------------------------------------------------*/

   virtual void MakeLambda1( cHpNum Tau ) = 0;

/* Writes in (the array poined by) Lambda1 (protected field of type LMRow) the
   result of Lambda + ( Tau / t ) * d, where (Lambda is like Lambda1 and) d is
   the optimal solution of the direction finding subproblem. 
   The "format" of Lambda and Lambda1 depends on PPar2, and it is analogous
   to that of Lam1 in Fi() [see above]: if PPar2 == 0, they are just "dense"
   NumVar-vectors of the Lagrangean Multipliers, L[ i ] being the multiplier
   relative to variable 'i'. Otherwise (PPar2 > 0), both Lambda and Lambda1
   are intended to be "restricted" to the set of Multipliers that are defined
   in the subproblem: the set is available in the protected fields LamBase
   (Index_Set) and LamDim (Index), that is, only the first LamDim entries of
   Lambda and Lambda1 are significative, and Lambda[ i ] is the Mutiplier
   relative to LamBase[ i ]. LamBase is guaranteed to be ordered in increasing
   sense and InINF-terminated. */

/*--------------------------------------------------------------------------*/

   virtual HpNum EpsilonD( void ) = 0;

/* The precision implemented by the subproblem solver:
   Lambda[ i ] <= EpsilonD() is considered as == 0. */

/*--------------------------------------------------------------------------*/

   virtual void SensitAnals( HpNum &lp , HpNum &cp ) = 0;

/* Gives a linear lower approximation on the optimal value of v (the predicted
   decrease from the Cutting Plane model) as a function of the parameter t in
   the subproblem. That is, returns two numbers lp and cp such that

     v( t ) >= t * lp + cp. */

/*--------------------------------------------------------------------------*/

#if( DO_AGGR )

   virtual void AggregateZ( cIndex where ) = 0;

/* Substitute item 'i' with the *optimal aggregate subgradient* Z from the
   latest *succesfull* call to SolveSubP() [see AggrPrimalSol() above].

   If DO_AGGR == 1, the method is only called after a succesfull SolveSubP(),
   and prior to the elimination of any of the "basic" items.

   If DO_AGGR == 2 instead, this method is (also) called after a failure
   reported by SolveSubP(): hence, the Z that is inserted must be that of the
   latest "good" solution of the subproblem. This may require saving some
   information (Z itself and/or the optimal base). However, it is still
   guaranteed that all the items in the latest optimal base will still be
   there when AggregateZ() is called.

   Note that i may be the name of an existing item, that must be deleted.

   The method has to *call AggrPrimalSol()*, passing it the optimal primal
   multipliers used for the aggregation: if ADD_CNST > 0, this may *not* be
   the full optimal primal solution of the subproblem since the multipliers
   relative to constraints must not be there (Z is an aggregate subgradient).
   */

#endif

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

  Index MaxIter;       // Maximum number of iterations

  Index BPar1;         // control related "Int" parameters
  Index BPar2;

  HpNum Incr;          // increase/decrease parameters
  HpNum Decr;

  HpNum m1;            // control related "Float" parameters
  HpNum m2;
  HpNum m3;

  HpNum tStar;         // optimality related parameters
  HpNum EpsLin;

  HpNum tMaior;        // step size (max, min and initial)
  HpNum tMinor;
  HpNum tInit;

  HpNum MPar1;         // step selection related parameters
  HpNum MPar2; 
  HpNum MPar3;

  Index PPar1;         // pricing related parameters
  Index PPar2;
  Index PPar3;

  Index MaxNumVar;     // maximum number of variables
  Index NumVar;        // actual number of variables (<= MaxNumVar)
  Index NNVars;        // actual number of NN variables (<= NumVar)
  Index FlrNme;        // the position of the "floor" in the Bundle

  LMRow Lambda;        // the current point
  LMRow Lambda1;       // the tentative point

  Index_Set LamBase;   // the set of indices of Lambda
  Index_Set Lam1Bse;   // the set of indices of Lambda1
  Index LamDim;        // dimension of LamBase

  BStatus Result;     // result of the latest call to Solve()

  BOOL BHasChgd;      // TRUE if LamBase has changed during the latest
                      // pricing (never set to TRUE if PPar2 == 0, unless
                      // at the first iteration)

  HpNum t;            // the (tremendous) t parameter

  HpNum FiLambda;     // Fi( Lambda )
  HpNum FiLambda1;    // Fi( Lambda1 )

  HpNum LowerBound;   // lower bound on the optimal value (may be - HpINF)
  char  LBMode;       // how LowerBound is used

  #if( SPRBL_FI )
   Index FiK;         // number of "components" of Fi()
   HpRow FiLambdak;   // vector of [ Fi[ k ]( Lambda ) ]
   HpRow FiLambda1k;  // vector of [ Fi[ k ]( Lambda1 ) ]

   Index whtsG1;      // for some methods in the interface with the subproblem
                      // solver, this tells to which component G1 belongs
  #endif

  #if( LOG_BND )
   ostream *BLog;      // the output stream object
   char BLLvl;         // the "level of verbosity"
  #endif

  #if( TIMERS_B )
   OPTtimers *TOTt;    // total time
   OPTtimers *FIt;     // Fi() time
  #endif

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

  inline void Delete( cIndex i );

/*--------------------------------------------------------------------------*/

  Index FindAPlace( void );

/*--------------------------------------------------------------------------*/

  inline HpNum Heuristic1( void );

/*--------------------------------------------------------------------------*/

  inline HpNum Heuristic2( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

  Index FiEvaltns;     // total number of Fi() calls
  Index GiEvaltns;     // total number of Gi() calls
  Index SCalls;        // nuber of calls to Solve() (the current included)
  Index ParIter;       // nuber of iterations in this run

  LMNum ScPr1;         // ScalarProduct( dir , G1 )
  HpNum Alfa1;         // linearization error of G1 w.r.t. Lambda

  HpNum dNorm;         // || dir ||^2
  HpNum Sigma;         // convex combination of the Alfa's
  HpNum v;             // t * dNorm + Sigma: the "predicted increase"
  HpNum vStar;         // tStar * dNorm + Sigma: the "maximum expected
                       // increase" used in the stopping criteria

  HpNum DeltaFi;       // FiLambda - FiLambda1
  HpNum FiBest;        // best value of Fi() found so far
  HpNum RfrncFi;       // the value of Fi() where the zero of the Cutting
                       // Plane model is fixed: it is usually == FiLambda,
                       // but it is always finite even if FiLambda == HpINF
  HpNum EpsU;          // precison required by the long-term t-strategy

  Index FreDim;        // number of free positions in the Bundle
  Index SgMxNm;        // max "name" in the Bundle (+ 1)
  Index_Set FreList;   // list (heap) of free positions in the Bundle
  #if( DO_AGGR )
   Index whisZ;        // name of the current "aggregate subgradient" Z
  #endif

  SIndex_Set OOBase;   // Out-Of-Base counters
  Index_Set InctvCtr;  // "out of base" counter for variables
  Index_Set nBase;     // temporary

  BOOL dChanges;       // TRUE if the subgradients returned by GetGi() change
                       // the solution of the subproblem
  BOOL StrctNZ;        // TRUE if small negative multipliers are not allowed
  BOOL SSDone;         // TRUE if the laste step was a SS

/*--------------------------------------------------------------------------*/

 };  // end( class Bundle )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void Bundle::SetStrctNZ( cBOOL SNZ )
{
 StrctNZ = SNZ;
 }

/*--------------------------------------------------------------------------*/

inline Bundle::BStatus Bundle::GetBStatus( void )
{
 return( Result );
 }

/*--------------------------------------------------------------------------*/

inline cLMRow Bundle::ReadSol( cIndex_Set &I , Index &D )
{
 I = LamBase;
 D = LamDim;

 return( Lambda );
 }

/*--------------------------------------------------------------------------*/

inline void Bundle::ReadSol( register LMRow L )
{
 if( LamBase )
 {
  register cLMRow tL = Lambda;
  register cIndex_Set tLB = LamBase;
  register Index h = *(tLB++);
  for( register Index i = 0 ; i < NumVar ; )
   if( i++ == h )
   {
    *(L++) = *(tL++);
    h = *(tLB++);
    }
   else
    *(L++) = 0;
  }
 else
 {
  register cLMRow tL = Lambda + NumVar;
  for( L += NumVar ; tL > Lambda ; )
   *(--L) = *(--tL);
  }
 }  // end( ReadSol( L ) )

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::ReadFiVal( void )
{
 return( FiLambda );
 }

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::ReadBestFiVal( void )
{
 return( FiBest );
 }

/*--------------------------------------------------------------------------*/

inline BOOL Bundle::CurrentIsLast( void )
{
 return( SSDone );
 }

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::ReaddNorm( void )
{
 return( dNorm );
 }

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::ReadSigma( void )
{
 return( Sigma );
 }

/*--------------------------------------------------------------------------*/

inline Index Bundle::FiEval( void )
{
 return( FiEvaltns );
 }

/*--------------------------------------------------------------------------*/

inline Index Bundle::GiEval( void )
{
 return( GiEvaltns );
 }

/*--------------------------------------------------------------------------*/

inline Index Bundle::NrCalls( void )
{
 return( SCalls );
 }

/*--------------------------------------------------------------------------*/

inline Index Bundle::NrIters( void )
{
 return( ParIter );
 }

/*--------------------------------------------------------------------------*/

inline Index Bundle::GetNumVar( void )
{
 return( NumVar );
 }

/*--------------------------------------------------------------------------*/

#if( TIMERS_B )

inline void Bundle::TimeBAll( double &t_us , double &t_ss )
{
 t_us = t_ss = 0;
 TOTt->Read( t_us , t_ss );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline void Bundle::TimeBFi( double &t_us , double &t_ss )
{
 t_us = t_ss = 0;
 FIt->Read( t_us , t_ss );
 }

#endif

/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
 };  // end( namespace Bundle_di_unipi_it )
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* Bundle.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File Bundle.h ------------------------------*/
/*--------------------------------------------------------------------------*/
