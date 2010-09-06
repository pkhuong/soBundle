/*--------------------------------------------------------------------------*/
/*----------------------------- File Bundle.C ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--    Minimization of a convex NonDifferentiable (possibly nonexact)    --*/
/*--       function over a polyhedral set via a "Bundle" algorithm.       --*/
/*--                                                                      --*/
/*--                            VERSION 2.12                              --*/
/*--                           12 - 10 - 2004                             --*/
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
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Bundle.h"

#include "OPTvect.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

// used in FiAndGi()- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if( SPRBL_FI )
 #define DeltaFik() FiLambda1k[ whtsG1 ] - FiLambdak[ whtsG1 ]
#else
 #define DeltaFik() - DeltaFi
#endif

// used in RemoveItem() - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if( ! ADD_CNST )
 #define BCSize() 0
#endif

// log-related macros - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if( LOG_BND )
 #define BLOG( l , x ) if( BLLvl > l ) *BLog << x

 #define BLOG2( l , c , x ) if( c ) BLOG( l , x )

 #define BLOG3( l , c , x , y ) if( BLLvl > l ) \
                                if( c ) *BLog << x; else *BLog << y
#else
 #define BLOG( l , x )

 #define BLOG2( l , c , x )

 #define BLOG3( l , c , x , y )
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
using namespace Bundle_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

const HpNum Nearly  = 1.01;
const HpNum Nearly2 = 1.02;

/*--------------------------------------------------------------------------*/
/*---------------------- IMPLEMENTATION OF Bundle --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

Bundle::Bundle( istream *iStrm , cIndex NV , cIndex ENV )
{
 // initialize algorithmic parameters - - - - - - - - - - - - - - - - - - - -

 DfltdSfInpt( iStrm , MaxIter , InINF );

 DfltdSfInpt( iStrm , BPar1 , Index( 10 ) );
 DfltdSfInpt( iStrm , BPar2 , Index( 100 ) );

 DfltdSfInpt( iStrm , Incr , HpNum( 10 ) );
 DfltdSfInpt( iStrm , Decr , HpNum( .1 ) );

 DfltdSfInpt( iStrm , m1 , HpNum( .1 ) );
 DfltdSfInpt( iStrm , m2 , HpNum( .3 ) );
 DfltdSfInpt( iStrm , m3 , HpNum( 3 ) );

 DfltdSfInpt( iStrm , tStar  , HpNum( 100 ) );
 DfltdSfInpt( iStrm , EpsLin , HpNum( 1e-6 ) );

 DfltdSfInpt( iStrm , tMaior , HpNum( 1e+6 ) );
 DfltdSfInpt( iStrm , tMinor , HpNum( 1e-6 ) );
 DfltdSfInpt( iStrm , tInit  , HpNum( 1 ) );

 DfltdSfInpt( iStrm , MPar1 , HpNum( 0 ) );
 DfltdSfInpt( iStrm , MPar2 , HpNum( 0 ) );
 DfltdSfInpt( iStrm , MPar3 , HpNum( 8 ) );

 DfltdSfInpt( iStrm , PPar1 , Index( 30 ) );
 DfltdSfInpt( iStrm , PPar2 , Index( 10 ) );
 DfltdSfInpt( iStrm , PPar3 , Index( 5 ) );

 // other initializations - - - - - - - - - - - - - - - - - - - - - - - - - -

 MaxNumVar = NV;
 NNVars = NumVar = min( NV , ENV );

 if( PPar2 )
 {
  LamBase = new Index[ MaxNumVar + 2 ];
  nBase   = new Index[ MaxNumVar + 2 ];
  *(LamBase++) = InINF;  // set a "floor" to the bases: this is used in
  *(nBase++) = InINF;    // FormD() when they are read high -> low

  LamBase[ LamDim = 0 ] = InINF;

  Lam1Bse = new Index[ MaxNumVar + 1 ];
  *Lam1Bse = InINF;

  if( PPar3 )
   InctvCtr = new Index[ MaxNumVar ];
  else
  {
   InctvCtr = NULL;
   PPar1 = InINF;
   }
  }
 else
 {
  LamBase = nBase = Lam1Bse = NULL;
  LamDim = NumVar;
  InctvCtr = NULL;
  }

 Lambda = new LMNum[ MaxNumVar ];
 VectAssign( Lambda , LMNum( 0 ) , NumVar );  // the default starting point

 Lambda1 = new LMNum[ MaxNumVar ];

 OOBase = new SIndex[ BPar2 ];
 VectAssign( OOBase , SInINF , BPar2 );

 FreList = new Index[ BPar2 ];

 FreDim = FiEvaltns = GiEvaltns = SCalls = SgMxNm = 0;
 FiLambda = FiLambda1 = FiBest = HpINF;
 LowerBound = - HpINF;
 RfrncFi = 0;
 LBMode = 0;
 BHasChgd = FALSE;
 FlrNme = InINF;
 #if( DO_AGGR )
  whisZ = InINF;
 #endif
 SSDone = StrctNZ = FALSE;
 Result = kError;

 #if( LOG_BND )
  BLog = NULL;
  BLLvl = 0;
 #endif

 #if( TIMERS_B )
  TOTt = new OPTtimers();
  FIt = new OPTtimers();
 #endif

 }  // end( Bundle )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

#if( SPRBL_FI )

 Bundle::void SetK( cIndex k )
 {
  FiK = k;

  FiLambdak  = new HpNum[ k ];
  FiLambda1k = new HpNum[ k ];

  for( register Index i = k ; i-- ; )
   FiLambdak[ i ] = RfrncFi;

  whtsG1 = InINF;
  }

#endif

/*--------------------------------------------------------------------------*/

void Bundle::SetUC( cIndex Strt , cIndex Stp , cIndex_Set Which )
{
 NNVars = NumVar - ( Stp - Strt );

 if( Which )
  while( *(Which++) < InINF )
   NNVars--;

 }  // end( Bundle:: SetUC )

/*--------------------------------------------------------------------------*/

void Bundle::SetLambda( register cLMRow tLambda )
{
 // first, compute Lambda - NewLambda for updating all the things - - - - - -
 // (Alfa[], the bounds ...) that depends on Lambda - - - - - - - - - - - - -

 if( PPar2 )
 {
  register LMRow tL = Lambda;        // re-construct the old point Lambda,
  register Index_Set tLB = LamBase;  // "dense", in Lambda1
  for( register Index h = 0 ; h < NumVar ; )
   if( h == *tLB )
   {
    Lambda1[ h++ ] = *(tL++);
    tLB++;
    }
   else
    Lambda1[ h++ ] = 0;
  }
 else
  Swap( Lambda , Lambda1 );  // remember the previous point in Lambda1

 if( tLambda )
  VectSubtract( Lambda1 , tLambda , NumVar );

 // now change the Current Point in the subproblem solver: pretend that the -
 // value of Fi() in the *new* Lambda is *the same as in the old Lambda*, - -
 // since the correction will be done as soon as Fi() is computed - - - - - -

 ChangeCurrPoint( Lambda1 , 0 );

 FiLambda = FiBest = HpINF;  // Fi( Lambda ) is not known

 // now, the actual setting of Lambda - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LMNum teD = EpsilonD();

 if( PPar2 )  // record the point in Lambda, "sparse" - - - - - - - - - - - -
 {
  LamDim = 0;

  if( tLambda )
   for( register Index i = 0 ; i < NumVar ; i++ )
    if( ABS( Lambda[ LamDim ] = *(tLambda++) ) > teD )
    {
     LamBase[ LamDim++ ] = i;
     if( PPar3 )
      InctvCtr[ i ] = 0;
     }

  LamBase[ LamDim ] = InINF;
  BHasChgd = TRUE;

  SetVars();  // set the current set of "active" variables
  }
 else  // record the point in Lambda, "dense" - - - - - - - - - - - - - - - -
  if( tLambda )
   VectAssign( Lambda , tLambda , NumVar );
  else
   VectAssign( Lambda , LMNum( 0 ) , NumVar );

 SSDone = FALSE;

 }  // end( SetLambda )

/*--------------------------------------------------------------------------*/

Index Bundle::SetLowerBound( cHpNum LwBnd )
{
 if( LBMode == 2 )       // set the "floor"
  if( LwBnd > - HpINF )  // set the bound - - - - - - - - - - - - - - - - - -
  {
   if( FlrNme < InINF )  // the "floor" is already there
    ChgAlfa( FlrNme , LowerBound - LwBnd );
   else
    if( ( FlrNme = FindAPlace() ) < InINF )
    {
     OOBase[ FlrNme ] = - SInINF;
     GetRHS( GetItem( FlrNme ) );

     HpNum fScPr;
     HpNum fAlfa = RfrncFi - LwBnd;
     SetItem( 0 , 0 , fAlfa , fScPr , NULL );
     }
   }
  else  // remove the bound - - - - - - - - - - - - - - - - - - - - - - - - -
   if( FlrNme < InINF )
   {
    Delete( FlrNme );
    FlrNme = InINF;
    }

 LowerBound = LwBnd;
 return( FlrNme );

 }  // end( SetLowerBound )

/*--------------------------------------------------------------------------*/

void Bundle::SetLBMode( const char LBMd )
{
 LBMode = LBMd;

 if( ( LBMd < 2 ) && ( FlrNme < InINF ) )
 {
  Delete( FlrNme );
  FlrNme = InINF;
  }
 }

/*--------------------------------------------------------------------------*/

#if( LOG_BND )

void Bundle::SetBLog( ostream *log , const char lvl )
{
 if( ( BLog = log ) )
  BLLvl = lvl;
 else
  BLLvl = 0;

 if( BLLvl > 1 )
  *BLog << endl << "Vars = " << NumVar << " (" << MaxNumVar
	<< ") ~ Max # = " << BPar2 << " ~ Rfrsh = " << BPar1 << " ~ t* = "
	<< tStar << " ~ Eps = " << EpsLin << endl
        << "t in [" << tMinor << ", " << tMaior << "] (tInit = " << tInit
        << ") ~ fFactors = (" << Incr << "/" << Decr << ")" << endl
        << "m1 = " << m1 << " ~ m2 = " << m2 << " ~ m3 = " << m3
        << " ~ Pricing: " << PPar1 << " - " << PPar2 << " - " << PPar3
	<< endl;
 }

#endif

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

Bundle::BStatus Bundle::Solve( void )
{
 #if( TIMERS_B )
  TOTt->Start();
 #endif

 Result = kOK;
 ParIter = 0;
 t = tInit;
 SCalls++;

 if( MPar1 )
  EpsU = MPar3 * EpsLin;

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 do {
  // construct the direction d- - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FormD();

  if( Result )  // problems in the subproblem solver
   break;

  // maintainance operations- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UpdtCntrs();  // update out-of-base counters and eliminate outdated info

  BLOG( 1 , endl << "{" << SCalls << "-" << ParIter << "-" << BSize() << "-"
	    << ReadBDim() << "} t = " << t << " ~ || d || = " << dNorm
	    << " ~ Sigma = " << Sigma << " ~ Fi = " );

  BLOG3( 1 , FiLambda == HpINF , " - INF" , - FiLambda );
  BLOG( 1 , endl << "           " );

  BLOG2( 1 , PPar2 , "LamDim = " << LamDim << " ~ " );
  BLOG2( 1 , MPar1 , "eU = " << EpsU << " ~ " );

  // hook for derived classes - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  EIStatus rs = EveryIteration();

  if( rs == kEIAbort )
  {
   BLOG( 1 , " EveryIteration():STOP" << endl );

   Result = kAbort;
   break;
   }

  if( rs == kEILoopNow )
  {
   BLOG( 1 , " EveryIteration():loop" << endl );

   continue;
   }

  // check for optimality - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( rs == kEINorm ) && ( FiLambda < HpINF ) )
   if( vStar <= EpsLin * max( ABS( FiLambda ) , HpNum( 1 ) ) )
    break;

  ParIter++;  // iterations where Fi() is not evaluated do not count

  // calculate Lambda1- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FormLambda1( t );

  // calculate Fi( Lambda1 )- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for(;;)
  {
   FiAndGi();

   if( ( ! dChanges ) && ( DeltaFi <= m1 * v ) )
    if( LooseSubgradients() )  // signal "danger of cycling", and, on demand
     continue;                 // of the derived class, re-compute Fi()

   break;
   }

  BLOG2( 1 , LowerBound > - HpINF , "UB = " << - LowerBound << " ~ " );
  BLOG( 1 , "Fi1 = " );

  if( FiLambda1 == - HpINF )
  {
   BLOG( 1 , "+ INF => STOP." << endl );

   Result = kUnbounded;
   break;
   }

  BLOG3( 1 , FiLambda1 == HpINF , " - INF ~ R.H.S. = " << Alfa1 << endl ,
	    - FiLambda1 << " ~ Alfa1 = " << Alfa1 << " ~ Gi1xd = "
	    << - ScPr1 << endl );

  // check the Lower Bound- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( FiLambda < HpINF )  // .. but only if Fi( Lambda ) is defined
  {
   if( FiBest <= LowerBound )
   {
    Result = kUnbounded;
    break;
    }

   if( LBMode == 1 )
    if( FiBest - EpsLin * ABS( FiBest ) <= LowerBound )
     break;
   }

  #if( ADD_CNST )
   // avoid the t-changing phase if Lambda1 is unfeasible - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // note: one possible alternative t-strategy would be to set t to the
   // largest value that would have produced a feasible point, i.e.
   // t := Alfa1 / ( - ScPr1 ) [if ScPr1 were calculated]

   if( FiLambda1 == HpINF )
    continue;
   else
  #endif
    if( FiLambda == HpINF )  // if reached feasibility- - - - - - - - - - - -
    {                        // - - - - - - - - - - - - - - - - - - - - - - -
     GotoLambda1();          // go to the feasible point
     continue;               // and start the actual minimization of Fi()
     }

  // the NS / SS decision - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( DeltaFi > m1 * v )  // SS( i )- - - - - - - - - - - - - - - - - - - - -
  {
   HpNum tt = Heuristic1();

   BLOG( 1 , endl << " SS: DeltaFi (" << DeltaFi << ") > m1 * v (" << m1 * v
	     << ") ~ Ht = " << tt << endl );

   GotoLambda1();

   t = min( min( tMaior , t * Incr ) , max( tt , t ) );
   }
  else  // NS( i )- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  {
   BLOG( 1 , endl << " NS( i ): DeltaFi (" << DeltaFi << ") <= m1 * v ("
	     << m1 * v << ")" );

   if( MPar1 )  // some Long-Term t-strategy active - - - - - - - - - - - - -
   {
    HpNum eU;
    cHpNum AFL = max( ABS( FiLambda ) , HpNum( 1 ) );

    if( MPar1 == 1 )           // Soft Long-Term t-strategy
    {
     if( max( vStar , v ) < EpsU * AFL )
      EpsU = max( max( vStar , v ) / AFL , EpsLin );

     eU = m2 * EpsU;
     }
    else                   // Hard Long-Term t-strategy
     eU = Nearly2 * EpsU;  // Nearly2 is used to avoid numerical problems

    if( v <= eU * AFL )  // v is "small" => t is
    {
     BLOG( 1 , " and small v => NS" << endl );

     continue;                       // do not allow decreasing t
     }
    }

   if( Alfa1 > m3 * Sigma )  // ! NS( ii )- - - - - - - - - - - - - - - - - -
   {
    HpNum tt = Heuristic2();

    BLOG( 1 , endl << " ! NS( ii ): Alfa1 (" << Alfa1 << ") > m3 * Sigma ("
              << m3 * Sigma << ") ~ Ht = " << tt << endl );

    if( MPar1 == 1 )
    {
     // don't let t become too small: EpsU * | FiLambda | is how much
     // the function is expected to increase in the future, so (at least
     // promised) increases of that order of magnitude are required

     HpNum vlin , vconst;
     SensitAnals( vlin , vconst );

     // v( t ) ~= - t * vlin - vconst with vlin and vconst <= 0; if
     // v( tNew ) must be >= m2 * EpsU * | FiLambda |, it must be
     // t >= ( m2 * EpsU * | FiLambda | + vconst ) / ( - vlin )

     if( - vlin < HpeM )  // v( t ) is [almost] constant =>
      tt = t;             // the CP model is [almost] bounded
     else
      tt = max( tt ,
		( m2 * EpsU * max( ABS( FiLambda ) , HpNum( 1 ) ) + vconst )
		/ ( - vlin ) );
     }

    t = max( max( tMinor , t * Decr ) , min( tt , t ) );
    }
   #if( LOG_BND )
    else
     BLOG( 1 , " and NS( ii ) => NS" << endl );
   #endif
   }   // end else( SS( i ) )

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  } while( ParIter < MaxIter );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ParIter == MaxIter )
  Result = kMaxIter;

 #if( TIMERS_B )
  TOTt->Stop();
 #endif

 return( Result );

 }  // end( Bundle::Solve )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

#if( ADD_CNST )

 Index Bundle::AddConstraint( cSgRow C , HpNum V , cSIndex Cntr ,
			      cIndex_Set CBse )
 {
  Index wh = FindAPlace();

  if( wh < InINF )
  {
   if( PPar2 )
    if( CBse )
     V -= ScalarProduct( Lambda , LamBase , C , CBse );
    else
     V -= ScalarProduct( Lambda , C , LamBase );
   else
    if( CBse )
    {
     register HpNum tV = 0;
     register cSgRow tC = C;
     register cIndex_Set tB = CBse;
     for( register Index h ; ( h = *(tB++) ) < InINF ; )
      tV += *(tC++) * Lambda[ h ];

     V -= tV;
     }
    else
     V -= ScalarProduct( Lambda , C , NumVar );

   VectAssign( GetItem( wh ) , C , NumVar );

   HpNum fScPr;
   SetItem( HpINF , 0 , V , fScPr , CBse );

   OOBase[ wh ] = - Cntr;

   if( V < - EpsilonD() )  // Lambda is unfeasible
    FiLambda = FiBest = HpINF;
   }

  return( wh );

  }  // end( Bundle::AddConstraint )

#endif

/*--------------------------------------------------------------------------*/

void Bundle::RemoveItem( cIndex Name )
{
 if( OOBase[ Name ] == SInINF )  // no such item exists ...
  return;                        // silently return

 Delete( Name );                 // delete it

 if( FlrNme == Name )            // if it was the "floor" ...
  FlrNme = InINF;

 if( ! SgMxNm )                  // it was the last item in the Bundle,
 {                               // which is now empty
  FiLambda = FiBest = HpINF;
  FreDim = 0;
  }
 }  // end( RemoveItem )

/*--------------------------------------------------------------------------*/

void Bundle::RemoveItems( void )
{
 RmvItems();

 FiLambda = FiBest = HpINF;
 FreDim = SgMxNm = 0;
 FlrNme = InINF;
 }

/*--------------------------------------------------------------------------*/

Index Bundle::AddVariable( SgRow NwVar , cLMNum IV , cBOOL NNV )
{
 if( NumVar >= MaxNumVar )
  return( InINF );

 Lambda[ LamDim ] = IV;

 if( PPar2 )  // L.V.G. "on"- - - - - - - - - - - - - - - - - - - - - - - - -
 {
  LamBase[ LamDim++ ] = NumVar;
  LamBase[ LamDim ] = InINF;

  if( PPar3 )
   InctvCtr[ NumVar ] = 0;
  }
 else  // L.V.G. "off - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  LamDim++;  // keep LamDim == NumVar

 if( ! NNV )
  AddVar( NumVar , NwVar , - LMINF );
 else
 {
  AddVar( NumVar , NwVar , - IV );
  NNVars++;
  }

 if( ABS( IV ) > EpsilonD() )
  FiLambda = FiBest = HpINF;

 return( NumVar++ );

 }  // end( AddVariable )

/*--------------------------------------------------------------------------*/

Index Bundle::RemoveVariable( cIndex i )
{
 if( i >= NumVar )  // there's no such thing as a variable 'i'
  return( InINF );  // quietly return

 if( NNVars )
  if( GetNN()[ i ] & 1 )
   NNVars--;

 NumVar--;

 SubstVar( i );  // do it in the MP solver: this has to be done *before*
                 // changing LamBase[] (if any) because it could be used
                 // by the MP solver

 if( i == NumVar )  // eliminating the last variable- - - - - - - - - - - - -
 {                  //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( PPar2 )  // L.V.G. "on" - - - - - - - - - - - - - - - - - - - - - - - -
  {
   if( LamBase[ LamDim - 1 ] == i )        // i is in base
   {
    if( ABS( Lambda[ --LamDim ] ) > EpsilonD() )
     FiLambda = FiBest = HpINF;

    LamBase[ LamDim ] = InINF;
    BHasChgd = TRUE;
    }
   }
  else  // L.V.G. "off" - - - - - - - - - - - - - - - - - - - - - - - - - - -
  {
   if( ABS( Lambda[ i ] ) > EpsilonD() )
    FiLambda = FiBest = HpINF;

   LamDim--;                        // keep LamDim == NumVar
   }
  }
 else  // rename NumVar as i- - - - - - - - - - - - - - - - - - - - - - - - -
 {     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( PPar2 )  // L.V.G. "on" - - - - - - - - - - - - - - - - - - - - - - - -
  {
   cIndex h = BinSearch( LamBase , LamDim , i );

   if( LamBase[ h ] == i )                 // i is in base
   {
    if( ABS( Lambda[ h ] ) > EpsilonD() )
     FiLambda = FiBest = HpINF;

    if( LamBase[ LamDim - 1 ] == NumVar )  // NumVar is in base
    {
     Lambda[ h ] = Lambda[ --LamDim ];
     LamBase[ LamDim ] = InINF;
     }
    else                                   // NumVar is *not* in base
    {
     ShiftVect( LamBase + h , (LamDim--) - h );
     ShiftVect( Lambda + h , LamDim - h );
     }

    BHasChgd = TRUE;
    }
   else                                    // i is *not* in base
    if( LamBase[ LamDim - 1 ] == NumVar )  // NumVar is in base
    {
     RotateRVect( Lambda + h , LamDim - h - 1 );
     ShiftRVect( LamBase + h , LamDim - h - 1 );
     LamBase[ h ] = i;

     BHasChgd = TRUE;
     }

   if( PPar3 )
    InctvCtr[ i ] = InctvCtr[ NumVar ];
   }
  else  // L.V.G. "off" - - - - - - - - - - - - - - - - - - - - - - - - - - -
  {
   if( ABS( Lambda[ i ] ) > EpsilonD() )
    FiLambda = FiBest = HpINF;

   Lambda[ i ] = Lambda[ NumVar ];
   LamDim--;                        // keep LamDim == NumVar
   }
  }  // end( else( i != NumVar ) )- - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( NumVar );

 }  // end( RemoveVariable )

/*--------------------------------------------------------------------------*/

void Bundle::ChangeRHSi( cIndex i , cSgNum DeltaRHSi )
{
 if( i >= NumVar )  // name of an undefined variable
  return;           // quietly return

 if( ! DeltaRHSi )  // actually no change is done
  return;           // quietly return

 ChgRHSi( i , DeltaRHSi );

 LMNum Lami;  // current value of Lambda[ i ]

 if( PPar2 )
 {
  Index h = BinSearch( LamBase , LamDim , i );

  if( LamBase[ h ] != i )  // the variable does not exist ...
  {
   ShiftRVect( Lambda + h , LamDim - h );
   Lambda[ h ] = Lami = 0;

   ShiftRVect( LamBase + h , (++LamDim) - h );  // ... so create it
   LamBase[ h ] = i;

   if( PPar3 )
    InctvCtr[ i ] = 0;
   BHasChgd = TRUE;
   }
  else
   Lami = Lambda[ h ];
  }
 else
  Lami = Lambda[ i ];

 if( Lami )
 {
  FiLambda += DeltaRHSi * Lami;
  RfrncFi = FiBest = FiLambda;
  }
 }  // end( ChangeRHSi )

/*--------------------------------------------------------------------------*/

void Bundle::TranslateSubgradients( cHpRow Deltas )
{
 if( FlrNme < InINF )
  LowerBound += Deltas[ FlrNme ];

 ChgAlfa( Deltas , TRUE );

 if( FiLambda < HpINF )
 {
  HpNum mA = CheckAlfa();

  if( mA < 0 )     // if there are negative Alfas
   RfrncFi -= mA;  // change the reference value to eliminate them

  FiLambda = FiBest = HpINF;
  }
 }  // end( TranslateSubgradients )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

Bundle::~Bundle()
{
 // output times and statistics - - - - - - - - - - - - - - - - - - - - - - -

  BLOG( 1 , endl << "Total Fi() evaluations = " << FiEvaltns << endl );

 #if( TIMERS_B )
  double u, s;
  TOTt->Read( u , s );

  BLOG( 1 , "Tot. time (s): user " << u << ", system " << s << endl );

  FIt->Read( u , s );

  BLOG( 1 , "Fi() time (s): user " << u << ", system " << s << endl );
 #endif

 // memory deallocation - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( SPRBL_FI )
  delete[] FiLambda1k;
  delete[] FiLambdak;
 #endif

 #if( TIMERS_B )
  delete FIt;
  delete TOTt;
 #endif

 delete[] FreList;
 delete[] OOBase;

 delete[] Lambda1;
 delete[] Lambda;

 delete[] InctvCtr;

 if( PPar2 )
 {
  delete[] Lam1Bse;
  delete[] --nBase;
  delete[] --LamBase;
  }
 }  // end( ~Bundle )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------- HOOKS FOR DERIVED CLASSES -------------------------*/
/*--------------------------------------------------------------------------*/

Bundle::EIStatus Bundle::EveryIteration( void )
{
 // Hard Long-Term t-strategy - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( MPar1 > 1 ) && ( FiLambda < HpINF ) )
 {
  cHpNum AFL = max( ABS( FiLambda ) , HpNum( 1 ) );

  if( vStar <= EpsU * AFL )
   if( EpsU <= EpsLin )
    return( kEINorm );
   else
   {
    BLOG( 1 , "Convergence, increasing accuracy." << endl << "           " );

    EpsU = max( vStar / ( MPar1 * AFL ) , EpsLin );
    }

  if( v <= EpsU * AFL )
  {
   BLOG( 1 , "small v => Augment( t )." << endl << "           " );

   HpNum vlin , vconst;
   SensitAnals( vlin , vconst );

   if( - vlin < HpeM )  // v( t ) is [almost] constant => it exists a t s.t.
    t = tStar;          // || d( t ) || [~]= 0 => the CP model is [~] bounded
   else
   {
    // v( tNew ) <= - vconst - tNew * vlin <= eps * |Fi|
    // => tNew = ( eps * |Fi| + vconst ) / ( - vlin )

    t = ( EpsU * AFL * Nearly + vconst ) / ( - vlin );
    t = min( tStar , t );
    }

   t = max( tMinor , t );  // safeguard
   return( kEILoopNow );

   }  // end if( Useless Step )
  }  // end if( Hard t-strategy )

 return( kEINorm );

 }  // end( Bundle::EveryIteration )

/*--------------------------------------------------------------------------*/
/*----------------------- OTHER PROTECTED METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

void Bundle::FormD( void )
{
 for(;;)
 {
  // first, solve the subproblem- - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  while( ( Result = SolveSubP( t , dNorm , Sigma ) ) == kError )  //- - - - -
  {
   BLOG( 1 , endl << " Warning: failure in subproblem solver [FormD]." );

   if( BSize() <= 1 )  // only one item in the Bundle
   {
    BLOG( 0 , endl << " ERROR: unable to recover from failure [FormD]." );
    return;
    }

   Index j = ReadBDim();
   Index i = InINF;

   // the last *removable* item in Base is eliminated - - - - - - - - - - - -

   for( register cIndex_Set tB = ReadBase() ; j-- ; )
    if( OOBase[ tB[ j ] ] >= 0 )
    {
     i = tB[ j ];
     break;
     }

   if( i == InINF )            // there are no *removable* items in Base- - -
    for( j = SgMxNm ; j-- ; )  // pick any removable item anywhere- - - - - -
     if( ( OOBase[ j ] >= 0 ) && ( OOBase[ j ] < SInINF ) )
     {
      i = j;
      break;
      }

   if( i == InINF )  // there are no removable items at all - - - - - - - - -
   {
    BLOG( 0 , endl << " ERROR: unable to recover from failure [FormD]." );
    return;
    }

   #if( DO_AGGR == 2 )
    if( whisZ == InINF )       // if Z is not in the Bundle
    {
     AggregateZ( whisZ = i );  // substitute i with Z
     OOBase[ whisZ ] = -1;     // and ensure it won't be eliminated
     }
    else                       // but if Z is already in the Bundle
   #endif
     Delete( i );              // ... just delete i

   }  // end while( error ) - - - - - - - - - - - - - - - - - - - - - - - - -

  vStar = dNorm * tStar + Sigma;
  v = dNorm * t + Sigma;
  #if( DO_AGGR )
   whisZ = InINF;  // Z is changed
  #endif

  if( ! PPar2 )  // no L.V.G.
   return;       // nothing else to do

  if( LamDim == NumVar )  // all variables are there
   break;                 // no "price in" to do (but possibly "price out")

  // LVG: price in- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // do pricing only for the first PPar1 iterations and then once every PPar2
  // iterations: however, do it also if convergence is detected or if the
  // subproblem is dual unfeasible

  #if( ADD_CNST )
   if( Result == kOK )
  #endif
    if( ( ParIter > PPar1 ) && ( ( ParIter - PPar1 ) % PPar2 ) &&
	( vStar > EpsLin * max( ABS( FiLambda ) , HpNum( 1 ) ) ) )
      return;

  cLMRow tdir = Price();  // ask for *all* the entries of d[]
  LMNum epsDir = EpsilonD() * t;
  register Index nBD = 0;
  register Index oBD = 0;

  if( ! NNVars )  // there are only UC variables- - - - - - - - - - - - - - -
  {
   for( register Index k = 0 ; k < NumVar ; k++ )
    if( LamBase[ oBD ] == k )
    {
     oBD++;
     nBase[ nBD++ ] = k;
     }
    else
     if( ABS( tdir[ k ] ) > epsDir )
     {
      nBase[ nBD++ ] = k;
      if( PPar3 )
       InctvCtr[ k ] = 0;

      AddVar( k , NULL , 0 );
      BLOG( 4 , endl << " Created variable " << k );
      }
   }
  else
   if( NNVars == NumVar )  // there are only NN variables - - - - - - - - - -
   {
    for( register Index k = 0 ; k < NumVar ; k++ )
     if( LamBase[ oBD ] == k )
     {
      oBD++;
      nBase[ nBD++ ] = k;
      }
     else
      if( tdir[ k ] < - epsDir )
      {
       nBase[ nBD++ ] = k;
       if( PPar3 )
        InctvCtr[ k ] = 0;

       AddVar( k , NULL , 0 );
       BLOG( 4 , endl << " Created variable " << k << " (>= 0)" );
       }
    }
   else  // there are both NN and UC variables- - - - - - - - - - - - - - - -
   {
    register char *IsNN = GetNN();
    for( register Index k = 0 ; k < NumVar ; k++ )
     if( LamBase[ oBD ] == k )
     {
      oBD++;
      nBase[ nBD++ ] = k;
      }
     else
      if( ( tdir[ k ] < - epsDir ) ||
	  ( ( tdir[ k ] > epsDir ) && ( ! ( IsNN[ k ] & 1 ) ) ) )
      {
       nBase[ nBD++ ] = k;
       if( PPar3 )
        InctvCtr[ k ] = 0;

       AddVar( k , NULL , 0 );
       BLOG( 4 , endl << " Created variable " << k );
       BLOG2( 4 , IsNN[ k ] & 1 , " (>= 0)" );
       }
    }  // end else( there are both NN and UC variables )- - - - - - - - - - -

  if( nBD == oBD )  // no changes in LamBase- - - - - - - - - - - - - - - - -
   break;
  else  // LamBase has changed- - - - - - - - - - - - - - - - - - - - - - - -
  {
   BHasChgd = TRUE;
   LamDim = nBD;

   nBase--; LamBase--; Lambda--;

   for( ; nBD > oBD ; nBD-- )
    if( LamBase[ oBD ] == nBase[ nBD ] )
     Lambda[ nBD ] = Lambda[ oBD-- ];
    else
     Lambda[ nBD ] = 0;

   nBase++; LamBase++; Lambda++;

   Swap( LamBase , nBase );
   LamBase[ LamDim ] = InINF;
   }

  // end price in - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  }  // end( for(;;) )

 if( ParIter > PPar1 )  // skip the first PPar1 iterations (all of them if- -
 {                      // PPar3 == 0 => PPar1 == InINF ) - - - - - - - - - -
  // LVG: price out - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cLMRow tdir = Price();
  LMNum epsDir = EpsilonD() * t;
  register Index nBD = 0;
  register Index oBD = 0;

  if( ! NNVars )  // there are only UC variables- - - - - - - - - - - - - - -
   for( register Index k ; ( k = LamBase[ oBD ] ) < InINF ; oBD++ )
    if( ( ABS( Lambda[ oBD ] ) <= epsDir ) &&
	( ABS( tdir[ k ] ) <= epsDir )     &&
        ( ++InctvCtr[ k ] >= PPar3 ) )
    {
     RemoveVar( k );
     BLOG( 4 , endl << " Eliminated variable " << k );
     }
    else
    {
     Lambda[ nBD ] = Lambda[ oBD ];
     InctvCtr[ LamBase[ nBD++ ] = k ] = 0;
     }
  else
   if( NNVars == NumVar )  // there are only NN variables - - - - - - - - - -
    for( register Index k ; ( k = LamBase[ oBD ] ) < InINF ; oBD++ )
     if( ( Lambda[ oBD ] <= epsDir ) &&
	 ( tdir[ k ] >= - epsDir )   &&
         ( ++InctvCtr[ k ] >= PPar3 ) )
     {
      RemoveVar( k );
      BLOG( 4 , endl << " Eliminated variable " << k << " (>= 0)" );
      }
     else
     {
      Lambda[ nBD ] = Lambda[ oBD ];
      InctvCtr[ LamBase[ nBD++ ] = k ] = 0;
      }
   else  // there are both NN and UC variables- - - - - - - - - - - - - - - -
   {
    register char *IsNN = GetNN();
    for( register Index k ; ( k = LamBase[ oBD ] ) < InINF ; oBD++ )
     if( ( ABS( Lambda[ oBD ] ) <= epsDir )                   &&
	 ( ( tdir[ k ] >= - epsDir ) &&
	   ( ( tdir[ k ] <= epsDir ) || ( IsNN[ k ] & 1 ) ) ) &&
         ( ++InctvCtr[ k ] >= PPar3 ) )
     {
      RemoveVar( k );
      BLOG( 4 , endl << " Eliminated variable " << k << " (>= 0)" );
      BLOG2( 4 , IsNN[ k ] & 1 , " (>= 0)" );
      }
     else
     {
      Lambda[ nBD ] = Lambda[ oBD ];
      InctvCtr[ LamBase[ nBD++ ] = k ] = 0;
      }
    }  // end else( there are both NN and UC variables )- - - - - - - - - - -

  if( nBD < oBD )  // some variables have been eliminated - - - - - - - - - -
  {
   LamBase[ LamDim = nBD ] = InINF;
   BHasChgd = TRUE;
   }
  }  // end if( skip the first PPar1 iterations ) - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 }  // end( FormD )

/*--------------------------------------------------------------------------*/

void Bundle::FormLambda1( cHpNum Tau )
{
 MakeLambda1( Tau );

 if( StrctNZ && NNVars )  // eliminate all the negative components- - - - - -
 {
  register LMRow tL1 = Lambda1;

  if( NNVars < NumVar )
  {
   register char *IsNN = GetNN();

   if( PPar2 )
   {
    register Index_Set tLB = LamBase;
    for( register Index h ; ( h = *(tLB++) ) < InINF ; tL1++ )
     if( ( IsNN[ h ] & 1 ) && ( *tL1 < 0 ) )
      *tL1 = 0;
    }
   else
    for( tL1 += LamDim , IsNN += LamDim ; tL1-- > Lambda1 ; )
     if( ( *(--IsNN) & 1 ) && ( *tL1 < 0 ) )
      *tL1 = 0;
   }
  else
   for(  tL1 += LamDim ; tL1-- > Lambda1 ; )
    if( *tL1 < 0 )
     *tL1 = 0;
  }

 if( BHasChgd && LamBase )
  VectAssign( Lam1Bse , LamBase , LamDim + 1 );

 }  // end( FormLambda1 )

/*--------------------------------------------------------------------------*/

void Bundle::FiAndGi( void )
{
 // find an available slot in the Bundle (if any) - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index wh = FindAPlace();
 if( wh == InINF )  // no space found ...
 {
  BLOG( 0 , endl << " ERROR: No space in the Bundle [FiAndGi]." << endl );

  FiLambda1 = -HpINF;
  return;
  }

 SgRow G1 = GetItem( wh );

 // call SetGi() the first time - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SetGi( G1 , wh );

 // call Fi() - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( TIMERS_B )
  FIt->Start();
 #endif

 if( LamDim < NumVar )
  FiLambda1 = Fi( Lambda1 , Lam1Bse , LamDim );
 else
  FiLambda1 = Fi( Lambda1 , NULL , NumVar );

 SSDone = FALSE;

 #if( TIMERS_B )
  FIt->Stop();
 #endif

 FiEvaltns++;

 if( FiLambda1 == - HpINF )
  return;

 HpNum Tau = t;

 #if( ADD_CNST )
  if( FiLambda1 == HpINF )  // Fi() is not defined in Lambda1
  {
   DeltaFi = HpINF;
   Tau = 0;
   }
  else
 #endif
   DeltaFi = RfrncFi - FiLambda1;

 if( FiLambda1 < FiBest )
  FiBest = FiLambda1;

 // then the (possibly many) call(s) to Gi()- - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 dChanges = FALSE;

 for(;;)
 {
  // call Gi()- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if( TIMERS_B )
   FIt->Start();
  #endif

  cIndex_Set SGBse;
  SIndex GoOn = GetGi( Alfa1 , SGBse );

  GiEvaltns++;

  #if( TIMERS_B )
   FIt->Stop();
  #endif

  #if( SPRBL_FI )
   whtsG1 = ABS( GoOn );
  #endif

  // calculate ScPr1 and Alfa1- - - - - - - - - - - - - - - - - - - - - - - -

  #if( ADD_CNST )
   if( FiLambda1 == HpINF )  // G1 is a constraint: "scale" the R.H.S.
   {
    if( PPar2 )
     if( SGBse )
      Alfa1 -= ScalarProduct( Lambda , LamBase , G1 , SGBse );
     else
      Alfa1 -= ScalarProduct( Lambda , G1 , LamBase );
    else
     if( SGBse )
     {
      register HpNum tA1 = 0;
      register cSgRow tG1 = G1;
      register cIndex_Set tB = SGBse;
      for( register Index h ; ( h = *(tB++) ) < InINF ; )
       tA1 += *(tG1++) * Lambda[ h ];

      Alfa1 -= tA1;
      }
     else
      Alfa1 -= ScalarProduct( Lambda , G1 , NumVar );

    if( FiLambda == HpINF )     // if Lambda is *feasible* ...
     Alfa1 = max( Alfa1 , 0 );  // Alfa1 *must* be >= 0

    OOBase[ wh ] = - ABS( GoOn );
    }
  #endif

  dChanges |= SetItem( DeltaFik() , Tau , Alfa1 , ScPr1 , SGBse );

  if( GoOn >= 0 )
   break;
  else
  {
   if( ( wh = FindAPlace() ) == InINF )
    break;

   G1 = GetItem( wh );
   SetGi( G1 , wh );
   }
  }  // end of G1-collecting loop - - - - - - - - - - - - - - - - - - - - - -

 // now check for negative Alfa1 and evenctually correct them - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLambda < HpINF )  // ... but only if Fi( Lambda ) is defined
 {
  HpNum mA = CheckAlfa();

  if( mA < 0 )
  {
   Alfa1 -= mA;
   FiLambda -= mA;               // Fi( Lambda ) was incorrect ...
   DeltaFi -= mA;                // ... hence DeltaFi was
   RfrncFi = FiBest = FiLambda;  // FiBest *may* have been set to a wrong
                                 // value: the only safe thing is to reset it

   BLOG( 1 , "Fi += Alfa1 = " << mA << " ~ " );
   }
  }
 }  // end( FiAndGi )

/*--------------------------------------------------------------------------*/

void Bundle::GotoLambda1( void )
{
 ChangeCurrPoint( t , FiLambda1 - RfrncFi );

 #if( LOG_BND )
  HpNum mA = CheckAlfa( TRUE );
  if( ( BLLvl > 1 ) && ( mA < 0 ) && ( FiLambda < HpINF ) )
   *BLog << endl << " Warning: some Alfai = " << mA << " < 0 [GotoLambda1].";
 #else
  CheckAlfa( TRUE );  // in theory, one should diminish FiLambda of the min.
                      // negative Alfa reported by CheckAlfa(); however,
                      // FiLambda is of no interest since we are moving
 #endif

 RfrncFi = FiLambda = FiLambda1;
 Swap( Lambda , Lambda1 );
 #if( SPRBL_FI )
  Swap( FiLambdak , FiLambda1k );
 #endif
 SSDone = TRUE;

 }  // end( GotoLambda1 )

/*--------------------------------------------------------------------------*/

void Bundle::UpdtCntrs( void )
{
 // increase all the OOBase[] counters but those == +/-SInINF - - - - - - - -
 // items whose OOBase[] becomes 0 (e.g. the newly entered items, which have
 // OOBase[] == -1) are set to +1, in such a way that all and only the items
 // in the optimal base has OOBase[] == 0

 register SIndex_Set tOO = OOBase + SgMxNm;
 for( ; tOO-- > OOBase ; )
  if( ( *tOO < SInINF ) && ( *tOO > -SInINF ) )
  {
   (*tOO)++;
   if( ! *tOO )
    (*tOO)++;
   }    

 // set to 0 the OOBase[] counter for items in base (if not < 0)- - - - - - -

 register cIndex_Set Bt = ReadBase();
 for( register Index i ; ( i = *(Bt++) ) < InINF ; )
  if( OOBase[ i ] > 0 )
   OOBase[ i ] = 0;

 // now remove outdated items - - - - - - - - - - - - - - - - - - - - - - - -

 for( tOO = OOBase + SgMxNm ; tOO-- > OOBase ; )
  if( ( *tOO < SInINF ) && ( *tOO > SIndex( BPar1 ) ) )
  {
   cIndex h = tOO - OOBase;
   Delete( h );
   DeletedI( h );
   }
 }  // end( UpdtCntrs )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

inline void Bundle::Delete( cIndex i )
{
 // deletes the item with name 'i'

 RmvItem( i );

 HeapIns( FreList , i , FreDim++ );
 OOBase[ i ] = SInINF;

 while( SgMxNm && ( OOBase[ SgMxNm - 1 ] == SInINF ) )
  SgMxNm--;
 }

/*--------------------------------------------------------------------------*/

Index Bundle::FindAPlace( void )
{
 // This method is used to return the index of an available position in the
 // Bundle to store a new item. The method implements the B-strategies of the
 // code, i.e. it can decide to discard some "old" subgradient if the Bundle
 // gets full, hence it must be used with care. The positions occupied by
 // newely obtained subgradients are "locked": if there are no possible
 // positions left, then InINF is returned.

 Index wh;

 if( FreDim )  // there are deleted items - - - - - - - - - - - - - - - - - -
 {
  if( ( wh = HeapDel( FreList , --FreDim ) ) >= SgMxNm )
   SgMxNm = wh + 1;
  }
 else
  if( SgMxNm < BPar2 )  // there is still space - - - - - - - - - - - - - - -
   wh = SgMxNm++;
  else             // there's no space: find an item to be eliminated - - - -
  {
   register cHpRow tA = ReadAlfa();  // find the removable item with largest
                                     // OOBase[] && smallest Alfa[]

   for( register Index i = wh = 0 ; ++i < BPar2 ; )
    if( ( OOBase[ i ] > OOBase[ wh ] ) ||
        ( ( OOBase[ i ] == OOBase[ wh ] ) && ( tA[ i ] > tA[ wh ] ) ) )
     wh = i;

   if( OOBase[ wh ] < 0 )  // exit if there is not any removable item
    return( InINF );

   if( ! OOBase[ wh ] )    // it is a basic item: re-select the exiting item
   {                       // as the one with the smallest value of Mult[]
    register Index h;
    register cIndex_Set tB = ReadBase();
    register cHpRow tM = ReadMult();
    register HpNum tMin = HpINF;

    for( wh = InINF ; ( h = *(tB++) ) < InINF ; tM++ )
     if( ( *tM < tMin ) && ( OOBase[ h ] >= 0 ) )
     {
      wh = h;
      tMin = *tM;
      }

    #if( DO_AGGR )
     if( whisZ == InINF )  // Z is not already in
     {
      Index wh2 = InINF;   // select item with the second smallest Mult[]

      for( tB = ReadBase() , tM = ReadMult() , tMin = HpINF ;
	   ( h = *(tB++) ) < InINF ; tM++ )
       if( ( h != wh ) && ( *tM < tMin ) && ( OOBase[ h ] >= 0 ) )
       {
        wh2 = h;
        tMin = *tM;
        }

      if( wh2 == InINF )   // exit if there is no space for Z
       return( InINF );
      else
      {
       AggregateZ( whisZ = wh2 );  // (construct Z and) put it in wh2
       OOBase[ whisZ ] = -1;       // Z won't be removed in this iteration
       }
      }
    #else
     BLOG( 1 , endl << " Warning: basic subgradient discarded. [FindAPlace]"
               << endl );
    #endif

    }  // end if( a Basic item )

   RmvItem( wh );

   }  // end else( no space ) - - - - - - - - - - - - - - - - - - - - - - - -

 OOBase[ wh ] = -1;

 return( wh );

 }  // end( FindAPlace )

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::Heuristic1( void )
{
 if( Alfa1 < HpeM )
 {
  if( DeltaFi > HpeM )
   return( tMaior );
  else
   return( tMinor );
  }
 else
  return( t * ( ( DeltaFi + Alfa1 ) / ( 2 * Alfa1 ) ) );
 }

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::Heuristic2( void )
{
 if( ABS( v - DeltaFi ) < HpeM )
  return( tMaior );
 else
  return( t * ( v / ( 2 * ( v - DeltaFi ) ) ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Bundle.C -----------------------------*/
/*--------------------------------------------------------------------------*/
