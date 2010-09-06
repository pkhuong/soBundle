/*--------------------------------------------------------------------------*/
/*--------------------------- File BMinQuad.C ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                							  --*/
/*-- Implementation of the BTT algorithm for solving the box-constrained  --*/
/*-- Quadratic Problems arising as descent direction finding subproblems  --*/
/*-- within (box)Constrained Bundle algorithms. Uses as a subroutine the  --*/
/*-- TT algorithm for linearly constrained problems without bounds.       --*/
/*--                							  --*/
/*--                            VERSION 3.20       			  --*/
/*--                           18 - 10 - 2004                             --*/
/*--                							  --*/
/*-- 		     Original Idea and Implementation by:		  --*/
/*--                                                                      --*/
/*--			      Antonio Frangioni        			  --*/
/*--                                                                      --*/
/*--   			   Operations Research Group			  --*/
/*--			  Dipartimento di Informatica			  --*/
/*--   			     Universita' di Pisa			  --*/
/*--                                                                      --*/
/*--                         Copyright 1992 - 2004                        --*/
/*--                							  --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BMinQuad.h"

#include "OPTvect.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--                                                                      --*/
/*--      Some small macro definitions, used throughout the code.         --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

#if( TWOSIDED )
 #define lb( h ) lb[ h ]
#else
 #define lb( h ) bounds[ h ]
#endif

#define size( x ) eDir * max( ABS( x ) , LMNum( 1 ) )

#if( ! SIGNAL_B2CHG )
 #define B2HasChgd()
#endif

#if( ! SIGNAL_MBCHG )
 #define MBHasChgd()
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace MinQuad_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/

static cHpNum MaxEpsD = 1e-3;   // Maximum, ...
static cHpNum IntEpsD = 1e-6;   // and Initial value for eD

static cHpNum MFactor = 10;     // eD is increased by multiplying
                                // it by this factor

/* Meaning of the bits in the char describing a variable:

   bit 0 ==> 0 if it is UC, 1 if it is NN;
   bit 1 ==> 0 if it is non-existing, 1 if it is existing;
   bit 2 ==> 0 if it is in MBase2, 1 if it is in Base2;
   bit 3 ==> 0 if it is in Base2 at the LB, 1 if it is in Base2 at the UB.
   bit 4 ==> 0 if it is "normal", 1 if it is "taboo". */

static const char kNNN   =  1;  // a Non-existing NN variable
static const char kEUC   =  2;  // an Existing UC variable (in MBase2)
static const char kENN   =  3;  // an Existing NN variable in MBase2
static const char kLBC   =  7;  // an Existing NN variable in Base2 at its LB
#if( TWOSIDED )
 static const char kUBC  = 15;  // an Existing NN variable in Base2 at its UB
#endif

static const char kIsNN  =  1;  // ( i &  1 ) <=> i is a NN variable
static const char kIsIn  =  2;  // ( i &  2 ) <=> i is exisiting
static const char kInB2  =  4;  // ( i &  4 ) <=> i is in Base2
#if( TWOSIDED )
 static const char kIsUB =  8;  // ( i &  8 ) <=> i is in Base2 at the UB
#endif
static const char kIsTb  = 16;  // ( i & 16 ) <=> i is "taboo"

/*--------------------------------------------------------------------------*/
/*-------------------- IMPLEMENTATION OF BMinQuad  -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR -- -----------------------------*/
/*--------------------------------------------------------------------------*/

BMinQuad::BMinQuad( void ) : MinQuad()
{
 #if( LOG_BMQ )
  BMQLog = &clog;

  #if( LOG_BMQ > 1 ) 
   BCalls = BSccss = 0;
   SumAverages = 0;
  #endif
 #endif

 #if( TIMERS_BMQ )
  BMQt = new OPTtimers();
 #endif

 }  // end( BMinQuad )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BMinQuad::SetMaxDim( Index m , Index n , Index SDim )
{
 if( MaxBDim )   // this is not the first call: - - - - - - - - - - - - - - -
  MemDealloc();  // deallocate everything - - - - - - - - - - - - - - - - - -

 MinQuad::SetMaxDim( m , n , SDim );  // "throw" the method of the base class

 if( MaxBDim )  // m != 0: allocate everything- - - - - - - - - - - - - - - -
 {              //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SpaceDim = SDim;

  RealAlfa = new HpNum[ m ];

  GS       = new char[ SpaceDim ];
  Base2    = new Index[ SpaceDim + 2 ];
  MBase2   = Base2 + SpaceDim + 1;

  di       = new LMNum[ SpaceDim ];
  #if( ! CNDVD_TMP )
   tmpdi   = new LMNum[ SpaceDim ];
  #endif
  bounds   = new LMNum[ SpaceDim ];

  #if( TWOSIDED )
   lb      = new LMNum[ SpaceDim ];
   ub      = new LMNum[ SpaceDim ];
  #endif

  for( register Index i = SpaceDim ; i-- ; )
  {
   GS[ i ] = kNNN;    // by default, all variables are NN
   lb( i ) = 0;       // with 0 LB
   #if( TWOSIDED )
    ub[ i ] = LMINF;  // and +INF UB
   #endif
   }

  *MBase2 = *Base2 = InINF;

  }  // end( if( m != 0 ) )- - - - - - - - - - - - - - - - - - - - - - - - - -

 // variables initialization - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NNStop = TLDim = MB2Dim = B2Dim = 0;
 MaxVarAdd = MaxVarRmv = InINF;
 eD = IntEpsD;
 Bf = HpINF;

 #if( ! BEXACT )
  bNorm = 0;
 #endif

 }  // end( BMinQuad::SetMaxDim )

/*--------------------------------------------------------------------------*/

void BMinQuad::SetEpsilonD( HpNum NeweD )
{
 eD = NeweD ? NeweD : IntEpsD;
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR ADDING / REMOVING / CHANGING DATA -------------*/
/*--------------------------------------------------------------------------*/

void BMinQuad::AddSubGrad( cIndex n , cHpNum alfan )
{
 Bf = HpINF;  // the objective function is allowed to increase

 RealAlfa[ n ] = alfan;
 MinQuad::AddSubGrad( n , alfan - GiTLB( n , bounds ) );

 }  // end( BMinQuad::AddSubGrad )

/*--------------------------------------------------------------------------*/

void BMinQuad::AddConstr( cIndex n , cHpNum alfan )
{
 Bf = HpINF;  // the objective function is allowed to increase

 RealAlfa[ n ] = alfan;
 MinQuad::AddConstr( n , alfan - GiTLB( n , bounds ) );

 }  // end( BMinQuad::AddConstr )

/*--------------------------------------------------------------------------*/

void BMinQuad::ChangeAlfa( cHpNum DeltaAlfa )
{
 for( register Index i = 0 ; i < NxtBIdx ; i++ )
  if( IsASubG( i ) )
   RealAlfa[ i ] += DeltaAlfa;

 MinQuad::ChangeAlfa( DeltaAlfa );

 }  // end( ChangeAlfa( HpNum , HpRow ) )

/*--------------------------------------------------------------------------*/

void BMinQuad::MoveAlongD( cHpNum Tau , cHpNum DeltaFi )
{
 // take care of "active" constraints - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 register cIndex_Set o = Base2;
 register Index h;

 if( Tau == PrvsTi )  //- - - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  for( ; ( h = *(o++) ) < InINF ; )
  {
   #if( TWOSIDED )
    register cLMNum dh = di[ h ] * PrvsTi;

    if( lb[ h ] > - LMINF )
     lb[ h ] += dh;

    if( ub[ h ] < + LMINF )
     ub[ h ] += dh;
   #endif

   bounds[ h ] = 0;
   }

  #if( ! BEXACT )
   bNorm = 0;
  #endif
  }
 else  // Tau != PrvsTi - - - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  register cHpNum tauR = 1 - Tau / PrvsTi; 

  for( ; ( h = *(o++) ) < InINF ; )
  {
   #if( TWOSIDED )
    register cLMNum dh = di[ h ] * Tau;

    if( lb[ h ] > - LMINF )
     lb[ h ] += dh;

    if( ub[ h ] < + LMINF )
     ub[ h ] += dh;
   #endif

   bounds[ h ] *= tauR;
   }

  #if( ! BEXACT )
   bNorm *= ( tauR * tauR );
  #endif

  }  // end else( Tau != PrvsTi ) - - - - - - - - - - - - - - - - - - - - - -

 // take care of "inactive" constraints - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( o = MBase2 ; ( h = *(o++) ) < InINF ; )
  #if( TWOSIDED )
  {
   register cLMNum dh = di[ h ] * Tau;

   if( ub[ h ] < + LMINF )
    ub[ h ] += dh;

   if( lb[ h ] > - LMINF )
    lb[ h ] += dh;
   }
  #else
   if( GS[ h ] & kIsNN )
    bounds[ h ] += di[ h ] * Tau;
  #endif

 // take care of Alfa[] - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Tau != PrvsTi )
 {
  VectSubtract( Alfa , RealAlfa , NxtBIdx );
  VectScale( Alfa , 1 - Tau / PrvsTi , NxtBIdx );
  }

 Swap( Alfa , RealAlfa );
 MinQuad::MoveAlongD( Tau , DeltaFi );
 Swap( Alfa , RealAlfa );

 if( Tau == PrvsTi )
  VectAssign( Alfa , RealAlfa , NxtBIdx );
 else
  VectSum( Alfa , RealAlfa , NxtBIdx );

 Bf = HpINF;  // the objective function is allowed to increase

 }  // end( BMinQuad::MoveAlongD )

/*--------------------------------------------------------------------------*/

void BMinQuad::ReadAlfa( HpRow NewAlfa )
{
 VectAssign( NewAlfa , RealAlfa , NxtBIdx );

 }  // end( BMinQuad::ReadAlfa( HpRow ) )

/*--------------------------------------------------------------------------*/

void BMinQuad::ChangeBounds( void )
{
 bNorm = 0;
 register Index i = 0;

 for( ; i < B2Dim ; i++ )
 {
  register cIndex h = Base2[ i ];

  #if( TWOSIDED )
   if( GS[ h ] & kIsUB )
    bounds[ h ] = ub[ h ];
   else
    bounds[ h ] = lb[ h ];

   if( ABS( bounds[ h ] ) == LMINF )
    CutOffConstr( i , h );
   else
  #endif
    bNorm += bounds[ h ] * bounds[ h ];
  }

 for( i = 0 ; i < ActBDim ; i++ )
  if( IsThere( i ) )
   Alfa[ i ] = RealAlfa[ i ] - GiTLB( i , bounds );

 AlfaChanged();

 Bf = HpINF;
 B2HasChgd();

 #if( LOG_BMQ > 3 )
  *BMQLog << "Modified bounds." << endl;
 #endif

 }  // end( ChangeBounds )

/*--------------------------------------------------------------------------*/

void BMinQuad::InitialSetUp( void )
{
 // set Base2[] vars- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( Alfa , RealAlfa , NxtBIdx );

 bNorm = 0;
 register const char *tGS = GS;
 register Index i = B2Dim = MB2Dim = 0;
 for( ; i < SpaceDim ; i++ )
 {
  register const char GSi = *(tGS++);
  if( ( GSi & kIsIn ) && ( GSi & kInB2 ) )
  {
   Base2[ B2Dim++ ] = i;

   #if( TWOSIDED )
    if( tGS[ i ] & kIsUB )
     bounds[ i ] = ub[ i ];
    else
     bounds[ i ] = lb[ i ];
   #endif

   register cHpNum lbi = bounds[ i ];
   if( ABS( lbi ) > HpeM )
   {
    bNorm += lbi * lbi;

    if( ActBDim )
    {
     register cSgRow NewDim = GiTilde( i );
     for( register Index j = 0 ; j < NxtBIdx ; j++ )
      if( IsThere( j ) )
       Alfa[ j ] -= NewDim[ j ] * lbi;
     }
    }
   }
  }

 Base2[ B2Dim ] = InINF;

 // set MBase2[] vars - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MBase2 = Base2 + SpaceDim + 1;
 *MBase2 = InINF;

 for( ; i-- ; )
 {
  register const char GSi = *(--tGS);
  if( ( GSi & kIsIn ) && ( ! ( GSi & kInB2 ) ) )
  {
   *(--MBase2) = i;
   MB2Dim++;
   }
  }

 // set NNStop- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( tGS = GS + ( i = SpaceDim ) ; i ; i-- )
  if( ( *(--tGS) & kENN ) == kENN )
   break;

 NNStop = i;

 // final operations- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 B2HasChgd();

 if( ActBDim )    // if there are items in the Bundle
 {
  AlfaChanged();  // signal that Alfa[] has changed
  ChangeQ();      // update Q and L in response to the change
  }

 Bf = HpINF;

 #if( LOG_BMQ > 3 )
  *BMQLog << "Created " << MB2Dim + B2Dim << " variables;" << endl
          << "Initial base of " << B2Dim << " variables built." << endl;
 #endif

 }  // end( InitialSetUp() )

/*--------------------------------------------------------------------------*/

void BMinQuad::AddVar( cIndex i )
{
 GS[ i ] |= kIsIn;
 if( ( GS[ i ] & kIsNN ) && ( i >= NNStop ) )
  NNStop = i + 1;

 if( GS[ i ] & kInB2 )
 {
  register cIndex h = BinSearch2( Base2 , B2Dim , i );
  ShiftRVect( Base2 + h , (++B2Dim) - h );
  Base2[ h ] = i;

  #if( TWOSIDED )
   if( GS[ i ] & kIsUB )
    bounds[ i ] = ub[ i ];
   else
    bounds[ i ] = lb[ i ];
  #endif

  #if( ! BEXACT )
   bNorm += bounds[ i ] * bounds[ i ];
  #endif

  AddBasicVariable( GiTilde( i ) , bounds[ i ] );
  }
 else
 {
  register cIndex h = BinSearch2( MBase2 , MB2Dim++ , i );
  ShiftVect( --MBase2 , h );
  MBase2[ h ] = i;

  AddSGSpaceDim( GiTilde( i ) );
  }

 B2HasChgd();
 Bf = HpINF;

 #if( LOG_BMQ > 3 )
  if( GS[ i ] & kIsNN )
  {
   #if( TWOSIDED )
    *BMQLog << "Created constrained variable " << i << " with bounds [";

    if( lb[ i ] > - LMINF )
     *BMQLog << lb[ i ] << ", "<< endl;
    else
     *BMQLog << "-INF, " << endl;

    if( ub[ i ] < LMINF )
     *BMQLog << ub[ i ] << "]." << endl;
    else
     *BMQLog << "+INF]." << endl;
   #else
    *BMQLog << "Created constrained variable " << i << " with lower bound "
            << bounds[ i ] << "." << endl;
   #endif
   }
  else
   *BMQLog << "Created unconstrained variable " << i << endl;
 #endif

 }  // end( AddVar )

/*--------------------------------------------------------------------------*/

void BMinQuad::RemoveVar( cIndex i )
{
 if( ! ( GS[ i ] & kIsIn ) )
  return;

 if( GS[ i ] & kInB2 )
 {
  register cIndex h = BinSearch1( Base2 , B2Dim , i );
  ShiftVect( Base2 + h , (B2Dim--) - h );

  #if( ! BEXACT )
   bNorm -= bounds[ i ] * bounds[ i ];
  #endif

  RemoveBasicVariable( GiTilde( i ) , bounds[ i ] );
  }
 else
 {
  register cIndex h = BinSearch1( MBase2 , MB2Dim-- , i );
  ShiftRVect( MBase2++ , h );

  CutSGSpaceDim( GiTilde( i ) );
  }

 B2HasChgd();
 Bf = HpINF;

 GS[ i ] &= kIsNN;  // keep the first bit
 while( NNStop && ( ( GS[ NNStop - 1 ] & kENN ) != kENN ) )
  NNStop--;

 #if( LOG_BMQ > 3 )
  *BMQLog << "Eliminated " << ( GS[ i ] & kIsNN ? "" : "un" )
          << "constrained variable " << i << "." << endl;
 #endif

 }  // end( RemoveVar )

/*--------------------------------------------------------------------------*/

void BMinQuad::MoveVar( cIndex i , cIndex j , cBOOL iIsLst )
{
 if( ( GS[ j ] & kIsIn ) || ( ! ( GS[ i ] & kIsIn ) ) )
  return;

 if( GS[ i ] & kInB2 )
 {
  register Index h = iIsLst ? B2Dim - 1 : BinSearch1( Base2 , B2Dim , i );
  register Index k = BinSearch2( Base2 , B2Dim , j );

  if( h < k )
   ShiftVect( Base2 + h , (--k) - h );
  else
   ShiftRVect( Base2 + k , h - k );

  Base2[ k ] = j;
  #if( TWOSIDED )
   bounds[ j ] = bounds[ i ];
  #endif
  }
 else
 {
  register Index h = iIsLst ? MB2Dim - 1 : BinSearch1( MBase2 , MB2Dim , i );
  register Index k = BinSearch2( MBase2 , MB2Dim , j );

  if( h < k )
   ShiftVect( MBase2 + h , (--k) - h );
  else
   ShiftRVect( MBase2 + k , h - k );

  MBase2[ k ] = j;
  }

 di[ j ] = di[ i ];
 lb( j ) = lb( i );
 #if( TWOSIDED )
  ub[ j ] = ub[ i ];
 #endif
 char tgsj = GS[ j ] & kIsNN;  // keep the first bit of j
 GS[ j ] = GS[ i ];
 GS[ i ] = tgsj;

 if( GS[ j ] & kIsNN )
  if( j >= NNStop )
   NNStop = j + 1;
  else
   while( NNStop && ( ( GS[ NNStop - 1 ] & kENN ) != kENN ) )
    NNStop--;

 B2HasChgd();

 #if( LOG_BMQ > 3 )
  *BMQLog << "Variable " << i << " is now renamed as " << j << "." << endl;
 #endif

 }  // end( MoveVar )

/*--------------------------------------------------------------------------*/

void BMinQuad::RenameVars( cIndex_Set whch )
{
 if( *whch == InINF )  // no variables to remove
  return;              // nothing to do

 register Index i = *whch;  // first variable to remove

 // find the position in Base2[] and MBase2[] - - - - - - - - - - - - - - - -

 register Index_Set tB2 = Base2;
 while( *tB2 < i )
  tB2++;

 register Index_Set tMB2 = MBase2;
 while( *tMB2 < i )
  tMB2++;

 // rename the variables- - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( register Index j = i ; j < SpaceDim ; j++ )
  if( *whch == j )
   whch++;
  else
  {
   register const char GSi = GS[ j ];
   GS[ i ] = GSi;

   if( GSi & kIsIn )
    if( GSi & kInB2 )
     *(tB2++) = i;
    else
     *(tMB2++) = i;

   di[ i ] = di[ j ];
   #if( TWOSIDED )
    lb[ i ] = lb[ j ];
    ub[ i ] = ub[ j ];
   #endif
   bounds[ i++ ] = bounds[ j ];
   }

 // now clean the bottom variables- - - - - - - - - - - - - - - - - - - - - -

 for( ; i < SpaceDim ; i++ )
 {
  GS[ i ] = kNNN;
  di[ i ] = 0;
  lb( i ) = 0;
  #if( TWOSIDED )
   ub[ i ] = LMINF;
  #endif
  }

 // update NNStop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 while( NNStop && ( ( GS[ NNStop - 1 ] & kENN ) != kENN ) )
  NNStop--;

 B2HasChgd();

 #if( LOG_BMQ > 3 )
  *BMQLog << "Variables renamed." << endl;
 #endif

 }  // end( RenameVars )

/*--------------------------------------------------------------------------*/

void BMinQuad::MakeVarCnstr( cIndex i )
{
 GS[ i ] |= kIsNN;
 if( i >= NNStop )
  NNStop = i + 1;

 if( GS[ i ] & kInB2 )
 {
  register Index h = BinSearch1( MBase2 , MB2Dim-- , i );
  ShiftRVect( MBase2++ , h );

  #if( TWOSIDED )
   if( GS[ i ] & kIsUB )
    bounds[ i ] = ub[ i ];
   else
    bounds[ i ] = lb[ i ];
  #endif

  #if( ! BEXACT )
   bNorm += bounds[ i ] * bounds[ i ];
  #endif

  CutSGSpaceDim( GiTilde( i ) , bounds[ i ] );

  h = BinSearch2( Base2 , B2Dim , i );
  ShiftRVect( Base2 + h , (++B2Dim) - h );
  Base2[ h ] = i;

  B2HasChgd();
  Bf = HpINF;
  }

 #if( LOG_BMQ > 3 )
  *BMQLog << "Added constraint ";

  #if( TWOSIDED )
   if( lb[ i ] > - LMINF )
    *BMQLog << "[" << lb[ i ] << ", ";
   else
    *BMQLog << "[-INF, ";

   if( ub[ i ] < LMINF )
    *BMQLog << ub[ i ] << "]";
   else
    *BMQLog << "+INF]";
  #else
   *BMQLog << " >= " << bounds[ i ];
  #endif

  *BMQLog << " on variable " << i << endl;
 #endif

 }  // end( MakeVarCnstr )

/*--------------------------------------------------------------------------*/

void BMinQuad::MakeVarUnCnstr( cIndex i )
{
 if( GS[ i ] & kInB2 )
 {
  register Index h = BinSearch1( Base2 , B2Dim , i );
  ShiftVect( Base2 + h , (B2Dim--) - h );

  h = BinSearch2( MBase2 , MB2Dim++ , i );
  ShiftVect( --MBase2 , h );
  MBase2[ h ] = i;

  AddSGSpaceDim( GiTilde( i ) , bounds[ i ] );

  GS[ i ] = kEUC;

  B2HasChgd();

  #if( ! BEXACT )
   bNorm -= bounds[ i ] * bounds[ i ];
  #endif
  }
 else
  GS[ i ] &= ~kIsNN;

 while( NNStop && ( ( GS[ NNStop - 1 ] & kENN ) != kENN ) )
  NNStop--;

 #if( LOG_BMQ > 3 )
  *BMQLog << "Variable " << i << " is now unconstrained." << endl;
 #endif

 }  // end( MakeVarUnCnstr )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

MinQuad::MQError BMinQuad::SolveQP( HpNum ti )
{
 #if( TIMERS_BMQ )
  BMQt->Start();
 #endif

 HpNum LastBf = HpINF;
 BOOL AugEps = FALSE;

 for(;;)  // problems handling loop - - - - - - - - - - - - - - - - - - - - -
 {
  CalcOptDir( ti );

  if( ( ! QPStatus ) || ( QPStatus == kFatal ) ||
      ( QPStatus == kQPPrimUnbndd ) )
   break;

  ClearTabooList();

  if( Bf < HpINF )
   if( LastBf > Bf + eD * ( ABS( Bf ) / SpaceDim ) )
   {
    LastBf = Bf;     // reset the counters if there have been
    AugEps = FALSE;  // at least one strictly decreasing step
    }

  if( ! AugEps )     // the first time (after a decrease)
  {
   ChangeQ();        // refresh the data structures
   AugEps = TRUE;    // signal that it has already been done
   continue;         // and retry;
   }

  eD *= MFactor;     // eD increase (decrease precision)

  if( eD > MaxEpsD )
  { 
   #if( LOG_BMQ )
    *BMQLog << endl << "ERROR: eD (" << eD << ") out of limits. [Solve[B]QP]" 
            << endl;
   #endif

   QPStatus = kFatal;
   break;
   }
  }  // end while( problems ) - - - - - - - - - - - - - - - - - - - - - - - -

 #if( TIMERS_BMQ )
  BMQt->Stop();
 #endif

 return( QPStatus );
 
 }  // end( BMinQuad::SolveQP )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR READING RESULTS ---------------------*/
/*--------------------------------------------------------------------------*/

void BMinQuad::ReadZ( register LMRow g )
{
 register cIndex_Set o = Base2;
 register Index h;

 for( ; ( h = *(o++) ) < InINF ; )
  g[ h ] = di[ h ];

 for( o = MBase2 ; ( h = *(o++) ) < InINF ; )
  g[ h ] = di[ h ];

 }  // end( ReadZ )

/*--------------------------------------------------------------------------*/

void BMinQuad::ReadZSprs( register LMRow g )
{
 register cIndex_Set tB2 = Base2;
 register cIndex_Set tMB2 = MBase2;
 register Index h = *tB2;
 register Index k = *tMB2;

 while( h != k )
  if( h < k )
  {
   *(g++) = di[ h ];
   h = *(++tB2);
   }
  else
  {
   *(g++) = di[ k ];
   k = *(++tMB2);
   }
 }  // end( ReadZSprs )

/*--------------------------------------------------------------------------*/

Index BMinQuad::ReadVNames( register Index_Set VNames )
{
 register cIndex_Set tB2 = Base2;
 register cIndex_Set tMB2 = MBase2;
 register Index h = *tB2;
 register Index k = *tMB2;

 while( h != k )
  if( h < k )
  {
   *(VNames++) = h;
   h = *(++tB2);
   }
  else
  {
   *(VNames++) = k;
   k = *(++tMB2);
   }

 return( B2Dim + MB2Dim );

 }  // end( ReadVNames )

/*--------------------------------------------------------------------------*/

void BMinQuad::ReadD( register LMRow d ,  cIndex CpyFrst )
{
 register cIndex_Set o = Base2;
 register Index h;

 if( CpyFrst )
 {
  h = *(o++);
  for( register Index i = 0 ; i < CpyFrst ; i++ )
   if( i == h )
   {
    d[ i ] = bounds[ h ];
    h = *(o++);
    }
   else
    d[ i ] = - PrvsTi * di[ i ];
  }
 else
 {
  for( ; ( h = *(o++) ) < InINF ; )
   d[ h ] = bounds[ h ];

  for( o = MBase2 ; ( h = *(o++) ) < InINF ; )
   d[ h ] = - PrvsTi * di[ h ];
  }
 }  // end( ReadD )

/*--------------------------------------------------------------------------*/

void BMinQuad::ReadDSprs( register LMRow d )
{
 register cIndex_Set tB2 = Base2;
 register cIndex_Set tMB2 = MBase2;
 register Index h = *tB2;
 register Index k = *tMB2;

 while( h != k )
  if( h < k )
  {
   *(d++) = bounds[ h ];
   h = *(++tB2);
   }
  else
  {
   *(d++) = - PrvsTi * di[ k ];
   k = *(++tMB2);
   }
 }  // end( ReadDSprs )

/*--------------------------------------------------------------------------*/

HpNum BMinQuad::DPerG( register cSgRow g )
{
 register HpNum scpr = 0;
 register cIndex_Set o = Base2;
 register Index h;

 for( ; ( h = *(o++) ) < InINF ; )
  scpr -= g[ h ] * bounds[ h ];

 scpr /= PrvsTi;

 for( o = MBase2 ; ( h = *(o++) ) < InINF ; )
  scpr += g[ h ] * di[ h ];

 return( scpr );

 }  // end( DPerG )

/*--------------------------------------------------------------------------*/

void BMinQuad::AddD( register LMRow L1 , register cLMRow L2 , cHpNum Tau )
{
 register cIndex_Set o = Base2;
 register Index h;

 if( Tau == PrvsTi )
  for( ; ( h = *(o++) ) < InINF ; )
   L1[ h ] = L2[ h ] + bounds[ h ];
 else
 {
  register HpNum TauR = Tau / PrvsTi;

  for( ; ( h = *(o++) ) < InINF ; )
   L1[ h ] = L2[ h ] + TauR * bounds[ h ];
  }

 for( o = MBase2 ; ( h = *(o++) ) < InINF ; )
  L1[ h ] = L2[ h ] - Tau * di[ h ];

 }  // end( AddD )

/*--------------------------------------------------------------------------*/

void BMinQuad::AddDSprs( register LMRow L1 , register cLMRow L2 , cHpNum Tau )
{
 register cIndex_Set tB2 = Base2;
 register cIndex_Set tMB2 = MBase2;
 register Index h = *tB2;
 register Index k = *tMB2;

 if( Tau == PrvsTi )
 {
  while( h != k )
   if( h < k )
   {
    *(L1++) = *(L2++) + bounds[ h ];
    h = *(++tB2);
    }
   else
   {
    *(L1++) = *(L2++) - Tau * di[ k ];
    k = *(++tMB2);
    }
   }
 else
 {
  register HpNum TauR = Tau / PrvsTi;

  while( h != k )
   if( h < k )
   {
    *(L1++) = *(L2++) + TauR * bounds[ h ];
    h = *(++tB2);
    }
   else
   {
    *(L1++) = *(L2++) - Tau * di[ k ];
    k = *(++tMB2);
    }
  }
 }  // end( AddDSprs )

/*--------------------------------------------------------------------------*/

void BMinQuad::SensitAnals1( HpNum &v1 , HpNum &v2 , HpNum &v3 )
{
 MinQuad::SensitAnals1( v1 , v2 , v3 );
 v3 += bNorm;
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

BMinQuad::~BMinQuad( void )
{
 // First: print out information - - - - - - - - - - - - - - - - - - - - - -

 #if( LOG_BMQ > 2 )
  *BMQLog << endl << "BCalls " << BCalls << " (faults " << BCalls - BSccss
          << ") ~ " << SpaceDim << " variables, " << SumAverages / BCalls
	  << " avg. constr.";
  #if( TIMERS_BMQ )
   *BMQLog << endl << "Solve[B]QP() time: user  " << BMQt->u << ", system "
           << BMQt->s << " seconds." << endl;
  #endif
 #elif( LOG_BMQ > 1 )
  *BMQLog << SpaceDim << "\t" << BCalls << "\t" << SumAverages / BCalls
	  << "\t";
  #if( TIMERS_BMQ )
   *BMQLog  << ( BMQt->u + BMQt->s ) << "\t";
  #endif
 #endif

 // Second: memory deallocation - - - - - - - - - - - - - - - - - - - - - - -

 #if( TIMERS_BMQ )
  delete BMQt;
 #endif

 if( MaxBDim )
  MemDealloc();

 }  // end( ~BMinQuad )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

inline void BMinQuad::CalcOptDir( HpNum ti )

/* Same as in the base class, but taking into account the box constraints (in
   a "smart" way): calls MinQuad::SolveQP() as a subroutine. */
{
 #if( LOG_BMQ > 1 )
  BCalls++;
  float TmpSumAverages = 0;
  unsigned long int BStep = 0;
 #endif

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /*- - - - - - - - - - - - - - Main Loop - - - - - - - - - - - - - - - - - -*/
 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 HpNum BfBot = HpINF;  // "security level", used to switch between the "quick
                       // and dirty" pricing and the "slow but clean" one
 tmpMVarAdd = MaxVarAdd;  // used to switch between the "relaxed" pricing in
                          // (many variables can be priced in in the same
                          // iteration) and the "strict" one (only one at a
                          // time) that is used with the "taboo" mechanism to
                          // avoid loops
 Index Moved;

 do
 {
  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /*- - - - - - - - - - - - - - Inner Loop  - - - - - - - - - - - - - - - -*/
  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Within this inner loop, non-violated box constraints are removed (in
    other words, negative Xsi[]'s are eliminated).
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  do
  {
   #if( LOG_BMQ > 1 )
    #if( LOG_BMQ > 3 )
     *BMQLog << endl << "[" << BCalls << "/" << BStep++ << "-" << B2Dim
	     << "]:";

     #if( LOG_BMQ > 4 )
      CheckDS();  // debug checking of data structures
     #endif
    #endif

    TmpSumAverages += B2Dim;
   #endif

   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // the "usual" (QP) is solved- - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( MinQuad::SolveQP( ti ) )
   {
    #if( LOG_BMQ > 1 )
     SumAverages += TmpSumAverages / BStep;
    #endif

    return;  // errors are just passed up
    }

   MBHasChgd();

   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // feasibility check- - - - - -- - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Moved = 0;

   if( BfBot < HpINF )  // "secure" pricing - - - - - - - - - - - - - - - - -
   {                    //- - - - - - - - - - - - - - - - - - - - - - - - - -
    #if( LAZY_D == 0 )
     CalculateZ( tmpdi );
    #elif( LAZY_D == 1 )
     CalculateZ( Base2 , tmpdi );
    #endif

    // check if tmpdi is feasible: if not, find the - - - - - - - - - - - - -
    // maximum feasible step along tmpd - di

    register LMNum step = LMINF;
    register cHpNum eDir = eD * max( BDim , Index( 1 ) );

    register cIndex_Set tB2 = Base2;
    for( register Index h ; ( h = *(tB2++) ) < InINF ; )
    {
     #if( LAZY_D == 2 )
      tmpdi[ h ] = CalculateZ( h );
     #endif

     register cLMNum d1h = - ti * tmpdi[ h ];

     #if( TWOSIDED )
      if( GS[ h ] & kIsUB )
      {
       if( d1h + size( d1h ) > ub[ h ] )
	continue;
       }
      else
     #endif
       if( d1h - size( d1h ) < lb( h ) )
	continue;

     register cLMNum dh = ti * di[ h ];
     step = min( step , ( dh + bounds[ h ] ) / ( dh + d1h ) );
     }

    if( step == LMINF )  // tmpdi is feasible - - - - - - - - - - - - - - - -
     Swap( di , tmpdi );
    else                 // move along tmpd - di of step- - - - - - - - - - -
    {
     register cIndex_Set tB2 = Base2;
     for( register Index h ; ( h = *tB2 ) < InINF ; )
     {
      di[ h ] += step * ( tmpdi[ h ] - di[ h ] );

      if( ABS( ti * di[ h ] + bounds[ h ] ) > size( bounds[ h ] ) )
       tB2++;
      else
      {
       #if( LOG_BMQ > 3 )
        *BMQLog << " Exiting constraint, d[ " << h << " ] = "
		<< - ti * di[ h ] << endl;
       #endif

       CutOffConstr( tB2 - Base2 , h );

       #if( ! BEXACT )
        bNorm -= bounds[ h ] * bounds[ h ];
       #endif

       Moved++;
       }
      }  // end for( tB2 )

     if( ! Moved )  // haven't found any constraint to eliminate- - - - - - -
     {
      #if( LOG_BMQ )
       *BMQLog << "Fault: unable to reach the boundary. [BCalcOptDir]"
	       << endl;
      #endif

      QPStatus = kLoop;
      return;
      }
     }   // end else( move along tmpdi - di )
    }    // end if( "secure" pricing )
   else  // "quick and dirty" pricing - - - - - - - - - - - - - - - - - - - -
   {     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #if( LAZY_D == 0 )
     CalculateZ( di );
    #elif( LAZY_D == 1 )
     CalculateZ( Base2 , di );
    #endif

    register cHpNum eDir = eD * max( BDim , Index( 1 ) );
    register cIndex_Set tB2 = Base2;

    for( register Index h ; ( h = *tB2 ) < InINF ; )
    {
     #if( LAZY_D == 2 )
      di[ h ] = CalculateZ( h );
     #endif

     register cLMNum dh = - ti * di[ h ];

     #if( TWOSIDED )
      if( GS[ h ] & kIsUB )
      {
       if( dh + size( dh ) > ub[ h ] )
       {
	tB2++;
        continue;
        }
       }
      else
     #endif
       if( dh - size( dh ) < lb( h ) )
       {
	tB2++;
        continue;
        }

     #if( LOG_BMQ > 3 )
      #if( TWOSIDED )
       if( GS[ h ] & kIsUB )
        *BMQLog << " Exiting UB constraint, d[ " << h << " ] = " << dh
	        << " < " << ub[ h ] << endl;
       else
      #else
        *BMQLog << " Exiting LB constraint, d[ " << h << " ] = " << dh
	        << " > " << lb( h ) << endl;
      #endif
     #endif
 
     #if( ! BEXACT )
      bNorm -= bounds[ h ] * bounds[ h ];
     #endif

     CutOffConstr( tB2 - Base2 , h );

     if( ++Moved >= MaxVarRmv )
      break;
     }
    }  // end else( "quick and dirty" pricing ) - - - - - - - - - - - - - - -
       // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   #if( SIGNAL_B2CHG )
    if( Moved )
     B2HasChgd();
   #endif

   } while( Moved );

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /*- - - - - - - - - - - - - End Inner Loop  - - - - - - - - - - - - - - -*/
  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  #if( BEXACT )
   bNorm = Norm( bounds , Base2 );
  #else
   if( ! B2Dim )
    bNorm = 0;   // from time to time, clean up bNorm
  #endif

  #if( LOG_BMQ > 3 )
   *BMQLog << " bNorm = " << bNorm << " ~";
  #endif

  HpNum newBf = ( Quad - bNorm / ( ti * ti ) ) / 2 + Lin / ti;

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    The actual base is feasible, so calculate Bf and test for decrease.
    If the "quick and dirty" rule is used, the decrease is *not guaranteed*
    and this control may fail: in this case, the "secure" rule will be used
    afterwards until a strict decrease (w.r.t. the "old" Bf at *this*
    iteration) is obtained.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  if( Bf < HpINF )  // skipping the tests if Bf is undefined
  {
   // Test the relative decrease of B by checking ( Bf - newBf ) / Bf against
   // eD / SpaceDim; the idea is that each of the SpaceDim variables can make
   // Bf to decrease of at most 1 / SpaceDim of its current value

   if( Bf - newBf < - eD * ( ABS( Bf ) / SpaceDim ) )  // an increasing step
    BfBot = min( BfBot , Bf );  // just start using the "secure" pricing

   if( Bf - newBf <= eD * ( ABS( Bf ) / SpaceDim ) )   // a "short" step
    if( tmpMVarAdd > 1 )  // first time:
     tmpMVarAdd = 1;      // just start using the "strict" pricing
    else                  // any subsequent time: mark Entrd as "taboo"
    {
     GS[ Entrd ] |= kIsTb;
     TLDim++;
     }
   else                                                // a "good step"
    if( newBf + eD * ( ABS( newBf ) / SpaceDim ) < BfBot )
    {
     BfBot = HpINF;           // it is safe to disable the "secure" rule ...
     tmpMVarAdd = MaxVarAdd;  // ... and also the "strict" pricing
     ClearTabooList();        // ... and also to clear the taboo list
     }
   }  // end if( Bf < HpINF )

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    If Bf is "small enough" we can detect optimality without the (costly)
    check of constraints violation: however, the entries of d[] in MBase2
    have to be "artificially" set to 0.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  if( ( Bf = newBf ) <= MinFVal )
  {
   register cIndex_Set o = MBase2;
   for( register Index h ; ( h = *(o++) ) < InINF ; )
    di[ h ] = 0;

   break;
   }

  #if( LOG_BMQ > 3 )
   *BMQLog << " Bf = " << Bf << endl;
  #endif

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Now check for the existence of violated (primal) constraints to be put in
    the "constraints base" to obtain a decrease of Bf; in other words, check
    for any Xsi[ h ] == 0 that, if let free to become > 0, can diminish the
    value of Bf.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  #if( LAZY_D == 1 )
   CalculateZ( MBase2 , di );
  #endif

  register cHpNum eDir = eD * max( BDim , Index( 1 ) );
  register cIndex_Set tMB2 = MBase2;

  for( register Index h ; ( h = *tMB2 ) < NNStop ; tMB2++ )
   if( GS[ h ] == kENN )  // only existing NN variables in MBase2 which are
   {                      // not "taboo" can enter Base2
    #if( LAZY_D == 2 )
     di[ h ] = CalculateZ( h );
    #endif

    register cLMNum dh = - ti * di[ h ];
    if( dh + size( dh ) < lb( h ) )
    {
     #if( LOG_BMQ > 3 )
      *BMQLog << " Entering violated LB constraint, d[ " << h << " ] = "
	      << dh << " < " << lb( h ) << endl;
     #endif

     GS[ h ] = kLBC;
    #if( TWOSIDED )
     bounds[ h ] = lb[ h ];
     }
    else
     if( dh - size( dh ) > ub[ h ] )
     {
      #if( LOG_BMQ > 3 )
       *BMQLog << " Entering violated UB constraint, d[ " << h << " ] = "
	       << dh << " >= " << ub[ h ] << endl;
      #endif

      GS[ h ] = kUBC;
      bounds[ h ] = ub[ h ];
    #endif
      }
     else
      continue;

    Index i = tMB2 - MBase2;
    ShiftRVect( MBase2++ , i );
    MB2Dim--;

    i = BinSearch2( Base2 , B2Dim , h );
    ShiftRVect( Base2 + i , (++B2Dim) - i );
    Base2[ i ] = Entrd = h;

    #if( ! BEXACT )
     bNorm += bounds[ h ] * bounds[ h ];
    #endif

    CutSGSpaceDim( GiTilde( h ) , bounds[ h ] );

    if( ++Moved >= tmpMVarAdd )
     break;

    }  // end for( h )

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  #if( SIGNAL_B2CHG )
   if( Moved )
    B2HasChgd();
  #endif

  if( ( ! Moved ) && ( Bf - eD * ( ABS( Bf ) / SpaceDim ) > newBf ) )
  {
   // if no variables enters in Base2 (some because they are "taboo"),
   // but Bf() is not (about) as good as BfBot, then kLoop is returned

   #if( LOG_BMQ > 2 )
    *BMQLog << endl << "Fault: Bf() can't reach the min. value " << BfBot
	    << " [BCalcOptDir]";
   #endif

   QPStatus = kLoop;
   return;
   }
  } while( Moved );

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /*- - - - - - - - - - - - - End Main Loop - - - - - - - - - - - - - - - - -*/
 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ClearTabooList();  // clear the taboo list

 for( register Index i = 0 ; i < NxtBIdx ; i++ )
  GTd[ i ] -= ( RealAlfa[ i ] - Alfa[ i ] ) / ti;  // correct GTd[]

 #if( LOG_BMQ > 1 )
  SumAverages += TmpSumAverages / BStep;
  BSccss++;

  #if( LOG_BMQ > 3 )
   *BMQLog << " Stop: BOptimum reached ~ Bf = " << Bf << endl << endl;
  #endif
 #endif

 }  // end( BMinQuad::CalcOptDir )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::CutOffConstr( Index i , Index h )

/* Removes the basic constraints on the h = Base2[ i ] variable.
   Doesn't handle bNorm. */
{
 ShiftVect( Base2 + i , (B2Dim--) - i );

 i = BinSearch2( MBase2 , MB2Dim++ , h );
 ShiftVect( --MBase2 , i );
 MBase2[ i ] = h;

 GS[ h ] = kENN;

 AddSGSpaceDim( GiTilde( h ) , bounds[ h ] );
 }

/*--------------------------------------------------------------------------*/

inline void BMinQuad::AddBasicVariable( cSgRow NewDim , cHpNum lbh )

/* This is called when a basic variable has to be created: in this case,
   nothing has to be done to Q and L, but Alfa (not RealAlfa) has to be
   changed accordingly to incorporate the corresponding row of GXsiLXsi. */
{
 if( ABS( lbh ) > HpeM )
  for( register Index i = 0 ; i < NxtBIdx ; i++ )
   if( IsThere( i ) )
    Alfa[ i ] -= NewDim[ i ] * lbh;

 AlfaChanged();  // this has to be called even if lbh == 0

 }  // end( AddBasicVariable )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::RemoveBasicVariable( cSgRow OldDim , cHpNum lbh )

/* This is called when a basic variable has to be eliminated: in this case,
   nothing has to be done to Q and L, but Alfa (not RealAlfa) has to be
   changed accordingly to remove the corresponding row of GXsiLXsi. */
{
 if( ABS( lbh ) > HpeM )
  for( register Index i = 0 ; i < NxtBIdx ; i++ )
   if( IsThere( i ) )
    Alfa[ i ] += OldDim[ i ] * lbh;

 AlfaChanged();  // this has to be called even if lbh == 0

 }  // end( RemoveBasicVariable )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::ClearTabooList( void )
{
 if( TLDim )
 {
  register cIndex_Set tB = MBase2;
  register Index h;

  for( ; ( h = *(tB++) ) < InINF ; )
   if( GS[ h ] & kIsTb )
   {
    GS[ h ] &= ~kIsTb;
    if( ! --TLDim )
     break;
    }

  if( TLDim )
   for( tB = Base2 ; ( h = *(tB++) ) < InINF ; )
    if( GS[ h ] & kIsTb )
    {
     GS[ h ] &= ~kIsTb;
     if( ! --TLDim )
      break;
     }
  }
 }  // end( ClearTabooList )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::MemDealloc( void )
{
 #if( TWOSIDED )
  delete[] ub;
  delete[] lb;
 #endif

 delete[] bounds;
 #if( ! CNDVD_TMP )
  delete[] tmpdi;
 #endif
 delete[] di;

 delete[] Base2;
 delete[] GS;

 delete[] RealAlfa;

 }  // end( MemDealloc )

/*--------------------------------------------------------------------------*/

#if( LOG_BMQ > 4 )

inline void BMinQuad::CheckDS( void )
{
 // check Base2[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 assert( Base2[ B2Dim ] == InINF );

 if( B2Dim )
 {
  register cIndex_Set tB2 = Base2;
  register Index dim = B2Dim - 1;
  register Index i = *(tB2++);
  for( register Index h ; ( h = *(tB2++) ) < InINF ; )
  {
   assert( i < h );
   i = h;
   dim--;
   }

  assert( ! dim );
  }

 // check MBase2[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 assert( MBase2[ MB2Dim ] == InINF );

 if( MB2Dim )
 {
  register cIndex_Set tMB2 = MBase2;
  register Index dim = MB2Dim - 1;
  register Index i = *(tMB2++);
  for( register Index h ; ( h = *(tMB2++) ) < InINF ; )
  {
   assert( i < h );
   i = h;
   dim--;
   }

  assert( ! dim );
  }

 // check GS[], Base2[] and MBase2[] - - - - - - - - - - - - - - - - - - - - -

 register cIndex_Set tB2 = Base2;
 register cIndex_Set tMB2 = MBase2;
 for( register Index i = 0 ; i < SpaceDim ; i++ )
  if( GS[ i ] & kIsIn )         // an existing variable
   if( GS[ i ] & kInB2 )        // in Base2[]
    assert( i == *(tB2++) );
   else
    assert( i == *(tMB2++) );  // in MBase2[]
  else                         // a non-existing variable
   assert( ( *tB2 != i ) && ( *tMB2 != i ) );

 // check RealAlfa[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( register Index i = 0 ; i < NxtBIdx ; i++ )
  if( IsThere( i ) )
  {
   cHpNum xy = RealAlfa[ i ];
   cHpNum xx = xy - Alfa[ i ] - GiTLB( i , bounds );

   if( CMP( xx , xy ) )
    *BMQLog << endl << "Test failed: RealAlfa[" << i << "] = " << xx;
   }

 }  // end( CheckDS )

#endif

/*--------------------------------------------------------------------------*/
/*---------------------- End File BMinQuad.C -------------------------------*/
/*--------------------------------------------------------------------------*/
