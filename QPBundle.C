/*--------------------------------------------------------------------------*/
/*--------------------------- File QPBundle.C ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--   Implementation of the (QP) direction finding subproblem for the    --*/
/*-- Bundle class using the specialized (QP) solvers implemented in the   --*/
/*-- [B]MinQuad class (refere to [B]MinQuad.h and references therein).    --*/
/*--                                                                      --*/
/*--                            VERSION 1.00                              --*/
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
/*--                         Copyright 1994 - 2004                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "QPBundle.h"

#include "OPTvect.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( ( G_IMPLM < 1 ) || ( G_IMPLM > 2 ) )

 #define HAVE_BTH 0

 // HAVE_BTH == 1 <=> both the row-wise and the column-wise representation
 // of G are available

 #define WHAT_4_Q G_IMPLM

 // WHAT_4_Q == 0 => the column-wise representation of G is used for
 // calculating the entries of Q[][], that is the (restricted) scalar products
 // between items; if WHAT_4_Q > 0 the row-wise implementation is used instead

 #define WHAT_4_D G_IMPLM

 // WHAT_4_D == 0 => the column-wise representation of G is used for
 // calculating the entries of d[], that is the convex combinations of entries
 // of items; if WHAT_4_Q > 0 the row-wise implementation is used instead

#else

 #define HAVE_BTH 1

 #if( G_IMPLM == 1 )
  #define WHAT_4_Q 0
  #define WHAT_4_D 1
 #else
  #define WHAT_4_Q 1
  #define WHAT_4_D 0
 #endif

#endif

#if( G_IMPLM > 0 )
 #define TSGi( i ) TSubG[ i ]
#endif

// used in AggregateZ()- - - - - - - - - - - - - - - - - - - - - - - - - - -

#if( ( DO_AGGR == 1 ) && ( ! ADD_CNST ) )
 #define ZBase MinQuad::ReadBase()
 #define ZMult MinQuad::ReadMult()
 #define ZBDim MinQuad::ReadBDim()
#endif

// used in SolveSubP() - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if( ! HV_NNVAR )
 #define SettmpD( d )
 #define GettmpD( d )
#else
 #if( ! CNDVD_TMP )
  #define SettmpD( d )
  #define GettmpD( d )
 #else
  #define GettmpD( d ) d = GettmpD()
 #endif
#endif

// used in GiTG[j]()- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if( HV_NNVAR )
 #define SprsV PPar2 || NNVars
#else
 #define MBase2 LamBase
 #define SprsV PPar2
#endif

// used in SetItem()- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define size( x ) eDir * max( ABS( x ) , LMNum( 1 ) )

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
using namespace Bundle_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  These functions are not implemented as methods of the class, since  --*/
/*--  they don't use directly its data structures.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

static inline QuNum ScalarProductBB( register cSgRow g1 , register cSgRow g2 ,
				     register cIndex_Set B )
{
 // ScalarProduct( g1{B} , g2{B} )

 register QuNum t = 0;
 for( register Index h ; ( h = *(B++) ) < InINF ; )
  t += g1[ h ] * g2[ h ];

 return( t );
 }

/*--------------------------------------------------------------------------*/

static inline void VectAssignBB( register LMRow g1 , register cSgRow g2 ,
				register cHpNum x , register cIndex_Set B )
{
 // g1{B} = x * g2{B} (element-wise)

 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] = x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/

static inline void VectSumBB( register LMRow g1 , register cSgRow g2 , 
			      register cHpNum x , register cIndex_Set B )
{
 // g1{B} += x * g2{B} (element-wise)

 for( register Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] += x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF QPBundle -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

QPBundle::QPBundle( istream *iStrm , cIndex NV , cIndex ENV )
          :
          Bundle( iStrm , NV , ENV ) ,
          #if( HV_NNVAR )
	   BMinQuad()
          #else
           MinQuad()
          #endif
{
 // memory allocation - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SetMaxDim( BPar2 , min( SltSize , BPar2 ) , NV );

 #if( ! HV_NNVAR )
  dir = new LMNum[ MaxNumVar ];
 #endif

 #if( G_IMPLM < 1 )
  TSGk = new SgNum[ BPar2 ];
 #endif
 
 #if( G_IMPLM > 0 )
  TSubG = new SgRow[ MaxNumVar ];
  VectAssign( TSubG , SgRow( NULL ) , MaxNumVar );

  #if( HAVE_BTH )
   BFCntr = 0;
   Buffer = NULL;

   if( PPar2 )
   {
    Buffer = new SgRow[ MaxNumVar ];  // allocate the buffer of rows

    for( Index i = BFCntr = min( GTSltSize , MaxNumVar ) ; i-- ; )
     Buffer[ i ] = new SgNum[ BPar2 ];

    NrAllRws = 0;
    }
   else
  #endif
   {
    // allocate the minimum necessary number of "chunks" of GTSltSize rows,
    // i.e. GTSltSize * ceil( NumVar / GTSltSize )

    NrAllRws = min( MaxNumVar , Index( GTSltSize * ( NumVar / GTSltSize ) +
				( NumVar % GTSltSize ? GTSltSize : 0 ) ) );

    for( Index i = 0 ; i < NrAllRws ; )
     TSubG[ i++ ] = new SgNum[ BPar2 ];
    }
 #endif

 DimMinQuad = min( BPar2 , SltSize );

 #if( G_IMPLM >= 3 )
  tmpG1 = new SgNum[ MaxNumVar ];
 #else
  SubG = new SgRow[ BPar2 ];

  for( Index i = 0 ; i < DimMinQuad ; )
   SubG[ i++ ] = new SgNum[ MaxNumVar ];
 #endif

 #if( DO_AGGR && ( ADD_CNST || ( DO_AGGR == 2 ) ) )
  ZBase = new Index[ BPar2 + 1 ];
  ZMult = new HpNum[ BPar2 ];
  ZBDim = 0;
 #endif

 // set up [B]MinQuad data- - - - - - - - - - - - - - - - - - - - - - - - - -

 SetPricing( ( MPar2 < 0 ? HpINF : MPar2 ) );

 if( ! PPar2 )
 {
  PPar1 = PPar3 = 0;

  #if( HV_NNVAR )
   register char *vs = GetVars();
   register const char vsi = NNVar() | IsVar() | AcVar();
   for( Index i = NumVar ; i-- ; )
    *(vs++) = vsi;

   InitialSetUp();
  #endif
  }

 // set up scalar fields- - - - - - - - - - - - - - - - - - - - - - - - - - -

 RHSp = NULL;
 ZCmptd = 0;
 MinNewAlfa = HpINF;

 #if( DO_AGGR && ADD_CNST )
  CBZDim = 0;
 #endif

 }  // end( QPBundle )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void QPBundle::SetUC( cIndex Strt , cIndex Stp , cIndex_Set Which )
{
 Bundle::SetUC( Strt , Stp , Which );  // "throw" the method of the base class

 #if( HV_NNVAR )
  register const char nnv = NNVar() | AcVar();
  register char *const vs = GetVars();
  register const char isv = IsVar();
  register Index i = 0;
  for( ; i < Strt ; )
  {
   vs[ i ] &= isv;
   vs[ i++ ] |= nnv;
   }

  for( ; i < Stp ; )
   vs[ i++ ] &= isv;

  for( ; i < NumVar ; )
  {
   vs[ i ] &= isv;
   vs[ i++ ] |= nnv;
   }

  if( Which )
   for( register Index h ; ( h = *(Which++) ) < InINF ; )
    vs[ h ] &= isv;

  if( ! PPar2 )
   InitialSetUp();
 #else
  assert( ! NNVars );
 #endif

 }  // end( QPBundle::SetUC )

/*--------------------------------------------------------------------------*/

void QPBundle::SetRHS( cSgRow tRHS )
{
 if( tRHS )
 {
  RHSp = new SgNum[ NumVar ];
  VectAssign( RHSp , tRHS , NumVar );
  }
 else
 {
  delete[] RHSp;
  RHSp = NULL;
  }
 }

/*--------------------------------------------------------------------------*/

void QPBundle::GetRHS( SgRow tRHS )
{
 if( RHSp )
  VectAssign( tRHS , RHSp , NumVar );
 else
  VectAssign( tRHS , SgNum( 0 ) , NumVar );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

QPBundle::~QPBundle()
{
 delete[] RHSp;

 #if( DO_AGGR && ( ADD_CNST || ( DO_AGGR == 2 ) ) )
  delete[] ZMult;
  delete[] ZBase;
 #endif

 #if( G_IMPLM < 3 )
  for( ; DimMinQuad-- ; )
   delete[] SubG[ DimMinQuad ];

  delete[] SubG;
 #else
  delete[] tmpG1;
 #endif

 #if( G_IMPLM > 0 )
  #if( G_IMPLM < 3 )
   for( ; BFCntr-- ; )
    delete[] Buffer[ BFCntr ];

   delete[] Buffer;
  #endif

  for( Index i = MaxNumVar ; i-- ; )
   delete[] TSubG[ i ];

  delete[] TSubG;
 #endif

 #if( G_IMPLM < 1 )
  delete[] TSGk;
 #endif

 #if( ! HV_NNVAR )
  delete[] dir;
 #endif

 }  // end( ~QPBundle )

/*--------------------------------------------------------------------------*/
/*-------------------- INTERFACE WITH THE BUNDLE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

SgRow QPBundle::GetItem( cIndex i )
{
 Insrtd = i;

 if( i >= DimMinQuad )  // allocate new stuff
 {
  #if( G_IMPLM < 3 )
   Index j = DimMinQuad;
  #endif

  do
   DimMinQuad = min( Index( DimMinQuad + SltSize ) , BPar2 );
  while( i >= DimMinQuad );

  SetCrrBDim( DimMinQuad );

  #if( G_IMPLM < 3 )
   for( ; j < DimMinQuad ; )
    SubG[ j++ ] = new SgNum[ MaxNumVar ];
  #endif
  }

 #if( G_IMPLM < 3 )
  return( SubG[ i ] );
 #else
  return( tmpG1 );
 #endif
 }

/*--------------------------------------------------------------------------*/

BOOL QPBundle::SetItem( cHpNum DFi , cHpNum Tau , HpNum &Ai , HpNum &ScPri ,
			cIndex_Set SGBse )
{
 #if( G_IMPLM < 3 )
  register SgRow tGi = SubG[ Insrtd ];
 #else
  register SgRow tGi = tmpG1;
 #endif

 if( SGBse )  // "densify" the item if it is "sparse" - - - - - - - - - - - -
 {
  register Index j = 0;
  while( SGBse[ j ] < InINF )
   j++;

  if( j )
  {
   register Index i = NumVar;
   if( i )
    for( j-- ; --i > j ; )
     if( i == SGBse[ j ] )
      tGi[ i ] = tGi[ j-- ];
     else
      tGi[ i ] = 0;
   }
  else
   VectAssign( tGi , SgNum( 0 ) , NumVar );
  }

 if( RHSp )  // consider the RHS (if present) - - - - - - - - - - - - - - - -
 {
  tGi += NumVar;
  for( register SgRow tRHS = RHSp + NumVar ; tRHS > RHSp ; )
  {
   register SgNum ttGi = *(--tGi);
   *tGi = ttGi - *(--tRHS);
   }
  }

 // update the row-wise data structure if necessary - - - - - - - - - - - - -

 #if( G_IMPLM > 0 )
  #if( HAVE_BTH )
   if( PPar2 )
   {
    register Index_Set LBt = LamBase;
    for( register Index h ; ( h = *(LBt++) ) < InINF ; )
     TSubG[ h ][ Insrtd ] = tGi[ h ];
    }
   else
 #endif
   {
    register SgMat tTSG = TSubG + NumVar;
    for( tGi += NumVar ; tTSG > TSubG ; )
     (*(--tTSG))[ Insrtd ] = *(--tGi);
    }
 #endif

 #if( HV_NNVAR )
  #if( ADD_CNST )
   if( DFi < HpINF )  // it is a subgradient
  #endif
    if( ( ! ActBDim ) && NNVars )  // the first subgradient ever- - - - - - -
    {
     register char *tV = GetVars();
     register cLMRow tLB = LowerBounds();
     register cHpNum eDir = EpsilonD() * max( BDim , Index( 1 ) );
     for( register Index cnt = NNVars ; cnt ; tV++ , tLB++ , tGi++ )
      if( *tV & NNVar() )
      {
       cnt--;
       if( *tV & IsVar() )
       {
        register cLMNum Gii = - t * (*tGi);
        if( Gii + size( Gii ) < *tLB )
         *tV |= AcVar();
        else
         *tV &= ~AcVar();
        }
       }

     InitialSetUp();

     #if( G_IMPLM < 3 )
      tGi = SubG[ Insrtd ];
     #else
      tGi = tmpG1;
     #endif
     }
 #endif

 if( Tau )  // calculate the scalar product - - - - - - - - - - - - - - - - -
  #if( HV_NNVAR )
   ScPri = DPerG( tGi );
  #else
   if( PPar2 )
    ScPri = ScalarProductBB( dir , tGi , LamBase );
   else
    ScPri = ScalarProduct( dir , tGi , NumVar );
  #endif
 else
  ScPri = 0;

 SetGTd( Insrtd , ScPri );

 // pass the item to the (QP) solver & check violation of dual constraint - -

 #if( ADD_CNST )
  if( DFi == HpINF )
  {
   AddConstr( Insrtd , Ai );
   return( TRUE );  // constraints always change the optimal solution
   }
  else
 #endif
  {
   if( Tau )
    Ai -= DFi + ScPri * Tau;

   if( Ai < MinNewAlfa )
    MinNewAlfa = Ai;

   AddSubGrad( Insrtd , Ai );

   if( Tau && ( Ai > 0 ) )  // a subgradient with Alfa <= 0 is always OK
   {
    /* This subgradient corresponds to the constraint

         v >= - PrvsTi * G1 * dir - Ai

       in the subproblem, where "- PrvsTi" because the direction is
       actually - t * dir: check if the constraint is violated by the
       previous optimal solution dir, i.e. if dir will change. */

    HpNum ro = PrvsTi * ScPri + Ai + Readv();
    HpNum aro = ABS( PrvsTi * ScPri ) + Ai - Readv();

    if( ro >= - eR * BDim * aro )
     return( FALSE );
    }

   return( TRUE );

   }  // end( else( is a subgradient ) )
 }  // end( QPBundle::SetItem )

/*--------------------------------------------------------------------------*/

void QPBundle::SetVars( void )
{
 assert( PPar2 );

 #if( HV_NNVAR )
  register cLMRow tL = Lambda;
  register char *vs = GetVars();
  register LMRow bs = LowerBounds();
  register Index_Set Actv = LamBase;
  for( register Index i = 0 ; i < NumVar ; i++ , bs++ )
   if( *Actv == i )
   {
    Actv++;
    *bs = - *(tL++);
    *(vs++) |= IsVar();  // the variable is in (whichever type it is)

    #if( HAVE_BTH )
     if( ! TSubG[ i ] )
     {
      GetNewGTRow( i );
      SetGTRow( i );
      }
    #endif
    }
   else
   {
    *(vs++) &= NNVar();  // the variable is out (but keep its type)

    #if( HAVE_BTH )
     if( TSubG[ i ] )
      PutGTRow( i );
    #endif
    }

  InitialSetUp();
 #else
  #if( HAVE_BTH )
   register Index_Set Actv = LamBase;
   for( register Index i = 0 ; i < NumVar ; i++ )
    if( *Actv == i )
    {
     Actv++;

     if( ! TSubG[ i ] )
     {
      GetNewGTRow( i );
      SetGTRow( i );
      }
     }
    else
     if( TSubG[ i ] )
      PutGTRow( i );
  #endif

  ChangeQ();  // update Q[][] in response to the change
 #endif

 }  // end( QPBundle::SetVars( void ) )

/*--------------------------------------------------------------------------*/

void QPBundle::AddVar( cIndex i , SgRow NwVar , cLMNum lb )
{
 // if necessary, allocate the memory for the new row- - - - - - - - - - - -

 #if( G_IMPLM > 0 )
  #if( HAVE_BTH )
   if( PPar2 )
    GetNewGTRow( i );
   else
  #endif
    while( i >= NrAllRws )
     for( Index j = min( Index( MaxNumVar - NrAllRws ) , GTSltSize ) ; j-- ; )
      TSubG[ NrAllRws++ ] = new SgNum[ BPar2 ];
 #endif

 if( ! NwVar )  // generated by the LVG strategy- - - - - - - - - - - - - - -
 {
  NwVar = TSGi( i );
  #if( HAVE_BTH )
   SetGTRow( i );
  #endif
  }
 else   // a "real" new variable- - - - - - - - - - - - - - - - - - - - - - -
 {
  SgNum RHSi = 0;  // take care of the R.H.S. and the floor

  if( RHSp )
   RHSi = RHSp[ i ] = NwVar[ BPar2 ];

  if( FlrNme < InINF )
   NwVar[ FlrNme ] = RHSi;

  // copy the new information into the internal data structures

  #if( G_IMPLM > 0 )
   VectAssign( TSubG[ i ] , NwVar , NxtBIdx );
  #endif

  #if( G_IMPLM < 3 )
   register SgMat sgt = SubG;
   register SgRow tTSGk = NwVar;
   for( register Index j = NxtBIdx ; j-- ; )
    (*(sgt++))[ i ] = *(tTSGk++);
  #endif

  #if( HV_NNVAR )
   if( lb == - LMINF )  // set the sign of the variable
    GetVars()[ i ] = 0;
   else
    GetVars()[ i ] = NNVar();
  #else
   assert( lb == - LMINF );
  #endif
  }

 // now add the variable to the subproblem- - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  LowerBounds()[ i ] = lb;
  BMinQuad::AddVar( i );
 #else
  AddSGSpaceDim( NwVar );
 #endif

 }  // end( QPBundle::AddVar )

/*--------------------------------------------------------------------------*/

void QPBundle::RemoveVar( cIndex i )
{
 assert( PPar2 );

 #if( HAVE_BTH )
  PutGTRow( i );
 #endif

 #if( HV_NNVAR )
  BMinQuad::RemoveVar( i );
 #else
  CutSGSpaceDim( TSGi( i ) );
 #endif

 }  // end( QPBundle::RemoveVar )

/*--------------------------------------------------------------------------*/

void QPBundle::SubstVar( cIndex i )
{
 // remove `i'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  BMinQuad::RemoveVar( i );
 #else
  if( ( ! PPar2 ) || 
      ( ( i == NumVar ) && ( LamBase[ LamDim - 1 ] == NumVar ) ) ||
      ( ( i < NumVar ) &&
	( LamBase[ BinSearch( LamBase , LamDim , i ) ] == i ) ) )
   CutSGSpaceDim( TSGi( i ) );
 #endif

 // update the Bundle (rows and/or columns and/or RHS)- - - - - - - - - - - -

 if( i < NumVar )  // it is not the last variable
 {
  #if( G_IMPLM > 0 )
   Swap( TSubG[ i ] , TSubG[ NumVar ] );
  #endif

  #if( G_IMPLM < 3 )
   register SgMat sgt = SubG;
   for( register Index j = NxtBIdx ; j-- ; )
   {
    register SgRow tSG = *(sgt++);
    tSG[ i ] = tSG[ NumVar ];
    }
  #endif

  if( RHSp )
   RHSp[ i ] = RHSp[ NumVar ];
  }

 // give back the row - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HAVE_BTH )
  if( PPar2 && TSubG[ NumVar ] )  // give back the row to the Buffer[]
  {
   Buffer[ BFCntr++ ] = TSubG[ NumVar ];
   TSubG[ NumVar ] = NULL;
   NrAllRws--;
   }
 #endif

 #if( G_IMPLM >= 4 )
  delete[] TSubG[ NumVar ];
  TSubG[ NumVar ] = NULL;
  NrAllRws--;
 #endif

 // give NumVar the name `i'- - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  if( i < NumVar )         // if it is not the last variable already
   MoveVar( NumVar , i );
 #endif

 }  // end( QPBundle::SubstVar )

/*--------------------------------------------------------------------------*/

Bundle::BStatus QPBundle::SolveSubP( cHpNum tt , HpNum &dNrm , HpNum &Sgm )
{
 // give a good stopping value to the (QP) solver - - - - - - - - - - - - - -

 if( FiLambda < HpINF )
  SetMinFVal( ( EpsLin * ABS( FiLambda ) ) / ( 2 * max( tt , tStar ) ) );

 // solve the (QP)- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ZCmptd = 0;
 SettmpD( Lambda1 );

 MQError qpres = SolveQP( tt );

 GettmpD( Lambda1 );

 if( qpres == MinQuad::kFatal )
 {
  SetEpsilonR();   // restore initial EpsilonR

  #if( HV_NNVAR )
   BMinQuad::SetEpsilonD();  // restore initial EpsilonD
  #endif

  return( kError );
  }

 if( qpres == MinQuad::kQPPrimUnbndd )
  return( kUnfeasible );

 // copy the optimal base to ZBase- - - - - - - - - - - - - - - - - - - - - -

 dNrm = MinQuad::ReaddNorm();
 Sgm = MinQuad::ReadSigma();

 #if( DO_AGGR )
  ZxD = dNrm;

  #if( HV_NNVAR )
  ZSigma = BMinQuad::ReadSigma( FALSE );
   ZxD += Sgm - ZSigma;
  #else
   ZSigma = Sgm;
  #endif

  #if( ADD_CNST )
   if( CBZDim = ReadCBDim() )
   {
    register cIndex_Set tB = MinQuad::ReadBase();
    register cHpRow tM = MinQuad::ReadMult();
    register Index i;

    for( ZBDim = 0 ; ( i = ZBase[ ZBDim ] = *(tB++) ) < InINF ; tM++ )
     if( ! IsAConst( i ) )
      ZMult[ ZBDim++ ] = *tM;
    }
   else
  #endif
  #if( ADD_CNST || ( DO_AGGR == 2 ) )
   {
    VectAssign( ZMult , MinQuad::ReadMult() , ZBDim = MinQuad::ReadBDim() );
    VectAssign( ZBase , MinQuad::ReadBase() , ZBDim + 1 );
    }
  #endif
 #endif

 // set ZCmptd- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  #if( LAZY_D == 2 )
   if( NNVars < NumVar )
    ZCmptd = 2;
   else
  #endif
    if( PPar2 )
     ZCmptd = 3;
    else
     ZCmptd = 4;
 #else
  ZCmptd = 1;
 #endif

 return( Bundle::kOK );

 }  // end( QPBundle::SolveSubP )

/*--------------------------------------------------------------------------*/

cLMRow QPBundle::Price( void )
{
 return( CompleteZ( MinQuad::ReadBase() , MinQuad::ReadMult() ) );

 }  // end( QPBundle::Price )

/*--------------------------------------------------------------------------*/

HpNum QPBundle::CheckAlfa( cBOOL All )
{
 register HpNum mA = MinNewAlfa;
 MinNewAlfa = HpINF;

 if( All )
 {
  register cHpRow tA = MinQuad::ReadAlfa();

  for( register Index i = 0 ; i < NxtBIdx ; tA++ )
   if( IsASubG( i++ ) && ( mA > *tA ) )
    mA = *tA;
  }

 if( mA < - max( ABS( FiLambda1 ) , HpNum( 1 ) ) * eR )
  ChangeAlfa( - mA );
 else
  mA = 0;

 return( mA );

 }  // end( QPBundle::CheckAlfa )

/*--------------------------------------------------------------------------*/

void QPBundle::ChgAlfa( cIndex i , cHpNum DeltaAlfai )
{
 ChangeAlfa( i , DeltaAlfai );
 }

/*--------------------------------------------------------------------------*/

void QPBundle::ChgAlfa( cHpRow DeltaAlfa , BOOL Mde )
{
 if( Mde )
  for( register Index i = 0 ; i < DimMinQuad ; )
   ChangeAlfa( i++ , *(DeltaAlfa++) );
 else
  ChangeAlfa( *DeltaAlfa );

 }  // end( QPBundle::ChgAlfa( cHpRow , BOOL ) )

/*--------------------------------------------------------------------------*/

void QPBundle::ChangeCurrPoint( cLMRow DLambda , cHpNum DFi )
{
 // first, update the bounds - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  VectSum( LowerBounds() , DLambda , NumVar );

  ChangeBounds();
 #endif

 // now change the Alfa of each item still in the Bundle- - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ActBDim )
 {
  #if( WHAT_4_Q )  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   register Index i = 0;
   for( ; i < NxtBIdx ; i++ )
    if( IsASubG( i ) )
     tmpa[ i ] = DFi;
    else
     tmpa[ i ] = 0;

   register cLMRow tDL = DLambda;
   register LMNum tDLi;

   if( PPar2 )
   {
    for( register Index_Set tB = LamBase ; ( i = *(tB++) ) < InINF ; )
     if( ( tDLi = *(tDL++) ) )
     {
      register SgRow tTSGi = TSubG[ i ] + NxtBIdx;
      for( register HpRow ttmpa = tmpa + NxtBIdx ; ttmpa > tmpa ; )
       *(--ttmpa) += *(--tTSGi) * tDLi;
      }
    }
   else
    for( i = 0 ; i < NumVar ; i++ )
     if( ( tDLi = *(tDL++) ) )
     {
      register SgRow tTSGi = TSubG[ i ] + NxtBIdx;
      for( register HpRow ttmpa = tmpa + NxtBIdx ; ttmpa > tmpa ; )
       *(--ttmpa) += *(--tTSGi) * tDLi;
      }

   for( register Index h = NxtBIdx ; h-- ; )
    ChangeAlfa( h , tmpa[ h ] );

  #else   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   for( register Index i = 0 ; i < NxtBIdx ; i++ )
    if( IsThere( i ) )
    {
     register HpNum res = IsASubG( i ) ? DFi : 0;
     register cSgRow tGh = SubG[ i ];

     if( PPar2 )
     {
      register Index_Set tB = LamBase;
      for( register Index h ; ( h = *(tB++) ) < InINF ; )
       res += tGh[ h ] * DLambda[ h ];
      }
     else
     {
      register cLMRow tDL = DLambda + NumVar;
      for( tGh += NumVar  ; tDL > DLambda ; )
       res += *(--tDL) * (*(--tGh));
      }

     ChangeAlfa( i , res );
     }

  #endif   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  }  // end if( ActBDim )
 }  // end( ChangeCurrPoint( DLambda , DFi )

/*--------------------------------------------------------------------------*/

void QPBundle::ChangeCurrPoint( cHpNum Tau , cHpNum DFi )
{
 MoveAlongD( Tau , DFi );

 }  // end( QPBundle::ChangeCurrPoint( Tau , DFi ) )

/*--------------------------------------------------------------------------*/

cHpRow QPBundle::ReadAlfa( void )
{
 return( MinQuad::ReadAlfa() );
 }

/*--------------------------------------------------------------------------*/

void QPBundle::ChgRHSi( cIndex i , cSgNum DeltaRHSi )
{
 // if necessary, allocate the i-th row - - - - - - - - - - - - - - - - - - -

 #if( HAVE_BTH )
  if( PPar3 && ( ! TSubG[ i ] ) )
  {
   GetNewGTRow( i );
   SetGTRow( i );
   }
 #endif

 // if the variable exists, remove it - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  if( PPar2 && ( ! ( GetVar( i ) & IsVar() ) ) )
   LowerBounds()[ i ] = 0;
  else
   RemoveVar( i );
 #endif

 register SgRow SGi = TSGi( i );

 #if( ! HV_NNVAR )
  if( ( ! PPar2 ) || ( LamBase[ BinSearch( LamBase , LamDim , i ) ] == i ) )
   CutSGSpaceDim( SGi );
 #endif

 // put there the new modified entries- - - - - - - - - - - - - - - - - - - -

 for( register Index j = 0 ; j < NxtBIdx ; j++ )
 {
  SGi[ j ] += DeltaRHSi;

  #if( G_IMPLM < 3 )
   SubG[ j ][ i ] = SGi[ j ];
  #endif
  }

 // now (re-)create variable i with new entries - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  BMinQuad::AddVar( i );
 #else
  AddSGSpaceDim( SGi );
 #endif

 }  // end( QPBundle::ChangeRHSi )

/*--------------------------------------------------------------------------*/

void QPBundle::MakeLambda1( cHpNum Tau )
{
 MakeZ();  // ensure that everything is ready

 #if( HV_NNVAR )
  if( PPar2 )
   AddDSprs( Lambda1 , Lambda , Tau );
  else
   AddD( Lambda1 , Lambda , Tau );
 #else
  if( PPar2 )
  {
   register LMRow tL = Lambda;
   register LMRow tL1 = Lambda1;
   register Index_Set tLB = LamBase;
   for( register Index h ; ( h = *(tLB++) ) < InINF ; )
    *(tL1++) = *(tL++) - Tau * dir[ h ];
   }
  else
   VectAdd( Lambda1 , Lambda , dir , LMNum( - Tau ) , NumVar );
 #endif

 }  // end( QPBundle::MakeLambda1 )

/*--------------------------------------------------------------------------*/

HpNum QPBundle::EpsilonD( void )
{
 #if( HV_NNVAR )
  return( BMinQuad::EpsilonD() );
 #else
  return( EpsD * ( BDim ? BDim : 1 ) );
 #endif
 }

/*--------------------------------------------------------------------------*/

void QPBundle::SensitAnals( HpNum &lp , HpNum &cp )
{
 HpNum foo;
 SensitAnals1( lp , cp , foo );
 }

/*--------------------------------------------------------------------------*/

#if( DO_AGGR )

 void QPBundle::AggregateZ( cIndex where )
 {
  AggrPrimalSol( ZMult , ZBase , ZBDim , where );

  #if( G_IMPLM < 3 )
   SgRow tmpG1 = SubG[ where ];
  #endif

  #if( ADD_CNST )  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   if( CBZDim )    // the optimal Base contains constraints
   {
    #if( ( ! WHAT_4_D ) || HAVE_BTH )
     // if FullZ() may use SubG[], it is necessary to check whether or not
     // where is in ZBase and, in case, deal with it explicitly

     register Index i;
     #if( WHAT_4_D && HAVE_BTH )
      if( ! PPar2 )  // the rows are used if both are there and not L.V.G.
       i = ZBDim;    // which implies that the columns are not used
      else
     #endif
       for( i = 0 ; i < ZBDim ; i++ )
        if( ZBase[ i ] == where )
         break;

     if( i < ZBDim )  // where is in ZBase
     {
      Swap( ZBase[ i ] , *ZBase );
      Swap( ZMult[ i ] , *ZMult );

      ScalarMult( tmpG1 , *ZMult , NumVar );

      for( register cIndex_Set tZB = ZBase ; *(++tZB) < InINF ; )
       VectSum( tmpG1 , SubG[ *tZB ] , NumVar );
      }
     else
    #endif
      FullZ( tmpG1 , ZBase , ZMult );

    if( ZCmptd )  // compute the scalar product
    {
     MakeZ();  // to do that, the direction must be calculated

     #if( HV_NNVAR )
      ZxD = DPerG( tmpG1 );
     #else
      if( PPar2 )
       ZxD = ScalarProductBB( dir , tmpG1 , LamBase );
      else
       ZxD = ScalarProduct( tmpG1 , dir , NumVar );
     #endif
     }
    }  // end if( CBZDim )
   else
  #endif  /* ! ADD_CNST- - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   {
    VectAssign( tmpG1 , CompleteZ( ZBase , ZMult ) , NumVar );

    }  // end else( ! CBZDim )

  // if necessary, copy Z in the row-wise data structure- - - - - - - - - - -

  #if( G_IMPLM > 0 )
   #if( HAVE_BTH )
    if( PPar2 )
    {
     register SgRow tZ = tmpG1;
     register Index_Set LBt = LamBase;
     for( register Index h ; ( h = *(LBt++) ) < InINF ; )
      TSubG[ h ][ where ] = tZ[ h ];
     }
    else
   #endif
    {
     register SgRow tz = tmpG1 + NumVar;
     for( register SgMat tTSG = TSubG + NumVar ; tTSG > TSubG ; )
      (*(--tTSG))[ where ] = *(--tz);
     }
  #endif

  // delete the old item & insert Z - - - - - - - - - - - - - - - - - - - - -

  ResetBundle( where );

  SetGTd( Insrtd = where , ZxD );

  AddSubGrad( where , ZSigma );

  }  // end( AggregateZ )

#endif

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

#if( LAZY_Q )

 QuNum QPBundle::GiTGj( cIndex i , cIndex j )
 {
  if( SprsV )
   if( i == j )
    return( Norm( SubG[ i ] , MBase2 ) );
   else
    return( ScalarProductBB( SubG[ i ] , SubG[ j ] , MBase2 ) );
  else
   if( i == j )
    return( Norm( SubG[ i ] , NumVar ) );
   else
    return( ScalarProduct( SubG[ i ] , SubG[ j ] , NumVar ) );

  }  // end( QPBundle::GiTGj )

#else  /*-------------------------------------------------------------------*/

 void QPBundle::GiTG( cIndex i , register QuRow Qi , cIndex iMax )
 {
  #if( WHAT_4_Q )  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   VectAssign( Qi , QuNum( 0 ) , iMax );

   if( SprsV )
   {
    register cIndex_Set tMB2 = MBase2;
    for( register Index h ; ( h = *(tMB2++) ) < InINF ; )
     VectSum( Qi , TSubG[ h ] , TSubG[ h ][ i ] , iMax );
    }
   else
   {
    register SgMat tTSG = TSubG;
    for( register Index h = NumVar ; h-- ; tTSG++ )
     VectSum( Qi , *tTSG , (*tTSG)[ i ] , iMax );
    }

  #else  /* WHAT_4_Q == 0- - - - - - - - - - - - - - - - - - - - - - - - - -*/

   register SgMat tSG = SubG;
   register SgRow SGi = tSG[ i ];
   register Index n = 0;

   if( SprsV )
   {
    for( ; n < i ; n++ , Qi++ , tSG++ )
     if( IsThere( n ) )
      *Qi = ScalarProductBB( SGi , *tSG , MBase2 );

    *(Qi++) = Norm( *(tSG++) , MBase2 );

    for( ; ++n < iMax ; Qi++ , tSG++ )
     if( IsThere( n ) )
      *Qi = ScalarProductBB( SGi , *tSG , MBase2 );
    }
   else
   {
    for( ; n < i ; n++ , Qi++ , tSG++ )
     if( IsThere( n ) )
      *Qi = ScalarProduct( SGi , *tSG , NumVar );

    *(Qi++) = Norm( *(tSG++) , NumVar );

    for( ; ++n < iMax ; Qi++ , tSG++ )
     if( IsThere( n ) )
      *Qi = ScalarProduct( SGi , *tSG , NumVar );
    }

  #endif   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  }  // end( QPBundle::GiTG )

#endif

/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )

 cSgRow QPBundle::GiTilde( cIndex i )
 {
  #if( G_IMPLM > 0 )
   return( TSubG[ i ] );
  #else
   return( TSGi( i ) );
  #endif

  }  // end( QPBundle::GiTilde )

#endif

/*--------------------------------------------------------------------------*/

LMNum QPBundle::CalculateZ( cIndex h )
{
 #if( WHAT_4_D )
  return( ScalarProduct( Mult , TSubG[ h ] , Base ) );
 #else
  register LMNum zh = 0;
  register HpRow Mlt = Mult;
  register Index_Set Bse = Base;
  for( register Index i ; ( i = *(Bse++) ) < InINF ; )
   zh += *(Mlt++) * SubG[ i ][ h ];

  return( zh );
 #endif

 }  // end( QPBundle::CalculateZ( < one > ) )

/*--------------------------------------------------------------------------*/

void QPBundle::CalculateZ( register cIndex_Set Wh , register LMRow z )
{
 #if( WHAT_4_D )

  for( register Index h ; ( h = *(Wh++) ) < InINF ; )
   z[ h ] = ScalarProduct( Mult , TSubG[ h ] , Base );

 #else

  if( BDim )
  {
   register cHpRow tM = Mult;
   register cIndex_Set tB = Base;
   VectAssignBB( z , SubG[ *tB ] , *tM , Wh );

   for( register Index h ; ( h = *(++tB) ) < InINF ; )
    VectSumBB( z , SubG[ h ] , *(++tM) , Wh );
   }
  else
   VectAssignB( z , LMNum( 0 ) , Wh );

 #endif

 }  // end( QPBundle::CalculateZ( < a set > ) )

/*--------------------------------------------------------------------------*/

void QPBundle::CalculateZ( register LMRow z )
{
 if( PPar2 )
  QPBundle::CalculateZ( LamBase , z );
 else
  FullZ( z , Base , Mult );

 }  // end( QPBundle::CalculateZ( < all > ) )

/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )

 LMNum QPBundle::GiTLB( cIndex h , cLMRow l )
 {
  #if( G_IMPLM >= 3 )
   if( h != Insrtd )
   {
    register LMNum res = 0;
    register Index_Set Bse = Base2;
    for( register Index i ; ( i = *(Bse++) ) < InINF ; )
     res += TSubG[ i ][ h ] * l[ i ];

    return( res );
    }
   else
    return( ScalarProductBB( l , tmpG1 , Base2 ) );
  #else
    return( ScalarProductBB( l , SubG[ h ] , Base2 ) );
  #endif
  }

#endif

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

#if( HAVE_BTH )

 inline void QPBundle::GetNewGTRow( cIndex i )
 {
  if( ! BFCntr )
   for( Index j = BFCntr = min( MaxNumVar - NrAllRws , GTSltSize ) ; j-- ; )
    Buffer[ j ] = new SgNum[ BPar2 ];

  TSubG[ i ] = Buffer[ --BFCntr ];

  NrAllRws++;
  }

/*--------------------------------------------------------------------------*/

 inline void QPBundle::SetGTRow( cIndex i )
 {
  // copy the entries from the columns into TSubG[ i ]

  register SgMat sgt = SubG;     
  register SgRow tTSGk = TSubG[ i ];
  for( register Index j = NxtBIdx ; j-- ; )
   *(tTSGk++) = (*(sgt++))[ i ];
  }

/*--------------------------------------------------------------------------*/

 inline void QPBundle::PutGTRow( cIndex i )
 {
  Buffer[ BFCntr++ ] = TSubG[ i ];
  TSubG[ i ] = NULL;
  NrAllRws--;
  }

#endif

/*--------------------------------------------------------------------------*/

#if( G_IMPLM < 1 )

 inline SgRow QPBundle::TSGi( cIndex i )
 {
  register SgMat sgt = SubG;
  register SgRow tTSGk = TSGk;
  for( register Index j = NxtBIdx ; j-- ; )
   *(tTSGk++) = (*(sgt++))[ i ];

  return( TSGk );
  }

#endif

/*--------------------------------------------------------------------------*/

inline void QPBundle::FullZ( register LMRow z , cIndex_Set tB , cHpRow tM )
{
 // computes all the entries of Z, comprised those not corresponding to
 // Multipliers that are defined in the QP subproblem; if HAVE_BTH and the
 // row-wise representation is preferred, however, it can be used only if
 // no L.V.G. is being used, since in this case not all the entries are in
 // principle available row-wise

 #if( WHAT_4_D && HAVE_BTH )   /*- - - - - - - - - - - - - - - - - - - - -*/
  if( ! PPar2 )
 #endif
  #if( WHAT_4_D )
  {
   register SgMat tTSG = TSubG + NumVar;
   for( z += NumVar ; tTSG > TSubG ; )
    *(--z) = ScalarProduct( tM , *(--tTSG) , tB );
   }
  #endif
 #if( WHAT_4_D && HAVE_BTH )
  else
 #endif
  #if( ( ! WHAT_4_D ) || HAVE_BTH )  /*- - - - - - - - - - - - - - - - - -*/
   if( BDim )
   {
    VectAssign( z , SubG[ *tB ] , *tM , NumVar );

    for( register Index h ; ( h = *(++tB) ) < InINF ; )
     VectSum( z , SubG[ h ] , *(++tM) , NumVar );
    }
   else
    VectAssign( z , LMNum( 0 ) , NumVar );

  #endif   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 }  // end( FullZ )

/*--------------------------------------------------------------------------*/

inline void QPBundle::MakeZ( void )
{
 #if( HV_NNVAR )  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  #if( LAZY_D == 2 )
   if( ZCmptd == 2 )
    #if( WHAT_4_D )
     if( NNVars )  // compute the entries of z corresponding to UC variables,
     {             // taking into account that those corresponding to NN
                   // variables are already there
      register cIndex_Set tMB2 = InActiveVars();
      for( register Index h ; ( h = *(tMB2++) ) < InINF ; )
       if( ! ( GetVar( h ) & NNVar() ) )
        SetD( h , ScalarProduct( Mult , TSubG[ h ] , Base ) );

      if( PPar2 )   // if not all variables are defined
       ZCmptd = 3;  // then there are still entries to compute
      else
       ZCmptd = 4;  // else all the entries are already there
      }
     else         // there are no NN variables => just compute all z
    #endif
     {
      FullZ( SetD() , MinQuad::ReadBase() , MinQuad::ReadMult() );
      ZCmptd = 4;
      }
  #endif
 #else  /* ! HV_NNVAR - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  if( ZCmptd == 1 )
  {
   CalculateZ( dir );

   if( PPar2 )
    ZCmptd = 3;
   else
    ZCmptd = 4;
   }

 #endif   /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 }  // end( MakeZ )

/*--------------------------------------------------------------------------*/

cLMRow QPBundle::CompleteZ( cIndex_Set tB , cHpRow tM )
{
 #if( HV_NNVAR )
  LMRow dir = SetD();
 #endif

 if( ZCmptd < 4 )  // compute the entries of d[] for non-defined variables
 {
  #if( G_IMPLM >= 3 )
   #if( HV_NNVAR && ( LAZY_D == 2 ) )
    if( ZCmptd == 2 )
    {
     register char *vs = GetVars();          // skip variables that are both
     register char NNP = NNVar() | IsVar();  // NN and defined in the QP
     for( register Index h = 0 ; h < NumVar ; h++ )
      if( ( *(vs++) & NNP ) != NNP )
       dir[ h ] = ScalarProduct( tM , TSubG[ h ] , tB );
     }
    else
   #endif
     if( ZCmptd == 3 )
     {
      register Index_Set tLB = LamBase;
      for( register Index h = 0 ; h < NumVar ; h++ )
       if( *tLB == h )
        tLB++;
       else
        dir[ h ] = ScalarProduct( tM , TSubG[ h ] , tB );
      }
     else  // ZCmptd == 1
  #endif
      FullZ( dir , tB , tM );

  ZCmptd = 4;
  }

 return( dir );

 }  // end( CompleteZ )

/*--------------------------------------------------------------------------*/
/*-------------------------- End File QPBundle.C ---------------------------*/
/*--------------------------------------------------------------------------*/
