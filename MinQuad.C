/*--------------------------------------------------------------------------*/
/*---------------------------- File MinQuad.C ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Implementation of the TT (Tall and Thin) algorithm for solving the   --*/
/*-- Quadratic Problems arising as tentative descent direction finding    --*/
/*-- subproblem within Bundle algorithms for the (linearly constrained)   --*/
/*-- minimization of NonDifferentiable convex functions.                  --*/
/*--                                                                      --*/
/*--                            VERSION 4.10      			  --*/
/*--                	       18 - 10 - 2004			       	  --*/
/*--                                                                      --*/
/*--                    Original Idea and Implementation by:              --*/
/*--                                                                      --*/
/*--                            Antonio Frangioni                         --*/
/*--                                                                      --*/
/*--                         Operations Research Group                    --*/
/*--                        Dipartimento di Informatica                   --*/
/*--                           Universita' di Pisa                        --*/
/*--                                                                      --*/
/*--                         Copyright 1992 - 2004                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MinQuad.h"

#include "OPTvect.h"

/*--------------------------------------------------------------------------*/
/*------------------------- UNDOCUMENTED SWITCHES --------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------- IFELEMEQZERO ---------------------------------*/

#define IFELEMEQZERO 0

/* In the [Add/Cut]SGSpaceDim() methods, it is possible to control if a given
   element of the [new/old] dimension is zero in order to avoid doing the
   relative sweep of [complex] Givens matirces. Surprisingly enough, if there
   are not many zeroes this has a rather bad effect on running times of codes
   heavily relying on these procedures. This switch is given to turn on/off
   this feature, so that it can be used for problems with many nonzeroes in
   the [new/old] dimensions. */

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--                                                                      --*/
/*--      Some small macro definitions, used throughout the code.         --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

#define stdsize() eR * Lmu

#define size( x ) stdsize() * max( ABS( x ) , HpNum( 1 ) )
#define gsize( x ) stdsize() * max( x , HpNum( 1 ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( VARCOEFF )
 #define e( i , x ) x * cf[ i ]
#else
 #define e( i , x ) ( G[ i ] & IsACnst ? 0 : x )
#endif

#if( VARCOEFF )
 #define e0( i , x ) x * cf[ i ]
#else
 #define e0( i , x ) x
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( IFELEMEQZERO )
 #define cond_if( x ) if( x )
#else
 #define cond_if( x )
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( LAZY_Q )
 #define Qif( x ) if( x )
#else
 #define Qif( x )
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
using namespace MinQuad_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static cHpNum MaxEpsR = 1e-4;   // Maximum, ...
static cHpNum IntEpsR = 1e-10;  // ... Initial and ...
static cHpNum MinEpsR = 1e-16;  // ... Minimum value for eR

static cHpNum MFactor = 10;     // eR is [in/de]creased by [multiply/divid]ing
                                // it by this factor

static const char Void    = 0;  // an empty slot in the Bundle
static const char SubG    = 1;  // a (non basic) subgradient
static const char Cnst    = 3;  // a (non basic) constraint
static const char SubGinB = 5;  // a basic subgradient
static const char CnstinB = 7;  // a basic constraint

static const char IsACnst = 2;  // ( i & 2 ) == 1 <=> i is a constraint,
                                // possibly a basic one
static const char IsInBse = 4;  // ( i & 4 ) == 1 <=> i is a basic item,
                                // regardless to its type

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  These functions are not implemented as methods of the class, since  --*/
/*--  they don't use directly its data structures.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--  Functions for triangular matrices.                        	  --*/
/*----------------------------------------------------------------------------

static inline void LPerV( register const cHpRow *L , cHpRow V ,
                          register HpRow Res , cIndex n )
{
 // Res := L * V, where L is an n-by-n lower triangular matrix and V,
 // Res are n-vectors

 for( register Index i = 0 ; i++ < n ; L++ )
 {
  register HpNum t = 0;
  register cHpRow v1 = V;

  for( register Index j = 0 ; j < i ; )
   t += (*L)[ j++ ] * (*(v1++));

  *(Res++) = t;
  }
 }

------------------------------------------------------------------------------

static inline void LTPerV( const cHpRow *L , cHpRow V , register HpRow Res ,
                           cIndex n )
{
 // Res := L{T} * V; as above, but now L{T} is upper triangular

 for( register Index i = 0 ; i < n ; )
 {
  register HpNum t = 0;
  register cHpRow v1 = (V++);

  for( register Index j = i++ ; j < n ; )
   t += L[ j++ ][ i ] * (*v1);

  *(Res++) = t;
  }
 }

----------------------------------------------------------------------------*/

static inline void SolveSystem( register const cHpRow *L , register cHpRow b ,
				HpRow Res , cIndex Dim )
{
 // solves the system L * Result = b where L is a Dim-by-Dim lower triangular
 // matrix and b, Result are Dim-vectors

 for( register Index i = 0 ; i < Dim ; )
 {
  register HpNum t = 0;
  register HpRow r = Res;
  register cHpRow l = *(L++);

  for( register Index j = i++ ; j-- ; )
   t += *(l++) * (*(r++));

  (*r) = ( *(b++) - t ) / (*l);
  }
 }  // end( SolveSystem )

/*--------------------------------------------------------------------------*/

static inline void SolveSystemT( const cHpRow *L , register cHpRow b ,
				 HpRow Res , Index Dim )
{
 register Index i = Dim--;  // solves the system L{T} * Result = b
 b += Dim;                  // as above, but now L{T} is upper triangular
 L += Dim;
 Res += Dim;

 for( ; i ; )
 {
  register HpNum t = 0;
  register HpRow r = Res;
  register const cHpRow *l = L;

  for( register Index j = Dim - (--i) ; j-- ; )
   t += (*(l--))[ i ] * (*(r--));

  (*r) = ( *(b--) - t ) / (*l)[ i ];
  }
 }  // end( SolveSystemT )

/*--------------------------------------------------------------------------*/

static inline void SolvePerSys( register const cHpRow *L ,
				register cHpRow b , register HpRow R ,
				register cIndex_Set B , cIndex D )

/* Solve the triangular system L * R = b{B}, L being a D x D lower triangular
   matrix and b{B} the D-vector having as i-th element b[ B[ i ] ]. */
{
 for( register Index i = 0 ; i < D ; )
 {
  register HpNum t = 0;
  register HpRow r = R;
  register cHpRow l = *(L++);

  for( register Index j = i++ ; j-- ; )
   t += *(l++) * (*(r++));

  (*r) = ( b[ *(B++) ] - t ) / (*l);
  }
 }

/*--------------------------------------------------------------------------*/
/*--  Functions for Givens matrices.                                      --*/
/*----------------------------------------------------------------------------

static inline void ApplyGivensMToRow( register cHpNum c , register cHpNum s ,
                                      HpNum &ti , HpNum &tj )

   Apply the Givens matrix G to the two numbers ti and tj, i.e. the i-th
   and j-th element of a row. The Givens matrix is G( i , j , c , s )
      _____________________
     | 1                   |
     |   \                 |
     |     c  ...   -s ... |  i       where c^2 + s^2 = 1
     |     :         :     |
     |     s  ...    c ... |  j
     |     :         : \   |
     |_____:_________:___1_|
           i         j
 
   hence, being B = G( i , j , c , s ) * A,

            / A[ h ]                     if h != i && h != j
   B[ h ] = | c * A[ i ] - s * A[ j ]    if h == i
            \ c * A[ j ] + s * A[ i ]    if h == j

   Note that Givens matrices are not simmetric, but obviuosly

    G( i , j , c , s ){T} = G( i , j , c , -s )
{
 register HpNum t1 = ti;
 register HpNum t2 = tj;

 ti = c * t1 - s * t2;
 tj = s * t1 + c * t2;
 }

----------------------------------------------------------------------------*/

static inline void ApplyGivensMToCol( register cHpNum c , register cHpNum s ,
                                      HpNum &ti , HpNum &tj )

/* Same as before, but now G is applied to a column. */
{
 register HpNum t1 = ti;
 register HpNum t2 = tj;

 ti = c * t1 + s * t2;
 tj = c * t2 - s * t1;
 }

/*--------------------------------------------------------------------------*/

static inline void ApplyCompGivensMToCol( register cHpNum c ,
                                          register cHpNum s ,
                                          HpNum &ti , HpNum &tj )

/* Same as before, but now G is a *complex* Givens matrix, i.e.
      _____________________
     | 1                   |
     |   \                 |
     |     c  ...  -Is ... |  i       where c^2 - s^2 = 1
     |     :         :     |
     |    Is  ...    c ... |  j       and I * I = -1
     |     :         : \   |
     |_____:_________:___1_|
           i         j        */
{
 register HpNum t1 = ti;
 register HpNum t2 = tj;

 ti = c * t1 + s * t2;
 tj = c * t2 + s * t1;
 }

/*--------------------------------------------------------------------------*/
/*--  Other functions.                                                    --*/
/*--------------------------------------------------------------------------*/

static inline void AdjusttBounds( HpNum M1j , HpNum M2j , HpNum eM1j ,
                                  HpNum eM2j , HpNum &tMin , HpNum &tMax )

/* Adjust the values of tMin and tMax, given that the constraint
   t * M1j <= M2j must be satisfied. eXXX is the ABSOLUTE error admitted
   on XXX, i.e. ABS( XXX ) < eXXX => XXX == 0. */
{
 if( ABS( M2j ) < eM2j )
  M2j = 0;

 if( M1j > eM1j )                  // M1j > 0 => t <= M2j / M1j
  tMax = min( tMax , M2j / M1j );
 else
  if( M1j < - eM1j )               // M1j < 0 => t >= M2j / M1j
   tMin = max( tMin , M2j / M1j );
   // else M1j = 0, and the constraint t * 0 <= M2j must be satisfied for
   // all t, since it is satisfied at least for one t.
 }
 
/*--------------------------------------------------------------------------*/

static inline void VectAdd( register HpRow g , register cHpRow g1 ,
			    register cHpRow g2 , register cHpNum x ,
			    register cHpNum y , register Index n )
/* g := x * g1 + y * g2 */
{
 for( ; n-- ; )
  *(g++) = x * (*(g1++)) + y * (*(g2++));
 }

/*--------------------------------------------------------------------------*/
/*--------------------- IMPLEMENTATION OF MinQuad  -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

MinQuad::MinQuad( void )
{
 MaxBDim = 0;  // signal that the real initialization has not been done yet

 #if( LOG_MQ )
  MQLog = &clog;

  #if( LOG_MQ > 1 )
   Calls = Step = Success = 0;
   SumBDim = SumActBDim = 0;
   SumGSBase = SumGSBundle = 0;
   Insert = Delete = 0;
   SGSIncr = SGSDecr = 0;
  #endif
 #endif

 #if( TIMERS_MQ )
  MQt = new OPTtimers();
 #endif

 }  // end( MinQuad )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void MinQuad::SetMaxDim( Index m , Index n , Index SDim )
{
 if( MaxBDim )   // this is not the first call: - - - - - - - - - - - - - - -
  MemDealloc();  // deallocate everything - - - - - - - - - - - - - - - - - -

 MaxBDim = m;

 if( MaxBDim  )  // m != 0: allocate everything - - - - - - - - - - - - - - -
 {               // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 // the max dim. of L is min( MaxBDim , SDim + 3 ) if SDim
                 // is known, and MaxBDim otherwise
  CrrBDim = min( n , m );
  SgSpDim = ( SDim ? min( MaxBDim , Index( SDim + 3 ) ) : MaxBDim );

  G     = new char[ MaxBDim ];
  Alfa  = new HpNum[ MaxBDim ];
  GTd   = new HpNum[ MaxBDim ];
  Order = new Index[ MaxBDim + 1 ];
  Q     = new QuRow[ MaxBDim ];
  L     = new HpRow[ SgSpDim ];
  Mult  = new HpNum[ MaxBDim + 1 ];
  tmpv  = new HpNum[ MaxBDim + 1 ];
  Base  = new Index[ SgSpDim + 1 ];
  z2    = new HpNum[ SgSpDim ];
  z1    = new HpNum[ SgSpDim ];

  #if( VARCOEFF )
   cf   = new HpNum[ MaxBDim ];
  #endif

  tmpa  = new HpNum[ MaxBDim ];

  #if( ! VASTE_MEM )
   tmpr = new HpNum[ MaxBDim ];
  #endif

  for( register Index i = 0 ; i < CrrBDim ; )
  {
   Q[ i ] = new QuNum[ i + 1 ];

   #if( LAZY_Q )
    register QuRow pq = Q[ i ];
    G[ i++ ] = Void;

    for( register Index j = i ; j-- ; )
     *(pq++) = -QuINF;
   #else
    G[ i++ ] = Void;
   #endif
   }

  #if( ! LAZY_Q )
   tempQ = new QuNum[ MaxBDim ];
  #endif

  for( register Index i = 0 ; i < min( SgSpDim , CrrBDim ) ; i++ )
   #if( VASTE_MEM < 2 )
    L[ i ] = new HpNum[ i + 1 ]
   #else
    L[ i ] = new HpNum[ SgSpDim ];
   #endif

  *Order = *Base = InINF;

  }  // end( if( m != 0 ) )- - - - - - - - - - - - - - - - - - - - - - - - - -

 // variables initialization - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ActBDim = ActCNum = BDim = CBDim = Dependent = LBIsInB = NxtBIdx = TLDim = 0;
 LwrBnd = Minf = HpINF;
 LBMult = EpsRo = 0;
 MinFVal = - HpINF;
 AddedOne = FALSE;
 eR = IntEpsR;
 PrvsTi = 1;

 }  // end( MinQuad::SetMaxDim )

/*--------------------------------------------------------------------------*/

void MinQuad::SetCrrBDim( cIndex Newn )
{
 if( ( Newn <= 0 ) || ( Newn > MaxBDim ) )  // invalid active dimension
  return;                                   // quietly return

 if( Newn > CrrBDim )  // increasing active Bundle dimension- - - - - - - - -
 {
  for( register Index i = CrrBDim ; i < Newn ; )
  {
   #if( LAZY_Q )
    register QuRow pq = Q[ i ] = new QuNum[ i + 1 ];
    G[ i++ ] = Void;

    for( register Index j = i ; j-- ; )
     *(pq++) = -QuINF;
   #else
    Q[ i ] = new QuNum[ i + 1 ];
    G[ i++ ] = Void;
   #endif
   }

  for( register Index i = CrrBDim ; i < min( SgSpDim , Newn ) ; i++ )
   #if( VASTE_MEM < 2 )
    L[ i ] = new HpNum[ i + 1 ];
   #else
    L[ i ] = new HpNum[ SgSpDim ];
   #endif
  }
 else  // decreasing active Bundle dimension- - - - - - - - - - - - - - - - -
 {
  for( register Index i = CrrBDim ; i-- > Newn ; )
   if( G[ i ] )
    ResetBundle( i );

  for( register Index i = CrrBDim ; i-- > Newn ; )
   delete[] Q[ i ];

  for( register Index i =  min( SgSpDim , CrrBDim ) ; i-- > Newn ; )
   delete[] L[ i ];
  }

 CrrBDim = Newn;
 
 }  // end( MinQuad::SetCrrBDim )
 
/*--------------------------------------------------------------------------*/

void MinQuad::SetEpsilonR( HpNum NeweR )
{
 eR = NeweR ? NeweR : IntEpsR;
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR ADDING / REMOVING / CHANGING DATA -------------*/
/*--------------------------------------------------------------------------*/

void MinQuad::AddSubGrad( cIndex n , cHpNum alfan )
{
 if( ( ActBDim >= CrrBDim ) || ( n >= CrrBDim ) || G[ n ] )
  return;

 Alfa[ n ] = alfan;
 G[ n ] = SubG;

 #if( VARCOEFF )
  cf[ n ] = 1;
 #endif

 // the name (index in G and Q) of the newly entered item is inserted in
 // Order[]: items are checked for optimality from the end of Order upwards

 Order[ ActBDim++ ] = n;
 Order[ ActBDim ] = InINF;

 if( n >= NxtBIdx )
  NxtBIdx = n + 1;

 #if( ! LAZY_Q )
  GiTG( n , tempQ , NxtBIdx );

  register QuRow tQn = Q[ n ];
  register QuRow ttQ = tempQ;
  register Index i = 0;

  for( ; i <= n ; i++ )
   *(tQn++) = *(ttQ++);

  register char *tG = G + i;
  register QuMat tQ = Q + i;

  for( ; i++ < NxtBIdx ; tQ++ , ttQ++ )
   if( *(tG++) )
    (*tQ)[ n ] = *ttQ;
 #endif

 #if( LOG_MQ > 3 )
  *MQLog << "Added new subgradient (" << n << "): #Bundle = " << ActBDim
         << endl;
 #endif

 AddedOne = TRUE;   // signal that a new item has been added

 }  // end( MinQuad::AddSubGrad )

/*--------------------------------------------------------------------------*/

void MinQuad::AddConstr( cIndex n, cHpNum alfan )
{
 if( ( ActBDim >= CrrBDim ) || ( n >= CrrBDim ) || G[ n ] )
  return;

 Alfa[ n ] = alfan;
 G[ n ] = Cnst;

 #if( VARCOEFF )
  cf[ n ] = 0;
 #endif

 Order[ ActBDim++ ] = n;
 Order[ ActBDim ] = InINF;
 ActCNum++;

 if( n >= NxtBIdx )
  NxtBIdx = n + 1;

 #if( ! LAZY_Q )
  GiTG( n , tempQ , NxtBIdx );

  register QuRow tQn = Q[ n ];
  register QuRow ttQ = tempQ;
  register Index i = 0;

  for( ; i <= n ; i++ )
   *(tQn++) = *(ttQ++);

  register char *tG = G + i;
  register QuMat tQ = Q + i;

  for( ; i++ < NxtBIdx ; tQ++ , ttQ++ )
   if( *(tG++) )
    (*tQ)[ n ] = *ttQ;
 #endif

 #if( LOG_MQ > 3 )
  *MQLog << "Added new constraint (" << n << "): #Bundle = " << ActBDim
         << endl;
 #endif

 AddedOne = TRUE;

 }  // end( MinQuad::AddConstr )

/*--------------------------------------------------------------------------*/

void MinQuad::ResetBundle( cIndex which )
{
 if( ( which >= CrrBDim ) || ( ! G[ which ] ) )  // invalid name
  return;                                        // quitely return

 if( G[ which ] & IsInBse )
 {
  register Index j = 0;
  while( Base[ j ] != which )  // find its position ..
   j++;

  CutSubGradFromBase( j );     // .. end eliminate it
  Minf = f = HpINF;            // Mult[] must be changed
  }

 register Index j = 0;
 while( Order[ j ] != which )  // find the position in Order
  j++;

 ShiftVect( Order + j , (ActBDim--) - j );

 #if( LAZY_Q )
  j = which;

  for( register QuRow pq = Q[ j++ ] ; j-- ; )
   *(pq++) = -QuINF;

  for( j = NxtBIdx ; j > which ; )
   Q[ --j ][ which ] = -QuINF;
 #endif

 #if( LOG_MQ > 3 )
  *MQLog << "Eliminated " << ( G[ which ] & IsACnst ? "constraint (" :
            "subgradient (" ) << which << "): #Bundle = " << ActBDim << endl;
 #endif

 if( G[ which ] & IsACnst )
  ActCNum--;

 G[ which ] = Void;

 while( NxtBIdx && ( G[ NxtBIdx - 1 ] == Void ) )
  NxtBIdx--;

 }  // end( ResetBundle( Index ) )

/*--------------------------------------------------------------------------*/

void MinQuad::ResetBundle( void )
{
 #if( LOG_MQ > 3 )
  *MQLog << "Complete reset of the bundle. [ResetBundle]" << endl;
 #endif

 for( register Index_Set o = Order ; ActBDim ; ActBDim-- )
 {
  #if( LAZY_Q )
   register Index j = *(o++);
   G[ j ] = Void;

   register QuRow pq = Q[ j++ ];

   for( ; j-- ; )
    *(pq++) = -QuINF;
  #else
   G[ *(o++) ] = Void;
  #endif
  }

 ActCNum = CBDim = BDim = Dependent = LBIsInB = NxtBIdx = 0;
 *Base = *Order = InINF;
 AddedOne = FALSE;

 }  // end( ResetBundle( void ) )

/*--------------------------------------------------------------------------*/

void MinQuad::ChangeAlfa( cHpNum DeltaAlfa )
{
 register Index_Set tO = Order;
 if( CBDim )
 {
  for( register Index n ; ( n = *(tO++) ) < InINF ; )
   if( ! ( G[ n ] & IsACnst ) )
    Alfa[ n ] += DeltaAlfa;
  }
 else
  for( register Index n ; ( n = *(tO++) ) < InINF ; )
   Alfa[ n ] += DeltaAlfa;

 VectSum( z2 , z1 , DeltaAlfa , BDim - Dependent );

 #if( ! EXACT )
  z1Tz2 += DeltaAlfa * z1Norm;
 #endif

 Lin = HpINF;  // Alfa is changed (but z2 is correct)

 #if( LOG_MQ > 3 )
  *MQLog << "Alfa shifted of " << DeltaAlfa << "." << endl;
 #endif

 }  // end( ChangeAlfa( HpNum , HpRow ) )

/*--------------------------------------------------------------------------*/

void MinQuad::MoveAlongD( cHpNum Tau , cHpNum DeltaFi )
{
 register Index_Set tO = Order;
 register HpNum dfR = DeltaFi + Tau * LastRo;

 if( Tau == PrvsTi )  //- - - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  for( register Index n ; ( n = *(tO++) ) < InINF ; )
   if( G[ n ] & IsInBse )
    Alfa[ n ] = e( n , dfR );
   else
    Alfa[ n ] += PrvsTi * GTd[ n ] + e( n , DeltaFi );

  if( LwrBnd < HpINF )
   if( LBIsInB )
    LwrBnd = dfR;
   else
    LwrBnd += DeltaFi;

  VectAssign( z2 , z1 , dfR , BDim - Dependent );

  #if( ! EXACT )
   z1Tz2 = dfR * z1Norm;
  #endif
  }
 else  // Tau != PrvsTi - - - - - - - - - - - - - - - - - - - - - - - - - - -
 {
  register HpNum tauR = 1 - Tau / PrvsTi;

  for( register Index n ; ( n = *(tO++) ) < InINF ; )
   if( G[ n ] & IsInBse )
    Alfa[ n ] = Alfa[ n ] * tauR + e( n , dfR );
   else
    Alfa[ n ] += Tau * GTd[ n ] + e( n , DeltaFi );

  if( LwrBnd < HpINF )
   if( LBIsInB )
    LwrBnd = LwrBnd * tauR + dfR;
   else
    LwrBnd += DeltaFi;

  VectScale( z2 , tauR , BDim - Dependent );
  VectSum( z2 , z1 , dfR , BDim - Dependent );

  #if( ! EXACT )
   z1Tz2 = tauR * z1Tz2 + dfR * z1Norm;
  #endif

  }  // end else( Tau != PrvsTi ) - - - - - - - - - - - - - - - - - - - - - -

 Lin = HpINF;  // Alfa is changed (but z2 is correct)

 #if( LOG_MQ > 3 )
  *MQLog << endl << "A step of " << Tau << " along d was taken.";
 #endif

 }  // end( MoveAlongD )

/*--------------------------------------------------------------------------*/

void MinQuad::ReadAlfa( HpRow NewAlfa )
{
 VectAssign( NewAlfa , Alfa , NxtBIdx );

 }  // end( MinQuad::ReadAlfa( HpRow ) )

/*--------------------------------------------------------------------------*/

void MinQuad::AddSGSpaceDim( cSgRow NewDim )
{
 // the temporary tmpv is used

 #if( LOG_MQ > 1 )
  SGSIncr++;
 #endif

 // update Q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 register Index i;
 for( register Index_Set tO = Order ; ( i = *(tO++) ) < InINF ; )
 {
  register SgNum NDi = NewDim[ i ];

  cond_if( NDi != 0 )
   for( register cIndex_Set tO2 = Order ; tO2 < tO ; )
   {
    register Index j = *(tO2++);
    if( i >= j )
    {
     Qif( Q[ i ][ j ] > - QuINF )
      Q[ i ][ j ] += NDi * NewDim[ j ];
     }
    else
     Qif( Q[ j ][ i ] > - QuINF )
      Q[ j ][ i ] += NDi * NewDim[ j ];
    }
  }

 #if( ! EXACT )
  if( Dependent > 1 )
  {
   LhNorm = LowQ( Base[ BDim - 2 ] );
   LjNorm = LowQ( Base[ BDim - 1 ] );
   }
  else
   if( Dependent )
    LhNorm = LowQ( Base[ BDim - 1 ] );
 #endif

 // update L and all the rest - - - - - - - - - - - - - - - - - - - - - - - -

 for( register Index i = BDim ; i-- ; )
  tmpv[ i ] = NewDim[ Base[ i ] ];
 
 HpNum Xsi = 0;
 HpNum Sigma = 0;
 cIndex BDMD = BDim - Dependent;
 for( register Index i = 0 ; i < BDMD ; i++ )
 {
  HpNum ai = tmpv[ i ];

  if( ABS( ai ) <= eR * BDim )  // if it's already 0, nothing to do
   continue;

  HpNum aj = L[ i ][ i ];
  HpNum d = L[ i ][ i ] = sqrt( ai * ai + aj * aj );

  register HpNum Gc = aj / d;    // Gc * ai + Gs * aj = 0
  register HpNum Gs = - ai / d;  // Gc ^ 2 + Gs ^ 2 = 1

  for( register Index j = i ; ++j < BDim ; )
   ApplyGivensMToCol( Gc , Gs , tmpv[ j ] , L[ j ][ i ] );

  ApplyGivensMToCol( Gc , Gs , Xsi   , z1[ i ] );
  ApplyGivensMToCol( Gc , Gs , Sigma , z2[ i ] );
  }

 #if( ! EXACT )
  z1Norm -= Xsi * Xsi;
  z1Tz2 -= Xsi * Sigma;
  // z2Norm -= Sigma * Sigma;
 #endif

 // check for gained 1 lin. independent - - - - - - - - - - - - - - - - - - -

 ChkNwLinInd();

 Quad = HpINF;  // f() must be changed: the "new" Quad is just
                // Quad + ( Sum{ i } NewDim[ Base[ i ] * Mult[ i ] )^2,
                // but it is not worth to keep it updated

 }  // end( MinQuad::AddSGSpaceDim( cSgRow ) )

/*--------------------------------------------------------------------------*/

void MinQuad::AddSGSpaceDim( cSgRow NewDim , cHpNum lbh )
{
 // the temporaries tmpv and tmpa are used

 if( ABS( lbh ) <= HpeM )
 {
  AddSGSpaceDim( NewDim );
  return;
  }

 #if( LOG_MQ > 1 )
  SGSIncr++;
 #endif

 // update Q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 register Index i;
 for( register cIndex_Set tO = Order ; ( i = *(tO++) ) < InINF ; )
 {
  register SgNum NDi = NewDim[ i ];

  cond_if( NDi != 0 )
   for( register cIndex_Set tO2 = Order ; tO2 < tO ; )
   {
    register Index j = *(tO2++);
    if( i >= j )
    {
     Qif( Q[ i ][ j ] > - QuINF )
      Q[ i ][ j ] += NDi * NewDim[ j ];
     }
    else
     Qif( Q[ j ][ i ] > - QuINF )
      Q[ j ][ i ] += NDi * NewDim[ j ];
    }
  }

 // update L and all the rest - - - - - - - - - - - - - - - - - - - - - - - -

 for( register Index i = BDim ; i-- ; )
  tmpv[ i ] = NewDim[ Base[ i ] ];

 HpNum Xsi = 0;
 HpNum Sigma = 0;
 HpNum AdjRow0 = 1;   
 cIndex BDMD = BDim - Dependent;
 for( register Index i = 0 ; i < BDMD ; i++ )
 {
  register HpNum ai = tmpv[ i ];

  if( ABS( ai ) <= eR * BDim )
   tmpa[ i ] = 0;
  else
  {
   HpNum aj = L[ i ][ i ];
   HpNum d = L[ i ][ i ] = sqrt( ai * ai + aj * aj );

   register HpNum Gc = aj / d;
   register HpNum Gs = - ai / d;

   for( register Index j = i ; ++j < BDim ; )
    ApplyGivensMToCol( Gc , Gs , tmpv[ j ] , L[ j ][ i ] );

   ApplyGivensMToCol( Gc , Gs , Xsi   , z1[ i ] );
   ApplyGivensMToCol( Gc , Gs , Sigma , z2[ i ] );

   tmpa[ i ] = - AdjRow0 * Gs;
   AdjRow0 *= Gc;
   }
  }

 #if( ! EXACT )
  Sigma += lbh * AdjRow0;

  z1Norm -= Xsi * Xsi;
  z1Tz2 -= Xsi * Sigma;
  // z2Norm -= Sigma * Sigma;
 #endif

 // modify Alfa: this MUST be done before checking for gained 1 lin.- - - - -
 // dependent - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 register Index h;
 for( register cIndex_Set tO = Order ; ( h = *(tO++) ) < InINF ; )
  Alfa[ h ] += NewDim[ h ] * lbh;

 VectSum( z2 , tmpa , lbh , BDim - Dependent );

 // check for gained 1 lin. independent - - - - - - - - - - - - - - - - - - -

 ChkNwLinInd();

 Quad = HpINF;  // f() must be changed

 }  // end( MinQuad::AddSGSpaceDim( cSgRow , HpNum ) )

/*--------------------------------------------------------------------------*/

void MinQuad::CutSGSpaceDim( cSgRow OldDim )
{
 // the temporary tmpv is used, and MoveSubgradToLastPos() may be invoked

 #if( LOG_MQ > 1 )
  SGSDecr++;
 #endif

 // update Q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 register Index i;
 for( register cIndex_Set tO = Order ; ( i = *(tO++) ) < InINF ; )
 {
  register cSgNum ODi = OldDim[ i ];

  cond_if( ODi != 0 )
   for( register cIndex_Set tO2 = Order ; tO2 < tO ; )
   {
    register cIndex j = *(tO2++);
    if( i >= j )
    {
     Qif( Q[ i ][ j ] > - QuINF )
      Q[ i ][ j ] -= ODi * OldDim[ j ];
     }
    else
     Qif( Q[ j ][ i ] > - QuINF )
      Q[ j ][ i ] -= ODi * OldDim[ j ];
    }
  }

 // update L and all the rest - - - - - - - - - - - - - - - - - - - - - - - -

 for( register Index i = BDim ; i-- ; )
  tmpv[ i ] = OldDim[ Base[ i ] ];

 HpNum Xsi = 0;
 HpNum Sigma = 0;

 for( register Index j = 0 ; j < BDim - Dependent ; )
 {
  HpNum ai = tmpv[ j ];

  if( ABS( ai ) <= eR * BDim )
   j++;
  else
  {
   HpNum aj = L[ j ][ j ];
   HpNum d = aj * aj - ai * ai;

   if( d <= size( LowQ( j ) ) )       // checking for lost 1 lin. indep.
   {                                  // (ignore Delta < 0)
    if( Dependent >= 2 - LBIsInB )    // this one would be in excess
     CutSubGradFromBase( BDim - 1 );  // very easy task

    register Index h = BDim - Dependent - 1;

    if( j < h )                       // not the last (independent) row
    {                                 // move it in last position, so that
     MoveSubgradToLastPos( j );       // it can be easily eliminated
     RotateVect( tmpv + j , h - j ); 
     }
    #if( ! EXACT )
     else
      h = j;

     z1Norm -= z1[ h ] * z1[ h ];
     z1Tz2 -= z1[ h ] * z2[ h ];
 
     if( Dependent )                  // already ...
      LjNorm = LowQ( Base[ BDim - 1 ] );

     LhNorm = LowQ( Base[ j ] );
    #endif

    Dependent++;
    }
   else  // d != 0
   {
    L[ j ][ j ] = d = sqrt( d );

    register HpNum Gc = aj / d;    // Gc * ai + Gs * aj = 0
    register HpNum Gs = - ai / d;  // Gc ^ 2 - Gs ^ 2 = 1

    for( register Index i = j ; ++i < BDim ; )
     ApplyCompGivensMToCol( Gc , Gs , tmpv[ i ] , L[ i ][ j ] );

    ApplyCompGivensMToCol( Gc , Gs , Xsi   , z1[ j ] );
    ApplyCompGivensMToCol( Gc , Gs , Sigma , z2[ j++ ] );

    }  // end else( d = 0 )
   }  // end if( ai != 0 )
  }  // end for( j )

 #if( ! EXACT )
  z1Norm += Xsi * Xsi;
  z1Tz2 += Xsi * Sigma;
  // z2Norm += Sigma * Sigma;

  if( Dependent > 1 )  
  {
   LhNorm = LowQ( Base[ BDim - 2 ] );
   LjNorm = LowQ( Base[ BDim - 1 ] );
   }
  else
   if( Dependent )
    LhNorm = LowQ( Base[ BDim - 1 ] );
 #endif

 // recalculate Lmu - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LiMax = LiMin = 0;
 UpdateLmu( 1 );

 Quad = HpINF;  // f() must be changed: the "new" Quad is just
                // Quad - ( Sum{ i } NewDim[ Base[ i ] * Mult[ i ] )^2,
                // but it is not worth to keep it updated

 }  // end( CutSGSpaceDim( cSgRow ) )

/*--------------------------------------------------------------------------*/

void MinQuad::CutSGSpaceDim( cSgRow OldDim , cHpNum lbh )
{
 // uses temporaries tmpv and tmpa and may call MoveSubgradToLastPos()

 if( ABS( lbh ) <= HpeM )
 {
  CutSGSpaceDim( OldDim );
  return;
  }

 #if( LOG_MQ > 1 )
  SGSDecr++;
 #endif

 // update Q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 register Index i;
 for( register cIndex_Set tO = Order ; ( i = *(tO++) ) < InINF ; )
 {
  register cSgNum ODi = OldDim[ i ];

  cond_if( ODi != 0 )
   for( register cIndex_Set tO2 = Order ; tO2 < tO ; )
   {
    register cIndex j = *(tO2++);
    if( i >= j )
    {
     Qif( Q[ i ][ j ] > - QuINF )
      Q[ i ][ j ] -= ODi * OldDim[ j ];
     }
    else
     Qif( Q[ j ][ i ] > - QuINF )
      Q[ j ][ i ] -= ODi * OldDim[ j ];
    }
  }

 #if( ! EXACT )
  if( Dependent > 1 )
  {
   LhNorm = LowQ( Base[ BDim - 2 ] );
   LjNorm = LowQ( Base[ BDim - 1 ] );
   }
  else
   if( Dependent )
    LhNorm = LowQ( Base[ BDim - 1 ] );
 #endif

 // update Alfa - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 register Index h;
 for( register cIndex_Set tO = Order ; ( h = *(tO++) ) < InINF ; )
  Alfa[ h ] -= OldDim[ h ] * lbh;

 // update L and all the rest - - - - - - - - - - - - - - - - - - - - - - - -

 for( register Index i = BDim ; i-- ; )
  tmpv[ i ] = OldDim[ Base[ i ] ];

 HpNum Xsi = 0;
 HpNum Sigma = 0;
 HpNum AdjRow0 = 1;   
 for( register Index j = 0 ; j < BDim - Dependent ; )
 {
  HpNum ai = tmpv[ j ];

  if( ABS( ai ) <= eR * BDim )
   tmpa[ j++ ] = 0;
  else
  {
   HpNum aj = L[ j ][ j ];
   HpNum d = aj * aj - ai * ai;

   if( d <= size( LowQ( j ) ) )       // checking for lost 1 lin. indep.
   {                                  // (ignore Delta < 0)
    if( Dependent >= 2 - LBIsInB )    // this one would be in excess
     CutSubGradFromBase( BDim - 1 );  // very easy

    register Index h = BDim - Dependent - 1;

    if( j < h )                       // not the last (independent) row
    {                                 // move it in last position, so that
     MoveSubgradToLastPos( j );       // it can be easily eliminated
     RotateVect( tmpv + j , h - j ); 
     RotateVect( tmpa + j , h - j ); 
     }
    #if( ! EXACT )
     else
      i = j;

     z1Norm -= z1[ h ] * z1[ h ];
     z1Tz2 -= z1[ h ] * z2[ h ];
 
     if( Dependent )                  // already ...
      LjNorm = LowQ( Base[ BDim - 1 ] );

     LhNorm = LowQ( Base[ j ] );
    #endif

    Dependent++;
    }
   else
   { 
    d = L[ j ][ j ] = sqrt( d );

    register HpNum Gc = aj / d;
    register HpNum Gs = - ai / d;

    for( register Index i = j ; ++i < BDim ; )
     ApplyCompGivensMToCol( Gc , Gs , tmpv[ i ] , L[ i ][ j ] );

    ApplyCompGivensMToCol( Gc , Gs , Xsi   , z1[ j ] );
    ApplyCompGivensMToCol( Gc , Gs , Sigma , z2[ j ] );

    tmpa[ j++ ] = AdjRow0 * Gs;
    AdjRow0 *= Gc;

    }  // end else( d = 0 )
   }  // end if( ai != 0 )
  }  // end for( j )

 #if( ! EXACT )
  Sigma += lbh * AdjRow0;

  z1Norm += Xsi * Xsi;
  z1Tz2 += Xsi * Sigma;
  // z2Norm += Sigma * Sigma;
 #endif

 // recalculate Lmu - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 LiMax = LiMin = 0;
 UpdateLmu( 1 );

 // modify z2 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectSum( z2 , tmpa , lbh , BDim - Dependent );

 Quad = HpINF;  // f() must be changed

 }  // end( CutSGSpaceDim( cSgRow , HpNum ) )

/*--------------------------------------------------------------------------*/

void MinQuad::ChangeQ( cBOOL KeepBase )
{
 #if( LAZY_Q )
  for( register Index i = NxtBIdx ; i ; )  // all the elements of Q[] must
  {                                        // be recomputed, if needed
   register Index j = i;
   for( register QuRow pq = Q[ --i ] ; j-- ; )
    *(pq++) = -QuINF;
   }
 #else
  for( register Index i = 0 ; i < NxtBIdx ; i++ )  // recompute *now* all
   GiTG( i , Q[ i ] , i );                         // the elements of Q[]
 #endif

 if( KeepBase )   // attempt a "warm start" anyway
 {
  RebuildBase();
  Quad = HpINF;
  }
 else             // resign to a "cold start"
 {
  register cIndex_Set tB = Base;
  for( register Index h ; ( h = *(tB++) ) < InINF ; )
   G[ h ] &= ~IsInBse;

  BDim = CBDim = Dependent = LBIsInB = 0;
  *Base = InINF;
  }
 }  // end( ChangeQ )

/*--------------------------------------------------------------------------*/

#if( VARCOEFF )

void MinQuad::ChangeCoefficient( cIndex i , cHpNum ci )
{
 if( ( i < CrrBDim ) && ( cf[ i ] != ci ) && G[ i ] )
  if( ci >= 0 )
  {
   if( ( ci == 0 ) && ( ! ( G[ i ] & IsACnst ) ) )
   {
    G[ i ] |= IsACnst;
    ActCNum++;

    if( G[ i ] & IsInBse )
     CBDim++;
    }

   if( ( ci > 0 ) && ( G[ i ] & IsACnst ) )
   {
    ActCNum--;

    if( G[ i ] & IsInBse )
     CBDim--;

    G[ i ] &= ~IsACnst;
    }

   cf[ i ] = ci;

   if( G[ i ] & IsInBse )
    f = z1Norm = HpINF;  // Mult and z1 must be changed
   }
 }  // end( ChangeCoefficient )

#endif

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

MinQuad::MQError MinQuad::SolveQP( HpNum ti )
{
 #if( TIMERS_MQ )
  MQt->Start();
 #endif

 HpNum Lastf = HpINF;
 BOOL First = TRUE;
 BOOL AugEps = FALSE;
 BOOL DecrEps = FALSE;

 for(;;)
 {
  // solve the problem- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  CalcOptDir( ti );

  if( ( ! QPStatus ) || ( QPStatus == kQPPrimUnbndd ) )
   break;

  if( f < HpINF )
   if( Lastf > f + max( eR * ABS( f ) , HpeM ) )
   {
    Lastf = f;     // reset the counters if there have been at least one
    First = TRUE;  // strictly decreasing step
    }

  if( First )  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  {            // rebuild Base[] without adjusting eR
   RebuildBase();
   First = FALSE;
   continue;
   }

  // try to tune eR - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( QPStatus == kLoop )
  {
   eR /= MFactor;  // eR decrease (increase precision)
   DecrEps = TRUE;

   #if( LOG_MQ > 2 )
    *MQLog << endl << "Decreasing eR (= " << eR << "). [SolveQP]";
   #endif
   }
  else
  {
   eR *= MFactor;  // eR increase (decrease precision)
   AugEps = TRUE;

   #if( LOG_MQ > 2 )
    *MQLog << endl << "Increasing eR (= " << eR << "). [SolveQP]";
   #endif
   }

  TLDim = 0;            // since eR is changed, the taboo list can be cleaned
  Minf = HpINF;         // ... and the old value of Minf shall not be trusted
  LocalOptima = FALSE;  // ... but CalcOptDir() must do something

  if( AugEps && DecrEps )
  {
   #if( LOG_MQ )
    *MQLog << endl << "ERROR: augment and decrease of eR. [SolveQP]";
   #endif

   QPStatus = kFatal;
   break;
   }

  if( ( eR > MaxEpsR ) || ( eR < MinEpsR ) )
  { 
   #if( LOG_MQ )
    *MQLog << endl << "ERROR: eR (" << eR << ") out of limits. [SolveQP]";
   #endif

   QPStatus = kFatal;
   break;
   }
  }  // end while( problems ) - - - - - - - - - - - - - - - - - - - - - - - -

 #if( TIMERS_MQ )
  MQt->Stop();
 #endif

 return( QPStatus );
 
 }  // end( SolveQP )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR READING RESULTS ---------------------*/
/*--------------------------------------------------------------------------*/

void MinQuad::ReadGTd( HpRow g )
{
 register Index_Set tO = Order;
 for( register Index n ; ( n = *(tO++) ) < InINF ; )
  g[ n ] = GTd[ n ];

 }  // end( MinQuad::ReadGTd )

/*--------------------------------------------------------------------------*/

void MinQuad::SensitAnals1( HpNum &v1 , HpNum &v2 , HpNum &v3 )
{
 if( Dependent == 0 )  // "standard" base
  if( LBIsInB )        // ... but with the LB in: actually a nonstandard base
  {
   v1 = 0;             // same formula as Dependent == 1 with LwrBnd == Alfa
   v2 = - LwrBnd;      // and LhTz1 == LhTz2 == 0 (as Lh == 0)
   v3 = v2 * v2 * z1Norm + 2 * v2 * z1Tz2 + Norm( z2 , BDim );
   }
  else                 // a "real" standard base
  {
   v1 = - 1 / z1Norm;
   v2 = v1 * z1Tz2;
   v3 = Norm( z2 , BDim ) + v2 * z1Tz2;
   }
 else                  // "nonstandard" base
 {
  v1 = 0;
  v2 = - ( Alfa[ Base[ BDim - 1 ] ] - LhTz2 )
           / ( e( Base[ BDim - 1 ] , 1 ) - LhTz1 );

  v3 = v2 * v2 * z1Norm + 2 * v2 * z1Tz2 + Norm( z2 , BDim - 1 );
  }
 }  // end( MinQuad::SensitAnals1 )

/*--------------------------------------------------------------------------*/

void MinQuad::SensitAnals2( HpNum v1 , HpNum v2 , HpRow x1 , HpRow x2 )
{
 if( ( Dependent == 0 ) && ( ! LBIsInB ) )  // "standard" base
 {
  /* x(t) = L{-1}{T} * ( ro(t) * z1 - z2(t) ) =
     = L{-1}{T} * [ ( ro1 + ro2 / t ) * z1 - ( 1 / t ) * z2 ) =
     = L{-1}{T} * [ ro1 * z1 + ( 1 / t ) * [ ro2 * z1 - z2 ] ]
     where ro{i} = - v{i}. */

  SolveSystemT( L , z1 , x1 , BDim );
  VectAssign( x2 , x1 , BDim );
  VectScale( x1 , -v1 , BDim );
  VectScale( x2 , -v2 , BDim );
  SolveSystemT( L , z2 , tmpa , BDim );
  VectSubtract( x2 , tmpa , BDim );
  }
 else                                       // "nonstandard" base
 {
  HpNum ex1 = 1;
  HpNum ex2 = z1Tz2 + v2 * z1Norm;

  if( ! LBIsInB )
  {
   ex1 /= e( Base[ BDim - 1 ] , 1 ) - LhTz1;
   ex2 /= e( Base[ BDim - 1 ] , 1 ) - LhTz1;

   x1[ BDim - 1 ] = ex1;
   x2[ BDim - 1 ] = ex2;
   }

  if( BDim )
  {
   SolveSystemT( L , L[ BDim - 1 ] , x1 , BDim - 1 );
   VectAssign( x2 , x1 , BDim  - 1 );
   VectScale( x1 , - ex1 , BDim - 1 );
   VectScale( x2 , - ex2 , BDim - 1 );
   SolveSystemT( L , z2 , tmpa , BDim  - 1 );
   VectSubtract( x2 , tmpa , BDim - 1 );
   SolveSystemT( L , z1 , tmpa , BDim - 1 );
   VectScale( x1 , -v2 , BDim - 1 );
   VectSum( x2 , tmpa , BDim - 1 );
   }
  }
 }  // end( MinQuad::SensitAnals2 )

/*--------------------------------------------------------------------------*/

void MinQuad::SensitAnals3( cHpNum v1 , cHpNum v2 , HpRow x1 , HpRow x2 ,
                            HpNum &tMin , HpNum &tMax )
{
 /* Since v( t ) = v1 + ( 1 / t ) * v2 = - ro( t ) = - ro1 - ro2 / t, the
    j-th constraint out of base (relative to a subgradient) is saisfied for
    all t such that

     t * ( ro1 - [ Qj ] * x1 ) <= ( Alfa[ j ] + [ Qj ] * x2 - ro2 )

                   ||                   ||

     t *          M1j          <=       M2j

    Therefore, M1j > 0 => t <= M2j / M1j for all t; if the converse holds,
    t >= M2j / M1j. Since a feasible t (the current one) must exist, it is
    impossible that M1j > 0 and M2j < 0, while it is possible that
    M1j <= 0 && M2j >= 0.
    If the j-th item (out of base) is a constraint, the formula is rather

     t * ( - [ Qj ] * x1 ) <= [ Qj ] * x2

                  ||                ||

     t *          M1j      <=       M2j    */

 tMin = 0;
 tMax = HpINF;
 cHpNum av1 = ABS( v1 );
 cHpNum av2 = ABS( v2 );

 register Index_Set tO = Order;
 for( register Index h ; ( h = *(tO++) ) < InINF ; )
  if( ! ( G[ h ] & IsInBse ) )
  {
   GjTGBk( h , Base , BDim , tmpa );

   HpNum tmp1 = ScalarProduct( x1 , tmpa , BDim );
   HpNum tmp2 = ScalarProduct( x2 , tmpa , BDim );

   if( G[ h ] & IsACnst )
    AdjusttBounds( -tmp1 , tmp2 , size( tmp1 ) , size( tmp2 ) , tMin , tMax );
   else
    AdjusttBounds( - v1 - tmp1 , Alfa[ h ] + tmp2 + v2 ,
		   gsize( av1 + ABS( tmp1 ) ) ,
                   gsize( ABS( Alfa[ h ] ) + ABS( tmp2 ) + av2 ) ,
                   tMin , tMax );
   }

 /* Do the same for the Lower Bound: tmpa == 0, which mplies .... */

 if( LBIsInB )
  AdjusttBounds( - v1 , LwrBnd + v2 , gsize( av1 ) ,
		 gsize( ABS( LwrBnd ) + av2 ) , tMin , tMax );

 /* Values of t that make some x[ i ] < 0 must also be avoided, i.e.

     x1[ i ] + ( 1 / t ) * x2[ i ] >= 0 ==> t * ( - x1[ i ] ) <= x2[ i ]

    but this is exactly the same situation as above.  */

 for( register Index i = BDim ; i-- ; )
  AdjusttBounds( - x1[ i ] , x2[ i ] , stdsize() , stdsize() , tMin , tMax );

 }  // end( MinQuad::SensitAnals3 )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

MinQuad::~MinQuad( void )
{
 // first: print out information - - - - - - - - - - - - - - - - - - - - - -

 #if( LOG_MQ > 2 )
  if( Success <= 0 )
   *MQLog << endl << "No succesfully completed calls." << endl;
  else
  {
   *MQLog << endl << "Calls " << Calls << " (faults " << Calls - Success
          << ") ~ Steps " << Step << endl
          << "Average dimensions: base " << SumBDim / Step << ", bundle "
          << SumActBDim / Calls << "." << endl
          << "% of constraints: base " << SumGSBase / Step * 100
          << ", bundle " << SumGSBundle / Calls * 100 << endl
          << "Total operations: " << Insert << " insertions, " << Delete
          << " deletions." << endl;

   if( SGSIncr || SGSDecr )
    *MQLog << "SGSpace modifications: " << SGSIncr << " augments, "
           << SGSDecr << " decreases." << endl;

   #if( TIMERS_MQ )
    *MQLog << "SolveQP() time: user  " << MQt->u << ", system " << MQt->s
           << " seconds." << endl;
   #endif
   }
 #elif( LOG_MQ > 1 )
  *MQLog << Calls << "\t" << SumBDim / Step << "\t" << SumActBDim / Calls
   #if( TIMERS_MQ )
         << "\t" << ( MQt->u + MQt->s )
   #endif
         << endl;
 #endif

 // second: memory deallocation - - - - - - - - - - - - - - - - - - - - - - -

 #if( TIMERS_MQ )
  delete MQt;
 #endif

 if( MaxBDim )
  MemDealloc();

 }  // end( ~MinQuad )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

inline void MinQuad::GjTGBk( cIndex j , register cIndex_Set Bk ,
			     register Index Dim , register HpRow v )

/* v = G_j{T} * G{Bk}, i.e the j-th row of matrix Q, but only those elements
   whose "names" are listed in Bk, and with the same order. */
{
 for( ; Dim-- ; )
  *(v++) = LowQ( j , *(Bk++) );
 }

/*--------------------------------------------------------------------------*/

void MinQuad::RebuildBase( void )

/* Emergency procedure: rebuilds a feasible Base and the relative data
   structures (L, Base, z1 ...) from scratch, eliminating the problems caused
   by the accumulation of rounding errors. "Negative Delta" errors are simply
   disregarded, but the corresponding items are discarded from Base[]. */
{
 // first reset the old Base- - - - - - - - - - - - - - - - - - - - - - - - -
 // note that LBIsInB is not touched: it cannot cause numerical problems

 Lmu = 1;
 LiMax = LiMin = BDim = CBDim = Dependent = 0;

 #if( ! EXACT )
  z1Norm = z1Tz2 = 0;
 #endif

 register Index_Set bp = Base;
 register Index i;
 for( ; ( i = *(bp++) ) < InINF ; )
  G[ i ] &= ~IsInBse;

 // and then reconstruct it - - - - - - - - - - - - - - - - - - - - - - - - -

 Index SGDltd = 0;
 #if( LOG_MQ > 4 )
  Index CNDltd = 0;
 #endif

 register HpRow mp = Mult;
 for( bp = Base ; ( i = *(bp++) ) < InINF ; )
 {
  Base[ BDim ] = i;
  Mult[ BDim ] = *(mp++);

  GjTGBk( i , Base , BDim - Dependent , tmpa );
  SolveSystem( L , tmpa , L[ BDim ] , BDim - Dependent );

  HpNum Delta = CalculateDelta( i );
  cIndex MxDep = 1 - LBIsInB;

  if( ( Delta > 0 ) ||
      ( ( ! Delta ) && ( ( Dependent < MxDep ) || ( *bp == InINF ) ) ) )
  {
   // if the Lower Bound is *not* in Base, an item is added to Base if:
   // - it's independent ( Delta > 0 ), or
   // - it's the 1st dependent ( Delta == 0 && Dependent < 1 ), or
   // - it's the 2nd dependent, but it's the last item in the old Base
   // if the Lower Bound is in Base, the "tollerance to dependent items" is
   // decreased by 1: a dependent item can only be the last one in Base

   BDim++;
   if( G[ i ] & IsACnst )
    CBDim++;

   G[ i ] |= IsInBse;

   if( Delta > 0 )
    AddSubGradToBase( BDim - 1 , Delta );  // Lmu, z1 ... are updated inside
   else
    Dependent++;
   }
  else  // Delta < 0, or the second dependent item has arrived "too early"
  {
   if( ! ( G[ i ] & IsACnst ) )
     SGDltd++;
   #if( LOG_MQ > 4 )
    else
     CNDltd++;
   #endif

   #if( ! EXACT )
    if( ! Delta )  // LjNorm has been "trashed" by CalculateDelta()
     LjNorm = Norm( L[ BDim - 1 ] , BDim - Dependent );
   #endif
   }
  }  // end for( i )

 Base[ BDim ] = InINF;

 if( SGDltd )
  f = HpINF;     // Mult[] are unfeasible
 else
  Quad = HpINF;  // just to be more exact, recompute Quad

 #if( LOG_MQ > 4 )
  *MQLog << endl << SGDltd << " subgradient(s)" << " + " << CNDltd
	 << " constraint(s)" << " eliminated. [RebuildBase]" << endl;
 #endif

 }  // end( RebuildBase )

/*--------------------------------------------------------------------------*/

inline void MinQuad::CalcOptDir( cHpNum ti )

/* Solves the problem with the current data and the given value of ti.
   Set QPStatus to:

    kOK             = Normal termination: this is reported also if the Bundle
                      is empty, the resulting "optimal" base being empty too.
    kQPPrimUnbndd   = The problem is primal unbounded.

   If the algorithm exits with one of the remaining codes, the reported
   solution is primal feasible, but in general not optimal, i.e. the resulting
   direction d is not (in principle) *neither decreasing nor (dual) feasible*.

    kNegativeDelta  = Delta < 0.
    kNegVarInBase   = A variabile <= 0 is in Base.
    kInvalidV       = v[ i ] >= 0 for each i in Base.

   These errors are usually caused by numerical problems in the data, or by
   the accumulation of rounding errors in the data structures: the caller
   should call RebuildBase() and/or increase eR (so to relax the precision
   requirements).

    kLoop           = Loop detected.

   This error is usually caused by a wrong choice of eR: simply, eR is too
   big, so the caller should decrease it. */
{
 QPStatus = kOK;

 #if( LOG_MQ > 5 )
  CheckQValidity();  // check Q
 #endif

 if( BDim || LBIsInB )  // warm start - - - - - - - - - - - - - - - - - - - -
 {                      //- - - - - - - - - - - - - - - - - - - - - - - - - -
                        // if none among Base, Alfa, ti, eR and the Bundle
                        // has changed, this is (still) a restricted optima
  #if( LOG_MQ > 3 )
   *MQLog << endl << "Warm start";
  #endif

  if( f == HpINF )  // Base[] is changed: Mult[] may not be feasible- - - - -
  {
   LocalOptima = FALSE;

   #if( LOG_MQ > 3 )
    *MQLog << " ~ diff. B";
   #endif

   if( BDim > CBDim )  // there are subgradients in Base[]
   {
    register HpNum ex = 0;   // calculate ex = Sum{ i } e( i ) * Mult[ i ]
    register HpRow tM = Mult + BDim;

    #if( VARCOEFF )
     for( register HpRow tcf = cf + BDim ; tcf > cf ; )
      ex += *(--tM) * (*(--tcf));
    #else
     if( CBDim )
     {
      for( register Index_Set tB = Base + BDim ; tM-- > Mult ; )
       if( ! ( G[ *(--tB) ] & IsACnst ) )
        ex += *tM;
      }
     else
      for( ; tM > Mult ; )
       ex += *(--tM);
    #endif

    if( LBIsInB )                           // the Lower Bound is in Base[]
     LBMult = max( 1 - ex , HpNum( 0 ) );   // it "absorbs" the change
    else                                    // the Lower Bound is *not* in
     if( ABS( ex - 1 ) > eR * BDim )        // ex != 1
      if( ABS( ex ) > eR * BDim )           // ex != 0: ( 1 / ex ) * Mult
       VectScale( Mult , 1 / ex , BDim );   // sum to 1
      else                                  // ex == 0: any value is fine
      {
       ex = 1 / HpNum( BDim - CBDim );

       #if( VARCOEFF )
        for( tcf = cf + BDim , tM = Mult + BDim ; tM-- > Mult ; )
         if( ABS( *(--tcf) ) > eR * BDim )
          *tM = ex / *tcf;
       #else
        VectAssign( Mult , ex , BDim );
       #endif
       }
    }
   else                // there are only constraints in Base[]
    if( LBIsInB )      // ... apart from the the Lower Bound
     LBMult = 1;       // ... which therefore must have multiplier 1

   #if( VARCOEFF )
    if( z1Norm == HpINF )  // it was because the coefficients have changed
    {
     Minf = HpINF;

     SolvePerSys( L , cf , z1 , Base , BDim - Dependent );

     #if( ! EXACT )
      z1Norm = Norm( z1 , BDim );

      if( z1Tz2 < HpINF )
       z1Tz2 = ScalarProduct( z1 , z2 , BDim - Dependent );
     #endif
     }
   #endif
   }

  if( Quad == HpINF )  // Q[] is changed- - - - - - - - - - - - - - - - - - -
  {
   Minf = f = HpINF;     // signal not to check for decrease of f(), but
   LocalOptima = FALSE;  // avoid to recalculate Quad now (a costly task)

   #if( LOG_MQ > 3 )
    *MQLog << " ~ diff. Q";
   #endif
   }

  if( Lin == HpINF )  // Alfa is changed- - - - - - - - - - - - - - - - - - -
  {
   Minf = HpINF;
   LocalOptima = FALSE;

   if( z1Tz2 == HpINF )  // and also z2 is changed
   {
    SolvePerSys( L , Alfa , z2 , Base , BDim - Dependent );

    #if( ! EXACT )
     z1Tz2 = ScalarProduct( z1 , z2 , BDim - Dependent );
    #endif
    }

   #if( LOG_MQ > 3 )
    *MQLog << " ~ diff. Alfa";
   #endif
   }

  if( ti != PrvsTi )  // ti is changed- - - - - - - - - - - - - - - - - - - -
  {
   Minf = HpINF;
   LocalOptima = FALSE;

   #if( LOG_MQ > 3 )
    *MQLog << " ~ diff. t (= " << ti << ")";
   #endif
   }

  #if( LOG_MQ > 3 )
   if( AddedOne )
    *MQLog << " ~ diff. Bundle";
  #endif

  if( LocalOptima && ( ! AddedOne ) )  // the point is still optimal- - - - -
  {
   #if( LOG_MQ > 3 )
    *MQLog << " ~ Stop: latest x[] still OPTIMAL." << endl;
   #endif

   return;
   }
  }
 else  // BDIm == 0: cold start - - - - - - - - - - - - - - - - - - - - - - -
 {     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Lmu = 1;
  LiMax = LiMin = 0;
  LocalOptima = TRUE;

  LastRo = - HpINF;
  f = Minf = HpINF;
  z1Norm = z1Tz2 = 0;

  #if( LOG_MQ > 3 )
   *MQLog << endl << "Cold start";
  #endif
  }  // end else( cold start )- - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 PrvsTi = ti;       // reset PrvsTi
 Ename = InINF;     // reset the entering item
 AddedOne = FALSE;  // reset AddedOne

 #if( LOG_MQ > 1 )
  Calls++;

  #if( LOG_MQ > 3 )
   Index TmpStep = 0;  // Step counter
  #endif
 #endif

 /*-------------------------------------------------------------------------*/
 /*--------------------------- Main Loop -----------------------------------*/
 /*-------------------------------------------------------------------------*/

 do
 {
  /*-----------------------------------------------------------------------*/
  /*-------------------------- Inner Loop ---------------------------------*/
  /*-----------------------------------------------------------------------*/

  while( ! LocalOptima )
  {
   /* The inner loop begins: the loop is repeated until a restricted optima
      is found, and items can only be deleted from (but never added to) Base.
      At each iteration, the Base is kept in "standard form".

      If Ename < InINF the solution may have the form [ x 0 ], if the new item
      has just entered the Base (with 0 multiplier): in this case not all the
      directions v s.t. e*v = 0 are feasible, it must also be v[ Ewher ] > 0.

      If Base contains at least a subgradient when entering the loop
      (=> LastRo < - HpINF) then it will contain at least a subgradient in
      every iteration. */

   #if( LOG_MQ > 1 )
    Step++;
    SumBDim += BDim;
    SumGSBase += float( CBDim ) / float( max( BDim , Index( 1 ) ) );

    #if( LOG_MQ > 3 )
     *MQLog << endl << "(" << Calls << "|" << TmpStep++ << "|" << BDim << "|"
            << Dependent << ")";

     #if( LOG_MQ > 5 )
      CheckBaseValidity();  // check L, z1, z2 ...
     #endif
    #endif
   #endif

   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /*- - - - - - Find the new optimal point / descent direction. - - - - - -*/
   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   switch( Dependent + 3 * LBIsInB ) {

    case( 0 ):  // Case 0: Q{B} is nonsingular- - - - - - - - - - - - - - - -
                //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    {
     if( BDim == CBDim )
     {
      z1Tz2 = z1Norm = 0;
      LastRo = - HpINF;

      VectAssign( tmpa , z2 , - 1 / ti , BDim );  // tmpa = - z2 / ti
      }
     else
     {
      #if( EXACT )
       z1Tz2 = ScalarProduct( z1 , z2 , BDim );     
       z1Norm = Norm( z1 , BDim );
      #endif

      LastRo = ( 1 + z1Tz2 / ti ) / z1Norm;

      // tmpa = ro * z1 - z2 / ti
      VectAdd( tmpa , z1 , z2 , LastRo , - 1 / ti , BDim );
      }

     SolveSystemT( L , tmpa , tmpv , BDim );  // tmpv = L{T}{-1} * tmpa

     // if v is feasible ( v >= 0 ) then it is optimal for the current
     // subproblem: return to check for global optimality

     LocalOptima = TRUE;
     break;

     }  // end( case( 0 ) ) - - - - - - - - - - - - - - - - - - - - - - - - -

    case( 1 ):  // Case 1: Q{B} is *singular* - - - - - - - - - - - - - - - -
                //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // being w the multipliers such that G{B} * w = G[ h ], check
                // e( h ) - e( S ) * w =
                // e( h ) - e( S ){T} * L{B}{T}{-1} * L{B}{-1} * Q[ h ] =
                //                     = e( h ) - L[ BDim - 1 ] * z1
    {
     // if z1 == 0, Case 1.0 can *never* happen

     HpNum ll = ( BDim == CBDim ? 0 :
                  ScalarProduct( L[ BDim - 1 ] , z1 , BDim - 1 ) );

     // note that, for the projected hessian being singular it must be
     // ll = e( Base[ last ] ): this means that size( ll ) ~= 1

     if( ABS( e( Base[ BDim - 1 ] , 1 ) - ll ) <= stdsize() )
     {
      // Case 1.0 : the projected hessian also is singular- - - - - - - - - -
      // find w that solves L{B}{T} * w = L[ BDim - 1 ]: [ w , -1 ] is an
      // infinite direction

      SolveSystemT( L , L[ BDim - 1 ] , tmpv , BDim - 1 );
      tmpv[ BDim - 1 ] = - 1;

      #if( LOG_MQ > 4 )
       *MQLog << " ~ i.d.d. (1)";
      #endif
      }
     else  // Case 1.1 : the projected hessian is nonsingular - - - - - - - -
     {
      LhTz1 = ll;
      LhTz2 = ScalarProduct( L[ BDim - 1 ] , z2 , BDim - 1 );

      LastRo = ( ( Alfa[ Base[ BDim - 1 ] ] - LhTz2 ) / ti )
                 / ( e( Base[ BDim - 1 ] , 1 ) - LhTz1 );

      #if( EXACT )
       z1Tz2 = ScalarProduct( z1 , z2 , BDim - 1 );      
       z1Norm = Norm( z1 , BDim - 1 );
      #endif

      tmpv[ BDim - 1 ] = ( 1 + z1Tz2 / ti - LastRo * z1Norm )
                                  / ( e( Base[ BDim - 1 ] , 1 ) - LhTz1 );

      // tmpa = ro * z1 - z2 / ti
      VectAdd( tmpa , z1 , z2 , LastRo , - 1 / ti , BDim - 1 );

      // tmpa += - tmpv[ BDim - 1 ] * L[ BDim - 1 ]
      VectSum( tmpa , L[ BDim - 1 ] , - tmpv[ BDim - 1 ] , BDim - 1 );

      SolveSystemT( L , tmpa , tmpv , BDim - 1 );

      LocalOptima = TRUE;  // as above, if v is feasible ( v >= 0 ) ...

      }  // end else( e( BDim - 1 ) - e( S ) * v == 0 ) - - - - - - - - - - -

     break;

     }  // end( case( 1 ) ) - - - - - - - - - - - - - - - - - - - - - - - - -

    case( 2 ):  // Case 2: 2 linearly dependent items in Base - - - - - - - - 
                //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // being v1 e v2 vectors such that
                // G{B} * v1 = G[ Base[ BDim - 2 ] ] = G[ h ] and
                // G{B} * v2 = G[ Base[ BDim - 1 ] ] = G[ j ]
                // an infinite direction for the current point can always be
                //  obtained; v1 and v2 are obtained by solving
                // L{B}{T} * v1 = L[ BDim - 2 ] = Lh
                // L{B}{T} * v2 = L[ BDim - 1 ] = Lj
    {
     HpNum ll = e( Base[ BDim - 1 ] , 1 );

     if( BDim > CBDim )  // if z1 == 0, Case 2.0 can *never* happen
      ll -= ScalarProduct( z1 , L[ BDim - 1 ] , BDim - 2 );

     if( ABS( ll ) <= stdsize() )  // size( ll ) ~= 1 as in Case 2
     {
      // Case 2.0 : [ v2 0 -1 ] is an infinite direction- - - - - - - - - - -
      // if e( S ) * v2 = e( j ), then [ v2 0 -1 ] is an infinite direction

      SolveSystemT( L , L[ BDim - 1 ] , tmpv , BDim - 2 );
      tmpv[ BDim - 2 ] =  0;
      tmpv[ BDim - 1 ] = -1;

      #if( LOG_MQ > 4 )
       *MQLog << " ~ i.d.d. (2)";
      #endif
      }
     else
     {
      // Case 2.1 : [ (v1 - ro * v2) - 1 ro ] is an infinite direction- - - -
      // the direction w = [ v1 - ro * v2 , - 1 , ro ], where
      // ro = ( e( h ) - e( S ) * v1 ) / ( e( j ) - e( S ) * v2 ) is an
      //  infinite direction, and v = v2, ll = e( j ) - e( S ) * v2 != 0 are
      // already at hand: for the actual calculation, exploit the fact that
      // v1 - ro * v2 = ( L{B}{T}{-1} * Lh ) - ro * ( L{B}{T}{-1} * Lj ) =
      // L{B}{T}{-1} ( Lh - ro * Lj )

      HpNum ro = ( e( Base[ BDim - 2 ] , 1 )
                   - ScalarProduct( z1 , L[ BDim - 2 ] , BDim - 2 ) ) / ll;

      // tmpa = Lh - ro * Lj
      VectAdd( tmpa , L[ BDim - 2 ] , L[ BDim - 1 ] , - ro , BDim - 2 );

      SolveSystemT( L , tmpa , tmpv , BDim - 2 );

      tmpv[ BDim - 2 ] = - 1;
      tmpv[ BDim - 1 ] = ro;

      #if( LOG_MQ > 4 )
       *MQLog << " ~ i.d.d. (3)";
      #endif

      }  // end else( e( S ) * v2 != e( j ) ) - - - - - - - - - - - - - - - -

     LocalOptima = FALSE;
     break;

     }  // end( case( 2 ) ) - - - - - - - - - - - - - - - - - - - - - - - - -

    case( 3 ):  // Corresponding to Dependent == 0 && LBIsInB == 1- - - - - -
                //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // it is a simplified version of Case 1, as the Lower Bound
                // in Base[] is an all-0 => linearly dependent item
    {
     LastRo = LwrBnd / ti;

     #if( EXACT )
      z1Tz2 = ScalarProduct( z1 , z2 , BDim - 1 );      
      z1Norm = Norm( z1 , BDim - 1 );
     #endif

     tmpv[ BDim ] = 1 + z1Tz2 / ti - LastRo * z1Norm;

     VectAdd( tmpa , z1 , z2 , LastRo , - 1 / ti , BDim );

     SolveSystemT( L , tmpa , tmpv , BDim );

     LocalOptima = TRUE;

     break;

     }  // end( case( 3 ) ) - - - - - - - - - - - - - - - - - - - - - - - - -

    case( 4 ):  // Corresponding to Dependent == 1 && LBIsInB == 1- - - - - -
                //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // it is a simplified version of Case 2, as the Lower Bound
                // in Base[] is an all-0 => linearly dependent item added to
                // the linearly dependent item already in Base[]
    {
     SolveSystemT( L , L[ BDim - 1 ] , tmpv , BDim - 1 );

     tmpv[ BDim - 1 ] = - 1;

     tmpv[ BDim ] = e( Base[ BDim - 1 ] , 1 )
                    - ScalarProduct( z1 , L[ BDim - 1 ] , BDim - 1 );

     #if( LOG_MQ > 4 )
      *MQLog << " ~ i.d.d. (4)";
     #endif

     LocalOptima = FALSE;

     }  // end( case( 4 ) ) - - - - - - - - - - - - - - - - - - - - - - - - -

    }  // end switch( Dependent ) - - - - - - - - - - - - - - - - - - - - - -
       // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /*- - - - - - - - - - - Movement to the new point - - - - - - - - - - - -*/
   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /* Here LocalOptima == TRUE if tmpv is optimal for the restricted problem
      *withouth the >= 0 constraints*: otherwise tmpv is an "infinite
      direction". */

   register cIndex RBDim = BDim + LBIsInB;  // the "real" BDim

   if( LocalOptima )       // tmpv contains a point- - - - - - - - - - - - -
    if( Feasible() )       // if it is also feasible
     Swap( Mult , tmpv );  // Mult = tmpv
    else                   // if it is unfeasible
    {
     VectSubtract( tmpv , Mult , RBDim );  // tmpv = tmpv - Mult
     LocalOptima = FALSE;                  // go "towards" v
     }
   else                    // tmpv contains a direction - - - - - - - - - - -
    if( ( Ename < InINF ) && ( Mult[ Ewher ] == 0 ) )  // Mult[] has the form
    {                                                  // [ x, 0 ]
     if( tmpv[ Ewher ] < 0 )
      for( register HpRow tv = tmpv + RBDim ; tv-- > tmpv ; )
       *tv = - (*tv);
     }
    else                         // there is no i such that Mult[ i ] == 0
     if( vPerAlfa( tmpv ) > 0 )  // check tmpv * Alfa[] to decide the sense
      for( register HpRow tv = tmpv + RBDim ; tv-- > tmpv ; )
       *tv = - (*tv);

   if( ! LocalOptima )  // tmpv is a direction, possibly infinite - - - - - -
   {                    //- - - - - - - - - - - - - - - - - - - - - - - - - -
                        // Mult = Mult + Step * v : the max feasible Step is
                        // min{ - Mult[ i ] / v[ i ] : v[ i ] < 0 }
    if( LBIsInB )
     Mult[ BDim ] = LBMult;

    register HpNum eps = - eR * RBDim;
    register HpRow tv = tmpv + RBDim;
    register HpNum Stp = HpINF;

    for( register HpRow tM = Mult + RBDim ; tM-- > Mult ; )
     if( ( *(--tv) < eps ) && ( - ( *tM / *tv ) < Stp ) )
      Stp = - ( *tM / *tv );

    if( Stp < 0 )
    {
     #if( LOG_MQ > 2 )
      *MQLog << endl << "Fault: negative step ( = " << Stp << ", Lmu = "
             << Lmu << "), x[ i ] < 0 [CalcOptDir]";
     #endif

     QPStatus = kNegVarInBase;
     return;
     }

    if( Stp == HpINF )
    {
     if( CBDim )  // check that tmpv really is an "infinite-infinite"
     {            // direction, that is an infinite direction along which an
                  // infinite step is possible): tmpv[ i ] *must be zero* for
                  // each subgradient i, while there *must* be at least a
                  // constraint j s.t. tmpv[ j ] is *strictly positive*
      tv = tmpv;
      eps = - eps;
      Stp = 0;

      register Index_Set tB = Base;
      for( register Index h ; ( h = *(tB++) ) < InINF ; )
       if( G[ h ] & IsACnst )
        Stp += *(tv++);
       else
        if( ABS( *(tv++) ) > eps )
        {
         Stp = 0;
         break;
         }

      if( Stp && LBIsInB && ABS( *tv ) > eps )
       Stp = 0;

      if( Stp > eps )
      {
       VectAssign( Mult , tmpv , BDim );
       f = HpINF;

       #if( LOG_MQ > 3 )
        *MQLog << endl << " ~ Stop: problem is UNBOUNDED." << endl;
       #endif

       QPStatus = kQPPrimUnbndd;
       return;
       }
      }  // end( if( CBDim ) )

     #if( LOG_MQ > 2 )
      *MQLog << endl << "Fault: invalid direction v (Lmu = " << Lmu
             << "). [CalcOptDir]";
     #endif

     QPStatus = kInvalidV;
     return;
     }

    VectSum( Mult , tmpv , Stp , RBDim );  // Mult += Step * v

    #if( LOG_MQ > 3 )
     *MQLog << " ~ move (step = " << Stp << ")";
    #endif

    }  // end if( ! LocalOptima )

   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /*- - - - - - - - - - - - - Items elimination - - - - - - - - - - - - - -*/
   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /* All items with 0 Multiplier are eliminated from the Base: the deletion
      goes from the last item to the first one, because deleting an item in
      position q costs O( BDim - Dependent - q ) (eliminating dependent items
      costs nothing). If LocalOptima == FALSE, *at least one* item *must* be
      eliminated. */

   BOOL Mvd = LocalOptima;

   if( LBIsInB )                          // special check for the multiplier
    if( ( LwrBnd = Mult[ BDim ] ) < eR )  // corresponding to the Lower Bound
    {
     Mvd = TRUE;

     #if( LOG_MQ > 3 )
      *MQLog << " ~ out LB = " << LwrBnd;
     #endif
     
     if( Ename == MaxBDim )
     {
      #if( LOG_MQ > 2 )
       *MQLog << endl
	      << "Fault: LB in & out in the same iteration. [CalcOptDir]";
      #endif

      QPStatus = kLoop;
      return;
      }
     }

   for( register HpRow tM = Mult + BDim ; tM > Mult ; )
    if( *(--tM) <= eR )  // Mult[ i ] is assumed to be <= 1, even though
    {                    // it may not be true for constraints
     Mvd = TRUE;
     Index i = tM - Mult;

     #if( LOG_MQ > 3 )
      *MQLog << " ~ out x[" << Base[ i ] << " = B[" << i << "]] = " << *tM;
     #endif

     // the item that has entered in the latest pricing (Ename) should not
     // exit: if this happens, mark it as "taboo" by putting it in the first
     // TLDim positions of Order[]

     if( Base[ i ] == Ename )
     {
      #if( LOG_MQ > 4 )
       *MQLog << " (was in => blocked)";
      #endif

      Swap( Order[ Eordr ] , Order[ TLDim++ ] );
      Ename = InINF;
      }

     CutSubGradFromBase( i );
     }

   if( ! Mvd )  // LocalOptima == FALSE but no items deleted
   {
    #if( LOG_MQ > 2 )
     *MQLog << endl << "Fault: can't find a local optima => loop (Lmu = "
            << Lmu << ") [CalcOptDir]";
    #endif

    QPStatus = kLoop;
    return;
    }

   }  // end( while( ! LocalOptima ) )

  /*------------------------------------------------------------------------*/
  /*------------------------- End Inner Loop -------------------------------*/
  /*------------------------------------------------------------------------*/

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /*- - - - - - - - Computation of f() and decrease test- - - - - - - - - -*/
  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /* Each time the algorithm visits a restricted optima, f() is calculated to
     test for strict decrease: if it has not been obtained, the item that is
     entered last (if any) is marked as "taboo". That item will no more be
     allowed to enter the Base until a strict decrease is obtained (each time
     it happens, the whole "taboo list" is cleaned). Also, the previous value
     of f() is kept (in Minf), since it must at least be reached - the Base
     that existed when the "bad" item entered was "ok", and it can still be
     reached. Note that there are (possibly long) sequences of steps where
     f() is not computed, and therefore the test is not performed: however,
     a restricted optimal point must be reached within BDim steps, since *at
     least one* item *must* be removed from Base at each step, while no items
     are added. */

  Lin = ScalarProduct( Mult , Alfa , Base );
  if( LBIsInB )
   Lin += LwrBnd * LBMult;

  HpNum Newf;  // the new value of f()

  if( LastRo == - HpINF )
  {
   Quad = - Lin / ti;
   Newf = - Quad / 2;
   }
  else
  {
   Quad = LastRo - Lin / ti;
   Newf = ( LastRo + Lin / ti ) / 2;
   }

  if( ( Ename < MaxBDim ) && ( f < HpINF ) )
  {
   /* Test for *strict*  decrease of f(), i.e. Newf < f - eR * | f |: skip the
      test if Ename >= MaxBDim, i.e. just after a reoptimization with changes
      in in the data (where f() need not necessarly decrease), or if Ename has
      already been put in the taboo list because it has been eliminated, or if
      Ename == MaxBDim, that is the entered item is the Lower Bound.
      The test is also skipped when f == HpINF: this can happen during a
      "phase transition", i.e. when the first subgradient enters in Base. */

   if( Newf >= f - max( eR * ABS( f ) , HpeM ) )
   {
    // the latest entered item promised a decrease but in fact it gave an
    // increase: "mark" it as "taboo" (the taboo list will be cleaned only
    // when a better value of f() is obtained)

    #if( LOG_MQ > 4 )
     *MQLog << " ~ Df = " << Newf - f << " => " << Ename << " blocked";
    #endif

    Swap( Order[ Eordr ] , Order[ TLDim++ ] );
    }
   }

  f = Newf;
  #if( LOG_MQ > 3 )
   *MQLog << " ~ f = " << f;
  #endif

  if( f + max( eR * ABS( f ) , HpeM ) <= Minf )
  {
   TLDim = 0;  // reset the "taboo list"
   Minf = f;
   }

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /*- - - - - - - - - - - - - Optimality Test - - - - - - - - - - - - - - -*/
  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  if( LastRo > - HpINF )  // if Mult[] is feasible
   if( f <= MinFVal )     // and f() is small enough
    break;                // skip the tests

  /* In a *restricted optima*, the dual constraints in Base are satisfied to
     equality (this is *not* true in general): therefore, they don't need to
     be checked. For all other items
                                                   / v   i subgradient
      ro[ i ] = Alfa[ i ] + G[ i ] * G{B} * Mult + |
                                                   \ 0   i constraint

     is the derivative in Beta = 0 of f( Mult + Beta * [ 0 .. 1 .. 0 ] ), the
     "1" being in position i: hence, a decrease should be obtained by adding
     item 'i' to the Base (the i-th dual constraint is violated).

     The test for decrease assumes that Mult[ i ] will be <= 1 / BDim: hence,
     the maximum decrease in f() is estimated as about ro[ i ] * ( 1 / BDim ),
     so that the item is inserted when ro[ i ] * ( 1 / BDim ) < - eR * | f |
     => ro[ i ] < - eR * BDim * | f |. */

  Ename = InINF;
  HpNum minro = - max( eR * BDim * ABS( f ) , HpeM );
  cHpNum brkro = - ( EpsRo == HpINF ? EpsRo : EpsRo * BDim * ABS( f ) );

  register Index_Set tO = Order + ActBDim;
  for( register Index h ; tO > Order + TLDim; )   // pricing loop - - - - - -
   if( ! ( G[ h = *(--tO) ] & IsInBse ) )         // a nonbasic item- - - - -
   {                                              //- - - - - - - - - - - - -
    register HpNum ro = 0;
    register HpRow tM = Mult;
    register Index_Set bj = Base;
    for( register Index j ; ( j = *(bj++) ) < InINF ; )
     ro += *(tM++) * LowQ( h , j );

    GTd[ h ] = ro;
    ro += Alfa[ h ] / ti;

    if( ! ( G[ h ] & IsACnst ) )
     if( LastRo > - HpINF )
      ro = ro - e0( h , LastRo );
     else
     {
      ro = - HpINF;      // the first subgradient entering in Base
      Minf = f = HpINF;  // switch to the "true" f() (hence reset Minf)
      }

    if( ro < minro )
    {
     Ename = h;           // the entering item
     Eordr = tO - Order;  // its position in Order

     if( ( minro = ro ) <= brkro )
      break;
     }
    }  // end if( nonbasic item ) - - - - - - - - - - - - - - - - - - - - - -
       // end for( pricing loop ) - - - - - - - - - - - - - - - - - - - - - -
       // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // a specialized pricing for the Lower Bound- - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( Ename < InINF ) && ( LwrBnd < HpINF ) && ( ! LBIsInB ) )
  {
   register HpNum ro = LwrBnd / ti;

   if( LastRo > - HpINF )
    ro = ro - LastRo;
   else
    ro = - ( Minf = f = HpINF );

   if( ro < minro )
   {
    LocalOptima = FALSE;

    #if( LOG_MQ > 3 )
     *MQLog << " ~ in LB, ro = " << ro;
    #endif

    LBMult = 0;
    Ewher = BDim;
    Ename = MaxBDim;

    continue;  // go directly to the beginning of the next outer iteration
    }
   }  // end( special pricing for the Lower Bound ) - - - - - - - - - - - - -

  if( Ename < InINF )  // a dual constraint is violated - - - - - - - - - - -
  {                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
   LocalOptima = FALSE;

   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /*- - - - - - - - - - Inserting the violated constraint - - - - - - - -*/
   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /* Reminders on data structures organization: the lower trapezoidal L{B}
      is the factorization of the submatrix Q{B} made of the rows and columns
      of Q whose indices are in B. Dependent <= 2 tells how many rows of L{B}
      have a 0 diagonal element (they will always be the last ones): that is,
      B = B' + B", where | B" | = Dependent (may be 0) and L{B'} is a lower
      triangular matrix with all nonzero diagonal elements. The two vectors
      z1 and z1 solve L{B'} * zj = v{B'}, with v respectively the vector of
      coefficients of the (only) constraint and Alfa. The zj's are updated
      while L{B} is, i.e.

      - when a G[ i ] is added, but it is linearly dependent (i goes to B")
        the zj's don't change;

      - when a linearly independent G[ i ] is added, i.e. i goes to B', and
          _________
         | L{B'} 0 |   is the new L{B'}
         |__l____d_|

        then [ zj , ( v[ i ] - l{T} * zj ) / d ] is the new zj;

      - when a G[ i ] (in B' at position q) is eliminated, the same k - q - 1
        Givens matrices that must be applied to L{B} to restore triangularity
        can be applied to the zj's to obtain the new ones.

      Hence, the Base is kept in "standard" form, i.e.

      - L{B} = | L{B'} |  <- full rank, lower triangular factor
               |  Lh   |  \  rows corresponding respectively to
               |  Lj   |  |  Gh = G[ h ] = G[ Base[ Dim - 2 ] ] and
                          |  Gj = G[ j ] = G[ Base[ Dim - 1 ] ], i.e.
                          |  L{B'} * Lh = G{B'}{T} * G[ h ] and
                          |  L{B'} * Lj = G{B'}{T} * G[ j ] */
   #if( LOG_MQ > 3 )
    *MQLog << " ~ in " << Ename << ", ro = " << minro;
   #endif

   Mult[ BDim ] = 0;
   Base[ Ewher = BDim ] = Ename;

   // the system L{B} * L[ BDim ] = tmpa = G{B} * G[ Ename ] is solved to
   // build the new row of L

   GjTGBk( Ename , Base , BDim - Dependent , tmpa );
   SolveSystem( L , tmpa , L[ Ewher ] , BDim - Dependent );

   HpNum Delta = CalculateDelta( Ename );

   // Delta = 0 <=> G[ Ename ] is linearly dependent from the first
   // BDim - Dependent items in base

   Base[ ++BDim ] = InINF;
   G[ Ename ] |= IsInBse;
   if( G[ Ename ] & IsACnst )
    CBDim++;

   if( Delta < 0 )
   {
    Dependent++;  // the item is "at least" dependent

    #if( LOG_MQ > 2 )
     *MQLog << endl << "Fault: Delta = " << Delta << " < 0, (Lmu = " << Lmu
            << ") [CalcOptDir]";
    #endif

    QPStatus = kNegativeDelta;
    return;
    }

   if( Delta )  // the new item is linearly independent from those in B'
    AddSubGradToBase( BDim - 1 , Delta );
   else
    Dependent++;
   }
  else  // no dual constraint is violated - - - - - - - - - - - - - - - - - -
        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( f - max( eR * ABS( f ) , HpeM ) > Minf )
   {
    // if no more items can be "priced in" (some because they are "taboo"),
    // but f() is not (about) as good as Minf, then kLoop is returned

    #if( LOG_MQ > 2 )
     *MQLog << endl << "Fault: f() can't reach the min. value " << Minf
            << " (Lmu = " << Lmu << ") [CalcOptDir]";
    #endif

    QPStatus = kLoop;
    return;
    }

  } while( ! LocalOptima );

 /*-------------------------------------------------------------------------*/
 /*-------------------------- End Main Loop --------------------------------*/
 /*-------------------------------------------------------------------------*/

 // calculate the entries of GTd[] relative to "basic" items

 register Index_Set tB = Base;
 if( CBDim )
  for( register Index i ; ( i = *(tB++) ) < InINF ; )
   GTd[ i ] = e( i , LastRo ) - Alfa[ i ] / ti;
 else
  for( register Index i ; ( i = *(tB++) ) < InINF ; )
   GTd[ i ] = LastRo - Alfa[ i ] / ti;

 // calculate the entries of GTd[] relative to "taboo" items

 for( tB = Order + TLDim ; tB > Order ; )
 {
  register Index h = *(--tB);
  register HpNum ro = 0;
  register HpRow tM = Mult;
  register Index_Set bj = Base;
  for( register Index j ; ( j = *(bj++) ) < InINF ; )
   ro += *(tM++) * LowQ( h , j );

  GTd[ h ] = ro;
  }

 TLDim = 0;  // clean up taboo list

 #if( LOG_MQ > 3 )
  *MQLog << " ~ Stop: OPTIMAL solution." << endl;
 #endif

 #if( LOG_MQ > 1 )
  Success++;
  SumActBDim += ActBDim;
  if( ActBDim )
   SumGSBundle += float( ActCNum ) / float( ActBDim );
 #endif

 }  // end( CalcOptDir )

/*--------------------------------------------------------------------------*/

inline BOOL MinQuad::Feasible( void )
{
 register Index n = BDim + LBIsInB;
 register cHpNum e = eR * n;

 for( register cHpRow v = tmpv ; n && ( *(v++) >= -e ) ; )
  n--;

 return( n == 0 );  // TRUE <=> tmpv[ i ] >= - e for each i
 }

/*--------------------------------------------------------------------------*/

inline HpNum MinQuad::vPerAlfa( register cHpRow v )
{
 HpNum res = ScalarProduct( v , Alfa , Base );
 if( LBIsInB )
  res += LwrBnd * v[ BDim ];

 return( res );
 }

/*--------------------------------------------------------------------------*/

inline HpNum MinQuad::CalculateDelta( cIndex i )

/* Delta = LowQ( i , i ) - Norm( L[ BDim ] ), plus the check if it is zero:
   being e = eR * Lmu * Norm( ... ), it returns

   - sqrt( Delta ) if Delta > e         (OK, nonzero)
   - 0             if -e <= Delta <= e  (OK, zero)
   - Delta         if Delta < -e < 0,   (problems: negative Delta). */
{
 HpNum nrm = Norm( L[ BDim ] , BDim - Dependent );
 HpNum Delta = LowQ( i ) - nrm;
 HpNum eps = gsize( nrm );

 if( Delta > eps )
  Delta = sqrt( Delta );
 else
  if( Delta >= - eps )
  {
   #if( ! EXACT )
    if( Dependent )
     LjNorm = nrm;
    else
     LhNorm = nrm;
   #endif
   
   Delta = 0;
   }

 return( Delta );

 }  // end( CalculateDelta )

/*--------------------------------------------------------------------------*/

inline void MinQuad::AddSubGradToBase( Index BD , cHpNum Delta )

/* Insert a new linearly independent item in Base. The item is put in row BD
   (either BDim - 1 or BDim - 2), so if Dependent == 1 and BD == BDim - 1 a
   row swap is performed first (Ewher is updated). L[ BD ] must contain the
   correct first BD - Dependent entries of the row. */
{
 #if( LOG_MQ > 1 )
  Insert++;
 #endif

 if( Dependent && ( BD == BDim - 1 ) )  // swap the last two rows of L
 {
  #if( VASTE_MEM < 2 )
   VectXcg( L[ BD ] , L[ BD - 1 ] , BD - 1 );
  #else
   Swap( L[ BD ] , L[ BD - 1 ] );
  #endif

  #if( LOG_MQ > 4 )
   *MQLog << " ~ s.l.2";
  #endif

  L[ BD ][ BD - 1 ] = ( LowQ( Base[ BD ] , Base[ BD - 1 ] ) -
                        ScalarProduct( L[ BD ] , L[ BD - 1 ] , BD - 1 ) )
                      / Delta;

  // it is known that lh * lj == G_h * G_j , i.e. that a zero must be placed
  // there (note that it implies that LhNorm don't change): however, this may
  // lead to imprecision, so the actual calculation is better

  Swap( Base[ BD ] , Base[ BD - 1 ] );
  Swap( Mult[ BD ] , Mult[ BD - 1 ] );

  if( Ewher == BD-- )  // this is happening to the newly entered item
   Ewher--;
  }

 L[ BD ][ BD ] = Delta;

 z2[ BD ] = ( Alfa[ Base[ BD ] ] - ScalarProduct( z2 , L[ BD ] , BD ) )
              / Delta;

 z1[ BD ] = ( e( Base[ BD ] , 1 ) - ScalarProduct( z1 , L[ BD ] , BD ) )
              / Delta;

 #if( ! EXACT )
  z1Norm += z1[ BD ] * z1[ BD ];
  z1Tz2 += z1[ BD ] * z2[ BD ];
 #endif

 if( Delta < L[ LiMin ][ LiMin ] )
 {
  LiMin = BD;
  Lmu = L[ LiMax ][ LiMax ] / Delta;
  }

 if( Delta > L[ LiMax ][ LiMax ] )
 {
  LiMax = BD;
  Lmu = Delta / L[ LiMin ][ LiMin ];
  }
 }  // end( AddSubGradToBase )

/*--------------------------------------------------------------------------*/

inline void MinQuad::ChkNwLinInd( void )

/* Is analogous to AddSubGradToBase() above, but is is called from
   AddSGSpaceDim() rather than from CalcOptDir() or RebuildBase().
   Uses the information contained in the temporary tmpv[]. */
{
 if( Dependent )
 {
  Index BD = BDim - Dependent;
  HpNum Delta = L[ BD ][ BD ] = tmpv[ BD ];

  #if( ! EXACT )
   LhNorm = LowQ( Base[ BD ] );
  #endif

  if( Dependent > 1 )
  {
   L[ BDim - 1 ][ BD ] = tmpv[ BDim - 1 ];

   #if( ! EXACT )
    LjNorm = LowQ( Base[ BDim - 1 ] );
   #endif
   }

  if( ( Dependent == 2 ) &&
      ( Delta * Delta <= gsize( LowQ( Base[ BD ] ) ) ) )
  {
   // there are 2 lin. dep. items, and the first is *not* going to enter

   Delta = tmpv[ BDim - 1 ];

   if( Delta * Delta > gsize( LowQ( Base[ BDim - 1 ] ) ) )
   {
    // ... but the second *is* going to enter, so swap them

    #if( VASTE_MEM < 2 )
     VectXcg( L[ BD ] , L[ BDim - 1 ] , BD );
    #else
     Swap( L[ BD ] , L[ BDim - 1 ] );
    #endif

    Swap( Base[ BD ] , Base[ BDim - 1 ] );
    Swap( Mult[ BD ] , Mult[ BDim - 1 ] );

    #if( ! EXACT )
     Swap( LhNorm , LjNorm );
    #endif
    }
   }

  if( Delta * Delta > gsize( LowQ( Base[ BD ] ) ) )
  {
   Dependent--;

   if( Delta < 0 )
   {
    L[ BD ][ BD ] = Delta = - Delta;

    if( Dependent )
     L[ BDim - 1 ][ BD ] = - L[ BDim - 1 ][ BD ];
    }

   z2[ BD ] = ( Alfa[ Base[ BD ] ] - ScalarProduct( z2 , L[ BD ] , BD ) )
                / Delta;

   z1[ BD ] = ( e( Base[ BD ] , 1 ) - ScalarProduct( z1 , L[ BD ] , BD ) )
                / Delta;

   #if( ! EXACT )
    z1Norm += z1[ BD ] * z1[ BD ];
    z1Tz2 += z1[ BD ] * z2[ BD ];

    if( Dependent )
     LhNorm = LjNorm + L[ BDim - 1 ][ BD ] * L[ BDim - 1 ][ BD ];
   #endif
   }
  }

 LiMax = LiMin = 0;
 UpdateLmu( 1 );

 }  // end( ChkNwLinInd )

/*--------------------------------------------------------------------------*/

void MinQuad::CutSubGradFromBase( cIndex s )

/* Eliminates the item `i', in Base at position s (i == Base[ s ]): updates L,
   Base, G, BDim, Dependent, Ewher, Mult, z1, z2 and Lmu. */
{
 G[ Base[ s ] ] &= ~IsInBse;
 if( G[ Base[ s ] ] & IsACnst )
  CBDim--;

 if( s >= BDim - Dependent )  // eliminating a dependent row- - - - - - - - -
 {
  Dependent--;
  BDim--;

  if( s < BDim )  // eliminate the (BDim - 2)-th but keep the (BDim - 1)-th
  {
   #if( LOG_MQ > 4 )
    *MQLog << " ~ elim. Lh, keep Lj";
   #endif
   
   #if( VASTE_MEM == 0 )
    VectAssign( L[ s ] , L[ BDim ] , s );
   #elif( VASTE_MEM == 1 )
    delete[] L[ s ];
    L[ s ] = L[ BDim ];
    L[ BDim ] = new HpNum[ BDim + 1 ];
   #else
    Swap( L[ s ] , L[ BDim ] );
   #endif

   #if( ! EXACT )
    LhNorm = LjNorm;
   #endif

   Base[ s ] = Base[ BDim ];
   Mult[ s ] = Mult[ BDim ];

   if( Ewher == BDim )  // the last item in Base it the newly entered one
    Ewher--;
   }
  #if( LOG_MQ > 4 )
  else
   *MQLog << " ~ elim. L" << ( Dependent ? "j" : "h" );
  #endif

  Base[ BDim ] = InINF;
  return;
  }

 // update the factorization- - - - - - - - - - - - - - - - - - - - - - - - -

 #if( LOG_MQ > 1 )
  Delete++;
 #endif

 /* Currently, L{B} is the lower trapezoidal matrix
                   ___________
                / | L1  0  0  | s - 1
    BDim        | | l   d  0  | 1 (in position s)
                \ | Z   v  L2 | \ BDim - s - 1
    Dependent     |_X___y__K2_| /

    This procedure finds and apply to L{B} an orthogonal matrix G, made of
    BDim - Dependent - s - 1 Givens matrices, such that
       ___________     _____     ___________
      | L1  0  0  |   | I 0 |   | L1  0  0  |
      | l   d  0  | * |_0_G_| = | l   d  0  |   where
      |_Z___v__L2_|             |_Z___0__L2'|
       ________
      | L1  0  |  is the factorization of Q{B'}, B' = B / { s }.
      |_Z___L2'|

    The same Givens matrices are also applied to z1 and z2 to update them. */

 register Index j = s + 1;
 for( ; j < BDim - Dependent ; j++ )
 {
  // a Givens matrix Gsj is found that combines the s-th column (s fixed) with
  // the j-th one (j variable) in such a way that the j-th element of the s-th
  // column of L becomes zero. For any j, the s-th column will only have
  // BDim - Dependent - j nonzeroes, all and only in rows >= j, since all
  // those in rows < j have already been put to 0

  HpNum ai = L[ j ][ s ];

  if( ABS( ai ) <= eR * BDim )  // if it's already 0, nothing to do
   continue;

  HpNum aj = L[ j ][ j ];
  HpNum d = L[ j ][ j ] = sqrt( ai * ai + aj * aj );
  L[ j ][ s ] = 0;

  register HpNum Gc = aj / d;
  register HpNum Gs = - ai / d;

  for( register Index h = j ; ++h < BDim ; )
   ApplyGivensMToCol( Gc , Gs , L[ h ][ s ] , L[ h ][ j ] );

  ApplyGivensMToCol( Gc , Gs , z1[ s ] , z1[ j ] );
  ApplyGivensMToCol( Gc , Gs , z2[ s ] , z2[ j ] );

  }  // end for( j )

 /* If Dependent > 0 those numbers are useful, but they are deleted in
    the next phase, where all the rows under the s-th (even the linearly
    dependent ones, if any) are shifted up and all the columns right to
    the s-th are shifted left, as depicted in the figure below
     _________         _______
    | L1 0 0  |       | L1 0  |
    | l  d 0  |  =>   |_Z__L2_|
    |_Z__v_L2_|                */

 HpNum DeltaLh;
 HpNum DeltaLj;

 if( Dependent )
 {
  DeltaLh = L[ BDim - Dependent ][ s ];  // h = BDim - Dependent

  if( Dependent > 1 )                    // j = BDim - 1, h = BDim - 2
   DeltaLj = L[ BDim - 1 ][ s ];
  }

 BDim--;

 register HpMat lpp = L + s;

 #if( VASTE_MEM )
  #if( VASTE_MEM == 1 )
   delete[] L[ s ];
  #else
   HpRow temp1 = L[ s ]; 
  #endif

  ShiftVect( lpp , BDim - s );  // ShiftVect destroys the first position
                                // it 'd be BDim - s - 1, but we did BDim--
  #if( VASTE_MEM == 1 )
   L[ BDim ] = new HpNum[ BDim + 1 ];  // allocate new memory
  #else
   L[ BDim ] = temp1;                  // restore the pointer
  #endif

  for( j = 0 ; j < BDim - s ; lpp++ )
   ShiftVect( (*lpp) + s , ++j );
 #else
  for( j = 0 ; j < BDim - s ; )  // cannot exchange the pointers
  {                              // => copy element-vise
   register HpRow temp1 = *lpp;
   register HpRow temp2 = *(++lpp);
   register Index h = s;

   for( ; h-- ; )
    *(temp1++) = *(temp2++);  // copy the first s columns

   temp2++;                   // column s + 1 on the second row

   for( h = ++j ; h-- ; )
    *(temp1++) = *(temp2++);  // copy the remaining j + 1 columns
   }
 #endif

 ShiftVect( Mult + s , BDim - s );
 ShiftVect( Base + s , BDim - s );

 if( Ewher > s )  // the newly entered item is "moved"
  Ewher--;

 #if( ! EXACT )
  z1Norm -= z1[ s ] * z1[ s ];
  z1Tz2 -= z1[ s ] * z2[ s ];  
 #endif
 ShiftVect( z1 + s , BDim - s );
 ShiftVect( z2 + s , BDim - s );

 j = s;  // update Lmu

 if( LiMax >= s )
 {
  LiMax = 0;
  j = 1;
  }

 if( LiMin >= s )
 {
  LiMin = 0;
  j = 1;
  }

 UpdateLmu( j );

 if( Dependent )  // control if, after the elimination, some previously
 {                // dependent item have become independent
  j = BDim - 1;
  register Index h = BDim - Dependent;

  #if( EXACT )
   DeltaLh = ABS( HpNum( LowQ( Base[ h ] ) ) - Norm( L[ h ] , h ) );
  #else
   HpNum Signum = sgn( DeltaLh ) * sgn( DeltaLj );

   DeltaLh = ABS( HpNum( LowQ( Base[ h ] ) - LhNorm + DeltaLh * DeltaLh ) );

   // this apparently funny way of recalculating DeltaLh and DeltaLj makes
   // the LL{T} factorization much more exact than the "naive" one
  #endif

  if( DeltaLh > gsize( LowQ( Base[ h ] ) ) )
  {
   #if( LOG_MQ > 4 )
    *MQLog << " ~ Lh indep. (" << Dependent << ")";
   #endif

   Dependent--;
   AddSubGradToBase( h , DeltaLh = sqrt( DeltaLh ) );

   if( Dependent )
   {
    #if( ! EXACT )
     LjNorm -= DeltaLj * DeltaLj;
     LhNorm = LowQ( Base[ j ] );
    #endif

    #if( ! EXACT )
     DeltaLj = ABS( HpNum( LowQ( Base[ j ] ) ) - LjNorm );

     if( Signum < 0 )
      DeltaLj = - sqrt( DeltaLj );
     else
      DeltaLj = sqrt( DeltaLj );
    #else
     DeltaLj = ( LowQ( Base[ h ] , Base[ j ] ) -
                 ScalarProduct( L[ j ] , L[ h ] , h ) ) / DeltaLh;    
    #endif

    L[ j ][ h ] = DeltaLj;
    }
   }  // end if( DeltaLh > 0 )

  if( Dependent > 1 )  // ... still
  {
   #if( EXACT )
    DeltaLj = ABS( HpNum( LowQ( Base[ j ] ) - Norm( L[ j ] , h ) ) );
   #else
    DeltaLj = ABS( HpNum( LowQ( Base[ j ] ) - LjNorm + DeltaLj * DeltaLj ) );
    // note that LhNorm doesn't change
   #endif

   if( DeltaLj > gsize( LowQ( Base[ j ] ) ) )
   {
    #if( LOG_MQ > 4 )
     *MQLog << " ~ Lj indep.";
    #endif

    Dependent--;
    AddSubGradToBase( j , DeltaLj = sqrt( DeltaLj ) );

    #if( EXACT )
     DeltaLh = ( LowQ( Base[ h ] , Base[ j ] ) -
                 ScalarProduct( L[ j ] , L[ h ] , h ) ) / DeltaLj;    
    #else
     if( Signum < 0 )
      DeltaLh = - sqrt( DeltaLh );
     else
      DeltaLh = sqrt( DeltaLh );
    #endif

    L[ j ][ h ] = DeltaLh;

    }  // end if( DeltaLj > 0 )
   }  // end if( Dependent still > 1 )
  }  // end if( Dependent )

 Base[ BDim ] = InINF;

 }  // end( CutSubGradFromBase )

/*--------------------------------------------------------------------------*/

inline void MinQuad::MoveSubgradToLastPos( cIndex s )

/* Moves the item named i, in base at position s (Base[ s ] = i), to the
   last (Independent) position of Base (BDim - Dependent - 1), i.e.
                      _____________
                   / | L1  0  0  0 | s - 1
           BDim    | | l   Ls 0  0 | 1 (in position s)
                   | | Z   v  L2 0 | \ BDim - s - 1
      Dependent    \ |_x___y__k2_0_| /
                    _______________                       _______________
                   | L1  0   0   0 |                     | L1  0   0   0 |
           =>      | l   Ls' d   0 |           =>        | Z   L2' 0   0 |
                   | Z   0   L2' 0 |                     | l   d   Ls' 0 |
                   |_x___y'__k2'_0_|                     |_x___k2'_y'__0_|

   It is assumed that Dependent <= 1 and that s < BDim - 1 - Dependent (i.e.
   s is not a dependent row). Updates L, Base, Mult, z1 and z2; does *not*
   update Lmu and Dependent doesn't change. */
{
 #if( LOG_MQ > 1 )
  Delete++;  // although the operation is logically equivalent to Del + Ins,
             // its complexity is identical to a Del
 #endif

 #if( VASTE_MEM < 2 )
  #if( ! VASTE_MEM )
   HpRow a = tmpr;
  #else
   HpRow a = new HpNum[ BDim ];
  #endif

  VectAssign( a , L[ s ] , s + 1 );
 #else
  HpRow a = L[ s ];
 #endif

 register Index i = s + 1;
 for( ; i < BDim - Dependent ; )
  a[ i++ ] = 0;

 // first phase: column operations to create the zeroes - - - - - - - - - - -

 register Index j = s + 1;
 for( ; j < BDim - Dependent ; j++ )
 {
  HpNum ai = L[ j ][ s ];

  if( ABS( ai ) <= eR * BDim )        // if it's already 0, nothing to do
   continue;

  HpNum aj = L[ j ][ j ];
  HpNum d = L[ j ][ j ] = sqrt( ai * ai + aj * aj );

  register HpNum Gc = aj / d;    // Gc * ai + Gs * aj = 0
  register HpNum Gs = - ai / d;  // Gc ^ 2 + Gs ^ 2 = 1

  for( i = j ; ++i < BDim ; )
   ApplyGivensMToCol( Gc , Gs , L[ i ][ s ] , L[ i ][ j ] );

  ApplyGivensMToCol( Gc , Gs , z1[ s ] , z1[ j ] );
  ApplyGivensMToCol( Gc , Gs , z2[ s ] , z2[ j ] );

  ApplyGivensMToCol( Gc , Gs , a[ s ] , a[ j ] );

  }  // end for( j )

 // second phase: row operations to restore triangularity - - - - - - - - - -

 register HpMat lpp = L + s;

 RotateVect( a + s , BDim - s - 1 - Dependent );

 #if( VASTE_MEM )
  #if( VASTE_MEM == 1 )
   delete[] L[ s ];
  #endif

  ShiftVect( lpp , BDim - s - 1 - Dependent );

  for( j = 0 ; j < BDim - s - 1 - Dependent ; lpp++ )
   ShiftVect( (*lpp) + s , ++j );

  L[ BDim - 1 - Dependent ] = a;
 #else
  // cannot exchange the pointers => copy element-vise

  for( j = 0 ; j < BDim - s - 1 - Dependent ; )
  {
   register HpRow temp1 = *(lpp++);
   register HpRow temp2 = *lpp;

   for( i = s ; i-- ; )
    *(temp1++) = *(temp2++);  // copy the first s columns

   temp2++;                   // column s + 1 on the second row

   for( i = ++j ; i-- ; )
    *(temp1++) = *(temp2++);  // copy the remaining j + 1 columns
   }

  VectAssign( L[ BDim - 1 - Dependent ] , a , BDim - Dependent );
 #endif

 RotateVect( Mult + s , BDim - s - 1 - Dependent );
 RotateVect( Base + s , BDim - s - 1 - Dependent );
 RotateVect( z1 + s , BDim - s - 1 - Dependent );
 RotateVect( z2 + s , BDim - s - 1 - Dependent );

 if( Dependent )
  RotateVect( L[ BDim - 1 ] + s , BDim - s - Dependent );

 // z1Norm, z1Tz2, LhNorm etc. don't change

 }  // end( MoveSubgradToLastPos )

/*--------------------------------------------------------------------------*/

inline void MinQuad::UpdateLmu( cIndex s )

/* Updates the condition number of L{B'}, i.e. the linearly independent part
   of the LL{T} factorization of Q{B}, "under" row s. */
{
 if( BDim <= Dependent )  // L is empty ...
  Lmu = 1;                // ... by def.
 else
 {
  register HpNum Max = L[ LiMax ][ LiMax ];
  register HpNum min = L[ LiMin ][ LiMin ];
  register HpMat pL = L + s;

  for( register Index i = s ; i < BDim - Dependent ; i++ )
  {
   register HpNum Lii = (*(pL++))[ i ];

   if( Max < Lii )
   {
    LiMax = i;
    Max = Lii;
    }
   else
    if( min > Lii )
    {
     LiMin = i;
     min = Lii;
     }
   }

  Lmu = Max / min;
  }
 }  // end( UpdateLmu )

/*--------------------------------------------------------------------------*/

inline void MinQuad::MemDealloc( void )
{
 // memory deallocation - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( register Index i = min( SgSpDim , CrrBDim ) ; i-- ; )
  delete[] L[ i ];

 #if( ! LAZY_Q )
  delete[] tempQ;
 #endif

 for( register Index i = CrrBDim ; i-- ; )
  delete[] Q[ i ];

 #if( ! VASTE_MEM )
  delete[] tmpr;
 #endif

 delete[] tmpa;

 #if( VARCOEFF )
  delete[] cf;
 #endif

 delete[] z1;
 delete[] z2;
 delete[] Base;
 delete[] tmpv;
 delete[] Mult;
 delete[] L;
 delete[] Q;
 delete[] Order;
 delete[] GTd;
 delete[] Alfa;
 delete[] G;

 }  // end( MemDealloc )

/*--------------------------------------------------------------------------*/

#if( LOG_MQ > 5 )

#define CMP( x1 , x2 ) ABS( x1 ) > eR * Lmu * max( ABS( x2 ) , HpNum( 1 ) )

void MinQuad::CheckBaseValidity( void )
{
 if( Dependent )
 {
  L[ BDim - 1 ][ BDim - 1 ] = 0;

  if( Dependent > 1 )
  {
   L[ BDim - 2 ][ BDim - 2 ] = 0;
   L[ BDim - 1 ][ BDim - 2 ] = 0;
   }
  }

 register Index i = 0;

 for( ; i < BDim ; i++ )
  for( register Index j = 0 ; j <= i ; j++ )
  {
   HpNum xy = LowQ( Base[ i ] , Base[ j ] );
   HpNum xx = xy - ScalarProduct( L[ i ] , L[ j ] , j + 1 );

   if( CMP( xx , xy ) )
    *MQLog << endl << "Test failed: D[" << i << "][" << j << "] = " << xx;
   }

 // L * a = Alfa{Base}
 SolvePerSys( L , Alfa , tmpa , Base , BDim - Dependent );

 for( i = 0 ; i < BDim - Dependent ; i++ )
 {
  HpNum xx = z2[ i ] - tmpa[ i ];

  if( CMP( xx , tmpa[ i ] ) )
   *MQLog << endl << "Test failed: z2[" << i << "] = " << xx;
  }

 if( CBDim )
  for( i = 0 ; i < BDim - Dependent ; i++ )
   tmpv[ i ] = e( Base[ i ] , 1 );
 else
  VectAssign( tmpv , 1 , BDim - Dependent );

 SolveSystem( L , tmpv , tmpa , BDim - Dependent );  // L * a = e({Base})

 for( i = 0 ; i < BDim - Dependent ; i++ )
 {
  HpNum xx = z1[ i ] - tmpa[ i ];

  if( CMP( xx , tmpa[ i ] ) )
   *MQLog << endl << "Test failed: z1[" << i << "] = " << xx;
  }

 HpNum xy , xx;

 #if( ! EXACT )
  xx = ( xy = Norm( z1 , BDim - Dependent ) ) - z1Norm;

  if( CMP( xx , xy ) )
   *MQLog << endl << "Test failed: z1Norm = " << xx;

  xx = ( xy = ScalarProduct( z1 , z2 , BDim - Dependent ) ) - z1Tz2;

  if( CMP( xx , xy ) )
   *MQLog << endl << "Test failed: z1Tz2 = " << xx;

  if( Dependent )
  {
   xx = ( xy = Norm( L[ BDim - Dependent ] , BDim - Dependent ) ) - LhNorm;

   if( CMP( xx , xy ) )
    *MQLog << endl << "Test failed: LhNorm = " << xx;

   if( Dependent > 1 )
   {
    xx = ( xy = Norm( L[ BDim - 1 ] , BDim - Dependent ) ) - LjNorm;

    if( CMP( xx , xy ) )
     *MQLog << endl << "Test failed: LjNorm = " << xx;
    }
   }
 #endif

 register HpNum Max = ABS( **L );
 register HpNum min = Max;
 register Index h = 1;
 register HpMat pL = L + h;

 for( ; h < BDim - Dependent ; h++ , pL++ )
  if( Max < ABS( (*pL)[ h ] ) )
   Max = ABS( (*pL)[ h ] );
  else
   if( min > ABS( (*pL)[ h ] ) )
     min = ABS( (*pL)[ h ] );

 xx = ( xy = Max / min ) - Lmu;

 if( CMP( xx , HpNum( 1 ) ) )
  *MQLog << endl << "Test failed: Lmu = " << xx;

 }  // end( CheckBaseValidity )

/*--------------------------------------------------------------------------*/

void MinQuad::CheckQValidity( void )
{
 register Index cnt = 0;
 for( register Index i = 0 ; i < MaxBDim ; i++ )
  if( G[ i ] )
   cnt++;

 assert( cnt == ActBDim );

 register Index_Set o = Order;
 for( register Index i = ActBDim ; i-- ; )
  assert( G[ *(o++) ] );

 for( register Index i = ActBDim ; i ; )
 {
  register Index j = i;
  register Index h = Order[ --i ];
  register Index_Set jth = Order;

  #if( LAZY_Q )
   for( ; j ; j-- , jth++ )
    if( h >= *jth )
    {
     if( Q[ h ][ *jth ] > - QuINF )
      if( Q[ h ][ *jth ] != GiTGj( h , *jth ) )
       *MQLog << endl << "Test failed: Q[" << h << "][" << *jth << "]";
     }
    else
     if( Q[ *jth ][ h ] > - QuINF )
      if( Q[ *jth ][ h ] != GiTGj( h , *jth ) )
       *MQLog << endl << "Test failed: Q[" << *jth << "][" << h << "]";
  #else
   GiTG( h , tempQ , NxtBIdx );

   for( ; j ; j-- , jth++ )
    if( h >= *jth )
    {
     if( Q[ h ][ *jth ] != tempQ[ *jth ] )
      *MQLog << endl << "Test failed: Q[" << h << "][" << *jth << "]";
     }
    else
     if( Q[ *jth ][ h ] != tempQ[ *jth ] )
      *MQLog << endl << "Test failed: Q[" << *jth << "][" << h << "]";
  #endif
  }
 }  // end( CheckQValidity )

#endif

/*--------------------------------------------------------------------------*/
/*----------------------- End File MinQuad.C -------------------------------*/
/*--------------------------------------------------------------------------*/
