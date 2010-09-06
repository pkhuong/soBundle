/*--------------------------------------------------------------------------*/
/*----------------------- File MainAlg.C -----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                							  --*/
/*-- Example on how to use the [B]TT algorithm within a NDO solver        --*/
/*--                							  --*/
/*--                            VERSION 1.01       			  --*/
/*--                	       06 - 12 - 1996			       	  --*/
/*--                							  --*/
/*-- 		     Original Idea and Implementation by:		  --*/
/*--                                                                      --*/
/*--			      Antonio Frangioni        			  --*/
/*--                                                                      --*/
/*--   			   Operations Research Group			  --*/
/*--			  Dipartimento di Informatica			  --*/
/*--   			     Universita' di Pisa			  --*/
/*--                							  --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define NEEDCONSTR 1

/* This macro is defined for this example: if NEEDCONSTR != 0, it is assumed
   that the NDO solver needs box constraints on the direction. This is done
   in order to show how a large part of the interface is not affected by
   the presence of constraints. */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#if( NEEDCONSTR )
 #include "BMinQuad.h"
#else
 #include "MinQuad.h"
#endif

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

namespace MinQuad_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/

/* Define the class containing your favourite algorithm, where you want
   to use [B]TT. The class has to be a derived class from [B]MinQuad, in
   order for you to provide the proper implementation of some methods. */

class MainAlg :
               #if( NEEDCONSTR )
	        public MinQuad
               #else
	        public BMinQuad
               #endif
{
 public:

 < your public methods, constructor, destructor ... >

 private:

 Index LamDim;         // current length of the subgradients

 SgMat Subg;           // a typical "naive" implementation of the Bundle:
                       // Subg[ i ] = i-th subgradient, so that
                       // Subg[ i ][ j ] = j-th entry of the i-th subgrad.

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 #if( ! NEEDCONSTR )

/*--------------------------------------------------------------------------*/
/*--                    Interface with MinQuad                            --*/
/*--------------------------------------------------------------------------*/

  #if( LAZY_Q )

   QuNum GiTGj( cIndex i , cIndex j )
   {
    // just compute the scalar product between Subg[ i ] and Subg[ j ]

    QuNum t = 0;
    for( Index h = 0 ; h < LamDim ; h++ )
     t += Subg[ i ][ h ] * Subg[ j ][ h ];

    return( t );
    }

  #else

   void GiTG( cIndex i , register QuRow Qi , cIndex iMax )
   {
    // compute the scalar product between Subg[ i ] and Subg[ j ] for every
    // existing item j

    for( Index j = 0 ; j < iMax ; j++ )
     if( IsThere( j ) )
     {
      QuNum t = 0;
      for( Index h = 0 ; h < LamDim ; h++ )
       t += Subg[ i ][ h ] * Subg[ j ][ h ];

      Qi[ j ] = t;
      }
    }

  #endif

 #else

/*--------------------------------------------------------------------------*/
/*--                   Interface with BMinQuad                            --*/
/*--------------------------------------------------------------------------*/

  #if( LAZY_Q )

   QuNum GiTGj( cIndex i , cIndex j )
   {
    // compute the scalar product between Subg[ i ] and Subg[ j ] "restricted"
    // to the MB2Dim indices in the InINF-terminated set MBase2[]

    QuNum t = 0;
    cIndex_Set tMB2 = MBase2;
    for( Index h ; ( h = *(tMB2++) ) < InINF ; )
     t += Subg[ i ][ h ] * Subg[ j ][ h ];

    return( t );
    }

  #else

   void GiTG( cIndex i , register QuRow Qi , cIndex iMax )
   {
    // compute the scalar product between Subg[ i ] and Subg[ j ], "restricted"
    // to the MB2Dim indices in the InINF-terminated set MBase2[], for every
    // existing item j

    for( Index j = 0 ; j < iMax ; j++ )
     if( IsThere( j ) )
     {
      QuNum t = 0;
      cIndex_Set tMB2 = MBase2;
      for( Index h ; ( h = *(tMB2++) ) < InINF ; )
       t += Subg[ i ][ h ] * Subg[ j ][ h ];

      Qi[ j ] = t;
      }
    }

  #endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  void GiTilde( const Index i , register SgRow GTld )
  {
   for( Index h = 0 ; h < MaxItemN() ; h++ )
    GTld[ h ] = Subg[ h ][ i ];
   }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  #if( LAZY_D == 2 )

   LMNum CalculateZ( cIndex h )
   {
    LMNum Zh = 0;
    for( Index i = 0 ; i < BDim ; i++ )
     Zh += Subg[ Base[ i ] ][ h ] * Mult[ i ];

    return( Zh );
    }

  #elif( LAZY_D == 1 )

   void CalculateZ( register cIndex_Set Wh , register LMRow z )
   {
    for( Index h ; ( h = *(Wh++) ) < InINF ; )
    {
     LMNum Zh = 0;
     for( Index i = 0 ; i < BDim ; i++ )
      Zh += Subg[ Base[ i ] ][ h ] * Mult[ i ];

     z[ h ] = Zh;
     }
    }

  #else

   void CalculateZ( register LMRow z )
   {
    for( Index h = 0 ; h < LamDim ; )
     z[ h++ ] = 0;

    for( Index i = 0 ; i < BDim ; i++ )
     for( Index h = 0 ; h < LamDim ; h++ )
      z[ h++ ] += Subg[ Base[ i ] ][ h ] * Mult[ i ];
    }

  #endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  LMNum GiTLB( cIndex i , register cLMRow l )
  {
   LMNum res = 0;

   for( Index j = 0 ; j < B2Dim ; j++ )
    res += Subg[ Base2[ j ] ] * l[ Base2[ j ] ];

   return( res );
   }

/*--------------------------------------------------------------------------*/

#endif

/*--------------------------------------------------------------------------*/
/*--                          Actual Algorithm                            --*/
/*--------------------------------------------------------------------------*/

 void MyActualAlgorithm( void )
 {
  // to add a subgradient/constraint with "name" i (Subg[ i ]) to the bundle

  AddSubGrad( i , < alfa of i > );
  AddConstr( i , < alfa of i > );

  // to remove the subgradient with "name" i

  ResetBundle( i );

  // therefore, things like the "B-strategy" can be very easily handled

  // to change the linearization errors (because of a movment in the dual
  // space), a straightforward possibility is

  ChangeAlfa( < vector of new alfas > );
  ChangeAlfa( i , < new alfa of i > );

  // but if the movemnt has been along the latest direction d, it's better
  // to use

  MoveAlongD( Tau , DeltaFi );

  // and read [a read-only pointer to] the vector of new alfas with

  ReadAlfa( NewAlfa );  [ or NewAlfa = ReadAlfa(); ]

  // to solve the current problem with trust-region parameter ti

  MQError Status = SolveQP( ti );

  // to read the optimal solution and the corresponding value

  cHpRow M = ReadMult();
  cIndex_Set B = ReadBase();
  Index Bd = ReadBDim();

  HpNum fVal = (1 / 2) * ReaddNorm() + ReadSigma() / ti

  #if( NEEDCONSTR )
   // a read-only pointer to the vector z = - d / ti can be queried with

   cLMRow z = ReadZ();
  #endif

  // in case they are useful, [a read-only pointer to] the vector of scalar
  // products between the latest d and all the subgradients can be obtained as

  ReadGTd( gtd );  [ or gtd = ReadGTd(); ]

  // unconstrained dual variables can be added and removed with

  AddSGSpaceDim( NewDim );
  RemoveVariable( OldDim );

  #if( NEEDCONSTR )
   // and in the constrained case, another whole bunch of methods are
   // provided for handling constrained and unconstrained variables

   AddVar( i );
   RemoveVar( i );

   MakeVarCostr( i );
   MakeVarUnCnstr( i );

   // for giving "starting informations"

   InitialSetUp();
  #endif

  // other methods are available for minor purposes like obtaining
  // sensitivity analisys informations ...

  }   // end( MyActualAlgorithm )

 };  // end( class MainAlg )

/*--------------------------------------------------------------------------*/

 };  // end( namespace MinQuad_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------- End File MainAlg.C -----------------------------*/
/*--------------------------------------------------------------------------*/
