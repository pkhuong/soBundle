/*--------------------------------------------------------------------------*/
/*----------------------------- File Main.C --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  Simple main() for using the [QP]Bundle class.                       --*/
/*--                                                                      --*/
/*--                            VERSION 1.21                              --*/
/*--                           18 - 10 - 2004                             --*/
/*--                							  --*/
/*-- 		     Original Idea and Implementation by:		  --*/
/*--                							  --*/
/*--			       Antonio Frangioni       			  --*/
/*--                							  --*/
/*--   			   Operations Research Group			  --*/
/*--			  Dipartimento di Informatica			  --*/
/*--   			     Universita' di Pisa			  --*/
/*--                                                                      --*/
/*--                         Copyright 1996 - 2004                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "bundleXX.h"

#include <fstream>

/*--------------------------------------------------------------------------*/
/*------------------------------- USING ------------------------------------*/
/*--------------------------------------------------------------------------*/

#if( OPT_USE_NAMESPACES )
using namespace Bundle_di_unipi_it;
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 // construct the solver object - - - - - - - - - - - - - - - - - - - - - - -

 ifstream iStrm( "ParValue" );

 // extract the number of variables at the beginning of the stream- - - - - -
 // note how a unique stream can be a convenient place for putting all the
 // parameters

 Index NFiVar;
 DfltdSfInpt( &iStrm , NFiVar , Index( 10 ) );

 // construct the solver object - - - - - - - - - - - - - - - - - - - - - - -

 Bundle *b = new bundleXX( &iStrm , NFiVar );

 // if HV_NNVAR == 0 (see QPBundle.h), all the variables *must* be- - - - - -
 // unconstrained in sign: this is done by the next call- - - - - - - - - - -

 #if( ! HV_NNVAR )
  b->SetUC( 0 , NFiVar );
 #endif

 // if the log has to be done, it is done on clog - - - - - - - - - - - - - -

 #if( LOG_BND )
  b->SetBLog( &clog , 2 );
 #endif

 // minimize the function - - - - - - - - - - - - - - - - - - - - - - - - - -

 b->Solve();

 // print out results - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 clog << endl << "Fi^* = " << b->ReadFiVal() << endl;

 LMRow LStar = new LMNum[ b->GetNumVar() ];
 b->ReadSol( LStar );

 for( Index i = 0 ; i < b->GetNumVar() ; i++ )
  clog << LStar[ i ] << " ";

 clog << endl;

 delete[] LStar;

 // delete the solver object- - - - - - - - - - - - - - - - - - - - - - - - -

 delete( b );

 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Main.C -------------------------------*/
/*--------------------------------------------------------------------------*/
