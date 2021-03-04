////////////////////////////////////////////
//This is a function to calc b^(m-2)
//
//For WEST use only
//
//Author:John Lombardi
////////////////////////////////////////////

#include <math.h>
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo,RcppEigen)]]

//While loop is slightly faster
//Order of code is as follows:
//For <4 multiplications we can just do regular matrix multiplication
//We can maybe get rid of the if(mult < 4) but I figure we save some time
//By not going through the if/else. probably not a significant amount of time

//The first part of code uses * which is overloaded. We could use a different version but
//the documentation of rcppEigen is pretty annoying. You can use some sort of triCPP function
//but I had trouble implementing it
//triCPP vs * is maybe 50% faster? So worth to figure out eventually

//Second part of code is using matrix mult by squaring to save computations
//Useful when doing lots of multiplic: a*a*a*a*a*a
//If doing A ^ 13 - we do (A^8 * A^2) * A = 5 total multiplic
//vs A*A*A.... 12 time

// [[Rcpp::export]]

SEXP matrixCalcFastTwo(Eigen::MatrixXd & matrixOne, int multTimes,
                       Eigen::MatrixXd & matrixOfDiag)
{
  //First, do matrix by squaring if num multiplic > 4
  int i = 0;

  Eigen::MatrixXd  matOutput = matrixOne;

  //normal matrix mult
  if (multTimes < 4)
  {
    while (i < multTimes)
    {
      matOutput = matOutput * matrixOne;

      i = i + 1;
    }

    return Rcpp::wrap(matOutput);

  }
  else //matrix by squaring
  {

    //Figure out how many multiplications we will do
    int numTimes = multTimes;

    while (numTimes > 1)
    {
      //Even case
      if (numTimes % 2 == 0)
      {
        matrixOne = matrixOne * matrixOne;

        numTimes = numTimes / 2;

      }//end even case
      //odd case
      else {

        //Use matrixOfDiag

        //First save A ^ 1
        matrixOfDiag = matrixOne * matrixOfDiag;

        //Then Do A^2
        matrixOne = matrixOne * matrixOne;

        //Get to an even case
        numTimes = (numTimes - 1) / 2;


      }//end odd case

    }//end WHILE

    //Multiply here: either will be matrix * (1) or
    //matrix * A ^ 1
    matrixOne = matrixOne * matrixOfDiag;

    return Rcpp::wrap(matrixOne);
  }// end matrix by squaring



}// end foo
