#include "TMath.h"
#include <iostream>
#include <stdio.h>
#include <math.h>

/////////////////////////////////////////////
//// A program to detect when 2 ellipses touch: i.e. ellipse edge detection
//// Ellipses are determined to touch according to the criteria presented in the Etayo et al 2006 paper:
//// https://doi.org/10.1016/j.cagd.2006.01.002
//// -----------------------------------------------------------
//// Takes the ellipse parameters (5 for each ellipse), and outputs whether they overlap, or not
//// The mathematics has already been worked out: explicit matrix parameters, cubic coefficients from
//// the pencil polynomial, thus, this is a quick process to determine ellipse collision
//// Author: Luke Batten
////////////////////////////////////////////

// Variables: angle (rads)
// semiMajor, semiminor, h, k (a.u.)

// Ellipse 1:
Double_t angle1 = 20 * TMath::Pi()/180;
Double_t semiMajor1 = 5;
Double_t semiMinor1 = 4;
Double_t h1 = 2;
Double_t k1 = 10;

// Ellipse 2:
Double_t angle2 = 90 * TMath::Pi()/180;
Double_t semiMajor2 = 4;
Double_t semiMinor2 = 1;
Double_t h2 = 2;
Double_t k2 = 15;

// Generate the matrices from the ellipse parameters
Double_t** getExplicitMatrix(Double_t angle, Double_t semiMajor, Double_t semiMinor, Double_t h, Double_t k)
{
  Double_t** explicitMatrix = 0;
  Int_t dim = 3;
  
  explicitMatrix = new Double_t*[dim];

  for (Int_t he = 0; he < dim; he++)
    {
      explicitMatrix[he] = new Double_t[he];
    }

  // a11
  explicitMatrix[0][0] = pow((cos(angle)/semiMajor),2) + pow((sin(angle)/semiMinor),2);

  // a12
  explicitMatrix[0][1] = (sin(angle)*cos(angle))/pow(semiMajor,2) - (sin(angle)*cos(angle))/pow(semiMinor,2);
  // a21
  explicitMatrix[1][0] = explicitMatrix[0][1];

  // a22
  explicitMatrix[1][1] = pow((sin(angle)/semiMajor),2) + pow((cos(angle)/semiMinor),2); //

  // a13
  explicitMatrix[0][2] = (-cos(angle) * ( h * cos(angle) + k*sin(angle) ))/pow(semiMajor,2) + (sin(angle) * ( k * cos(angle) - h*sin(angle) ))/pow(semiMinor,2);
  // a31
  explicitMatrix[2][0] = explicitMatrix[0][2];

  // a23
  explicitMatrix[1][2] = (-sin(angle) * ( h * cos(angle) + k*sin(angle) ))/pow(semiMajor,2) + (cos(angle) * ( h * sin(angle) - k*cos(angle) ))/pow(semiMinor,2);
  // a32
  explicitMatrix[2][1] = explicitMatrix[1][2];

  //a33
  explicitMatrix[2][2] = pow(((h*cos(angle) + k*sin(angle)))/(semiMajor),2) + pow(((h*sin(angle) - k*cos(angle)))/(semiMinor1),2)  - 1; // 

  
  return explicitMatrix;
}

// Right up to here

// Attain the pencil characteristic polynomial coefficients
// PCPE = pencil characteristic polynomial coefficients
// We attain a cubic polynomial in lambda, just the coefficients are extracted
Double_t* getPCPE(Double_t **matrix1, Double_t **matrix2)
{
  const Int_t numCoefficients = 3; // There are technically 4, but we make one of them monic, are required
  static Double_t polynomialCoefficients[numCoefficients] = {};

  // The pencil determinant returns 48 terms, which make up the 4 coefficients
  // As there are many terms, try and make things a little neater
  // Matrix 1:
  Double_t a11 = matrix1[0][0]; Double_t a12 = matrix1[0][1]; Double_t a13 = matrix1[0][2];
  Double_t a21 = matrix1[1][0]; Double_t a22 = matrix1[1][1]; Double_t a23 = matrix1[1][2];
  Double_t a31 = matrix1[2][0]; Double_t a32 = matrix1[2][1]; Double_t a33 = matrix1[2][2];
  // Matrix 2:
  Double_t b11 = matrix2[0][0]; Double_t b12 = matrix2[0][1]; Double_t b13 = matrix2[0][2];
  Double_t b21 = matrix2[1][0]; Double_t b22 = matrix2[1][1]; Double_t b23 = matrix2[1][2];
  Double_t b31 = matrix2[2][0]; Double_t b32 = matrix2[2][1]; Double_t b33 = matrix2[2][2];

  //// 3rd order
  Double_t dNonMonic = -a13*a22*a31 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 + a11*a22*a33;
  
  //// 0th order
  Double_t cNonMonic = -b13*b22*b31 + b12*b23*b31 + b13*b21*b32 - b11*b23*b32 - b12*b21*b33 + b11*b22*b33;
  // c
  polynomialCoefficients[2] = cNonMonic/dNonMonic;
  
  //// 1st order
  Double_t bNonMonicTerms1_6 = -a33*b12*b21 + a32*b13*b21 + a33*b11*b22 - a31*b13*b22 - a32*b11*b23 + a31*b12*b23;
  Double_t bNonMonicTerms7_12 = a23*b12*b31 - a22*b13*b31 - a13*b22*b31 + a12*b23*b31 - a23*b11*b32 + a21*b13*b32;
  Double_t bNonMonicTerms13_18 = a13*b21*b32 - a11*b23*b32 + a22*b11*b33 - a21*b12*b33 - a12*b21*b33 + a11*b22*b33;
  //
  Double_t bNonMonic =  bNonMonicTerms1_6 + bNonMonicTerms7_12 + bNonMonicTerms13_18;
  // b
  polynomialCoefficients[1] = bNonMonic/dNonMonic;

  //// 2nd order
  Double_t aNonMonicTerms1_6 = -a23*a32*b11 + a22*a33*b11 + a23*a31*b12 - a21*a33*b12 - a22*a31*b13 + a21*a32*b13;
  Double_t aNonMonicTerms7_12 = a13*a32*b21 - a12*a33*b21 - a13*a31*b22 + a11*a33*b22 + a12*a31*b23 - a11*a32*b23;
  Double_t aNonMonicTerms13_18 = -a13*a22*b31 + a12*a23*b31 + a13*a21*b32 - a11*a23*b32 - a12*a21*b33 + a11*a22*b33;
  //
  Double_t aNonMonic = aNonMonicTerms1_6 + aNonMonicTerms7_12 + aNonMonicTerms13_18;
  // a
  polynomialCoefficients[0] = aNonMonic/dNonMonic;

  // Return the coefficients the characteristic polynomial
  return polynomialCoefficients;
}

// To use this function:
// Example which uses the parameters set at the beginning of the program: overlapEllipses()
// General usage by entering 10 parameters: overlapEllipses(0.6, 1, 2, 3, 4, 1.3, 5, 6, 7, 8)
// Note: angles are in RADIANS, NOT degrees
// Verbose option is off by default, if you want it on add 1 to the end of your function list

void overlapEllipses(Double_t angleA = angle1, Double_t semiMajorA = semiMajor1, Double_t semiMinorA = semiMinor1, Double_t hA = h1, Double_t kA = k1, Double_t angleB = angle2, Double_t semiMajorB = semiMajor2, Double_t semiMinorB = semiMinor2, Double_t hB = h2, Double_t kB = k2,  Bool_t verbose = 0)
{

  std::cout << "Determining if ellipses overlap..." << std::endl;

  // Find matrices used to find the characteristic polynomial
  Double_t **matrixA = getExplicitMatrix(angleA, semiMajorA, semiMinorA, hA, kA);
  Double_t **matrixB = getExplicitMatrix(angleA, semiMajorA, semiMinorA, hA, kA);

  if(verbose == 1)
    {
      
      for(Int_t i=0;i<3;i++)
	{
	  for(Int_t j=0;j<3;j++)
	    {
	      std::cout<<matrixA[i][j] << "   ";
	    }
	  std::cout << std::endl;
	}

      std::cout << std::endl;

      for(Int_t i=0;i<3;i++)
	{
	  for(Int_t j=0;j<3;j++)
	    {
	      std::cout<<matrixB[i][j] << "   ";
	    }
	  std::cout << std::endl;
	}
  
    }

  // Find the characteristic polynomial of the pencil:
  Double_t *f = getPCPE(matrixA, matrixB);
  
  Double_t a = f[0];
  Double_t b = f[1];
  Double_t c = f[2];

  if(verbose == 1)
    {
      std::cout << "a = " << a << std::endl;
      std::cout << "b = " << b << std::endl;
      std::cout << "c = " << c << std::endl;
    }

  TText *t21;
  Double_t textPositionX = 0.2;
  Double_t textPositionY = 0.8;
  
  // PCPEs have to pass either condition 1 or 2 to overlap
  if(a >= 0 && ( -3*b + pow(a,2) ) > 0 && ( 3*a*c + b*pow(a,2) - 4*pow(b,2) ) < 0 && ( -27*pow(c,2) + 18*c*b*a + pow(a,2)*pow(b,2) - 4*pow(a,3)*c - 4*pow(b,3)) > 0 )
    {
      std::cout << "Ellipses do NOT overlap!" << std::endl;
      t21 = new TText(textPositionX, textPositionY,"Ellipses do NOT overlap");
      t21->SetTextColor(kRed);
    }
  else if(a < 0 && ( -3*b + pow(a,2) ) > 0 && ( -27*pow(c,2) + 18*c*b*a + pow(a,2)*pow(b,2) - 4*pow(a,3)*c - 4*pow(b,3)) > 0)
    {
      std::cout << "Ellipses do NOT overlap!" << std::endl;
      t21 = new TText(textPositionX, textPositionY,"Ellipses do NOT overlap");
      t21->SetTextColor(kRed);
    }
  else
    {
      std::cout << "Ellipses overlap!" << std::endl;
      t21 = new TText(textPositionX, textPositionY,"Ellipses overlap");
      t21->SetTextColor(kGreen);
    }

  /// Sketch ellipses

  const Int_t numEllipses = 2;
  Double_t coordX[numEllipses] = {};
  Double_t coordY[numEllipses] = {};
  coordX[0] = hA;  coordX[1] = hA;  coordY[0] = kA;  coordY[1] = kA; 
  TGraph *grEllipsePoints1 = new TGraph(numEllipses,coordX,coordY);
  grEllipsePoints1->Draw("ap");

  // Must convert back to degrees
  TEllipse *e1 = new TEllipse(hA, kA, semiMajorA, semiMinorA, 0, 360, angleA * 180/TMath::Pi());
  TEllipse *e2 = new TEllipse(hB, kB, semiMajorB, semiMinorB, 0, 360, angleB * 180/TMath::Pi());
  e1->Draw("same");
  e1->SetFillColorAlpha(2, 0.5);
  e2->Draw("same");
  e2->SetFillColorAlpha(4, 0.5);

  t21->SetTextSize(0.04);
  t21->Draw();

  grEllipsePoints1->SetTitle("Ellipse plot");
  grEllipsePoints1->GetXaxis()->SetTitle("x");
  grEllipsePoints1->GetYaxis()->SetTitle("y");

  Int_t minX, maxX, minY, maxY = 0;

  // x range
  if( (hA-semiMajorA) < (hB-semiMajorB))
    {
      minX = hA - semiMajorA*2;
    }
  else
    {
      minX = hB - semiMajorB*2;
    }
  
  if( (hA+semiMajorA) > (hB+semiMajorB))
    {
      maxX = hA + semiMajorA*2;
    }
  else
    {
      maxX = hB + semiMajorB*2;
    }

  // y range
  if( (kA-semiMajorA) < (kB-semiMajorB))
    {
      minY = kA - semiMajorA*2;
    }
  else
    {
      minY = kB - semiMajorB*2;
    }
  
  if( (kA+semiMajorA) > (kB+semiMajorB))
    {
      maxY = kA + semiMajorA*2;
    }
  else
    {
      maxY = kB + semiMajorB*2;
    }
  
  grEllipsePoints1->GetXaxis()->SetLimits(minX,maxX);
  grEllipsePoints1->SetMaximum(maxY);
  grEllipsePoints1->SetMinimum(minY);
  
  return;
  
}
