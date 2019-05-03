#include "FunctionUtilities.h"
#include <math.h>
#include <string.h>

double sgn(double val) {
	return (double(0) < val) - (val < double(0));
}


/*
These four functions are taken from Matt's:
UCLA CAM 18-26, Matt Jacobs, Flavien Leger, Wuchen Li, and Stanley Osher, Solving Large-Scale Optimization Problems with a Convergence Rate Independent of Grid Size, May 2018
*/
double CubicSolve(double b, double c, double d) {
	constexpr double real3rdRoot1 = -.5;   // equals cos(2*M_PI/3);
	constexpr double im3rdRoot1 = 0.86602540378;   //equals sin(2*M_PI/3);
	constexpr double real3rdRoot2 = -.5;  //  equals cos(4*M_PI/3)=real3rdRoot1;
	constexpr double im3rdRoot2 = -0.86602540378;  //equals sin(4*M_PI/3)=-im3rdRoot1;

	double b3over3 = (b / 3)*(b / 3)*(b / 3);

	double p = c - b*(b / 3);
	double q = d + 2 * b3over3 - b*(c / 3);
	double solution = 0;

	if (p == 0) {

		solution = -sgn(q)*exp(log(fabs(q)) / 3.0);

	}
	else {
		double discrim = (q / 2)*(q / 2) + (p / 3)*(p / 3)*(p / 3);

		double s = sqrt(fabs(discrim));

		if (discrim<0) {

			double theta = atan2(s, -q / 2);

			double x = s*s + q*q / 4;
			double rc = exp(log(x) / 6);

			double thetac = theta / 3;

			double real = rc*cos(thetac);
			double im = rc*sin(thetac);

			double solution1 = 2 * real;


			double solution2 = 2 * (real*real3rdRoot1 - im*im3rdRoot1);
			double solution3 = 2 * (real*real3rdRoot2 - im*im3rdRoot2);

			solution = fmax(solution1, fmax(solution2, solution3));


		}
		else if (discrim>0) {

			double u3 = -q / 2 + s;
			double v3 = -q / 2 - s;

			double u = sgn(u3)*exp(log(fabs(u3)) / 3);
			double v = sgn(v3)*exp(log(fabs(v3)) / 3);

			solution = u + v;

		}
		else {
			solution = fmax(3 * q / p, -3 * q / (2 * p));

		}
	}

	return solution - b / 3;

}

double GetRealRootsOfCubic(double a, double b, double c, double d) {
	return CubicSolve(b / a, c / a, d / a);
}


/* Calculate the divergence of a vector field

THIS FUNCTION IS VERY DELICATE MODIFY AT YOUR OWN RISK

div grad u == laplacian u otherwise code will fail

*/
double calc_divergence(double *ux, double *uy, double *divergence, int n1, int n2) {
  double divTot = 0;
  for (int i = 0; i<n2 - 1; i++) {
    for (int j = 0; j<n1 - 1; j++) {

      int xp = (j + 1);
      int yp = (i + 1);

      divergence[i*n1 + j] = n1*(ux[i*n1 + xp] - ux[i*n1 + j]) + n2*(uy[yp*n1 + j] - uy[i*n1 + j]);
      divTot += pow(divergence[i*n1 + j], 2);
    }
    divergence[i*n1 + n1 - 1] = -n1*ux[i*n1 + n1 - 1] + n2*(uy[(i + 1)*n1 + n1 - 1] - uy[i*n1 + n1 - 1]);
    divTot += pow(divergence[i*n1 + n1 - 1], 2);
  }
  for (int j = 0; j<n1 - 1; j++) {
    divergence[(n2 - 1)*n1 + j] = n1*(ux[(n2 - 1)*n1 + j + 1] - ux[(n2 - 1)*n1 + j]) - n2*uy[(n2 - 1)*n1 + j];
    divTot += pow(divergence[(n2 - 1)*n1 + j], 2);

  }
  divergence[n2*n1 - 1] = -(n1*ux[n1*n2 - 1] + n2*uy[n1*n2 - 1]);
  divTot += pow(divergence[n2*n1 - 1], 2);

  divTot /= n1*n2;
  divTot = sqrt(divTot);
  return divTot;
}

/* Calculate the gradient of a function

THIS FUNCTION IS VERY DELICATE MODIFY AT YOUR OWN RISK

div grad u == laplacian u otherwise code will fail

*/

void calc_gradient(double *ux, double *uy, double *potential, int n1, int n2) {

  ux[0] = 0;
  memset(uy, 0, n1 * sizeof(double));
  for (int j = 1; j<n1; j++) {
    ux[j] = n1*(potential[j] - potential[j - 1]);
  }


  for (int i = 1; i<n2; i++) {
    ux[n1*i] = 0;
    uy[n1*i] = n2*(potential[i*n1] - potential[(i - 1)*n1]);
    for (int j = 1; j<n1; j++) {
      int xm = j - 1;
      int ym = i - 1;
      ux[i*n1 + j] = n1*(potential[i*n1 + j] - potential[i*n1 + xm]);
      uy[i*n1 + j] = n2*(potential[i*n1 + j] - potential[ym*n1 + j]);
    }
  }

}