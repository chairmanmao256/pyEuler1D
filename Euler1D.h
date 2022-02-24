// Finite volume method for 1-D Euler equation
#ifndef EULER1D_H_
#define EULER1D_H_

#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>

class Euler1D
{
private:
     std::string boundary_type;
     int ncells; // The number of cellls. The virtual cells are not included.
     double** sol = new double* [3];      // The cell-center solution vector.               ncells columns
     double** sol_old = new double* [3];  // Store the cell-center value at the last step   ncells columns
     double** rec_l = new double* [3];    // The reconstructed left value at interfaces.    ncells + 1 columns 
     double** rec_r = new double* [3];    // The reconstructed right value at interfaces.   ncells + 1 columns
     double** flux = new double* [3];     // The fluxes at interfaces.                      ncells + 1 columns
     double** avg = new double* [3];      // Averaged primitive variables at interfaces.    ncells + 1 columns
     double* mesh;   // The x mesh of the solution.
     double vL[3][1], vLL[3][1], vR[3][1], vRR[3][1]; // Virtual cells to impose boundary conditions. (Second order reconstruction is used)
     double* get_eig(double u, double a)
     {
         double eps;
         double* lambda = new double[3]; //store the eig value
         eps = 0.1*a;  // Threshold to invoke entropy modification
         if (abs(u-a) < eps)
         {
             lambda[0] = (pow(abs(u-a),2) + pow(eps,2))/(2.0*eps);
         }
         else
         {
             lambda[0] = abs(u-a);
         }
         if (abs(u) < eps)
         {
             lambda[1] = (pow(abs(u),2) + pow(eps,2))/(2.0*eps);
         }
         else
         {
             lambda[1] = abs(u);
         }
         if (abs(u + a) < eps)
         {
             lambda[2] = (pow(abs(u+a),2) + pow(eps,2))/(2.0*eps);
         }
         else
         {
             lambda[2] = abs(u+a);
         }
         return lambda;
     }
     void prim2con(int flag)
     {
         double temp0,temp1,temp2;
         if (flag == 0) // do it on cell center values
         {
             for (int i = 0; i < ncells; i++)
             {
                 temp0 = sol[0][i];
                 temp1 = sol[0][i]*sol[1][i];
                 temp2 = sol[0][i]*(0.5*sol[1][i]*sol[1][i]) + 1.0/(1.4-1.0)*sol[2][i];
                 sol[0][i] = temp0;
                 sol[1][i] = temp1;
                 sol[2][i] = temp2;
             }
             temp0 = vL[0][0];
             temp1 = vL[0][0]*vL[1][0];
             temp2 = vL[0][0]*(0.5*vL[1][0]*vL[1][0]) + 1.0/(1.4-1.0)*vL[2][0];
             vL[0][0] = temp0;
             vL[1][0] = temp1;
             vL[2][0] = temp2;

             temp0 = vLL[0][0];
             temp1 = vLL[0][0]*vLL[1][0];
             temp2 = vLL[0][0]*(0.5*vLL[1][0]*vLL[1][0]) + 1.0/(1.4-1.0)*vLL[2][0];
             vLL[0][0] = temp0;
             vLL[1][0] = temp1;
             vLL[2][0] = temp2;

             temp0 = vR[0][0];
             temp1 = vR[0][0]*vR[1][0];
             temp2 = vR[0][0]*(0.5*vR[1][0]*vR[1][0]) + 1.0/(1.4-1.0)*vR[2][0];
             vR[0][0] = temp0;
             vR[1][0] = temp1;
             vR[2][0] = temp2;

             temp0 = vRR[0][0];
             temp1 = vRR[0][0]*vRR[1][0];
             temp2 = vRR[0][0]*(0.5*vRR[1][0]*vRR[1][0]) + 1.0/(1.4-1.0)*vRR[2][0];
             vRR[0][0] = temp0;
             vRR[1][0] = temp1;
             vRR[2][0] = temp2;

         }
         else // do it on reconstructed values
         {
             for (int i = 0; i < ncells+1; i++)
             {
                 temp0 = rec_l[0][i];
                 temp1 = rec_l[0][i]*rec_l[1][i];
                 temp2 = rec_l[0][i]*(0.5*rec_l[1][i]*rec_l[1][i]) + 1.0/(1.4-1.0)*rec_l[2][i];
                 rec_l[0][i] = temp0;
                 rec_l[1][i] = temp1;
                 rec_l[2][i] = temp2;

                 temp0 = rec_r[0][i];
                 temp1 = rec_r[0][i]*rec_r[1][i];
                 temp2 = rec_r[0][i]*(0.5*rec_r[1][i]*rec_r[1][i]) + 1.0/(1.4-1.0)*rec_r[2][i];
                 rec_r[0][i] = temp0;
                 rec_r[1][i] = temp1;
                 rec_r[2][i] = temp2;
             }
         }
     }
     void con2prim(int flag)
     {
         double temp0,temp1,temp2;
         if (flag == 0) // do it on cell center values
         {
             for (int i = 0; i < ncells; i++)
             {
                 temp0 = sol[0][i];
                 temp1 = sol[1][i]/sol[0][i];
                 temp2 = (1.4-1.0)*sol[2][i]-(1.4-1.0)/2.0*sol[1][i]*sol[1][i]/sol[0][i];
                 sol[0][i] = temp0;
                 sol[1][i] = temp1;
                 sol[2][i] = temp2;
             }
             temp0 = vL[0][0];
             temp1 = vL[1][0]/vL[0][0];
             temp2 = (1.4-1.0)*vL[2][0]-(1.4-1.0)/2.0*vL[1][0]*vL[1][0]/vL[0][0];
             vL[0][0] = temp0;
             vL[1][0] = temp1;
             vL[2][0] = temp2;

             temp0 = vLL[0][0];
             temp1 = vLL[1][0]/vLL[0][0];
             temp2 = (1.4-1.0)*vLL[2][0]-(1.4-1.0)/2.0*vLL[1][0]*vLL[1][0]/vLL[0][0];
             vLL[0][0] = temp0;
             vLL[1][0] = temp1;
             vLL[2][0] = temp2;

             temp0 = vR[0][0];
             temp1 = vR[1][0]/vR[0][0];
             temp2 = (1.4-1.0)*vR[2][0]-(1.4-1.0)/2.0*vR[1][0]*vR[1][0]/vR[0][0];
             vR[0][0] = temp0;
             vR[1][0] = temp1;
             vR[2][0] = temp2;

             temp0 = vRR[0][0];
             temp1 = vRR[1][0]/vRR[0][0];
             temp2 = (1.4-1.0)*vRR[2][0]-(1.4-1.0)/2.0*vRR[1][0]*vRR[1][0]/vRR[0][0];
             vRR[0][0] = temp0;
             vRR[1][0] = temp1;
             vRR[2][0] = temp2;


         }
         else // do it on reconstructed values
         {
             for (int i = 0; i < ncells+1; i++)
             {
                 temp0 = rec_l[0][i];
                 temp1 = rec_l[1][i]/rec_l[0][i];
                 temp2 = (1.4-1.0)*rec_l[2][i]-(1.4-1.0)/2.0*rec_l[1][i]*rec_l[1][i]/rec_l[0][i];
                 rec_l[0][i] = temp0;
                 rec_l[1][i] = temp1;
                 rec_l[2][i] = temp2;

                 temp0 = rec_r[0][i];
                 temp1 = rec_r[1][i]/rec_r[0][i];
                 temp2 = (1.4-1.0)*rec_r[2][i]-(1.4-1.0)/2.0*rec_r[1][i]*rec_r[1][i]/rec_r[0][i];
                 rec_r[0][i] = temp0;
                 rec_r[1][i] = temp1;
                 rec_r[2][i] = temp2;
             }
         }
     }
     double minmod(double a, double b)
     {
         double val,fac;
         if (a*b>0)
         {
             if (a>0){fac = 1.0;}
             else{fac = -1.0;}
             val = fac*std::min(abs(a),abs(b));
         }
         else
         {
             val = 0.0;
         }
         return val;
     } 
     void setold()
     {
         for (int i = 0; i < 3; i++)
         {
             for (int j = 0; j < ncells; j++)
             {
                 sol_old[i][j] = sol[i][j];
             }
         }
         //std::cout<<"Old value stored!\n";
     }
public:
    Euler1D(std::string & bctype, int n = 100); //Constructor.
    ~Euler1D();                   //Destructor
    void initialize();
    void set_boundary();
    void reconstruction_0();      //Zero-th order reconstruction
    void reconstruction_TVD();    //Second order TVD reconstruction
    void avg_roe();               //Calculate Roe's average at interfaces
    void cal_roe_flux();
    void writefile(const std::string & filename);
    void timeadvancement(double dt, int oldset);
};

#endif