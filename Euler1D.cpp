// Finite volume method for 1-D Euler equation
#include <fstream>
#include "Euler1D.h"

// Constructor
Euler1D::Euler1D(std::string & bctype, int n)
{
    double dx;
    boundary_type = bctype;
    ncells = n;
    dx = 1.0/double(ncells); // The cell length

    // 2D arrays construction:
    for (int i = 0;i < 3;i++)
    {
        sol[i] = new double[ncells];
        sol_old[i] = new double[ncells];
        rec_l[i] = new double[ncells+1];
        rec_r[i] = new double[ncells+1];
        flux[i] = new double[ncells+1];
        avg[i] = new double[ncells+1];
        vL[i][1] = 0.0;
        vLL[i][1] = 0.0;
        vR[i][1] = 0.0;
        vRR[i][1] = 0.0;
    }

    // x mesh setup:
    mesh = new double[ncells];
    mesh[0] = 0.5*dx;
    for (int i = 1;i < ncells;i++)
    {
        mesh[i] = mesh[i-1] + dx;
    }
}

// Class Destructor
Euler1D::~Euler1D()
{
    std::cout<<"Clean up Euler1D object!";
    for (int i = 0;i < 3;i++)
    {
        delete [] sol[i];
        delete [] sol_old[i];
        delete [] rec_l[i];
        delete [] rec_r[i];
        delete [] flux[i];
        delete [] avg[i];
    }
    delete [] mesh;
}

// Other methods
void Euler1D::initialize()
{
    std::cout<<"Initializing...\n";
    for (int i = 0;i < ncells;i++)
    {
        if (mesh[i]<=0.6 && mesh[i]>=0.4)
        {
            sol[0][i] = 1.0;
            sol[1][i] = 0.75;
            sol[2][i] = 1.0;
        }
        else if (mesh[i] > 0.6)
        {
            sol[0][i] = 0.125;
            sol[1][i] = 0.0;
            sol[2][i] = 0.1;
        }
        else
        {
            sol[0][i] = 1.2;
            sol[1][i] = 0.5;
            sol[2][i] = 1.2;
        }
    }
    std::cout<<"Initilaztion completed...\n";
}

void Euler1D::set_boundary()
{
    if (boundary_type == "Inf")
    {
        // Left boundary
        vL[0][0] = sol[0][0];
        vL[1][0] = sol[1][0];
        vL[2][0] = sol[2][0];
        vLL[0][0] = sol[0][0];
        vLL[1][0] = sol[1][0];
        vLL[2][0] = sol[2][0];

        // Right boundary
        vR[0][0] = sol[0][ncells-1];
        vR[1][0] = sol[1][ncells-1];
        vR[2][0] = sol[2][ncells-1];
        vRR[0][0] = sol[0][ncells-1];
        vRR[1][0] = sol[1][ncells-1];
        vRR[2][0] = sol[2][ncells-1];
    }
    else if (boundary_type == "Periodic")
    {
        // Left boundary
        vL[0][0] = sol[0][ncells-1];
        vL[1][0] = sol[1][ncells-1];
        vL[2][0] = sol[2][ncells-1];
        vLL[0][0] = sol[0][ncells-2];
        vLL[1][0] = sol[1][ncells-2];
        vLL[2][0] = sol[2][ncells-2];

        // Right boundary
        vR[0][0] = sol[0][0];
        vR[1][0] = sol[1][0];
        vR[2][0] = sol[2][0];
        vRR[0][0] = sol[0][1];
        vRR[1][0] = sol[1][1];
        vRR[2][0] = sol[2][1];
    }
    else
    {
        std::cout<<"Boundary type: "<<boundary_type<<" is not available!";
    }
}

void Euler1D::reconstruction_0()
{
    // boundary faces reconstruction
    rec_l[0][0] = vL[0][0];
    rec_l[1][0] = vL[1][0];
    rec_l[2][0] = vL[2][0];
    rec_r[0][ncells] = vR[0][0];
    rec_r[1][ncells] = vR[1][0];
    rec_r[2][ncells] = vR[2][0];

    // interior faces
    for (int i = 1;i < ncells + 1;i++)
    {
        rec_l[0][i] = sol[0][i-1];
        rec_l[1][i] = sol[1][i-1];
        rec_l[2][i] = sol[2][i-1];

        rec_r[0][i-1] = sol[0][i-1];
        rec_r[1][i-1] = sol[1][i-1];
        rec_r[2][i-1] = sol[2][i-1];
    }
}

void Euler1D::reconstruction_TVD()
{
    //std::cout<<"This is TVD reconstruction. \n";
    double lefteig[3][3], RM[3][3];
    double R[3],RR[3],L[3],LL[3], Dr[3], Dl[3], templ[3], tempr[3]; // conservative variables
    double dx, rho, u, p, a, h, rdiff, ldiff;
    int j;
    dx = 1.0/double(ncells);
    prim2con(0);// convert sol to conservative variables

    for (int i = 0; i < ncells + 1; i++)
    {
        // get the right, right right, left and left left cell value
        if (i==0)
        {
            for(j = 0; j < 3; j++)
            {
                L[j] = vL[j][0];
                LL[j] = vLL[j][0];
                R[j] = sol[j][i];
                RR[j] = sol[j][i+1];
            }
        }
        else if (i==1)
        {
            for(j = 0; j < 3; j++)
            {
                L[j] = sol[j][0];
                LL[j] = vL[j][0];
                R[j] = sol[j][i];
                RR[j] = sol[j][i+1];
            }
        }
        else if (i==ncells - 1)
        {
            for(j = 0; j < 3; j++)
            {
                L[j] = sol[j][i-1];
                LL[j] = sol[j][i-2];
                R[j] = sol[j][i];
                RR[j] = vR[j][0];
            }
        }
        else if (i==ncells)
        {
            for(j = 0; j < 3; j++)
            {
                L[j] = sol[j][i-1];
                LL[j] = sol[j][i-2];
                R[j] = vR[j][0];
                RR[j] = vRR[j][0];
            }
        }
        else
        {
            for(j = 0; j < 3; j++)
            {
                L[j] = sol[j][i-1];
                LL[j] = sol[j][i-2];
                R[j] = sol[j][i];
                RR[j] = sol[j][i+1];
            }
        }

        // averaged values
        rho = avg[0][i];
        u = avg[1][i];
        p = avg[2][i];
        a = sqrt(1.4*p/rho);
        h = 1.4/(1.4-1.0)*p/rho;

        // Calculate the left eigen vector matrix
        lefteig[0][0] = (- pow(u,3) + a*u*u + 2.0*h*u)/(2.0*(- a*u*u + 2.0*a*h));
        lefteig[0][1] = -(- u*u + 2.0*a*u + 2.0*h)/(2.0*(- a*u*u + 2.0*a*h));
        lefteig[0][2] = 1.0/(- u*u + 2.0*h);
        lefteig[1][0] = (2.0*(- u*u + h))/(- u*u + 2.0*h);
        lefteig[1][1] = (2.0*u)/(- u*u + 2.0*h);
        lefteig[1][2] = -2.0/(- u*u + 2.0*h);
        lefteig[2][0] = (pow(u,3) + a*u*u - 2.0*h*u)/(2.0*(- a*u*u + 2.0*a*h));
        lefteig[2][1] = -(u*u + 2.0*a*u - 2.0*h)/(2.0*(- a*u*u + 2.0*a*h));
        lefteig[2][2] = 1.0/(- u*u + 2.0*h);

        /*std::cout<<"rho:  "<<rho<<" a:  "<<a<<" h:  "<<h<<" u:  "<<u<<" p:  "<<p<<"\n";
        std::cout<<"LEFT MATRIX: "<<lefteig[0][0]<<"  "<<lefteig[0][1]<<"  "<<lefteig[0][2]<<"\n";
        std::cout<<"             "<<lefteig[1][0]<<"  "<<lefteig[1][1]<<"  "<<lefteig[1][2]<<"\n";
        std::cout<<"             "<<lefteig[2][0]<<"  "<<lefteig[2][1]<<"  "<<lefteig[2][2]<<"\n";
        std::cout<<"\n";*/

        // Calculate the right matrix
        RM[0][0] = 1.0;
        RM[0][1] = 1.0;
        RM[0][2] = 1.0;
        RM[1][0] = u-a;
        RM[1][1] = u;
        RM[1][2] = u+a;
        RM[2][0] = h-u*a;
        RM[2][1] = 0.5*u*u;
        RM[2][2] = h+u*a;


        // Calculate the slope
        for (j = 0; j < 3; j++)
        {
            rdiff = (RR[0]-R[0])*lefteig[j][0] + (RR[1]-R[1])*lefteig[j][1] + (RR[2]-R[2])*lefteig[j][2];
            ldiff = (R[0]-L[0])*lefteig[j][0] + (R[1]-L[1])*lefteig[j][1] + (R[2]-L[2])*lefteig[j][2];
            Dr[j] = minmod(rdiff,ldiff);

            rdiff = (L[0]-LL[0])*lefteig[j][0] + (L[1]-LL[1])*lefteig[j][1] + (L[2]-LL[2])*lefteig[j][2];
            ldiff = (R[0]-L[0])*lefteig[j][0] + (R[1]-L[1])*lefteig[j][1] + (R[2]-L[2])*lefteig[j][2];
            Dl[j] = minmod(rdiff,ldiff);
        }

        // Reconstruct the interface value
        for (j = 0; j < 3; j++)
        {
            templ[j] = lefteig[j][0]*L[0] + lefteig[j][1]*L[1] + lefteig[j][2]*L[2] + Dl[j]*0.5;
            tempr[j] = lefteig[j][0]*R[0] + lefteig[j][1]*R[1] + lefteig[j][2]*R[2] - Dr[j]*0.5;
        }

        // Multiply R matrix (characteristic variables ---> ordinary variables)
        for (j = 0; j < 3; j++)
        {
            rec_l[j][i] = RM[j][0]*templ[0] + RM[j][1]*templ[1] + RM[j][2]*templ[2];
            rec_r[j][i] = RM[j][0]*tempr[0] + RM[j][1]*tempr[1] + RM[j][2]*tempr[2];
        }

    }
    con2prim(1); // Convert conservative variables to primitive variables for reconstructed values
    con2prim(0); // Convert conservative variables to primitive variables for cell center values
}

void Euler1D::avg_roe()
{
    //std::cout<<"This is avg_roe! \n";
    double HL,HR,Hav,a2;
    for (int i = 0; i < ncells+1; i++)
    {
        avg[0][i] = sqrt(rec_l[0][i]*rec_r[0][i]);
        avg[1][i] = (sqrt(rec_l[0][i])*rec_l[1][i] + sqrt(rec_r[0][i])*rec_r[1][i])/(sqrt(rec_l[0][i])+sqrt(rec_r[0][i]));
        HL = 1.4/(1.4-1.0)*rec_l[2][i]/rec_l[0][i];
        HR = 1.4/(1.4-1.0)*rec_r[2][i]/rec_r[0][i];
        Hav = (HR*sqrt(rec_r[0][i]) + HL*sqrt(rec_l[0][i]))/(sqrt(rec_l[0][i])+sqrt(rec_r[0][i]));
        a2 = (1.4-1.0)*(Hav - 0.5*avg[1][i]*avg[1][i]);
        avg[2][i] = avg[0][i]*a2/1.4;
    }
}

void Euler1D::cal_roe_flux()
{
    //std::cout<<"This is cal_roe_flux! \n";
    double del[3], f1[3], f2[3], f3[3], fr[3], fl[3],F[3],k[3]; // Vectors. Maybe replaced by objects in vecmath
    double rho, u, p, a, H; // Store the averaged variables
    double* lam;
    
    for (int i = 0; i < ncells+1; i++)
    {
        del[0] = rec_r[0][i]-rec_l[0][i];
        del[1] = rec_r[1][i]-rec_l[1][i];
        del[2] = rec_r[2][i]-rec_l[2][i];

        rho = avg[0][i];
        u = avg[1][i];
        p = avg[2][i];
        a = sqrt(1.4*p/rho);
        H = 1.4/(1.4-1.0)*p/rho;

        lam = get_eig(u,a);

        k[0] = lam[0]*(del[2]-rho*a*del[1])/(2.0*a*a);
        k[1] = lam[1]*(a*a*del[0]-del[2])/(a*a);
        k[2] = lam[2]*(del[2]+rho*a*del[1])/(2*a*a);

        delete [] lam;

        f1[0] = 1.0;
        f1[1] = u - a;
        f1[2] = H - u*a;

        f2[0] = 1.0;
        f2[1] = u;
        f2[2] = 0.5*u*u;

        f3[0] = 1.0;
        f3[1] = u + a;
        f3[2] = H + u*a;

        fr[0] = rec_r[0][i]*rec_r[1][i];
        fr[1] = rec_r[0][i]*pow(rec_r[1][i],2) + rec_r[2][i];
        fr[2] = rec_r[1][i]*(0.5*rec_r[0][i]*pow(rec_r[1][i],2) + 1.4*rec_r[2][i]/(1.4-1.0));

        fl[0] = rec_l[0][i]*rec_l[1][i];
        fl[1] = rec_l[0][i]*pow(rec_l[1][i],2) + rec_l[2][i];
        fl[2] = rec_l[1][i]*(0.5*rec_l[0][i]*pow(rec_l[1][i],2) + 1.4*rec_l[2][i]/(1.4-1.0));

        flux[0][i] = 0.5*(fr[0]+fl[0])-0.5*(k[0]*f1[0] + k[1]*f2[0] + k[2]*f3[0]);
        flux[1][i] = 0.5*(fr[1]+fl[1])-0.5*(k[0]*f1[1] + k[1]*f2[1] + k[2]*f3[1]);
        flux[2][i] = 0.5*(fr[2]+fl[2])-0.5*(k[0]*f1[2] + k[1]*f2[2] + k[2]*f3[2]);
    }
}

void Euler1D::writefile(const std::string & filename)
{
    std::cout<<"Opening the file...\n";
    std::ofstream outfile;
    outfile.open(filename);
    //outfile<<"variables=x,rho,u,p\n";
    //outfile<<" zone i="<<ncells<<",f=point\n";
    for (int i = 0; i < ncells; i++)
    {
        outfile<<mesh[i]<<"  "<<sol[0][i]<<"  "<<sol[1][i]<<"  "<<sol[2][i]<<"\n";
        //std::cout<<mesh[i]<<"  "<<sol[0][i]<<"  "<<sol[1][i]<<"  "<<sol[2][i]<<"\n";
    }
    outfile.close();
}

void Euler1D::timeadvancement(double dt, int oldset)
{
    double dx;
    dx = 1.0/ncells;

    prim2con(0);    //Convert primitive to conservative, for cell center.
    if (oldset > 0) //Store the current sol. Invoke at the first Runge-Kutta step.
    {
        setold();   //Store the old-conservative value!!
    }
    
    for(int i = 0; i < ncells; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sol[j][i] = sol_old[j][i] - dt/dx*(flux[j][i+1]-flux[j][i]);
        }
    }
    con2prim(0);   //Convert conservative to primitive, for cell center, after time advancement
}
