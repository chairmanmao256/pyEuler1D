//User interface to Solve 1-D Euler equations
#include <iostream>
#include <string>
#include "Euler1D.h"


int main()
{
    double CFL, dtt, dt, dx, t;
    int step, N, bc;
    std::string bctype;

    std::cout<<"***********************************************************************\n";
    std::cout<<"***********************Pysical problem setup***************************\n";
    std::cout<<"***********************************************************************\n";
    std::cout<<"\n";
    std::cout<<"Specify the boundary condtion. Current boundary condtions available:\n";
    std::cout<<"\n";
    std::cout<<"----------------------------------------\n";
    std::cout<<"| Infinity boundary.           ------1 |\n";
    std::cout<<"| Periodic boundary.           ------2 |\n";
    std::cout<<"----------------------------------------\n";
    std::cout<<"\n";
    std::cout<<"Enter the label following a boundary condtion to specify it.\n";
    std::cin>>bc;
    while(bc>2||bc<0)
    {
        std::cout<<"Please specify a boundary condtion from the list below:\n";
        std::cout<<"----------------------------------------\n";
        std::cout<<"| Infinity boundary.           ------1 |\n";
        std::cout<<"| Periodic boundary.           ------2 |\n";
        std::cout<<"----------------------------------------\n";
        std::cin>>bc;
    }

    if(bc == 1)
    {
        bctype = "Inf";
    }
    else
    {
        bctype = "Periodic";
    }

    std::cout<<"\n";
    

    std::cout<<"***********************************************************************\n";
    std::cout<<"********************Space time mesh specification**********************\n";
    std::cout<<"***********************************************************************\n";
    std::cout<<"\n";
    std::cout<<"Enter the number of cells in x direction: \n";
    std::cin>>N;
    dx = 1.0/double(N);
    std::cout<<"The cell length is: dx = "<<dx<<"\n";
    std::cout<<"Enter the CFL number. The length of time step is calculated by CFL*dx \n";
    std::cin>>CFL;
    dtt = CFL*dx;
    std::cout<<"The time step is: dt = "<<dtt<<"\n";
    std::cout<<"Enter the final time where the solution process ends: \n";
    std::cin>>t;
    if (t>0.25) {std::cout<<"Warning! The structure of solution may penetrate the boundary! \n";}
    std::cout<<"\n";
    step = int(t/dtt);
    
    std::cout<<"***********************************************************************\n";
    std::cout<<"*******************************Execution*******************************\n";
    std::cout<<"***********************************************************************\n";
    std::cout<<"\n";

    

    Euler1D test(bctype,N);
    test.initialize();
    test.set_boundary();

    for (int i = 0; i < step; i++)
    {
        for (int j = 0; j < 4; j++)
        {   
            test.reconstruction_0();
            test.avg_roe();
            test.reconstruction_TVD();
            test.avg_roe();
            test.cal_roe_flux();
            dt = dtt/double(4-j);
            test.timeadvancement(dt,1-j); //only wen j = 0 store the old value
            test.set_boundary();
        }
    }

    test.writefile("output.plt");

    return 0;
}