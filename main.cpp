/*
One dimensional Wavelet Adapted Mesh Refinement Code
Romit Maulik - Oklahoma State University - Computational Fluid Dynamics Laboratory
email: romit.maulik@okstate.edu

Systems to solve:
Burgers - Baeza, A., Martínez-Gavara, A., & Mulet, P. (2012). Adaptation based on
interpolation errors for high order mesh refinement methods applied to conservation laws.
Applied Numerical Mathematics, 62(4), 278-296.

Shallow Water Equations - Ambrosi, D. "Approximation of shallow water equations by Roe's Riemann solver."
International journal for numerical methods in fluids 20.2 (1995): 157-168.

Sod Shock Tube - http://www.csun.edu/~jb715473/examples/euler1d.htm

Brio-Wu Shock Tube - http://www.csun.edu/~jb715473/examples/mhd1d.htm

TRSSE Problem - A high-order cell-centered Lagrangian scheme for one-dimensional
elastic–plastic problems - Cheng, Toro, Computers and Fluids 2015 (122) 136-152

Nonlinear Elasticity - Titarev V, Romanski E, Int. J. Num. Methods. Engg. 2008, 73(7): 897-926

AMR Routine:
Roussel, O., Schneider, K., Tsigulin, A., & Bockhorn, H. (2003).
A conservative fully adaptive multiresolution algorithm for parabolic PDEs.
Journal of Computational Physics, 188(2), 493-523.

Finite Volume Central Flux Computations:
Hyman, James M., Robert J. Knapp, and James C. Scovel.
"High order finite volume approximations of differential operators on nonuniform grids."
Physica D: Nonlinear Phenomena 60.1 (1992): 112-138.

Shock Capturing options:

MUSCL - MinMod/MC
San, Omer, and Kursat Kara.
"Numerical assessments of high-order accurate shock capturing schemes:
Kelvin–Helmholtz type vortical structures in high-resolutions."
Computers & Fluids 89 (2014): 254-276.

WENO3/WENO5
San, Omer, and Kursat Kara. "Evaluation of Riemann flux solvers for WENO reconstruction schemes: Kelvin–Helmholtz instability."
Computers & Fluids 117 (2015): 24-41.

WENO6 - Self Derived

Riemann Solver:
Rusanov -- Maulik, Romit, and Omer San.
"Resolution and Energy Dissipation Characteristics of Implicit LES and Explicit Filtering Models for Compressible Turbulence."
Fluids 2.2 (2017): 14.

FORCE --
San, Omer, and Kursat Kara. "Evaluation of Riemann flux solvers for WENO reconstruction schemes: Kelvin–Helmholtz instability."
Computers & Fluids 117 (2015): 24-41.

*/

#include "Control_Parameters.h"
#include "Global.h"
#include "Function_Def.h"
#include "nle_routines.h"
#include "Init_Domain.h"
#include "Datastruct.h"
#include "Wavelet.h"
#include "Output.h"
#include "Face_Fluxes.h"
#include "Force_flux.h"
#include "Rusanov_Solver.h"
#include "Source_terms.h"
#include "Time_Int.h"
#include "Static_Eval.h"

int main()
{
    unordered_map<int,Cell*>Cellvect;
//    analyse_static(Cellvect);


    create_domain(Cellvect);

    cout<<"Enter number of timestep outputs "<<endl;
    cin>>snap_steps;
    evolution(Cellvect);


    return 0;
}
