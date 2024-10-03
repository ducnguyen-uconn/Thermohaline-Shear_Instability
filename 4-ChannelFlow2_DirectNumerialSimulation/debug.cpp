/**
 * This file is a part of channelflow version 2.0 https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include <iomanip>
#include <iostream>
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/utilfuncs.h"

using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        cout << "================================================================\n";
        cout << "This program integrates a plane Couette flow from a random\n";
        cout << "initial condition. Velocity fields are saved at intervals dT=1.0\n";
        cout << "in a data-couette/ directory, in channelflow's binary data file\n";
        cout << "format." << endl << endl;

        // Define gridsize
        const int Nx = 16; // Nx
        const int Ny = 15; // Nz ()
        const int Nz = 16; // Ny

        // Define box size
        const Real Lx = pi; // Lx
        const Real a = 1e-6; // min_Lz = a
        const Real b = 1.0; // max_Lz = b
        const Real Lz = pi; // Ly

        // Define flow parameters
        const Real Reynolds = 400.0; // Reynolds number
        const Real nu = 1.0 / Reynolds; // kinematic viscosity
        const Real dPdx = 0.0; // mean pressure gradient

        // Define integration parameters
        const int n = 32;         // take n steps between printouts
        const Real dt = 1.0 / n;  // integration timestep 
        const Real T = 30.0;      // integrate from t=0 to t=T

        fftw_loadwisdom();
        // Define DNS parameters
        DNSFlags flags;
        flags.baseflow = LaminarBase;

        // set time-stepping algorithm
        // CNFE1 or SBDF1: 1st-order Crank-Nicolson, Foward-Euler or 1st-order Semi-implicit Backward Differentiation Formula
        // CNAB2: 2nd-order Crank-Nicolson, Adams-Bashforth
        // CNRK2: 2nd-order semi-implicit Crank-Nicolson, Runge-Kutta algorithm
        // SMRK2: 2nd-order semi-implicit Runge-Kutta
        // SBDF2, SBDF3, SBDF4: 2nd, 3rd, and 4th-order Semi-implicit Backward Differentiation Formulae
        flags.timestepping = SBDF3; // CNFE1, CNAB2, CNRK2, SMRK2, SBDF1,SBDF2,SBDF3[default], SBDF4
        
        // set initialization timestepping algorithm
        // Some of the time-stepping algorithms listed above (SBDF in particular) require data from N previous time steps
        // This indicates that DNS class instead takes its first N steps with an initialization timestepping algorithm that requires no previous data
        flags.initstepping = CNRK2; // CNFE1, CNRK2[default], SMRK2

        // set form of nonlinear term of Navier-Stokes calculation
        flags.nonlinearity = Rotational; // Rotational, SkewSymmetric[default], Alternating, Linearized
        // Nonlinear terms are calculated with collocation methods
        flags.dealiasing = DealiasXZ; // DealiasXZ=2/3, NoDealiasing, DealiasY=3/2, DealiasXYZ
        // flags.nonlinearity = SkewSymmetric;
        // flags.dealiasing   = NoDealiasing;

        // boundary conditions
        flags.ulowerwall = -1.0; // boundary condition U(z=minLz=a)=-1.0
        flags.uupperwall = 1.0; // boundary condition U(z=maxLz=b)=1.0

        flags.taucorrection = true;
        flags.constraint = PressureGradient;  // enforce constant pressure gradient
        flags.dPdx = dPdx;
        flags.dt = dt;
        flags.nu = nu;

        // Construct data fields: 3d velocity and 1d pressure
        cout << "building velocity and pressure fields..." << flush;
        vector<FlowField> fields = {
            FlowField(Nx, Ny, Nz, 3, Lx, Lz, a, b), // define 3d velocity field
            FlowField(Nx, Ny, Nz, 1, Lx, Lz, a, b)  // define 1d pressure field
        };
        cout << "done" << endl;

        // Define size and smoothness of initial disturbance
        Real spectralDecay = 0.5;
        Real magnitude = 0.1;
        int kxmax = 3;
        int kzmax = 3;
        // Perturb velocity field
        fields[0].addPerturbations(kxmax, kzmax, 1.0, spectralDecay);
        fields[0] *= magnitude / L2Norm(fields[0]);

        // Construct Navier-Stoke integrator, set integration method
        cout << "building DNS..." << flush;
        DNS dns(fields, flags);
        cout << "done" << endl;

        mkdir("data");
        Real cfl = dns.CFL(fields[0]);
        for (Real t = 0; t <= T; t += n * dt) {
            cout << "         t == " << t << endl;
            cout << "       CFL == " << cfl << endl;
            cout << " L2Norm(u) == " << L2Norm(fields[0]) << endl;
            cout << "divNorm(u) == " << divNorm(fields[0]) << endl;
            cout << "      dPdx == " << dns.dPdx() << endl;
            cout << "     Ubulk == " << dns.Ubulk() << endl;

            // Write velocity and modified pressure fields to disk
            fields[0].save("data/u" + i2s(int(t)));
            fields[1].save("data/q" + i2s(int(t)));

            // Take n steps of length dt
            dns.advance(fields, n);
            cout << endl;
        }
    }
    cfMPI_Finalize();
}
