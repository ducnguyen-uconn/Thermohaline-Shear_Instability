/**
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include <fstream>
#include <iomanip>
#include <iostream>
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/symmetry.h"
#include "channelflow/utilfuncs.h"

using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        string purpose(
            "Construct a random field with zero divergence and Dirichlet BCs\n"
            "The field is u = sum_{jkl} a_{jkl} T_l(y) exp(2 pi i (jx/Lx + kz/Lz)) where\n"
            "a_{jkl} = (random # in [-1 1]) * smoothness^{|l| + |j| + |k|}, with corrections\n"
            "to assure u has zero divergence and no-slip BCs, and with rescaling\n"
            "u = magnitude/L2Norm(u)");

        ArgList args(argc, argv, purpose);

        const int Nx = args.getint("-Nx", "--Nx", "# x gridpoints");
        const int Ny = args.getint("-Ny", "--Ny", "# y gridpoints");
        const int Nz = args.getint("-Nz", "--Nz", "# z gridpoints");
        const Real alpha = args.getreal("-a", "--alpha", 0, "Lx = 2 pi/alpha");
        const Real gamma = args.getreal("-g", "--gamma", 0, "Lz = 2 pi/gamma");
        const Real lx = (alpha == 0.0) ? args.getreal("-lx", "--lx", 0.0, "Lx = 2 pi lx") : 1 / alpha;
        const Real lz = (gamma == 0.0) ? args.getreal("-lz", "--lz", 0.0, "Lz = 2 pi lz") : 1 / gamma;
        const Real Lx = (lx == 0.0) ? args.getreal("-Lx", "--Lx", "streamwise (x) box length") : 2 * pi * lx;
        const Real Lz = (lz == 0.0) ? args.getreal("-Lz", "--Lz", "spanwise   (z) box length") : 2 * pi * lz;
        const Real ymin = args.getreal("-ymin", "--ymin", -0.5, "lower wall height (y is wallnormal) ");
        const Real ymax = args.getreal("-ymax", "--ymax", +0.5, "upper wall height (y is wallnormal) ");
        const int seed = args.getint("-sd", "--seed", 1, "seed for random number generator");
        const Real smooth = args.getreal("-s", "--smoothness", 0.4, "smoothness of field, 0 < s < 1");
        const Real magn = args.getreal("-m", "--magnitude", 0.2, "magnitude  of field, 0 < m < 1");
        const bool meanfl = args.getflag("-mf", "--meanflow", "perturb the mean");

        const string symmstr = args.getstr("-symms", "--symmetries", "", "file of symmetries to satisfy");

        const string uname = args.getstr(3, "<fieldname>", "output file");
        const string tname = args.getstr(2, "<fieldname>", "output file");
        const string sname = args.getstr(1, "<fieldname>", "output file");
        args.check();
        args.save("./");

        srand48(seed);

        FlowField u(Nx, Ny, Nz, 3, Lx, Lz, ymin, ymax);
        u.addPerturbations(u.kxmaxDealiased(), u.kzmaxDealiased(), 1.0, 1 - smooth, meanfl);

        SymmetryList s;
        if (symmstr.length() > 0) {
            s = SymmetryList(symmstr);
            cout << "Restricting random field to invariant subspace generated by symmetries" << endl;
            cout << s << endl;
        }

        for (int i = 0; i < s.length(); ++i)
            u += s[i](u);

        u *= magn / L2Norm(u);
        u.setPadded(true);
        u.save(uname);

        FlowField temp(Nx, Ny, Nz, 1, Lx, Lz, ymin, ymax);
        temp.addPerturbations(u.kxmaxDealiased(), u.kzmaxDealiased(), 1.0, 1 - smooth, meanfl);

        for (int i = 0; i < s.length(); ++i)
            temp += s[i](temp);

        temp *= magn / L2Norm(temp);
        temp.setPadded(true);
        temp.save(tname);

        FlowField salt(Nx, Ny, Nz, 1, Lx, Lz, ymin, ymax);
        salt.addPerturbations(u.kxmaxDealiased(), u.kzmaxDealiased(), 1.0, 1 - smooth, meanfl);

        for (int i = 0; i < s.length(); ++i)
            salt += s[i](salt);

        salt *= magn / L2Norm(salt);
        salt.setPadded(true);
        salt.save(sname);
    }
    cfMPI_Finalize();
}