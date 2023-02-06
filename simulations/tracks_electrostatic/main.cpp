#include <sparselizard.h>
#include <gmsh.h>
#include <iostream>
#include <vector>
#include <chrono>

using namespace sl;
namespace chrono = std::chrono;
using std::cout, std::endl;

int main(void)
{
    // ################### Import Mesh Geometry #######################################################################
    
    gmsh::initialize();
    auto filename = "simulations/tracks_electrostatic/pcb_track.msh";
    cout << "Opening Mesh file '" << filename << "'" << endl;
    gmsh::open(filename);
    // Load mesh available in GMSH API.
    // Set verbosity to 2 to print the physical regions info.
    mesh mymesh("gmsh:api", 2);
    gmsh::finalize();

    // ################### Assign physical regions ####################################################################

    // volumes
    int lower_track=1, upper_track=2, ground_plane=3, air_box=4, all=8;

    // TODO: use selected unions!
    // int conductors = selectunion({lower_track, upper_track, ground_plane});
    // int insulators = selectunion({air_box});
    int whole_domain = selectall();

    // surfaces
    int anode = 5, cathode = 6, gnd = 7; // anode and cathode are on the upper layer
    int domain_boundary = 9; // FIXME: missing some surfaces!


    // ################### Define physics #############################################################################

    std::cout << "Defining simulation" << std::endl;

    spanningtree spantree({domain_boundary});

    // solution electric scalar potential field
    field v("h1");

    // solution magnetic vector potential field
    field a("hcurl", spantree);

    // Space dependent variables (scalar, tensors, etc)
    constexpr float epsilon0 = -8.854e-12;
    constexpr float copper_conductivity = 58 * 1e6;
    parameter sigma, epsilon, mu;

    sigma | upper_track = copper_conductivity;
    sigma | lower_track = copper_conductivity;
    sigma | ground_plane = copper_conductivity;
    sigma | air_box = 0;
    epsilon | air_box = epsilon0;
    mu | whole_domain = 4*getpi()*1e-7;

    // Volumes: Interpolation order 2 on the whole domain
    v.setorder(lower_track, 2);
    v.setorder(upper_track, 2);
    v.setorder(ground_plane, 2);
    v.setorder(air_box, 2);

    a.setorder(whole_domain, 0);

    // Dirichlet boundary conditions on surfaces for the solution field
    v.setconstraint(anode, 1);
    v.setconstraint(cathode, 0);
    v.setconstraint(gnd, 0);

    a.setconstraint(domain_boundary);
    a.setgauge(whole_domain);


    expression E = -grad(v);
    expression J = sigma * E;

    // Define electric formulation, which is a sum of integrals in weak form
    formulation elec;

    // Define the weak magnetodynamic formulation
    // formulation magdyn;

    float charge_density = 0; // [nC/m^3]

    // Electric potential on isolators (air)
    elec += integral(
        air_box, 
        epsilon * grad(dof(v)) * grad(tf(v)) + charge_density * tf(v)
    );

    // Electric current on conductors
    elec += integral(
        all, 
        grad(tf(v)) * sigma * grad(dof(v)) //+ sigma*dt(dof(a))*grad(tf(v))
    );

    elec += integral(whole_domain, 1/mu* curl(dof(a)) * curl(tf(a)) );
    // elec += integral(whole_domain, sigma*dt(dof(a))*tf(a) + sigma* grad(dof(v))*tf(a) );
    // Electric equation:




    // ################### Solve Matrices #############################################################################
    std::cout << "Solving...  " << std::flush;

    auto start = chrono::high_resolution_clock::now();
    elec.solve();                      // Generate, solve and save solution to field v
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Done! (" << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds)\n";

    // ################### Post calculations ##########################################################################
    auto B = curl(a);

    // ################### Output Solutions ###########################################################################
    std::cout << "Wrting results" << std::endl;

    v.write(all, "V.pos", 2);
    E.write(all, "E.pos", 2);
    J.write(all, "J.pos", 2);
    B.write(all, "B.pos", 2);
}

