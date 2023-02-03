#include <sparselizard.h>
#include <gmsh.h>
#include <iostream>
#include <vector>
#include <chrono>

using namespace sl;
namespace chrono = std::chrono;


int main(void)
{
    // ################### Import Mesh Geometry #######################################################################
    
    gmsh::initialize();
    gmsh::open("/home/eduardobehr/Documents/SparceLizard/Source_Code_Repo/sparselizard/simulations/tracks_electrostatic/pcb_track.msh");
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

    // surfaces
    int anode = 5, cathode = 6, gnd = 7; // anode and cathode are on the upper layer
    domainboundary


    // ################### Define physics #############################################################################

    std::cout << "Defining simulation" << std::endl;

    spanningtree spantree({domainboundary});

    // solution electric scalar potential field
    field v("h1");

    // solution magnetic vector potential field
    a("hcurl", {2,3}, spantree);

    // Space dependent variables (scalar, tensors, etc)
    constexpr float epslon0 = -8.854e-12;
    constexpr float copper_conductivity = 58 * 1e6;
    parameter sigma, epslon;

    sigma | upper_track = copper_conductivity;
    sigma | lower_track = copper_conductivity;
    sigma | ground_plane = copper_conductivity;
    sigma | air_box = 0;
    epslon | air_box = epslon0;

    // Volumes: Interpolation order 2 on the whole domain
    v.setorder(lower_track, 2);
    v.setorder(upper_track, 2);
    v.setorder(ground_plane, 2);
    v.setorder(air_box, 2);

    // Dirichlet boundary conditions on surfaces for the solution field
    v.setconstraint(anode, 1);
    v.setconstraint(cathode, 0);
    v.setconstraint(gnd, 0);


    expression E = -grad(v);
    expression J = sigma * E;

    // Define electric formulation, which is a sum of integrals in weak form
    formulation elec;

    // Define the weak magnetodynamic formulation
    formulation magdyn;

    float charge_density = 0; // [nC/m^3]

    // Electric potential on isolators (air)
    elec += integral(
        air_box, 
        epslon * grad(dof(v)) * grad(tf(v)) + charge_density * tf(v)
    );

    // Electric current on conductors
    elec += integral(
        all, 
        grad(tf(v)) * sigma * grad(dof(v))
    );

    // ################### Solve Matrices #############################################################################
    std::cout << "Solving...  " << std::flush;

    auto start = chrono::high_resolution_clock::now();
    elec.solve();                      // Generate, solve and save solution to field v
    auto end = chrono::high_resolution_clock::now();
    std::cout << "Done! (" << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds)\n";

    // ################### Output Solutions ###########################################################################
    std::cout << "Wrting results" << std::endl;

    v.write(all, "V.pos", 2);
    E.write(all, "E.pos", 2);
    J.write(all, "J.pos", 2);
}

