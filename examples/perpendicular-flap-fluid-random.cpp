/** @file perpendicular-flap-fluid-random.cpp

    @brief Fluid participant for the PreCICE example "perpendicular-flap" with random forces.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Li (2023 - ..., TU Delft)
*/

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEUtils.h>
#include <random>
#include <chrono>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Parse command line
    bool plot = false;
    index_t plotmod = 10;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    std::string precice_config;
    real_t forceAmplitude = 100.0;  // Maximum amplitude of random forces
    real_t forceMean = 0.0;         // Mean value of random forces
    
    gsCmdLine cmd("Fluid solver with random forces for PreCICE perpendicular flap.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    cmd.addReal("a", "amplitude", "Maximum amplitude of random forces", forceAmplitude);
    cmd.addReal("", "mean", "Mean value of random forces", forceMean);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Check config file
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    // Generate fluid domain (matching the solid domain)
    gsMultiPatch<> patches;
    
    gsKnotVector<> KV1 (0,1,2,2) ;
    gsKnotVector<> KV2 (0,1,24,2) ;
    
    gsTensorBSplineBasis<2> basis(KV1,KV2);
    gsMatrix<> coefs = basis.anchors();
    coefs.transposeInPlace();
    coefs.col(0).array() *= 0.1;
    coefs.col(0).array() -= 0.05;
    
    patches.addPatch(basis.makeGeometry(give(coefs)));

    // Create bases
    gsMultiBasis<> bases(patches);

    // Set degree
    bases.setDegree( bases.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        bases.uniformRefine();

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";

    // Set the interface for the precice coupling (same as solid)
    std::vector<patchSide> couplingInterfaces(3);
    couplingInterfaces[0] = patchSide(0,boundary::east);
    couplingInterfaces[1] = patchSide(0,boundary::north);
    couplingInterfaces[2] = patchSide(0,boundary::west);

    /*
     * Initialize the preCICE participant
     */
    std::string participantName = "Fluid";
    gsPreCICE<real_t> participant(participantName, precice_config);

    /*
     * Data initialization
     * 
     * - Meshes:
     *   + FluidMesh: This mesh contains the interface points
     * 
     * - Data:
     *   + Stress: Forces/stresses that the fluid applies to the solid
     *   + Displacement: Displacement from the solid (read data)
     */
    
    std::string FluidMesh        = "Fluid-Mesh";
    std::string StressData       = "Stress";
    std::string DisplacementData = "Displacement";

    // Get the quadrature nodes on the coupling interface
    gsOptionList quadOptions = gsAssembler<>::defaultOptions();
    
    // Get the quadrature points
    gsMatrix<> quad_uv = gsQuadrature::getAllNodes(bases.basis(0),quadOptions,couplingInterfaces);
    gsMatrix<> quad_xy = patches.patch(0).eval(quad_uv);
    gsVector<index_t> quad_xyIDs;
    participant.addMesh(FluidMesh,quad_xy,quad_xyIDs);

    // Define precice interface
    real_t precice_dt = participant.initialize();

    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<real_t> dist(forceMean, forceAmplitude);
    
    // Initialize force data
    gsMatrix<> quad_stress(2, quad_xy.cols());
    gsMatrix<> quad_displacement(2, quad_xy.cols());
    quad_stress.setZero();
    quad_displacement.setZero();

    // Define the solution collection for Paraview
    gsParaviewCollection collection("./output/fluid_forces");

    index_t timestep = 0;
    real_t time = 0;
    real_t dt = precice_dt;

    // Store force history for visualization
    std::vector<gsMatrix<>> forceHistory;

    // Time integration loop
    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
        {
            gsInfo << "Checkpoint written at timestep " << timestep << "\n";
        }

        // Read displacement data from solid
        participant.readData(FluidMesh, DisplacementData, quad_xyIDs, quad_displacement);
        
        // Generate random forces at each quadrature point
        gsInfo << "Generating random forces for timestep " << timestep << "...\n";
        
        for (index_t i = 0; i < quad_stress.cols(); ++i)
        {
            // Generate random forces in x and y directions
            // Using normal distribution with specified mean and amplitude
            quad_stress(0, i) = dist(gen);  // Force in x-direction
            quad_stress(1, i) = dist(gen);  // Force in y-direction
            
            // Apply spatial variation based on position
            real_t x = quad_xy(0, i);
            real_t y = quad_xy(1, i);
            
            // Add some spatial correlation to make forces more realistic
            real_t spatialFactor = 1.0 + 0.5 * std::sin(10.0 * x) * std::cos(10.0 * y);
            quad_stress(0, i) *= spatialFactor;
            quad_stress(1, i) *= spatialFactor;
            
            // Add time-varying component
            real_t timeFactor = 1.0 + 0.3 * std::sin(2.0 * M_PI * time);
            quad_stress(0, i) *= timeFactor;
            quad_stress(1, i) *= timeFactor;
        }
        
        // Store force data for history
        forceHistory.push_back(quad_stress);
        
        // Write forces to interface
        participant.writeData(FluidMesh, StressData, quad_xyIDs, quad_stress);
        
        // Log some statistics
        real_t maxForce = quad_stress.cwiseAbs().maxCoeff();
        real_t meanForceX = quad_stress.row(0).mean();
        real_t meanForceY = quad_stress.row(1).mean();
        gsInfo << "Force statistics - Max: " << maxForce 
               << ", Mean X: " << meanForceX 
               << ", Mean Y: " << meanForceY << "\n";

        // Do the coupling
        precice_dt = participant.advance(dt);
        
        if (participant.requiresReadingCheckpoint())
        {
            gsInfo << "Reading checkpoint at timestep " << timestep << "\n";
            // In a real fluid solver, we would restore the fluid state here
        }
        else
        {
            // Advance time
            time += dt;
            timestep++;
            
            // Visualization output
            if (plot && timestep % plotmod == 0)
            {
                // Create a field showing the forces
                gsMultiPatch<> forceField;
                gsMatrix<> forceCoefs(quad_stress.rows(), bases.basis(0).size());
                
                // Simple projection of forces to basis functions
                // In practice, you might want a more sophisticated projection
                for (index_t i = 0; i < bases.basis(0).size(); ++i)
                {
                    forceCoefs(0, i) = quad_stress(0, i % quad_stress.cols());
                    forceCoefs(1, i) = quad_stress(1, i % quad_stress.cols());
                }
                
                forceField.addPatch(bases.basis(0).makeGeometry(give(forceCoefs)));
                
                gsField<> field(patches, forceField);
                std::string fileName = "./output/fluid_forces_" + util::to_string(timestep);
                gsWriteParaview<>(field, fileName, 500);
                fileName = "fluid_forces_" + util::to_string(timestep) + "0";
                collection.addTimestep(fileName, time, ".vts");
            }
        }
        
        // Adjust timestep if needed
        dt = std::min(dt, precice_dt);
    }

    if (plot)
    {
        collection.save();
    }

    // Output some final statistics
    gsInfo << "\nSimulation completed after " << timestep << " timesteps\n";
    gsInfo << "Total simulation time: " << time << "\n";

    return EXIT_SUCCESS;
}