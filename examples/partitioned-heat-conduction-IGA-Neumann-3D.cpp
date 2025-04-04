/** @file partitioned-heat-conduction-IGA-Neumann-3D.cpp

    @brief 3D Heat equation participant for a double coupled heat equation (Neumann side)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Li (TU Delft, 2023-...)
*/



#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsPreCICEUtils.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = true;
    index_t plotmod = 1;
    index_t numRefine  = 2;
    index_t numElevate = 0;
    short_t side = 0;
    real_t alpha = 3;
    real_t beta  = 1.2;
    real_t gamma = 2.5;  // New coefficient for 3D exact solution
    real_t time  = 0;
    real_t k_temp = 1;

    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("3D Flow over heated plate for PreCICE (Neumann participant).");
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    cmd.addInt("r","numRefine", "Number of uniform refinements", numRefine);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
    //! [Read input file]
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    /*
     * Initialize the preCICE participant
     *
     * 
     * Neumann side:
     * Receive: flux vector q(u) from the other side as Neumann condition
     * Write: temperature as Dirichlet boundary condition to the other side
     *
     */
    std::string participantName = "Neumann";
    gsPreCICE<real_t> participant(participantName, precice_config);


    //! [Read input file]

    /*
     * Data initialization
     *
     * This participant provides mesh information (knot vector, control points) to the Dirichlet,
     * and it creates its own mesh for flux data.
     * And reads heat flux, writes temperature to the Dirichlet participant.
     * The follow meshes and data are made available:
     *
     * - Meshes:
     *   + GeometryKnotMesh         This mesh contains the knots as mesh vertices
     *   + GeometryControlPointMesh This mesh contains the control points as mesh vertices
     *   + FluxMesh                 This mesh contains the integration points as mesh vertices
     *
     * - Data:
     *   + TemperatureData:         This data is defined on the GeometryControlPointMesh and stores the temperature of the control points
     *   + FluxData:                This data is defined on the FluxMesh and stores heat flux
     */

    gsMultiPatch<> patches;
    // Create a 3D cube, from (0,0,0) to (1,1,1)
    patches.addPatch(gsNurbsCreator<>::BSplineCube());


    
    // Translate 1 unit along the x-axis
    gsMatrix<> translation(3,1);
    translation << 1.0, 0.0, 0.0;
    patches.patch(0).translate(translation);

    gsDebugVar(patches.patch(0).coefs());
    
    for (int r =0; r < numRefine; ++r)
        patches.uniformRefine();

    gsMultiBasis<> basesNeumann(patches);
    // ----------------------------------------------------------------------------------------------

    /// Create from patches and boundary/interface information

    //Mesh provided by Neumann
    std::string GeometryKnotMesh                = "Geometry-Knot-Mesh";
    std::string GeometryControlPointMesh        = "Geometry-Control-Point-Mesh";

    // Mesh provided by Dirichlet
    std::string FluxKnotMesh                   = "Flux-Knot-Mesh";
    std::string FluxControlPointMesh           = "Flux-Control-Point-Mesh";

    // Data provided by Neumann
    std::string TemperatureData                 = "Temperature-Data";

    // Data provided by Dirichlet
    std::string FluxControlPointData           = "Flux-Control-Point-Data";

    // Setup bounding box onto the force mesh - Neumann is the receiver
    gsMatrix<> bbox(3,2);  // 3D bounding box, as the coupling interface is 3D
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300, // Y dimension limits
            -1e300, 1e300; // Z dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());
    // Neumann participant needs to set access region for Flux mesh, as it will receive these meshes from Dirichlet
    participant.setMeshAccessRegion(FluxControlPointMesh, bbox);

    // ----------------------------------------------------------------------------------------------

    std::vector<patchSide> couplingInterface(1);
    // Coupling on the west side (x=1 plane), corresponding to the east side of Dirichlet
    couplingInterface[0] = patchSide(0,boundary::west);
    std::vector<gsGeometry<>::uPtr> boundaries(couplingInterface.size());

    std::vector<gsGeometry<>::uPtr> couplingBoundaries(couplingInterface.size());

    couplingBoundaries[0] = patches.patch(0).boundary(couplingInterface[0].side());
    
    // Add geometry knot mesh - for 3D problems, the boundary is now 2D
    gsDebugVar(couplingBoundaries[0]->basis().basis(0));
    gsMatrix<> geomKnotMatrix = knotsToMatrix(couplingBoundaries[0]->basis().basis(0));
    gsDebugVar(geomKnotMatrix.dim());  
    
    // Neumann provides Geometry-Knot-Mesh
    participant.addMesh(GeometryKnotMesh, geomKnotMatrix);

    // Add temperature control points mesh
    gsVector<index_t> temperatureControlPointsIDs;
    gsMatrix<> temperatureControlPoints = couplingBoundaries[0]->coefs();
    
    // Transpose the control points matrix
    temperatureControlPoints.transposeInPlace();
    gsDebugVar(temperatureControlPoints.dim());
    
    // Add to mesh using 3D points - this is the mesh provided by the Neumann participant
    participant.addMesh(GeometryControlPointMesh, temperatureControlPoints, temperatureControlPointsIDs);

    // Initialize PreCICE - Neumann participant creates two meshes: GeometryKnotMesh and GeometryControlPointMesh
    real_t precice_dt = participant.initialize();

    // Get mesh provided by Dirichlet participant from the API - do not add or modify these meshes
    gsVector<index_t> fluxKnotIDs;
    gsMatrix<> fluxKnots;
    participant.getMeshVertexIDsAndCoordinates(FluxKnotMesh, fluxKnotIDs, fluxKnots);
    gsDebugVar(fluxKnots);

    // Get control point mesh provided by Dirichlet participant from the API - do not add or modify these meshes
    gsVector<index_t> fluxControlPointIDs;
    gsMatrix<> fluxControlPoint;
    participant.getMeshVertexIDsAndCoordinates(FluxControlPointMesh, fluxControlPointIDs, fluxControlPoint);
    gsDebugVar(fluxControlPoint);
    // Initialize 2D flux data - fluxControlPoint stores 2D points (y,z)
    gsMatrix<> fluxData(fluxControlPoint.rows(), fluxControlPoint.cols());
    fluxData.setZero();

    gsMultiPatch<> fluxMesh;
    gsBasis<> * basis = knotMatrixToBasis<real_t>(fluxKnots).get();
    fluxControlPoint.transposeInPlace();
    fluxMesh.addPatch(give(basis->makeGeometry(fluxControlPoint)));
    gsWriteParaview(fluxMesh, "fluxMesh");
    gsDebugVar(fluxMesh);

    // Set external heat-flux to zero - Update for 3D
    // Using Gaussian function to create a local heat source at (2,1,1)
    real_t amplitude = 5.0; // Heat source intensity
    real_t sigma = 0.1;     // Heat source distribution width
    std::string sourceExpr = std::to_string(amplitude) + "*exp(-((x-2)^2+(y-1)^2+(z-1)^2)/(" + std::to_string(2*sigma*sigma) + "))";
    gsFunctionExpr<> f(sourceExpr, 3);
    
    // 3D exact solution - for translated coordinate system
    gsFunctionExpr<> u_ex("1+(z-1)^2+" + std::to_string(alpha) + "*y^2+" + std::to_string(gamma) + "*x^2 + " 
                         + std::to_string(beta) + "*" + std::to_string(time), 3);
    
    gsBoundaryConditions<> bcInfo;

    // Apply Neumann condition on the boundary (q·n = 0)
    gsConstantFunction<> zeroflux(0.0,3);
    // Set boundary conditions for 3D domain
    bcInfo.addCondition(0, boundary::east, condition_type::neumann, &fluxMesh.patch(0), 0, false, 0); // Coupling boundary
    bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet, &u_ex);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &u_ex);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &u_ex);
    bcInfo.addCondition(0, boundary::front, condition_type::dirichlet, &u_ex);
    bcInfo.addCondition(0, boundary::back,  condition_type::dirichlet, &u_ex);

    bcInfo.setGeoMap(patches);

    // ----------------------------------------------------------------------------------------------
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    //-----------------------Get discretization space---------------
    gsExprAssembler<> A(1,1);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    A.setIntegrationElements(basesNeumann);

    gsExprEvaluator<> ev(A);

    geometryMap G = A.getMap(patches);

    // Set the discretization space
    space u = A.getSpace(basesNeumann);

    //----------------------RHS of the heat equation \frac{\partial u}{\partial t} = \nabla \cdot (\kappa \nabla u) + f ---------------------
    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Set the solution and initial condition
    gsMatrix<> solVector, solVector_ini;
    solution u_sol = A.getSolution(u, solVector);
    solution u_ini = A.getSolution(u, solVector);

    // Initialize the system and assemble the mass matrix
    // $M = \int_{\Omega} \varphi_i \varphi_j , d\Omega$
    u.setup(bcInfo, dirichlet::homogeneous, 0);
    A.initSystem();
    gsDebugVar(A.numDofs());
    A.assemble( u * u.tr() * meas(G));
    gsSparseMatrix<> M = A.matrix();

    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner
    gsSparseSolver<>::CGDiagonal solver;

    real_t dt = 0.1;

    // Project u_wall as initial condition
    // $u_{\text{ini}}(x) = \int_{\Omega} \varphi_i u_{\text{ex}}(x) , d\Omega$
    auto uex = A.getCoeff(u_ex, G);
    // RHS of the projection
    u.setup(bcInfo, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble( u * u.tr() * meas(G), u * uex * meas(G) );
    solver.compute(A.matrix());
    solVector_ini = solVector = solver.solve(A.rhs());

    // Initialize temperature data - scalar values
    gsMatrix<> temperatureData(3, temperatureControlPoints.cols());
    // Initialize data
    if (participant.requiresWritingCheckpoint())
    {
        // Write temperature data
        for (index_t k=0; k!=temperatureControlPoints.cols(); k++)
        {
            gsMatrix<> temp;
            temp = ev.eval(u_sol, temperatureControlPoints.col(k));
            temperatureData(0,k) = temp.at(0);
            temperatureData(1,k) = temp.at(1);
            temperatureData(2,k) = temp.at(2);
        }
        gsDebugVar(temperatureData);
        participant.writeData(GeometryControlPointMesh, TemperatureData, temperatureControlPointsIDs, temperatureData);
    }

    // Read thermal flux data
    participant.readData(FluxControlPointMesh, FluxControlPointData, fluxControlPointIDs, fluxData);
    gsDebugVar(fluxData);
    
    // Update boundary conditions - using received thermal flux
    fluxMesh.patch(0).coefs() = fluxData;
    
    gsDebugVar(fluxMesh.patch(0).coefs());
    A.initSystem();

    // Assemble stiffness matrix $K = \int_{\Omega} \kappa \nabla \varphi_i \cdot \nabla \varphi_j , d\Omega$    
    A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G), u * uex * meas(G) );
    gsSparseMatrix<> K = A.matrix();
    // Discretize time, the system becomes: $(M + \Delta t K) u^{n+1} = M u^n + \Delta t f$

    gsVector<> F = dt*A.rhs() + M*solVector;
    gsVector<> F0 = F;
    gsVector<> F_checkpoint = F;
    gsMatrix<> solVector_checkpoint = solVector;

    // Create 3D visualization collections
    gsParaviewCollection collection("solution_neumann_3D");
    collection.options().setSwitch("plotElements", true);
    gsParaviewCollection exact_collection("exact_solution_neumann_3D");
    exact_collection.options().setSwitch("plotElements", true);
    gsParaviewCollection error_collection("error_solution_neumann_3D");
    error_collection.options().setSwitch("plotElements", true);
  
    index_t timestep = 0;
    index_t timestep_checkpoint = 0;
    real_t t_checkpoint = 0;
    if (plot)
    {
        std::string fileName = "solution_neumann_3D" + util::to_string(timestep);
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt("plot.npts", 1000);
        ev.writeParaview( u_sol, G, fileName);
        for (size_t p=0; p!=patches.nPatches(); p++)
        {
          fileName = "solution_neumann_3D" + util::to_string(timestep) + std::to_string(p);
          collection.addTimestep(fileName,time,".vts");
        //   collection
        }

        fileName = "exact_solution_neumann_3D" + util::to_string(timestep);
        ev.writeParaview( uex, G, fileName);
        for (size_t p=0; p!=patches.nPatches(); p++)
        {
          fileName = "exact_solution_neumann_3D" + util::to_string(timestep) + std::to_string(p);
          exact_collection.addTimestep(fileName,time,".vts");
        }
        
        auto error = A.getCoeff(u_ex, G) - u_sol;
        fileName = "error_solution_neumann_3D" + util::to_string(timestep);
        ev.writeParaview(error, G, fileName);
        for (size_t p=0; p!=patches.nPatches(); p++)
        {
          fileName = "error_solution_neumann_3D" + util::to_string(timestep) + std::to_string(p);
          error_collection.addTimestep(fileName,time,".vts");
        }
    }
    while (participant.isCouplingOngoing())
    {
        // Update 3D exact solution
        u_ex = gsFunctionExpr<>("1+(z-1)^2+" + std::to_string(alpha) + "*y^2+" + std::to_string(gamma) + "*x^2 + " 
                               + std::to_string(beta) + "*" + std::to_string(time), 3);
        u.setup(bcInfo, dirichlet::l2Projection, 0);
        A.initSystem();
        A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );
        K = A.matrix();
        F = dt*A.rhs() + M*solVector;
        
        // save checkpoint
        if (participant.requiresWritingCheckpoint())
        {
            gsDebugVar("Write checkpoint");
            F_checkpoint = F0;
            t_checkpoint = time;
            timestep_checkpoint = timestep;
            solVector_checkpoint = solVector;
        }

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt,precice_dt);

        // solve gismo timestep
        gsInfo << "Solving timestep " << timestep*dt << "...";
        solVector = solver.compute(M + dt*K).solve(F);
        gsInfo<<"Finished\n";

        // Calculate and write temperature data
        for (index_t k=0; k!=temperatureControlPoints.cols(); k++)
        {
            // Convert control point back to original format for temperature calculation
            gsMatrix<> controlPoint = temperatureControlPoints.col(k);
            gsMatrix<> temp;
            temp = ev.eval(u_sol, controlPoint);
            temperatureData(0,k) = temp.at(0);
            temperatureData(1,k) = temp.at(1);
            temperatureData(2,k) = temp.at(2);
        }
        gsDebugVar(temperatureData);
        participant.writeData(GeometryControlPointMesh, TemperatureData, temperatureControlPointsIDs, temperatureData);
        
        
        // Read thermal flux data
        participant.readData(FluxControlPointMesh, FluxControlPointData, fluxControlPointIDs, fluxData);
        fluxMesh.patch(0).coefs() = fluxData;
        gsDebugVar(fluxData);

        // do the coupling
        precice_dt = participant.advance(dt);

        // advance variables
        time += dt;
        timestep += 1;
        F0 = F;

        if (participant.requiresReadingCheckpoint())
        {
            gsDebugVar("Read checkpoint");
            F0 = F_checkpoint;
            time = t_checkpoint;
            timestep = timestep_checkpoint;
            solVector = solVector_checkpoint;
        }
        else
        {
            if (timestep % plotmod==0 && plot)
            {
                std::string fileName = "solution_neumann_3D" + util::to_string(timestep);
                ev.options().setSwitch("plot.elements", true);
                ev.options().setInt("plot.npts", 1000);
                ev.writeParaview( u_sol, G, fileName);
                for (size_t p=0; p!=patches.nPatches(); p++)
                {
                  fileName = "solution_neumann_3D" + util::to_string(timestep) + std::to_string(p);
                  collection.addTimestep(fileName,time,".vts");
                }

                fileName = "exact_solution_neumann_3D" + util::to_string(timestep);
                ev.writeParaview( uex, G, fileName);
                for (size_t p=0; p!=patches.nPatches(); p++)
                {
                  fileName = "exact_solution_neumann_3D" + util::to_string(timestep) + std::to_string(p);
                  exact_collection.addTimestep(fileName,time,".vts");
                }
                
                auto error = A.getCoeff(u_ex, G) - u_sol;
                fileName = "error_solution_neumann_3D" + util::to_string(timestep);
                ev.writeParaview(error, G, fileName);
                for (size_t p=0; p!=patches.nPatches(); p++)
                {
                  fileName = "error_solution_neumann_3D" + util::to_string(timestep) + std::to_string(p);
                  error_collection.addTimestep(fileName,time,".vts");
                }
            }
        }
    }

    if (plot)
    {
        collection.save();
        exact_collection.save();
        error_collection.save();
    }
    return  EXIT_SUCCESS;
} 