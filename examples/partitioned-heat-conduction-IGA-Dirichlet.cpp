/** @file heat-equation-coupling.cpp

    @brief Heat equation participant for a double coupled heat equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (University of Pavia),  J.Li (TU Delft, 2023-...)
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
    real_t time  = 0;
    real_t k_temp = 1;

    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    // cmd.addInt("l","loadCase", "Load case: 0=constant load, 1='spring' load", loadCase);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    //cmd.addSwitch("readTime", "Get the read time", get_readTime);
    //cmd.addSwitch("writeTime", "Get the write time", get_writeTime);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
    //! [Read input file]
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    /*
     * Initialize the preCICE participant
     *
     * 
     * Dirichlet side:
     * Receive: temperature change as Dirichlet boundary condition
     * Write: flux vector q(u) to the other side as Neumann condition
     *
     */
    std::string participantName = "Dirichlet";
    gsPreCICE<real_t> participant(participantName, precice_config);



    //! [Read input file]

     /*
     * Data initialization
     *
     * This participant receives mesh information (knot vector, control points) from the Neumann,
     * and it creates its own mesh based on the information. 
     * And writes heat flux, reads temperature from the Neumann participant.
     * The follow meshes and data are made available:
     *
     * - Meshes:
     *   + KnotMesh             This mesh contains the knots as mesh vertices
     *   + ControlPointMesh:    This mesh contains the control points as mesh vertices
     *   + FluxMesh:            This mesh contains the integration points as mesh vertices
     *
     * - Data:
     *   + ControlPointData:    This data is defined on the ControlPointMesh and stores the temperature of the control points
     *   + FluzData:            This data is defined on the ForceMesh and stores pressure/forces
     */

    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.0));
    for (int r =0; r < numRefine; ++r)
        patches.uniformRefine();

    gsMultiBasis<> basesDirichlet(patches);
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

    // Setup bounding box onto the force mesh
    gsMatrix<> bbox(2,2);
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300; // Y dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());
    participant.setMeshAccessRegion(GeometryControlPointMesh, bbox);


    // ----------------------------------------------------------------------------------------------

    std::vector<patchSide> couplingInterface(1);
    couplingInterface[0] = patchSide(0,boundary::north);
    std::vector<gsGeometry<>::uPtr> boundaries(couplingInterface.size());




    gsMultiBasis<> fluxBases;
    std::vector<gsGeometry<>::uPtr> couplingBoundaries(couplingInterface.size());

    couplingBoundaries[0] = patches.patch(0).boundary(couplingInterface[0].side());

    // for(index_t i= 0; i < couplingInterface.size();++i)
    // {
    //     couplingBoundaries[i] = patches.patch(0).boundary(couplingInterface[i].side()); // Add boundary coefficients to a vector
    //     auto& basis_temp = couplingBoundaries[i]->basis(); // Assuming this returns a pointer to gsBasis
    //     fluxBases.addBasis(&basis_temp);
    // }
    
    // Add flux knot mesh
    gsVector<> fluxKnotMatrix = knotsToVector(couplingBoundaries[0]->basis());
    gsDebugVar(fluxKnotMatrix);  
    participant.addMesh(FluxKnotMesh, fluxKnotMatrix.transpose());


    // Add flux control points mesh
    gsVector<index_t> fluxControlPointsIDs;
    gsMatrix<> fluxControlPoints = couplingBoundaries[0]->coefs();
    gsDebugVar(fluxControlPoints);
    participant.addMesh(FluxControlPointMesh,fluxControlPoints.transpose(), fluxControlPointsIDs);


    real_t precice_dt = participant.initialize();

    //Get the temperature mesh from direct-access="true"direct-access="true"the API
    gsVector<index_t> geometryKnotIDs;
    gsMatrix<> geometryKnots;
    participant.getMeshVertexIDsAndCoordinates(GeometryKnotMesh,geometryKnotIDs,geometryKnots);

    gsDebugVar(geometryKnots);

    //Get the temperature mesh from direct-access="true"direct-access="true"the API
    gsVector<index_t> geometryControlPointIDs;
    gsMatrix<> geometryControlPoint;
    participant.getMeshVertexIDsAndCoordinates(GeometryControlPointMesh,geometryControlPointIDs,geometryControlPoint);

    gsMatrix<> tempData(geometryControlPoint.rows(),geometryControlPoint.cols());
    tempData.setZero();

    gsMultiPatch<> tempMesh;
    gsBasis<> * basis = knotMatrixToBasis<real_t>(geometryKnots.row(0)).get();
    tempMesh.addPatch(give(basis->makeGeometry(tempData.row(0).transpose())));
    gsDebugVar(tempMesh);


    // Set external heat-flux to zero
    gsConstantFunction<> f(beta-2-2*alpha,2);
    gsFunctionExpr<> u_ex("1+x^2+" + std::to_string(alpha) + "*y^2 + " + std::to_string(beta) + "*" + std::to_string(time),2);
    
    gsBoundaryConditions<> bcInfo;

    // gsPreCICEFunction<real_t> g_CD(&interface,meshName,(side==0 ? tempName : fluxName),patches,1);
    gsFunction<> * g_C = &u_ex;
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, g_C, 0, false, 0); //For initialization, will be changed later
    bcInfo.addCondition(0, boundary::east,  condition_type::dirichlet, &u_ex);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &u_ex);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &u_ex);

    bcInfo.setGeoMap(patches);

    // ----------------------------------------------------------------------------------------------
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    //-----------------------Get discretization space---------------
    gsExprAssembler<> A(1,1);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    A.setIntegrationElements(basesDirichlet);

    gsExprEvaluator<> ev(A);

    geometryMap G = A.getMap(patches);

    // Set the discretization space
    space u = A.getSpace(basesDirichlet); // Use the Dirichlet basis function to discretize solution vector u
    // u_h(x) = \sum_{i=1}^{N}u_i\varphi_i(x), where \varphi_i(x) is the basis function and u_i is the corresponding node values 

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

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    // $u_{\text{ini}}(x) = \int_{\Omega} \varphi_i u_{\text{ex}}(x) , d\Omega$
    auto uex = A.getCoeff(u_ex, G);
    // RHS of the projection
    u.setup(bcInfo, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble( u * u.tr() * meas(G), u * uex * meas(G) );
    solver.compute(A.matrix());
    solVector_ini = solVector = solver.solve(A.rhs());

    gsMatrix<> fluxData(2,fluxControlPoints.transpose().cols()), tmp2;
    // Write initial data
    if (participant.requiresWritingCheckpoint())
    {
        for (index_t k=0; k!=fluxControlPoints.rows(); k++)
        {
            gsWarn<<"Write the flux here!!!\n";
            tmp2 = ev.eval( - jac(u_sol) * nv(G).normalized(),fluxControlPoints.transpose().col(k));
            gsInfo << "Got here\n";
            fluxData(0,k) = tmp2.at(0);
        }
        gsDebugVar(fluxData);
        // gsDebugVar(result);
        participant.writeData(FluxControlPointMesh, FluxControlPointData, fluxControlPointsIDs, fluxData);

    }
    participant.readData(GeometryControlPointMesh,TemperatureData,geometryControlPointIDs,tempData);
    tempMesh.patch(0).coefs()=tempData.row(0).transpose(); 
    gsDebugVar(tempMesh.patch(0).coefs());
    g_C = &tempMesh.patch(0); //Update the boundary condition for the east coupled boundary
    A.initSystem();

    // Assemble stiffness matrix $K = \int_{\Omega} \kappa \nabla \varphi_i \cdot \nabla \varphi_j , d\Omega$    
    A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G), u * uex * meas(G) );
    gsSparseMatrix<> K = A.matrix();
    // Discretize time, the system becomes: $(M + \Delta t K) u^{n+1} = M u^n + \Delta t f$

    gsVector<> F = dt*A.rhs() + M*solVector;
    gsVector<> F0 = F;
    gsVector<> F_checkpoint = F;
    gsMatrix<> solVector_checkpoint = solVector;

    gsParaviewCollection collection("solution_dirichlet");
    collection.options().setSwitch("plotElements",true);
    gsParaviewCollection exact_collection("exact_solution_dirichlet");
  

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;
    real_t t_checkpoint = 0;
    if (plot)
    {
        std::string fileName = "solution_dirichlet" + util::to_string(timestep);
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt("plot.npts", 1000);
        ev.writeParaview( u_sol   , G, fileName);
        for (size_t p=0; p!=patches.nPatches(); p++)
        {
          fileName = "solution_dirichlet" + util::to_string(timestep) + std::to_string(p);
          collection.addTimestep(fileName,time,".vts");
        }
        // fileName = "exact_solution_dirichlet" + util::to_string(timestep);
        // ev.writeParaview( uex   , G, fileName);
        // for (size_t p=0; p!=patches.nPatches(); p++)
        // {
        //   fileName = "exact_solution_dirichlet" + util::to_string(timestep) + std::to_string(p);
        //   exact_collection.addTimestep(fileName,time,".vts");
        // }

        // ev.writeParaview( u_sol   , G, "initial_solution_dirichlet");

    }
    while (participant.isCouplingOngoing())
    {

        u_ex = gsFunctionExpr<>("1+x^2+" + std::to_string(alpha) + "*y^2 + " + std::to_string(beta) + "*" + std::to_string(time),2);
        u.setup(bcInfo, dirichlet::l2Projection, 0); // NOTE:
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
        // write heat fluxes to interface
        // gsMatrix<> result(2,fluxControlPoints.transpose().cols()), tmp;

        for (index_t k=0; k!=fluxControlPoints.transpose().cols(); k++)
        {
            gsWarn<<"Write the flux here!!!\n";
            //Calculate the heat flux $q=-k\nabla u$, here we use jac(u_sol) to calculate $\nabla u$
            // Use normal vector to calculate the flux in the normal direction
            tmp2 = ev.eval(-jac(u_sol) * nv(G).normalized(),fluxControlPoints.transpose().col(k));
            gsInfo << "Got here\n";
            fluxData(0,k) = tmp2.at(0);
        }
        participant.writeData(FluxControlPointMesh, FluxControlPointData, fluxControlPointsIDs, fluxData);
        participant.readData(GeometryControlPointMesh,TemperatureData,geometryControlPointIDs,tempData);
        gsDebugVar(tempData);
        tempMesh.patch(0).coefs()=tempData.row(0).transpose(); 
        gsDebugVar(tempData.row(0).transpose());
        gsDebugVar(tempMesh.patch(0).coefs());
        g_C = &tempMesh.patch(0);
        // g_C = &tempMesh.patch(0); //Update the boundary condition for the east coupled boundary

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
                std::string fileName = "solution_dirichlet" + util::to_string(timestep);
                ev.options().setSwitch("plot.elements", true);
                ev.options().setInt("plot.npts", 1000);
                ev.writeParaview( u_sol   , G, fileName);
                for (size_t p=0; p!=patches.nPatches(); p++)
                {
                  fileName = "solution_dirichlet" + util::to_string(timestep) + std::to_string(p);
                  collection.addTimestep(fileName,time,".vts");
                }

            }
        }
    }

    if (plot)
    {
        collection.save();
        exact_collection.save();
    }
    return  EXIT_SUCCESS;
}
