/// This is the fluid-structure interaction benchmark FSI2 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006.
/// Author(s): H.M. Verhelst (University of Pavia),  J.Li (TU Delft, 2023-...)
///
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsALE.h>
#include <gsElasticity/gsPartitionedFSI.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>

#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsPreCICEFunction.h>
// #include <gsPreCICE/gsPreCICEVectorFunction.h>

using namespace gismo;


int main(int argc, char* argv[])
{
    gsInfo << "Testing the two-way fluid-structure interaction solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    // beam parameters
    bool plot = true;
    index_t plotmod = 1;
    index_t numUniRef  = 3;
    index_t numElevate = 0;
    short_t side = 0;
    real_t alpha = 3;
    real_t beta  = 1.2;
    real_t time  = 0;
    real_t k_temp = 1;


    std::string precice_config("../precice_config.xml");

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the two-way fluid-structure interaction solver in 2D.");
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    // cmd.addReal("t","time","Time span, sec",timeSpan);
    // cmd.addReal("s","step","Time step",timeStep);
    // cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    // cmd.addInt("v","verbosity","Amount of info printed to the prompt: 0 - none, 1 - crucial, 2 - all",verbosity);
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(1.0,0.0,2.0,1.0));

    // creating bases
    for (index_t i = 0; i < numUniRef; ++i)
        patches.uniformRefine();
    gsMultiBasis<> basisTemperature(patches);

    // Set external heat-flux to zero
    gsConstantFunction<> f(beta-2-2*alpha,2);
    gsFunctionExpr<> u_ex("1+x^2+" + std::to_string(alpha) + "*y^2 + " + std::to_string(beta) + "*" + std::to_string(time),2);


    //=============================================//
        // Define preCICE setup //
    //=============================================//

    std::string participantName = "Neumann";
    gsPreCICE<real_t> participant(participantName, precice_config);

    // Mesh provided by solid
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
    participant.setMeshAccessRegion(FluxControlPointMesh, bbox);

    // Setup geometry control point mesh
    gsVector<index_t> geometryControlPointIDs;
    gsMatrix<> geometryControlPoints = patches.patch(0).coefs().transpose();

    std::vector<patchSide> couplingInterface(1);
    couplingInterface[0] = patchSide(0,boundary::south);
    // std::vector<gsGeometry<>::uPtr> boundaries(couplingInterface.size());

    gsMultiBasis<> temperatureBases;
    std::vector<gsGeometry<>::uPtr> couplingBoundaries(couplingInterface.size());
    couplingBoundaries[0] = patches.patch(0).boundary(couplingInterface[0].side());

    gsMatrix<> temperatureControlPoints = couplingBoundaries[0]->coefs();
    participant.addMesh(GeometryControlPointMesh, temperatureControlPoints.transpose(), geometryControlPointIDs);

    // gsMatrix<> boundaryPoints(temperatureControlPoints.cols(),temperatureControlPoints.rows());
    // boundaryPoints.row(0) = temperatureControlPoints.col(0).transpose();

    // gsDebugVar(boundaryPoints);
    gsMatrix<> boundaryPoints = temperatureControlPoints;
    boundaryPoints.col(0).array() -= 1;
    gsDebugVar(boundaryPoints);

    gsVector<> temperatureKnotMatrix = knotsToVector(couplingBoundaries[0]->basis());
    gsDebugVar(temperatureKnotMatrix);  
    participant.addMesh(GeometryKnotMesh, temperatureKnotMatrix.transpose());



    // Setup the geometry knot mesh

    //participant.initialize();

    // Initialize preCICE (send solid mesh to API)
    real_t precice_dt = participant.initialize();


    //Get the force mesh from direct-access="true"direct-access="true"the API
    gsVector<index_t> fluxKnotIDs;
    gsMatrix<> fluxKnots;
    participant.getMeshVertexIDsAndCoordinates(FluxKnotMesh,fluxKnotIDs,fluxKnots);

    gsDebugVar(fluxKnots);

    // Gives a full tensor product basis

    gsBasis<> * basis = knotMatrixToBasis<real_t>(fluxKnots.row(0)).get();

    gsDebugVar(basis->size());


    gsMatrix<> fluxControlPoints(basis->size(),1);
    fluxControlPoints.setZero();


    // Gives the coefficients of the control points ONLY ON THE BOUNDARIES

    // IMPORTANT: THE control points now are all control points!!!
    gsVector<index_t> fluxControlPointIDs;
    participant.getMeshVertexIDsAndCoordinates(FluxControlPointMesh, fluxControlPointIDs,fluxControlPoints);
    gsDebugVar(fluxControlPoints);




    // // Step 2: Regenerate the geometry
    gsMultiPatch<> fluxMesh; //Geometry object belongs to gsFunctionSet
    fluxMesh.addPatch(give(basis->makeGeometry(fluxControlPoints.transpose())));

    gsDebugVar(fluxMesh.patch(0).coefs());


    // NOTE: forceMesh should have domainDim 2!!
    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//
    // Setup the boundary condition for the Neumann side
    gsBoundaryConditions<> bcInfo;
    gsConstantFunction<> g_D(0,2); // Dirichlet
    bcInfo.addCondition(0, boundary::south,  condition_type::neumann  , &fluxMesh.patch(0));
    bcInfo.addCondition(0, boundary::east,  condition_type::dirichlet  , &u_ex);
    bcInfo.addCondition(0, boundary::west,  condition_type::dirichlet  , &u_ex);
    bcInfo.addCondition(0, boundary::north,  condition_type::dirichlet  , &u_ex);

    bcInfo.setGeoMap(patches);
    gsDebugVar(bcInfo);

    // ----------------------------------------------------------------------------------------------

    // Generate system matrix and load vector
    // gsInfo << "Assembling mass and stiffness...\n";
    // 
    // ----------------------------------------------------------------------------------------------
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    A.setIntegrationElements(basisTemperature);

    gsExprEvaluator<> ev(A);


    // Set the geometry map
    geometryMap G = A.getMap(patches);

        // Set the discretization space
    space u = A.getSpace(basisTemperature);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Set the solution
    gsMatrix<> solVector, solVector_ini;
    solution u_sol = A.getSolution(u, solVector);
    solution u_ini = A.getSolution(u, solVector);

    u.setup(bcInfo, dirichlet::homogeneous, 0); // Add boundary condition
    A.initSystem();
    gsDebugVar(A.numDofs());
    A.assemble( u * u.tr() * meas(G));
    gsSparseMatrix<> M = A.matrix();

    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner
    gsSparseSolver<>::CGDiagonal solver;

    real_t dt = 0.1;

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    auto uex = A.getCoeff(u_ex, G);

    // RHS of the projection
    u.setup(bcInfo, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble( u * u.tr() * meas(G), u * uex * meas(G) );
    solver.compute(A.matrix());
    solVector_ini = solVector = solver.solve(A.rhs());

    gsMatrix<> results(2,boundaryPoints.transpose().cols()), tmp2;


    // Write initial data
    if (participant.requiresWritingCheckpoint())
    {
        
        for (index_t k=0; k!=boundaryPoints.transpose().cols(); k++)
        {
                gsInfo << "Write Temperature Here! \n";
                tmp2 = ev.eval(u_sol,boundaryPoints.transpose().col(k));
                gsDebugVar(tmp2);
                results(0,k) = tmp2.at(0);
        }
        gsDebugVar(results);
        participant.writeData(GeometryControlPointMesh, TemperatureData, geometryControlPointIDs, results);
    }
    participant.readData(FluxControlPointMesh,FluxControlPointData,fluxControlPointIDs,fluxControlPoints);
    fluxMesh.patch(0).coefs() = fluxControlPoints;
    gsDebugVar(fluxControlPoints);
    A.initSystem();
    A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G), u * uex * meas(G) );
    gsSparseMatrix<> K = A.matrix();

    gsVector<> F = dt*A.rhs() + M*solVector;
    gsVector<> F0 = F;
    gsVector<> F_checkpoint = F;
    gsMatrix<> solVector_checkpoint = solVector;

    gsParaviewCollection collection("solution_neumann");
    collection.options().setSwitch("plotElements",true);
    gsParaviewCollection exact_collection("exact_solution_neumann");

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;
    real_t t_checkpoint = 0;
    if (plot)
    {
        std::string fileName = "solution_neumann" + util::to_string(timestep);
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt("plot.npts", 1000);
        ev.writeParaview( u_sol   , G, fileName);
        for (size_t p=0; p!=patches.nPatches(); p++)
        {
          fileName = "solution_neumann" + util::to_string(timestep) + std::to_string(p);
          collection.addTimestep(fileName,time,".vts");
        }
        ev.writeParaview( u_sol   , G, "initial_solution_neumann");

    }

    while (participant.isCouplingOngoing())
    {

        u_ex = gsFunctionExpr<>("1+x^2+" + std::to_string(alpha) + "*y^2 + " + std::to_string(beta) + "*" + std::to_string(time),2);
        u.setup(bcInfo, dirichlet::l2Projection, 0); // NOTE:

        A.initSystem();
        A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );
        K = A.matrix();
        // auto g_Neumann = A.getBdrFunction(G);
        // A.assembleBdr(bcInfo.get("Neumann"), u * g_Neumann.val() * nv(G).norm() );
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
        // 
        // 
        gsDebugVar(fluxControlPoints);
        for (index_t k=0; k!=boundaryPoints.transpose().cols(); k++)
        {
                gsInfo << "Write Temperature Here! \n";
                tmp2 = ev.eval(u_sol,boundaryPoints.transpose().col(k).transpose());
                results(0,k) = tmp2.at(0);
        }
        gsDebugVar(results);
        participant.writeData(GeometryControlPointMesh, TemperatureData, geometryControlPointIDs, results);
        participant.readData(FluxControlPointMesh,FluxControlPointData,fluxControlPointIDs,fluxControlPoints);
        fluxMesh.patch(0).coefs() = fluxControlPoints;
        gsDebugVar(fluxControlPoints);
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
                std::string fileName = "solution_neumann" + util::to_string(timestep);
                ev.options().setSwitch("plot.elements", true);
                ev.options().setInt("plot.npts", 1000);
                ev.writeParaview( u_sol   , G, fileName);
                for (size_t p=0; p!=patches.nPatches(); p++)
                {
                  fileName = "solution_neumann" + util::to_string(timestep) + std::to_string(p);
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

    return EXIT_SUCCESS;
}
