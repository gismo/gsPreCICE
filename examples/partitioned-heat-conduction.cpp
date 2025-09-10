/** @file partitioned-heat-conduction.cpp

    @brief Heat equation participant for a double coupled heat equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (2019-2024, TUDelft), J. Li (2023-..., TUDelft)
*/

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEFunction.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t plotmod = 1;
    index_t numRefine  = 2; // Number of uniform refinement (increase the number also requires to modify the rbf radius in the precice_config.xml)
    index_t numElevate = 0;
    short_t side = 0;
    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Coupled heat equation using PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    cmd.addInt("s","side", "Patchside of interface: 0 (Dirichlet), 1 (Neumann)", side);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]

    gsMultiPatch<> patches;
    if (side==0) //left
        patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,1.0,1.0));
    else if (side==1) //right
        patches.addPatch(gsNurbsCreator<>::BSplineRectangle(1.0,0.0,2.0,1.0));
    else
        GISMO_ERROR("Side unknown");

    real_t alpha = 3;
    real_t beta  = 1.3;
    real_t time  = 0;
    real_t k_temp = 1;

    // Set external heat-flux to zero
    gsConstantFunction<> f(beta-2-2*alpha,2);
    gsFunctionExpr<> u_ex("1+x^2+" + std::to_string(alpha) + "*y^2 + " + std::to_string(beta) + "*" + std::to_string(time),2);

    gsMultiBasis<> bases(patches); // true: poly-splines (not NURBS)

    bases.setDegree( bases.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        bases.uniformRefine();
    numRefine = 0;

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";

// ----------------------------------------------------------------------------------------------
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    A.setIntegrationElements(bases);

    gsExprEvaluator<> ev(A);

// ----------------------------------------------------------------------------------------------
    boxSide couplingSide;
    if (side==0) //left
        couplingSide = boxSide(2);
    else if (side==1) //right
        couplingSide = boxSide(1);
    else
        GISMO_ERROR("Side unknown");

    patchSide couplingInterface(0,couplingSide);
    // Use STL-style domain iterators on the boundary side
    const gsBasis<real_t>& basisOnPatch = bases.basis(couplingInterface.patch);
    typename gsBasis<real_t>::domainIter domIt    = basisOnPatch.domain()->beginBdr(couplingInterface.side());
    typename gsBasis<real_t>::domainIter domItEnd = basisOnPatch.domain()->endBdr(couplingInterface.side());
    index_t rows = patches.targetDim();
    gsMatrix<> nodes;
    // Start iteration over elements
    gsVector<> tmp;
    index_t k=0;

    gsOptionList quadOptions = A.options();
    // Quadrature options can be customized via A.options()

    index_t quadSize = 0;
    typename gsQuadRule<real_t>::uPtr QuRule; // Quadrature rule  ---->OUT
    for (; domIt < domItEnd; ++domIt, ++k)
    {
        QuRule = gsQuadrature::getPtr(basisOnPatch, quadOptions, couplingInterface.side().direction());
        quadSize += QuRule->numNodes();
    }
    gsMatrix<> uv(rows,quadSize); // Coordinates of the quadrature points in parameter space
    gsMatrix<> xy(rows,quadSize); // Coordinates of the quadrature points in physical space

    index_t offset = 0;

    // Reset iterator and loop again to fill uv/xy
    domIt = basisOnPatch.domain()->beginBdr(couplingInterface.side());
    for (; domIt < domItEnd; ++domIt, ++k )
    {
        QuRule = gsQuadrature::getPtr(basisOnPatch, quadOptions, couplingInterface.side().direction());
        // Map the Quadrature rule to the element
        QuRule->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                       nodes, tmp);
        uv.block(0,offset,rows,QuRule->numNodes()) = nodes;

        gsMatrix<> tmp2;
        patches.patch(couplingInterface.patch).eval_into(nodes,tmp2);
        xy.block(0,offset,rows,QuRule->numNodes()) = patches.patch(couplingInterface.patch).eval(nodes);
        offset += QuRule->numNodes();
    }

    std::string membername;
    if (side==0) //left participant
        membername = "Dirichlet";
    else if (side==1) //right participant
        membername = "Neumann";
    else
        GISMO_ERROR("Side unknown");

    // PreCICE participants and meshes per provided config
    // Dirichlet provides Dirichlet-Mesh, reads Temperature on Dirichlet-Mesh, writes Heat-Flux on Neumann-Mesh
    // Neumann   provides Neumann-Mesh,   reads Heat-Flux on Neumann-Mesh,   writes Temperature on Dirichlet-Mesh
    const std::string localMeshName  = membername + "-Mesh";
    const std::string remoteMeshName = (membername=="Dirichlet") ? std::string("Neumann-Mesh") : std::string("Dirichlet-Mesh");
    const std::string tempName = "Temperature";
    const std::string fluxName = "Heat-Flux";

    gsPreCICE<real_t> interface(membername, precice_config);
    // Register local mesh coordinates
    interface.addMesh(localMeshName, xy);

    // Define access region for the received mesh in serial mode
    // Use a generous bounding box pattern seen in other examples
    {
        gsMatrix<> bbox(rows, 2);
        bbox.col(0).setConstant(-1e300);
        bbox.col(1).setConstant( 1e300);
        bbox.transposeInPlace();
        bbox.resize(1, bbox.rows() * bbox.cols());
        interface.setMeshAccessRegion(remoteMeshName, bbox);
    }

    interface.requiresInitialData();

    real_t precice_dt = interface.initialize();

    // Direct-access: fetch partner mesh vertex IDs and coordinates.
    // Evaluate quantities at partner coordinates (no reordering needed).
    gsVector<index_t> remoteIDs;
    gsMatrix<>        remoteCoords;
    interface.getMeshVertexIDsAndCoordinates(remoteMeshName, remoteIDs, remoteCoords);
    // Invert partner physical coords to param coords on our coupling patch
    gsMatrix<> uvRemote;
    patches.patch(couplingInterface.patch).invertPoints(remoteCoords, uvRemote, 1e-10);
    // Ensure boundary parameter is exact for evalBdr
    const int bdir = couplingInterface.side().direction();
    const real_t bpar = couplingInterface.side().parameter();
    for (index_t i = 0; i != uvRemote.cols(); ++i)
        uvRemote(bdir, i) = bpar;

// ----------------------------------------------------------------------------------------------

    gsBoundaryConditions<> bcInfo;
    // Read functions according to participant role and config
    // Dirichlet reads Temperature on Dirichlet-Mesh; Neumann reads Heat-Flux on Neumann-Mesh
    gsPreCICEFunction<real_t> g_CD(&interface, (membername=="Dirichlet" ? localMeshName  : remoteMeshName), tempName, patches, 1);
    gsPreCICEFunction<real_t> g_CN(&interface, (membername=="Neumann"   ? localMeshName  : remoteMeshName), fluxName, patches, 1);
    gsFunction<> * g_C = &u_ex;

    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &u_ex, 0, false, 0);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &u_ex, 0, false, 0);
    if (side==0)
    {
        bcInfo.addCondition(0, couplingSide,condition_type::dirichlet  , g_C, 0, false, 0);
        bcInfo.addCondition(0, couplingSide.opposite(),condition_type::dirichlet, &u_ex, 0, false, 0);
    }
    else
    {
        bcInfo.addCondition(0, couplingSide,condition_type::neumann  , &g_CN);
        bcInfo.addCondition(0, couplingSide.opposite(),condition_type::dirichlet, &u_ex, 0, false, 0);
    }

    bcInfo.setGeoMap(patches);
    // gsDebugVar(bcInfo);

// ----------------------------------------------------------------------------------------------

    // Generate system matrix and load vector
    // gsInfo << "Assembling mass and stiffness...\n";

    // Set the geometry map
    geometryMap G = A.getMap(patches);

    // Set the discretization space
    space u = A.getSpace(bases);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Set the solution
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    u.setup(bcInfo, dirichlet::homogeneous, 0);
    A.initSystem();
    A.assemble( u * u.tr() * meas(G));
    gsSparseMatrix<> M = A.matrix();


    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner 
    gsSparseSolver<>::CGDiagonal solver;

    real_t dt = 0.01;

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    auto uex = A.getCoeff(u_ex, G);
    // RHS of the projection
    u.setup(bcInfo, dirichlet::l2Projection, 0);
    A.initSystem();
    A.assemble( u * u.tr() * meas(G), u * uex * meas(G) );
    solver.compute(A.matrix());
    solVector = solver.solve(A.rhs());
    gsMatrix<> result(1,uv.cols()), tmp2;
    // Write initial data
    if (interface.requiresWritingCheckpoint())
    {
        for (index_t k=0; k!=uv.cols(); k++)
        {
            if (side==0)
            {
                // Normal heat flux q_n = -k * grad(u) · n on the coupling boundary
                tmp2 = ev.evalBdr( - k_temp * (igrad(u_sol, G) * nv(G).normalized()), uv.col(k), couplingInterface);
                result(0,k) = tmp2.at(0);
            }
            else
            {
                tmp2 = ev.eval(u_sol,uv.col(k));
                result(0,k) = tmp2.at(0);
            }
        }
        // Build values at partner coordinates (order matches remoteIDs)
        result.resize(1, uvRemote.cols());
        for (index_t k=0; k!=uvRemote.cols(); ++k)
        {
            if (side==0) // Dirichlet writes flux
            {
                tmp2 = ev.evalBdr( - k_temp * (igrad(u_sol, G) * nv(G).normalized()), uvRemote.col(k), couplingInterface);
                result(0,k) = tmp2.at(0);
            }
            else // Neumann writes temperature
            {
                tmp2 = ev.evalBdr(u_sol, uvRemote.col(k), couplingInterface);
                result(0,k) = tmp2.at(0);
            }
        }
        // Write on partner mesh using its vertex IDs (direct-access config)
        interface.writeData(remoteMeshName, (side==0 ? fluxName : tempName), remoteIDs, result);
    }

    // Initialize the RHS for assembly
    if (side==0)
        g_C = &g_CD;
    A.initSystem();
    A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G), u * uex * meas(G) );
    gsSparseMatrix<> K = A.matrix();

    // Assemble the RHS
    gsVector<> F = dt*A.rhs() + M*solVector;
    gsVector<> F0 = F;
    gsVector<> F_checkpoint = F;
    gsMatrix<> solVector_checkpoint = solVector;

    gsParaviewCollection collection(membername + "_solution");
    gsParaviewCollection exact_collection(membername + "_exact_solution");

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;
    real_t t_checkpoint = 0;
    if (plot)
    {
        std::string fileName = membername + "_solution_" + util::to_string(timestep);
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt("plot.npts", 1000);
        ev.writeParaview( u_sol   , G, fileName);
        for (size_t p=0; p!=patches.nPatches(); p++)
        {
          fileName = membername + "_solution_" + util::to_string(timestep) + std::to_string(p);
          collection.addTimestep(fileName,time,".vts");
        }

        // write exact solution at initial time
        fileName = membername + "_exact_solution_" + util::to_string(timestep);
        {
            auto uex_coeff = A.getCoeff(u_ex, G);
            ev.writeParaview( uex_coeff , G, fileName);
        }
        for (size_t p=0; p!=patches.nPatches(); p++)
        {
          std::string fileName2 = membername + "_exact_solution_" + util::to_string(timestep) + std::to_string(p);
          exact_collection.addTimestep(fileName2,time,".vts");
        }

        ev.writeParaview( u_sol   , G, membername + "_initial_solution");

    }

    // time += precice_dt;
    while (interface.isCouplingOngoing())
    {
        // read temperature from interface
            u_ex = gsFunctionExpr<>("1+x^2+" + std::to_string(alpha) + "*y^2 + " + std::to_string(beta) + "*" + std::to_string(time),2);
            u.setup(bcInfo, dirichlet::l2Projection, 0); // NOTE: dirichlet::l2Projection is used to project the Dirichlet BCs

            A.initSystem();
            A.assemble( k_temp * igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );
            K = A.matrix();
            auto g_Neumann = A.getBdrFunction(G);
            A.assembleBdr(bcInfo.get("Neumann"), u * g_Neumann.val() * nv(G).norm() );
            F = dt*A.rhs() + M*solVector;

        // save checkpoint
        if (interface.requiresWritingCheckpoint())
        {
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
        gsMatrix<> result(1,uv.cols()), tmp;
        for (index_t k=0; k!=uv.cols(); k++)
        {
            if (side==0)
            {
                // Normal heat flux q_n = -k * grad(u) · n on the coupling boundary
                tmp = ev.evalBdr( - k_temp * (igrad(u_sol, G) * nv(G).normalized()), uv.col(k), couplingInterface);
                result(0,k) = tmp.at(0);
            }
            else
            {
                tmp = ev.eval(u_sol,uv.col(k));
                result(0,k) = tmp.at(0);
            }
        }
        // Build values at partner coordinates (order matches remoteIDs)
        result.resize(1, uvRemote.cols());
        for (index_t k=0; k!=uvRemote.cols(); ++k)
        {
            if (side==0) // Dirichlet writes flux
            {
                tmp = ev.evalBdr( - k_temp * (igrad(u_sol, G) * nv(G).normalized()), uvRemote.col(k), couplingInterface);
                result(0,k) = tmp.at(0);
            }
            else // Neumann writes temperature
            {
                tmp = ev.evalBdr(u_sol, uvRemote.col(k), couplingInterface);
                result(0,k) = tmp.at(0);
            }
        }
        // Write on partner mesh using its vertex IDs (direct-access config)
        interface.writeData(remoteMeshName, (side==0 ? fluxName : tempName), remoteIDs, result);

        // do the coupling
        precice_dt = interface.advance(dt);

        // advance variables
        time += dt;
        timestep += 1;
        F0 = F;

        if (interface.requiresReadingCheckpoint())
        {
            F0 = F_checkpoint;
            time = t_checkpoint;
            timestep = timestep_checkpoint;
            solVector = solVector_checkpoint;
        }
        else
        {
            if (timestep % plotmod==0 && plot)
            {
                std::string fileName = membername + "_solution_" + util::to_string(timestep);
                ev.options().setSwitch("plot.elements", true);
                ev.options().setInt("plot.npts", 1000);
                ev.writeParaview( u_sol   , G, fileName);
                for (size_t p=0; p!=patches.nPatches(); p++)
                {
                  fileName = membername + "_solution_" + util::to_string(timestep) + std::to_string(p);
                  collection.addTimestep(fileName,time,".vts");
                }

                // write exact solution at current time
                std::string fileNameExact = membername + "_exact_solution_" + util::to_string(timestep);
                {
                    auto uex_coeff = A.getCoeff(u_ex, G);
                    ev.writeParaview( uex_coeff , G, fileNameExact);
                }
                for (size_t p=0; p!=patches.nPatches(); p++)
                {
                  std::string fileName2 = membername + "_exact_solution_" + util::to_string(timestep) + std::to_string(p);
                  exact_collection.addTimestep(fileName2,time,".vts");
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
