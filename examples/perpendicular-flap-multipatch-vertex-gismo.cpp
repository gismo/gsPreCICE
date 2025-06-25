/** @file solid-gismo-elasticity.cpp

    @brief Elasticity participant for the PreCICE example "perpendicular-flap".


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J.Li (2023 - ..., TU Delft), H.M. Verhelst (2019 - 2024, TU Delft)
*/

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsLookupFunction.h>

#include <gsElasticity/src/gsMassAssembler.h>
#include <gsElasticity/src/gsElasticityAssembler.h>

#ifdef gsStructuralAnalysis_ENABLED
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicExplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicImplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBathe.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicWilson.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicRK4.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>
#endif

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t plotmod = 10;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    std::string precice_config;
    bool nonlinear = false;
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6 RK4

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    cmd.addInt("M", "method","1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6: RK4",method);
    cmd.addSwitch("nonlinear", "Use nonlinear elasticity", nonlinear);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    // Generate domain
    gsMultiPatch<> patches;
    // get from XML

    gsKnotVector<> KV1 (0,1,2,2) ;
    gsKnotVector<> KV2 (0,1,12,2) ;

    gsTensorBSplineBasis<2> basis(KV1,KV2);
    gsMatrix<> coefs = basis.anchors();
    coefs.transposeInPlace();
    coefs.col(0).array() *= 0.1;
    coefs.col(0).array() -= 0.05;
    coefs.col(1).array() *= 0.5;


    patches.addPatch(basis.makeGeometry(coefs));
    coefs.col(1).array() += 0.5;
    patches.addPatch(basis.makeGeometry(coefs));
    patches.computeTopology();
    // Create bases
    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    // Set degree
    bases.setDegree( bases.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        bases.uniformRefine();
    numRefine = 0;

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";

    // get from XML
    real_t rho = 3000;
    real_t E = 4000000;
    real_t nu = 0.3;

    // get from XML ??
    // Set the interface for the precice coupling
    std::vector<std::vector<patchSide> > patchCouplingInterfaces(patches.nPatches());
    patchCouplingInterfaces[0].push_back(patchSide(0,boundary::east));
    patchCouplingInterfaces[0].push_back(patchSide(0,boundary::west));
    patchCouplingInterfaces[1].push_back(patchSide(1,boundary::east));
    patchCouplingInterfaces[1].push_back(patchSide(1,boundary::north));
    patchCouplingInterfaces[1].push_back(patchSide(1,boundary::west));

    /*
     * Initialize the preCICE participant
     *
     *
     */
    // get from XML
    std::string participantName = "Solid";
    gsPreCICE<real_t> participant(participantName, precice_config);

    /*
     * Data initialization
     *
     * This participant manages the geometry. The follow meshes and data are made available:
     *
     * - Meshes:
     *   + KnotMesh             This mesh contains the knots as mesh vertices
     *   + ControlPointMesh:    This mesh contains the control points as mesh vertices
     *   + ForceMesh:           This mesh contains the integration points as mesh vertices
     *
     * - Data:
     *   + ControlPointData:    This data is defined on the ControlPointMesh and stores the displacement of the control points
     *   + ForceData:           This data is defined on the ForceMesh and stores pressure/forces
     */

    std::string SolidMesh        = "Solid-Mesh";
    std::string StressData       = "Stress";
    std::string DisplacementData = "Displacement";

    // Get the quadrature nodes on the coupling interface
    gsOptionList quadOptions = gsAssembler<>::defaultOptions();

    // Get the quadrature points

    // Containers for uv and xy coordinates of the quadrature points per patch
    std::vector<gsMatrix<>> patchQuad_uv(patches.nPatches());
    std::vector<gsMatrix<>> patchQuad_xy(patches.nPatches());
    // Global containers for uv and xy coordinates of the quadrature points (to be sent to preCICE)
    gsMatrix<> quad_uv;
    gsMatrix<> quad_xy;

    // Offset per patch
    std::vector<index_t> patchOffsets(patches.nPatches()+1);
    patchOffsets[0] = 0; // first patch starts at 0
    for (size_t p=0; p!=patches.nPatches(); ++p)
    {
        // Get the quadrature points for the current patch
        patchQuad_uv[p] = gsQuadrature::getAllNodes(bases.basis(p),quadOptions,patchCouplingInterfaces[p]);
        patchQuad_xy[p] = patches.patch(p).eval(patchQuad_uv[p]);
        patchOffsets[p+1] = patchOffsets[p]+patchQuad_uv[p].cols(); // store the number of quadrature points per patch
    }

    quad_uv.resize(patches.domainDim(), patchOffsets.back());
    quad_xy.resize(patches.geoDim(), patchOffsets.back());
    for (size_t p=0; p!=patches.nPatches(); ++p)
    {
        // Copy the quadrature points of the current patch to the global containers
        quad_uv.block(0,patchOffsets[p],patches.domainDim(),patchQuad_uv[p].cols()) = patchQuad_uv[p];
        quad_xy.block(0,patchOffsets[p],patches.geoDim(),patchQuad_xy[p].cols()) = patchQuad_xy[p];
    }
    gsVector<index_t> quad_xyIDs; // needed for writing
    participant.addMesh(SolidMesh,quad_xy,quad_xyIDs);
    gsInfo << "Mesh added successfully\n";

    gsDebugVar(quad_xy);
    gsDebugVar(quad_uv);

    // Check if initial data is required (call before initialize)
    bool needsInitialData = participant.requiresInitialData();
    
    // Define precice interface
    real_t precice_dt = participant.initialize();
    
    // Write initial data after initialize if required
    if (needsInitialData)
    {
        gsInfo << "Writing initial displacement data...\n";
        gsMatrix<> initialDisplacement(patches.geoDim(), quad_xy.cols());
        initialDisplacement.setZero(); // Initial displacement is zero
        participant.writeData(SolidMesh, DisplacementData, quad_xyIDs, initialDisplacement);
    }

// ----------------------------------------------------------------------------------------------

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet side
    gsConstantFunction<> g_D(0,patches.geoDim());
    // Coupling side
    // gsFunctionExpr<> g_C("1","0",patches.geoDim());

    // get from XML ???
    // Initialize quad_stress matrix for storing stress data from PreCICE
    gsMatrix<> quad_stress(2,quad_xy.cols());
    quad_stress.setZero();
    
    // Create lookup function with quad_xy points and quad_stress data
    gsLookupFunction<real_t> g_L(quad_xy, quad_stress);

    // gsPreCICEFunction<real_t> g_C(&participant,SolidMesh,StressData,patches,patches.geoDim(),false);
    // Add all BCs
    // Coupling interface
    for (size_t p=0; p!=patchCouplingInterfaces.size(); ++p)
        for (size_t i=0; i!=patchCouplingInterfaces[p].size(); ++i)
            bcInfo.addCondition(patchCouplingInterfaces[p][i].patch, patchCouplingInterfaces[p][i].side(), condition_type::neumann, &g_L, -1, true);

    // Bottom side (prescribed temp)
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D, 1);
    // Assign geometry map
    bcInfo.setGeoMap(patches);

// ----------------------------------------------------------------------------------------------
    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    // source function, rhs
    gsConstantFunction<> gravity(0.,0.,2);

    // creating mass assembler
    gsMassAssembler<real_t> massAssembler(patches,bases,bcInfo,gravity);
    massAssembler.options().setReal("Density",rho);
    massAssembler.assemble();

    // creating stiffness assembler.
    gsElasticityAssembler<real_t> assembler(patches,bases,bcInfo,g);
    assembler.options().setReal("YoungsModulus",E);
    assembler.options().setReal("PoissonsRatio",nu);
    assembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    assembler.assemble();

    gsMatrix<real_t> Minv;
    gsSparseMatrix<> M = massAssembler.matrix();
    gsSparseMatrix<> K = assembler.matrix();
    gsSparseMatrix<> K_T;

    // Time step
    real_t dt = precice_dt;

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    // RHS of the projection
    gsMatrix<> solVector;
    solVector.setZero(assembler.numDofs(),1);

    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
    // Assemble the RHS
    gsVector<> U_checkpoint, V_checkpoint, A_checkpoint, U, V, A;

    U_checkpoint = U = gsVector<real_t>::Zero(assembler.numDofs(),1);
    V_checkpoint = V = gsVector<real_t>::Zero(assembler.numDofs(),1);
    A_checkpoint = A = gsVector<real_t>::Zero(assembler.numDofs(),1);

    // Define the solution collection for Paraview
    gsParaviewCollection collection("./output/solution");

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;

    // Function for the Jacobian
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&fixedDofs](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> & m) {
        // to do: add time dependency of forcing
        assembler.assemble(x, fixedDofs);
        m = assembler.matrix();
        return true;
    };

    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::TResidual_t Residual = [&assembler,&fixedDofs](gsMatrix<real_t> const &x, real_t t, gsVector<real_t> & result)
    {
        assembler.assemble(x,fixedDofs);
        result = assembler.rhs();
        return true;
    };

    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::TForce_t TForce = [&assembler](real_t, gsVector<real_t> & result)
    {
        assembler.assemble();
        result = assembler.rhs();
        return true;
    };


    gsSparseMatrix<> C = gsSparseMatrix<>(assembler.numDofs(),assembler.numDofs());
    gsStructuralAnalysisOps<real_t>::Damping_t Damping = [&C](const gsVector<real_t> &, gsSparseMatrix<real_t> & m) { m = C; return true; };
    gsStructuralAnalysisOps<real_t>::Mass_t    Mass    = [&M](                          gsSparseMatrix<real_t> & m) { m = M; return true; };
    gsStructuralAnalysisOps<real_t>::Stiffness_t Stiffness = [&K](                      gsSparseMatrix<real_t> & m) { m = K; return true; };

    // if (nonlinear)
    //     timeIntegrator = new gsDynamicNewmark<real_t,true>(Mass,Damping,Jacobian,Residual);
    // else
    //     timeIntegrator = new gsDynamicNewmark<real_t,false>(Mass,Damping,Stiffness,TForce);



    gsDynamicBase<real_t> * timeIntegrator;
    if (method==1)
        timeIntegrator = new gsDynamicExplicitEuler<real_t,true>(Mass,Damping,Jacobian,Residual);
    else if (method==2)
        timeIntegrator = new gsDynamicImplicitEuler<real_t,true>(Mass,Damping,Jacobian,Residual);
    else if (method==3 && nonlinear==false)
        timeIntegrator = new gsDynamicNewmark<real_t,false>(Mass,Damping,Stiffness,TForce);
    else if (method==3 && nonlinear==true)
            timeIntegrator = new gsDynamicNewmark<real_t,true>(Mass,Damping,Jacobian,Residual);
    else if (method==4)
        timeIntegrator = new gsDynamicBathe<real_t,true>(Mass,Damping,Jacobian,Residual);
    else if (method==5)
    {
        timeIntegrator = new gsDynamicWilson<real_t,true>(Mass,Damping,Jacobian,Residual);
        timeIntegrator->options().setReal("gamma",1.4);
    }
    else if (method==6)
        timeIntegrator = new gsDynamicRK4<real_t,true>(Mass,Damping,Jacobian,Residual);
    else
        GISMO_ERROR("Method "<<method<<" not known");

    timeIntegrator->options().setReal("DT",dt);
    timeIntegrator->options().setReal("TolU",1e-6);
    timeIntegrator->options().setSwitch("Verbose",true);

    real_t time = 0;

    // Plot initial solution
    if (plot)
    {
        gsMultiPatch<> solution;
        assembler.constructSolution(solVector,fixedDofs,solution);

        gsField<> solField(patches,solution);
        std::string fileName = "./output/solution" + util::to_string(timestep);
        gsWriteParaview<>(solField, fileName, 500);
        for (size_t p=0; p!=patches.nPatches(); ++p)
        {
            fileName = "solution" + util::to_string(timestep) + util::to_string(p);
            collection.addPart(fileName,time,"Solution",p);
        }
    }

    gsMatrix<> refPoints(2,1);
    refPoints.col(0)<<0.5,1;
    gsVector<index_t> refPatches(1);
    refPatches[0] = 1; // patch 0

    gsStructuralAnalysisOutput<real_t> writer("pointData.csv",refPoints);
    writer.init({"x","y"},{"time"}); // point1 - x, point1 - y, time

    gsMatrix<> pointDataMatrix;
    gsMatrix<> otherDataMatrix(1,1);
    
    // Time integration loop
    gsInfo << "Starting coupling time loop...\n";
    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
        {
            U_checkpoint = U;
            V_checkpoint = V;
            A_checkpoint = A;

            gsInfo<<"Checkpoint written:\n";
            gsInfo<<"\t ||U|| = "<<U.norm()<<"\n";
            gsInfo<<"\t ||V|| = "<<V.norm()<<"\n";
            gsInfo<<"\t ||A|| = "<<A.norm()<<"\n";


            timestep_checkpoint = timestep;
        }

        participant.readData(SolidMesh,StressData,quad_xyIDs,quad_stress);
        g_L.update();
        // solve gismo timestep
        gsInfo << "Solving timestep " << time << "...\n";
        timeIntegrator->step(time,dt,U,V,A);
        solVector = U;
        gsInfo<<"Finished\n";

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt,precice_dt);


        gsMultiPatch<> solution;
        assembler.constructSolution(solVector,fixedDofs,solution);
        gsMatrix<> result(patches.geoDim(),patchOffsets.back());
        gsMatrix<> tmpResult;
        for (size_t p=0; p!=patches.nPatches(); ++p)
        {
            // evaluate the solution at the quadrature points
            solution.patch(p).eval_into(patchQuad_uv[p],tmpResult);
            result.block(0,patchOffsets[p],patches.geoDim(),patchQuad_uv[p].cols()) = tmpResult;
        }
        participant.writeData(SolidMesh,DisplacementData,quad_xyIDs,result);

        // do the coupling
        precice_dt =participant.advance(dt);

        if (participant.requiresReadingCheckpoint())
        {
            U = U_checkpoint;
            V = V_checkpoint;
            A = A_checkpoint;
            timestep = timestep_checkpoint;
        }
        else
        {
            // gsTimeIntegrator advances
            // advance variables
            time += dt;
            timestep++;

            gsField<> solField(patches,solution);
            if (timestep % plotmod==0) // Generate Paraview output for visualization
            {
                if (plot)
                {
                    // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
                    std::string fileName = "./output/solution" + util::to_string(timestep);
                    gsWriteParaview<>(solField, fileName, 500);
                    for (size_t p=0; p!=patches.nPatches(); ++p)
                    {
                        fileName = "solution" + util::to_string(timestep) + util::to_string(p);
                        collection.addPart(fileName,time,"Solution",p);
                    }
                }

                gsMatrix<> tmpData;
                pointDataMatrix.resize(solution.targetDim(),refPoints.cols());
                for (size_t k=0; k!=refPoints.cols(); ++k)
                {
                    // evaluate the solution at the reference points
                    solution.patch(refPatches[k]).eval_into(refPoints.col(k),tmpData);
                    if (pointDataMatrix.rows()==0)
                        pointDataMatrix.resize(tmpData.rows(),refPoints.cols());
                    pointDataMatrix.col(k) = tmpData.col(0);
                }

                solution.patch(0).eval_into(refPoints,pointDataMatrix);
                otherDataMatrix<<time;
                writer.add(pointDataMatrix,otherDataMatrix);
            }
        }
    }

    if (plot)
    {
        collection.save();
    }


    return  EXIT_SUCCESS;
}
