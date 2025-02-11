#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsLookupFunction.h>

#include <gsAssembler/gsExprEvaluator.h>

// #include <gsPreCICE/gsPreCICEVectorFunction.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElasticityAssembler.h>
 
#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsKLShell/src/gsThinShellAssembler.h>
#include <gsKLShell/src/gsMaterialMatrixContainer.h>
#include <gsKLShell/src/gsMaterialMatrixEval.h>
#include <gsKLShell/src/gsMaterialMatrixIntegrate.h>
#endif

 
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
    bool write = false;
    bool get_readTime = false;
    bool get_writeTime = false;
    index_t plotmod = 1;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    std::string precice_config;

    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6 RK4
    std::string dirname = "./output";

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    //cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    cmd.addInt("m", "method","1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson",method);
    cmd.addSwitch("write", "Create a file with point data", write);
    cmd.addSwitch("readTime", "Get the read time", get_readTime);
    cmd.addSwitch("writeTime", "Get the write time", get_writeTime);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");
    // Generate domain
    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,0.2,1.0));// HUGO: Is this size correct?
    gsMatrix<> coefs(3,patches.patch(0).coefs().rows());

    coefs.row(0).setZero();
    coefs.row(1) = patches.patch(0).coefs().col(1).transpose();
    coefs.row(2) = patches.patch(0).coefs().col(0).transpose();
    patches.embed(3);
    patches.patch(0).coefs() = coefs.transpose();


    patches.addAutoBoundaries();
 
    gsDebugVar(patches.patch(0).coefs());

    gsWriteParaview(patches,"./output/geometry"); 

    // Embed the 2D geometry to 3D
    gsMultiPatch<> solutions;
    // patches.embed(3);
    // source function, rhs
    gsConstantFunction<> g(0.,0.,3);

    // source function, rhs
    gsConstantFunction<> gravity(0.,0.,3);
    // Create bases
    // p-refine
    if (numElevate!=0)
        patches.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        patches.uniformRefine();
    // Create bases
    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    //----------------------------------[! Setup expression assembler for normal vectors] ----------------------------------
    // Set the expression assembler
    gsExprAssembler<> A_expr(1,1);
    // Set the parameter mesh as the integration mesh
    A_expr.setIntegrationElements(bases);
    // Set the expression evaluator

    gsExprEvaluator<> ev(A_expr);

    // Define integrant variables
    auto G = ev.getMap(patches);
    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n"; 

    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;
    // real_t nu = 0.3;
    // real_t mu = E / (2.0 * (1.0 + nu));
    // real_t lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    /*
     * Initialize the preCICE participant
     *
     *
     */
    std::string participantName = "Solid";
    gsPreCICE<real_t> participant(participantName, precice_config);

    // ----------------------------------------------------------------------------------------------
    // typedef gsExprAssembler<>::geometryMap geometryMap;
    // typedef gsExprAssembler<>::space       space;
    // typedef gsExprAssembler<>::solution    solution;

    /*
     * Data initialization
     *
     * This participant creates its own mesh, and it writes and reads displacements and stresses on its mesh.
     * The follow meshes and data are made available:
     *
     * - Meshes:
     *   + SolidMesh:               This mesh contains the integration points as mesh vertices
     *
     * - Data:
     *   + DisplacementData:    This data is defined on the SolidMesh and stores the displacement at the integration points
     *   + StressData:           This data is defined on the SolidMesh and stores pressure/forces at the integration points
     */
    std::string SolidMesh        = "Solid-Mesh";
    std::string StressData       = "Stress";
    std::string DisplacementData = "Displacement";
    // std::string ForceMesh        = "Fluid-Mesh";


    // HUGO: These are incorrect. Question: should they be boundary::??
    // std::vector<patchSide> couplingInterfaces(3);
    // couplingInterfaces[0] = patchSide(0,boundary::west);
    // couplingInterfaces[1] = patchSide(0,boundary::north);
    // couplingInterfaces[2] = patchSide(0,boundary::east);
    // gsOptionList quadOptions = A_quad.options(); 

    // Step 1: SolidMesh
    // Get the quadrature nodes on the coupling interface
    gsOptionList quadOptions = gsExprAssembler<>::defaultOptions();
    // Get the quadrature points
    gsVector<index_t> quadPointIDs;

    //We only want the left boundary quadrature points
    gsMatrix<> quadPointsAll = gsQuadrature::getAllNodes(bases.basis(0), quadOptions);

 

    //Give a side for the function
    index_t numQuadPt = 0;
    index_t counter = 0;
    real_t target = quadPointsAll(0,0);
    for (index_t i = 0; i < quadPointsAll.cols(); ++i) {
        if (quadPointsAll(0, i) == target)
        {  // Check first row for target
            numQuadPt++;
        }
    }
    gsDebugVar(quadPointsAll.dim());
    gsMatrix<> quadPoints(2, numQuadPt);
    gsMatrix<> quadPointNormals(3, numQuadPt);
    gsMatrix<> tempNormal_param(2,1);
    gsMatrix<> tempNormal(3,1);
    // Second pass: Fill quadPoints with the target values
    for (index_t i = 0; i < quadPointsAll.cols(); ++i) {
        if (quadPointsAll(0, i) == target) {  // Check for target again
            quadPoints.col(counter) = quadPointsAll.col(i);  // Fill the first row
            // quadPoints(1, counter) = quadPointsAll(1, i);  // Fill the second row
            // tempNormal= quadPointsAll.data().normal(i);
            tempNormal = ev.eval( sn(G).normalized(), quadPointsAll.col(i) );
            quadPointNormals.col(counter) = tempNormal;
            counter++;
        }
    }
    gsDebugVar(quadPointNormals);
    // Evaluate the y coordinates
    gsDebugVar(quadPoints);
    gsMatrix<> phyQuads;

    patches.patch(0).eval_into(quadPoints, phyQuads);
    gsDebugVar(phyQuads);
    gsMatrix<> phyQuadsAll;
    patches.patch(0).eval_into(quadPointsAll, phyQuadsAll);

   
    gsDebugVar(phyQuadsAll);
    gsDebugVar(phyQuadsAll.col(0));
    gsMatrix<> comPtLeft(numQuadPt,2);
    comPtLeft.setZero();
    gsMatrix<> comPtRight(numQuadPt,2);
    comPtRight.setZero();

    comPtLeft.col(1).transpose() = phyQuads.row(1);
    comPtRight.col(1).transpose() = phyQuads.row(1);
    gsDebugVar(comPtLeft);

    comPtLeft.col(0) = comPtLeft.col(0).array() -0.05;
    comPtRight.col(0) = comPtRight.col(0).array() +0.05;


    gsMatrix<> comPt(2*numQuadPt,2);
    comPt << comPtLeft, comPtRight;

    gsDebugVar(comPt);
    participant.addMesh(SolidMesh,comPt.transpose(),quadPointIDs); //Set the vertices to be datapoints (quadpoints in physical domain)

    real_t precice_dt = participant.initialize();
    // Define boundary conditions

    gsBoundaryConditions<> bcInfo;

    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 1);
    bcInfo.addCondition(0, boundary::south, condition_type::clamped, 0, 0, false, -1);
    bcInfo.setGeoMap(patches);

    // gsMatrix<> comPointsData(2, comPt.rows());

    // gsMatrix<> comForceData(3, comPt.rows());

 

    // // gsLookupFunction<real_t> comForce(comPt, comForceData);

    // comForceData.setRandom();

 

    // // Create a new matrix for the mid-plane force

    gsMatrix<> averagePt(2, numQuadPt);


    // // Calculate the average for each corresponding element
    // //Is this really an average? (maybe need to consider normal vector)
    // //
    gsVector<index_t> comForceIDs;
    // for (index_t i = 0; i < numQuadPt; ++i)
    // {
    //     // Add the first half and second half values and divide by 2
    //     averagePt(0, i) = (comForceData(0, i) + comForceData(0, i + numQuadPt)) / 2.0; // First row average
    //     averagePt(1, i) = (comForceData(1, i) + comForceData(1, i + numQuadPt)) / 2.0; // Second row average
    // }
    // gsDebugVar(averagePt);

    index_t num_repeat = quadPointsAll.cols()/numQuadPt;
    gsMatrix<> quadPointDatay = gsMatrix<>::Zero(averagePt.cols(), num_repeat);
    gsMatrix<> quadPointDataz = gsMatrix<>::Zero(averagePt.cols(), num_repeat);

    gsDebugVar("Here");
    gsDebugVar(quadPointDatay);
 
    for (index_t i = 0; i < num_repeat; ++i)
    {
        quadPointDatay.col(i) = averagePt.row(0).transpose();
        quadPointDataz.col(i) = averagePt.row(1).transpose();
    }
    quadPointDatay.resize(1, quadPointsAll.cols());
    quadPointDataz.resize(1, quadPointsAll.cols());
    gsMatrix<> quadPointsData(3, quadPointsAll.cols());
    quadPointsData.setZero();
    gsDebugVar(quadPointsData.dim());
    gsDebugVar(quadPointsAll.dim());
    // quadPointsData.row(0) << quadPointDatax;
    // quadPointsData.row(1) << quadPointDatay;
    // gsDebugVar(quadPointsData);

      // Set up the material matrices
    gsFunctionExpr<> E_modulus(std::to_string(E),3);
    gsFunctionExpr<> PoissonRatio(std::to_string(nu),3);
    gsFunctionExpr<> Density(std::to_string(rho),3);

    //Define thickness
    real_t thickness = 0.1;
 

    gsFunctionExpr<> t(std::to_string(thickness),3);
    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E_modulus;
    parameters[1] = &PoissonRatio;

    gsOptionList options;

 

    // gsMaterialMatrixBase<real_t>* materialMatrix;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3, real_t>(patches, t, parameters, Density, options);


    // Question: how to map the force information onto quad points as a surface force for gsThinShellAssembler?
    gsLookupFunction<real_t> surfForce(quadPointsAll, quadPointsData);
    // gsMatrix<> displacementData = gsMatrix<>::Zero(3, comPt.rows());
    // displacementData.setRandom();

 
    gsThinShellAssembler<3, real_t, false> assembler(patches, bases, bcInfo, surfForce, materialMatrix);
    gsOptionList assemblerOptions = options.wrapIntoGroup("Assembler");

    assembler.assemble();
    assembler.setOptions(assemblerOptions);

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;


    // Compute the mass matrix (since it is constant over time)
    assembler.assembleMass();
    gsSparseMatrix<> M = assembler.massMatrix();
    assembler.assemble();
    gsSparseMatrix<> K = assembler.matrix();

 
    // Define the solution collection for Paraview
    gsFileManager::mkdir(dirname);
    gsParaviewCollection collection(dirname + "/solution");

    // Time step
    // real_t dt = 0.1;
    real_t t_read = 0;
    real_t t_write = 0;
    real_t dt = precice_dt;


    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler,&solutions](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> & m)
    {
        // to do: add time dependency of forcing
        // For the shell
        ThinShellAssemblerStatus status;
        assembler.constructSolution(x, solutions);
        status = assembler.assembleMatrix(solutions);
        m = assembler.matrix();
        return true;
    };


    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::TResidual_t Residual = [&assembler,&solutions](gsMatrix<real_t> const &x, real_t t, gsVector<real_t> & result)
    {
        //Add assemble vector JL
        ThinShellAssemblerStatus status;
        assembler.constructSolution(x,solutions);
        status = assembler.assembleVector(solutions);
        result = assembler.rhs();
        return true;
    };
    gsSparseMatrix<> C = gsSparseMatrix<>(assembler.numDofs(),assembler.numDofs());
    gsStructuralAnalysisOps<real_t>::Damping_t Damping = [&C](const gsVector<real_t> &, gsSparseMatrix<real_t> & m) { m = C; return true; };
    gsStructuralAnalysisOps<real_t>::Mass_t    Mass    = [&M](                          gsSparseMatrix<real_t> & m) { m = M; return true; };

    gsDynamicBase<real_t> * timeIntegrator;
    if (method==1)
        timeIntegrator = new gsDynamicExplicitEuler<real_t,true>(Mass,Damping,Jacobian,Residual);
    else if (method==2)
        timeIntegrator = new gsDynamicImplicitEuler<real_t,true>(Mass,Damping,Jacobian,Residual);
    else if (method==3)
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
    timeIntegrator->options().setReal("TolU",1e-3);
    timeIntegrator->options().setSwitch("Verbose",true);

 

    // Project u_wall as ini
    // tial condition (violates Dirichlet side on precice interface)
    // RHS of the projection
    gsMatrix<> solVector;
    solVector.setZero(assembler.numDofs(),1);

    // Assemble the RHS
    gsVector<> F = assembler.rhs();
    gsDebugVar(F);
    gsVector<> F_checkpoint, U_checkpoint, V_checkpoint, A_checkpoint, U, V, A;

    F_checkpoint = F;
    U_checkpoint = U = gsVector<real_t>::Zero(assembler.numDofs(),1);
    V_checkpoint = V = gsVector<real_t>::Zero(assembler.numDofs(),1);
    A_checkpoint = A = gsVector<real_t>::Zero(assembler.numDofs(),1);
 
    real_t time = 0;

    if (plot)
    {
        gsMultiPatch<> solution;
        gsVector<> displacements = U;
        solution = assembler.constructDisplacement(displacements);
        solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
        gsField<> solField(patches,solution);
        std::string fileName = dirname + "/solution" + util::to_string(timestep);
        gsWriteParaview<>(solField, fileName, 500);
        fileName = "solution" + util::to_string(timestep) + "0";
        collection.addTimestep(fileName,time,".vts");
    }

    gsMatrix<> points(2,1);
    points.col(0)<<0.5,1;

 

    gsStructuralAnalysisOutput<real_t> writer("./output/pointData.csv",points);
    writer.init({"x","y","z"},{"time"}); // point1 - x, point1 - y, point1 - z, time
 

    gsMatrix<> pointDataMatrix;
    gsMatrix<> otherDataMatrix(1,1);

    gsMatrix<>comForceData(3, comPt.rows());
    gsMatrix<>quadPointsNormalData(3, quadPointsAll.cols());
    // comForceData.setZero();
 

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

        // assembler.assemble();

        // F = assembler.rhs();

 

        gsDebugVar(bcInfo);
        participant.readData(SolidMesh,StressData,quadPointIDs,comForceData);

        gsDebugVar(comForceData);

 
        //Is this really an average? (maybe need to consider normal vector)

        for (index_t i = 0; i < numQuadPt; ++i)
        {
            // Add the first half and second half values and divide by 2
            averagePt.col(i) = (comForceData.col(i) * quadPointNormals(i) + comForceData.col(i + numQuadPt) * quadPointNormals(i + numQuadPt)) / 2.0;
        }

        quadPointDatay = gsMatrix<>::Zero(averagePt.cols(), num_repeat);
        quadPointDataz = gsMatrix<>::Zero(averagePt.cols(), num_repeat);


        for (index_t i = 0; i < num_repeat; ++i)
        {
            quadPointDatay.col(i) = averagePt.row(1).transpose();
            quadPointDataz.col(i) = averagePt.row(0).transpose();
        }

       quadPointDatay.resize(1, quadPointsAll.cols());
       quadPointDataz.resize(1, quadPointsAll.cols());
       quadPointsData.row(1) << quadPointDatay;
       quadPointsData.row(0) << quadPointDataz;

        gsDebugVar(quadPointsData);
        if (get_readTime)
            t_read += participant.readTime();
        // forceMesh.patch(0).coefs() = forceControlPoints.transpose();

        // forceMesh.embed(3);
        assembler.assemble();
        F = assembler.rhs();
        // solve gismo timestep
        gsInfo << "Solving timestep " << time << "...\n";
        timeIntegrator->step(time,dt,U,V,A);
        solVector = U;

        gsDebugVar(solVector);

        gsInfo<<"Finished\n";

 

        // potentially adjust non-matching timestep sizes

        dt = std::min(dt,precice_dt);

        gsMultiPatch<> solution;
        gsVector<> displacements = U;
        solution = assembler.constructDisplacement(displacements); 
        gsDebugVar(displacements);

        gsMatrix<> centralPointDisp = solution.patch(0).eval(quadPoints);
        gsMatrix<> dispLeft(2,numQuadPt);
        dispLeft.setZero();
        dispLeft.row(0) = centralPointDisp.row(0);
        dispLeft.row(1) = centralPointDisp.row(1);

        gsMatrix<> dispRight(2,numQuadPt);
        dispRight.setZero();
        dispRight.row(0) = centralPointDisp.row(0);
        dispRight.row(1) = centralPointDisp.row(1);


        gsMatrix<> dispPoints(2, dispLeft.cols() + dispRight.cols());
        dispPoints << dispLeft, dispRight;

 
        // write heat fluxes to interface

        participant.writeData(SolidMesh,DisplacementData,quadPointIDs,dispPoints);

        if (get_writeTime)
            t_write +=participant.writeTime();

 

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
            // gsTimeIntegrator advances the time step                  
            // advance variables
            time += dt;
            timestep++;

            gsField<> solField(patches,solution);
            if (timestep % plotmod==0 && plot)
            {
                // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
                std::string fileName = dirname + "/solution" + util::to_string(timestep);
                gsWriteParaview<>(solField, fileName, 500);
                fileName = "solution" + util::to_string(timestep) + "0";
                collection.addTimestep(fileName,time,".vts");
            }

            solution.patch(0).eval_into(points,pointDataMatrix);
            otherDataMatrix<<time;
            writer.add(pointDataMatrix,otherDataMatrix);
        }

    }

    if (get_readTime)
    {
        gsInfo << "Read time: " << t_read << "\n";
    }

    if (get_writeTime)
    {
        gsInfo << "Write time: " << t_write << "\n";
    }

    if (plot)

    {

        collection.save();

    }

 

    delete timeIntegrator;

    return  EXIT_SUCCESS;

 

}