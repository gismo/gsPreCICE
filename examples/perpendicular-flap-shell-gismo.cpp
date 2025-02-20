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
const double PI = 3.14159265;
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
    index_t n = 5;
    index_t m = 5;
    index_t degree = 3;
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6 RK4

    std::string output("");
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
    cmd.addInt   ("n", "dof1", "Number of basis function in one direction"  , n);

    cmd.addInt   ("d", "degree", "Degree of a surface", degree);
    cmd.addString("o", "output", "Name of the output file.", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]
    // GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

     // Adjust values to the minimum required
    degree = math::max( (index_t)(0), degree    );
    n      = math::max(n, degree + 1);
    m      = math::max(m, degree + 1);

    gsInfo << "----------------------\n\n"
              << "n: " << n << "\n\n"
              << "degree: " << degree << "\n\n"
              << "output: " << output << "\n\n"
              << "----------------------\n\n";

    // 1. construction of a knot vector for each direction
    gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
    gsKnotVector<> kv2(0, 1, m - degree - 1, degree + 1);

    // 2. construction of a basis
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);

    // 3. construction of a coefficients
    gsMatrix<> greville = basis.anchors();
    gsMatrix<> coefs (greville.cols(), 3);

    for (index_t col = 0; col != greville.cols(); col++)
    {
        real_t x = greville(0, col);
        real_t y = greville(1, col);

        coefs(col, 0) = x;
        coefs(col, 1) = y;
        coefs(col, 2) = 0;
    }

    // 4. putting basis and coefficients toghether
    gsTensorBSpline<2, real_t>  surface(basis, coefs);

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(basis);
    auto G = ev.getMap(surface);

    // Make the volumetric basis
    //----------------------------------Generating the initial volume----------------------------------
    gsKnotVector<> kv3(0, 1, 0, 3);
    gsTensorBSplineBasis<3, real_t> vbasis(kv1, kv2, kv3);
    gsDebugVar(vbasis);
    gsMatrix<> coefs3D(3*greville.cols(), 3);

    gsFunctionExpr<> thickness("0.1",2);
    gsMatrix<> N;
    gsMatrix<> T = thickness.eval(greville);
    gsMatrix<> dcoefs(greville.cols(), 3);
    dcoefs.setZero();//not needed
    for (index_t k=0; k!=greville.cols(); k++)
    {
        N = ev.eval(sn(G).normalized(),greville.col(k));
        dcoefs.row(k).transpose() = N*T(0,k);
    }
    coefs3D.block(0                ,0,greville.cols(),3) = coefs-dcoefs;
    coefs3D.block(greville.cols()  ,0,greville.cols(),3) = coefs;
    coefs3D.block(2*greville.cols(),0,greville.cols(),3) = coefs+dcoefs;

    coefs3D << coefs-dcoefs, coefs, coefs+dcoefs;
    gsGeometry<>::uPtr volume = vbasis.makeGeometry(coefs3D);
    gsWriteParaview(*volume, "volume",1000,true);
    //----------------------------------Generating the boundary volume----------------------------------

    //----------------------------------Generating the initial surface----------------------------------
    // direction to eliminate
    short_t dir = 2;
    // ASSUMES THICKNESS IS ALWAYS IN DIR2
    // for (index_t k=0; k!=volume.coefs().rows()/volume.basis().component(2).size(); k++)

    index_t componentSize = volume->coefs().rows()/3;
    gsMatrix<> coefsSurf(componentSize,3);
    coefsSurf.setZero();
    gsMatrix<> dcoefsSurf(componentSize,3);

    gsDebugVar(componentSize);

    gsDebugVar(coefsSurf.dim());
    gsDebugVar(coefs3D.dim());

    gsDebugVar(N);



    for (index_t k=0; k!=componentSize; k++)
    {
        coefsSurf.row(k) = (coefs3D.row(k)+coefs3D.row(k+2*componentSize))/2;
        dcoefsSurf.row(k) = (coefs3D.row(k+2*componentSize)-coefs3D.row(k))/2; // think about the normal
    }
    coefsSurf.transposeInPlace();
    dcoefsSurf.transposeInPlace();
    gsWriteParaviewPoints(coefsSurf, "coefsSurf");
    gsWriteParaviewPoints(dcoefsSurf, "dcoefsSurf");

    gsGeometry<>::uPtr surf_mid = basis.makeGeometry(coefsSurf.transpose());
    gsMultiPatch<> mid_surface_geom;
    mid_surface_geom.addPatch(std::move(surf_mid)); // Use std::move to transfer ownership. The coefs for mid_surface_geom is 3D.



    gsWriteParaview(mid_surface_geom, "surf_mid",1000,true);

    gsMatrix<> quad_shell_uv = gsQuadrature::getAllNodes(basis.basis(0),quadOptions);
    gsMatrix<> quad_shell_xy = mid_surface_geom.eval(quad_shell_uv); // Quadrature points in the physical domain



    // Embed the 2D geometry to 3D
    gsMultiPatch<> solutions;
    // patches.embed(3);
    // source function, rhs
    gsConstantFunction<> g(0.,0.,3);

    // source function, rhs
    gsConstantFunction<> gravity(0.,0.,3);

    //----------------------------------[! Setup expression assembler for normal vectors] ----------------------------------


    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;
    // real_t nu = 0.3;
    real_t mu = E / (2.0 * (1.0 + nu));
    real_t lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    /*
     * Initialize the preCICE participant
     *
     *
     */
    std::string participantName = "Solid";
    // gsPreCICE<real_t> participant(participantName, precice_config);

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
    std::string StressData       = "Force";
    std::string DisplacementData = "Displacement";

    std::vector<patchSide> couplingInterfaces(2);
    couplingInterfaces[0] = patchSide(0,boundary::front);
    // couplingInterfaces[1] = patchSide(0,boundary::north);
    couplingInterfaces[2] = patchSide(0,boundary::back);

    index_t numQuadPtFront = 64;
    // index_t numQuadPtNorth = 24;
    index_t numQuadPtBack = 64;


    // Step 1: SolidMesh
    // Get the quadrature nodes on the coupling interface
    gsOptionList quadOptions = gsExprAssembler<>::defaultOptions();


    //We only want the left boundary quadrature points
    gsMatrix<> quad_uv = gsQuadrature::getAllNodes(vbasis.basis(0),quadOptions,couplingInterfaces);
    gsDebugVar(quad_uv);

    gsMatrix<> quad_xy = volume->eval(quad_uv); // Quadrature points in the physical domain

    gsDebugVar(quad_xy.dim());

    gsWriteParaviewPoints(quad_xy, "quadPointsAll");

    gsVector<index_t> quadPointIDs;
    // participant.addMesh(SolidMesh,quad_xy,quadPointIDs); //Set the vertices to be datapoints (quadpoints in physical domain)

    // real_t precice_dt = participant.initialize();

    gsBoundaryConditions<> bcInfo;

    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 1);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 2);
          
    bcInfo.setGeoMap(mid_surface_geom);

    // Set up the material matrices
    gsFunctionExpr<> E_modulus(std::to_string(E),3);
    gsFunctionExpr<> PoissonRatio(std::to_string(nu),3);
    gsFunctionExpr<> Density(std::to_string(rho),3);

    gsFunctionExpr<> t("0.1",3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E_modulus;
    parameters[1] = &PoissonRatio;

    gsOptionList options;

    // gsMaterialMatrixBase<real_t>* materialMatrix;
    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3, real_t>(patches, t, parameters, Density, options);

    gsLookupFunction<real_t> surfForce(quad_shell_xy, quadPointsData);

    gsThinShellAssembler<3, real_t, false> assembler(patches, basis, bcInfo, surfForce, materialMatrix);
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

    // gsMatrix<> points(2,1);
    // points.col(0)<<0.5,1;

    // gsStructuralAnalysisOutput<real_t> writer("./output/pointData.csv",points);
    // writer.init({"x","y","z"},{"time"}); // point1 - x, point1 - y, point1 - z, time

    gsMatrix<> pointDataMatrix;
    gsMatrix<> otherDataMatrix(1,1);

    gsMatrix<> ForceData(3, quad_xy.cols());
    gsMatrix<> quadPointsNormalData(3, quadPointsAll.cols());
    // comForceData.setZero();

    // Evaluate the normal vector for the mid surface
    gsExprEvaluator<> ev_shell;
    ev_shell.setIntegrationElements(basis);
    auto G_shell = ev_shell.getMap(mid_surface_geom);
    

    gsMatrix<> N_shell(3, quad_xy.rows());

    gsMatrix<> T_shell = thickness.eval(quad_uv);

    gsMatrix<> deformed_thickness(greville.cols(), 3);

    for (index_t k=0; k!=quad_uv.rows(); k++)
    {
        N_shell.row(k) = ev.eval( sn(G_shell).normalized(), quad_uv.row(k) );
        deformed_thickness.row(k) = N_shell.row(k) * T_shell(0,k);
    }

    gsMatrix<> ForceDataShell(quad_shell_uv.rows(), 3);

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
        
        participant.readData(SolidMesh,StressData,quadPointIDs,ForceData);

        gsDebugVar(ForceData);

     
        //Is this really an average? (maybe need to consider normal vector)
        for (index_t i = 0; i < quad_uv.rows(); ++i) 
        {
            // Add the first half and second half values and divide by 2
            quadPointsData.row(i) = (ForceData.row(i) * N + ForceData.col(i + numQuadPtFront ) * N) / 2.0;
        }


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
        gsInfo<<"Finished\n";

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt,precice_dt);

        gsMultiPatch<> solution;
        gsVector<> displacements = U;
        solution = assembler.constructDisplacement(displacements);

        gsDebugVar(displacements);

        gsMatrix<> MidPointDisp = solution.patch(0).eval(quad_shell_uv);

        gsMatrix<> DisplacementData(quad_xy.rows(), 3);

        DisplacementData << MidPointDisp, MidPointDisp;

        // gsMatrix<> deformed_thickness(greville.cols(), 3);

        // for (index_t k=0; k!=quad_uv.rows(); k++)
        // {
        //     N_shell.row(k) = ev.eval( sn(G_shell).normalized(), quad_uv.row(k) );
        //     deformed_thickness.row(k) = N_shell.row(k) * T_shell(0,k);
        // }

        // gsMatrix<> deformedcoefs3D(3*greville.cols(), 3);

        // deformedcoefs3D.block(0                ,0,greville.cols(),3) = displacements - deformed_thickness;
        // deformedcoefs3D.block(greville.cols()  ,0,greville.cols(),3) = displacements;
        // deformedcoefs3D.block(2*greville.cols(),0,greville.cols(),3) = displacements + deformed_thickness;


        // gsMatrix<> dispPoints(2, dispLeft.cols() + dispRight.cols());
        // dispPoints << dispLeft, dispRight;


        // write heat fluxes to interface
        participant.writeData(SolidMesh,DisplacementData,quadPointIDs,DisplacementData);
        
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

