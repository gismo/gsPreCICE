/**
 * @file perpendicular-flap-shell-gismo.cpp
 * @brief Calculate mechanical response of perpendicular flap shell element
 * 
 * This program calculates the mechanical response of a perpendicular flap shell element,
 * including deformation, stress and strain.
 * It uses the gsPreCICE library for structural analysis and gsThinShellAssembler for shell element analysis.
 * 
 * @author: J. Li (2023 -..., TU Delft ), H. Verhelst (U. Pavia)
 * 
 */

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsLookupFunction.h>

#include <gsAssembler/gsExprEvaluator.h>

#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElasticityAssembler.h>


#ifdef gsKLShell_ENABLED
#include <gsKLShell/src/getMaterialMatrix.h>
#include <gsKLShell/src/gsThinShellAssembler.h>
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
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");
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
    // 1. Construction of knot vectors for each direction
    gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
    gsKnotVector<> kv2(0, 1, m - degree - 1, degree + 1);

    // 2. construction of a basis
    // 2. Construction of a tensor B-spline basis
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);

    // 3. construction of a coefficients
    // 3. Generate control points for vertical flap geometry
    gsMatrix<> greville = basis.anchors();
    gsMatrix<> coefs (greville.cols(), 3);

    for (index_t col = 0; col != greville.cols(); col++)
    {
        real_t x = greville(0, col);
        real_t y = greville(1, col);

        // Position the flap vertically in the YZ plane at x=0
        coefs(col, 0) = 0;         // Fixed at x=0
        coefs(col, 1) = y;         // Y-coordinate 
        coefs(col, 2) = x;         // Map parameter x to physical z for vertical orientation
    }

    // 4. putting basis and coefficients toghether
    // 4. Create the surface geometry
    gsTensorBSpline<2, real_t>  surface(basis, coefs);

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(basis);
    auto G = ev.getMap(surface);

    // Make the volumetric basis
    //----------------------------------Generating the initial volume----------------------------------
    // Create a volumetric shell with thickness
    gsKnotVector<> kv3(0, 1, 0, 3);
    gsTensorBSplineBasis<3, real_t> vbasis(kv1, kv2, kv3);
    gsMatrix<> coefs3D(3*greville.cols(), 3);

    // Define constant thickness for the shell
    gsFunctionExpr<> thickness("0.1",2);
    gsMatrix<> N;
    gsMatrix<> T = thickness.eval(greville);
    gsMatrix<> dcoefs(greville.cols(), 3);
    dcoefs.setZero();
    // Calculate normal vectors at each control point
    for (index_t k=0; k!=greville.cols(); k++)
    {
        N = ev.eval(sn(G).normalized(),greville.col(k));
        dcoefs.row(k).transpose() = N*T(0,k);
    }
    
    // 填充3D控制点（只保留一次，移除重复部分）
    // Fill 3D control points with bottom, middle and top layers
    coefs3D.block(0                ,0,greville.cols(),3) = coefs-dcoefs;  // Bottom layer
    coefs3D.block(greville.cols()  ,0,greville.cols(),3) = coefs;         // Middle layer
    coefs3D.block(2*greville.cols(),0,greville.cols(),3) = coefs+dcoefs;  // Top layer

    gsGeometry<>::uPtr volume = vbasis.makeGeometry(coefs3D);
    gsWriteParaview(*volume, "volume",1000,true);
    //----------------------------------Generating the boundary volume----------------------------------

    //----------------------------------Generating the initial surface----------------------------------
    // Extract the middle surface geometry for shell analysis
    short_t dir = 2;

    index_t componentSize = volume->coefs().rows()/3;
    gsMatrix<> coefsSurf(componentSize,3);
    coefsSurf.setZero();
    gsMatrix<> dcoefsSurf(componentSize,3);

    // Calculate mid-surface control points and thickness vectors
    for (index_t k=0; k!=componentSize; k++)
    {
        coefsSurf.row(k) = (coefs3D.row(k)+coefs3D.row(k+2*componentSize))/2;  // Mid-surface points
        dcoefsSurf.row(k) = (coefs3D.row(k+2*componentSize)-coefs3D.row(k))/2; // Thickness direction vectors
    }
    coefsSurf.transposeInPlace();
    dcoefsSurf.transposeInPlace();
    gsWriteParaviewPoints(coefsSurf, "coefsSurf");
    gsWriteParaviewPoints(dcoefsSurf, "dcoefsSurf");

    // Create mid-surface geometry for shell analysis
    gsGeometry<>::uPtr surf_mid = basis.makeGeometry(coefsSurf.transpose());
    gsMultiPatch<> mid_surface_geom;
    mid_surface_geom.addPatch(std::move(surf_mid));

    gsWriteParaview(mid_surface_geom, "surf_mid",1000,true);

    gsOptionList quadOptions = gsExprAssembler<>::defaultOptions();
    gsMatrix<> quad_shell_uv = gsQuadrature::getAllNodes(basis.basis(0),quadOptions);
    gsMatrix<> quad_shell_xy = mid_surface_geom.patch(0).eval(quad_shell_uv); // Quadrature points in the physical domain



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
    gsPreCICE<real_t> participant(participantName, precice_config);

    gsDebugVar("Got here 1");

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

    std::vector<patchSide> couplingInterfaces(2);
    couplingInterfaces[0] = patchSide(0,boundary::front);
    // couplingInterfaces[1] = patchSide(0,boundary::north);
    couplingInterfaces[1] = patchSide(0,boundary::back);

    index_t numQuadPtFront = 64;
    // index_t numQuadPtNorth = 24;
    index_t numQuadPtBack = 64;


    // Step 1: SolidMesh
    // Get the quadrature nodes on the coupling interface



    //We only want the left boundary quadrature points
    gsDebugVar("Got here 2");
    gsMatrix<> quad_uv = gsQuadrature::getAllNodes(vbasis.basis(0), quadOptions, couplingInterfaces);
    gsMatrix<> quad_xy = volume->eval(quad_uv); // Quadrature points in the physical domain

    gsDebugVar(quad_xy);

    gsWriteParaviewPoints(quad_xy, "quadPointsAll");

    gsVector<index_t> quadPointIDs;
    participant.addMesh(SolidMesh,quad_xy,quadPointIDs); //Set the vertices to be datapoints (quadpoints in physical domain)

    real_t precice_dt = participant.initialize();

    gsBoundaryConditions<> bcInfo;

    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 1);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 2);

    bcInfo.addCondition(0, boundary::south, condition_type::clamped, 0, 0, false, -1);
          
    bcInfo.setGeoMap(mid_surface_geom);

    // Set up the material matrices
    gsFunctionExpr<> E_modulus(std::to_string(E), 3);
    gsFunctionExpr<> PoissonRatio(std::to_string(nu), 3);
    gsFunctionExpr<> Density(std::to_string(rho), 3);

    gsFunctionExpr<> thickness_3d("0.1", 3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E_modulus;
    parameters[1] = &PoissonRatio;

    gsOptionList options;



    // Question: how to map the force information onto quad points as a surface force for gsThinShellAssembler?
    gsMatrix<> quadPointsData(quad_shell_xy.rows(), quad_shell_xy.cols()); // Ensure the dimensions match
    quadPointsData.setZero();

    gsDebugVar(quad_shell_xy);
    gsLookupFunction<real_t> surfForce(quad_shell_xy, quadPointsData);

    // gsMatrix<> displacementData = gsMatrix<>::Zero(3, comPt.rows());
    // displacementData.setRandom();

 

    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;
    options.addInt("Material","Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden",0);
    options.addInt("Implementation","Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral",1);
    materialMatrix = getMaterialMatrix<3, real_t>(mid_surface_geom, thickness_3d, parameters, Density, options);

    gsMultiBasis<> bases(basis);

    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true >(mid_surface_geom, bases, bcInfo, surfForce, materialMatrix);

    
    gsOptionList assemblerOptions = options.wrapIntoGroup("Assembler");

    assembler->assemble();
    assembler->setOptions(assemblerOptions);

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;

    // Compute the mass matrix (since it is constant over time)
    assembler->assembleMass();
    gsSparseMatrix<> M = assembler->massMatrix();
    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();

    gsFileManager::mkdir(dirname);
    gsParaviewCollection collection(dirname + "/solution");
    gsParaviewCollection collection_volume(dirname + "/deformed_volume");

    // Time step
    // real_t dt = 0.1;
    real_t t_read = 0;
    real_t t_write = 0;
    real_t dt = precice_dt;

    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = [&assembler, &solutions](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> &m)
    {
        assembler->constructSolution(x, solutions);
        ThinShellAssemblerStatus status;
        status = assembler->assembleMatrix(solutions);
        m = assembler->matrix();
        return true;
    };

    // Function for the Residual
    gsStructuralAnalysisOps<real_t>::TResidual_t Residual = [&assembler, &solutions](gsMatrix<real_t> const &x, real_t t, gsVector<real_t> &result)
    {
        assembler->constructSolution(x, solutions);
        ThinShellAssemblerStatus status;
        status = assembler->assembleVector(solutions);
        result = assembler->rhs();
        return true;
    };

    gsSparseMatrix<> C = gsSparseMatrix<>(assembler->numDofs(), assembler->numDofs());
    gsStructuralAnalysisOps<real_t>::Damping_t Damping = [&C](const gsVector<real_t> &, gsSparseMatrix<real_t> &m) { m = C; return true; };
    gsStructuralAnalysisOps<real_t>::Mass_t Mass = [&M](gsSparseMatrix<real_t> &m) { m = M; return true; };

    gsDynamicBase<real_t> *timeIntegrator;
    if (method == 1)
        timeIntegrator = new gsDynamicExplicitEuler<real_t, true>(Mass, Damping, Jacobian, Residual);
    else if (method == 2)
        timeIntegrator = new gsDynamicImplicitEuler<real_t, true>(Mass, Damping, Jacobian, Residual);
    else if (method == 3)
        timeIntegrator = new gsDynamicNewmark<real_t, true>(Mass, Damping, Jacobian, Residual);
    else if (method == 4)
        timeIntegrator = new gsDynamicBathe<real_t, true>(Mass, Damping, Jacobian, Residual);
    else if (method == 5)
    {
        timeIntegrator = new gsDynamicWilson<real_t, true>(Mass, Damping, Jacobian, Residual);
        timeIntegrator->options().setReal("gamma", 1.4);
    }
    else if (method == 6)
        timeIntegrator = new gsDynamicRK4<real_t, true>(Mass, Damping, Jacobian, Residual);
    else
        GISMO_ERROR("Method " << method << " not known");

    timeIntegrator->options().setReal("DT", dt);
    timeIntegrator->options().setReal("TolU", 1e-3);
    timeIntegrator->options().setSwitch("Verbose", true);

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    // RHS of the projection
    gsMatrix<> solVector;
    solVector.setZero(assembler->numDofs(), 1);

    // Assemble the RHS
    gsVector<> F = assembler->rhs();



    gsVector<> F_checkpoint, U_checkpoint, V_checkpoint, A_checkpoint, U, V, A;

    F_checkpoint = F;
    U_checkpoint = U = gsVector<real_t>::Zero(assembler->numDofs(), 1);
    V_checkpoint = V = gsVector<real_t>::Zero(assembler->numDofs(), 1);
    A_checkpoint = A = gsVector<real_t>::Zero(assembler->numDofs(), 1);

    real_t time = 0;

    if (plot)
    {
        gsMultiPatch<> solution;
        gsVector<> displacements = U;
        solution = assembler->constructDisplacement(displacements);
        solution.patch(0).coefs() -= mid_surface_geom.patch(0).coefs(); // assuming 1 patch here
        gsField<> solField(mid_surface_geom, solution);
        

        // Create consistent file naming format for initial step
        std::string fileName = dirname + "/solution" + util::to_string(timestep);
        gsWriteParaview<>(solField, fileName, 500);
        std::string fileNameSol = "solution" + util::to_string(timestep) + "0";
        collection.addTimestep(fileNameSol, time, ".vts");
        

        // Also create output for initial volume deformation
        std::string fileName2 = dirname + "/deformed_volume" + util::to_string(timestep);
        gsWriteParaview(*volume, fileName2, 1000, true);
        std::string fileNameVol = "deformed_volume" + util::to_string(timestep);
        collection_volume.addTimestep(fileNameVol, time, ".vts");
        

    }

    // gsMatrix<> points(2, 1);
    // points.col(0) << 0.5, 1;

    // gsStructuralAnalysisOutput<real_t> writer("./output/pointData.csv", points);
    // writer.init({"x", "y", "z"}, {"time"}); // point1 - x, point1 - y, point1 - z, time

    gsMatrix<> pointDataMatrix;
    gsMatrix<> otherDataMatrix(1, 1);

    gsMatrix<> ForceData(3, quad_xy.cols());
    // comForceData.setZero();

    // Evaluate the normal vector for the mid surface
    gsExprEvaluator<> ev_shell;
    ev_shell.setIntegrationElements(basis);
    auto G_shell = ev_shell.getMap(mid_surface_geom);

    gsMatrix<> N_shell(3, quad_xy.rows());

    gsMatrix<> T_shell = thickness_3d.eval(quad_uv);

    // gsMatrix<> deformed_thickness(greville.cols(), 3);

    // for (index_t k = 0; k != quad_uv.rows(); k++)
    // {
    //     N_shell.row(k) = ev.eval(sn(G_shell).normalized(), quad_uv.row(k));
    //     deformed_thickness.row(k) = N_shell.row(k) * T_shell(0, k);
    // }

    gsDebugVar("Got here 3");

    gsMatrix<> ForceDataShell(quad_shell_uv.rows(), 3);
    gsMatrix<> deformed_thickness(greville.cols(), 3);
    gsMatrix<> deformedcoefs3D(3 * greville.cols(), 3);

    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
        {
            U_checkpoint = U;
            V_checkpoint = V;
            A_checkpoint = A;

            gsInfo << "Checkpoint written:\n";
            gsInfo << "\t ||U|| = " << U.norm() << "\n";
            gsInfo << "\t ||V|| = " << V.norm() << "\n";
            gsInfo << "\t ||A|| = " << A.norm() << "\n";

            timestep_checkpoint = timestep;
        }

        // assembler.assemble();
        // F = assembler.rhs();

        gsDebugVar(bcInfo);

        participant.readData(SolidMesh, StressData, quadPointIDs, ForceData);



        // Is this really an average? (maybe need to consider normal vector)
        gsDebugVar(ForceData.rows());
        gsDebugVar(ForceData.cols());
        gsDebugVar(N.rows());
        gsDebugVar(N.cols());
        gsDebugVar(quadPointsData.rows());
        gsDebugVar(quadPointsData.cols());
        gsDebugVar(quad_uv.cols());
        
        for (index_t i = 0; i < quadPointsData.cols(); ++i)
        {
            // Add the first half and second half values and divide by 2
            quadPointsData.col(i) = ((ForceData.col(i)) + (ForceData.col(i + 64))) / 2.0;

        }

        gsDebugVar(quadPointsData);

        if (get_readTime)
            t_read += participant.readTime();
        // forceMesh.patch(0).coefs() = forceControlPoints.transpose();

        // forceMesh.embed(3);
        assembler->assemble();
        F = assembler->rhs();

        // solve gismo timestep
        // Solve the time step using the time integrator
        gsInfo << "Solving timestep " << time << "...\n";
        timeIntegrator->step(time, dt, U, V, A);
        solVector = U;
        gsInfo << "Finished\n";

        // 调整不匹配的时间步大小
        // Adjust time step size to match PreCICE requirements
        dt = std::min(dt, precice_dt);

        // Construct displacement solution
        gsMultiPatch<> solution;
        gsVector<> displacements = U;
        solution = assembler->constructDisplacement(displacements);

        // Evaluate displacements at quadrature points and control points
        gsMatrix<> MidPointDisp = solution.patch(0).eval(quad_shell_uv);
        gsMatrix<> MidPointDisp2 = solution.patch(0).eval(greville);

        // 计算可视化的变形厚度
        // Compute deformed thickness vectors for visualization
        T_shell = thickness_3d.eval(quad_uv);
        for (index_t k = 0; k != greville.cols(); k++)
        {
            N_shell = ev.eval(sn(G_shell).normalized(), greville.col(k)).transpose();
            deformed_thickness.row(k) = N_shell * T_shell(0,k);
        }


        // Get the shell control points
        gsMatrix<> shell_coefs = mid_surface_geom.patch(0).coefs();
        

        // Create deformed volume control points
        deformedcoefs3D.resize(3 * greville.cols(), 3);
        

        // For vertically positioned shell, normal vectors point in x direction
        // Bottom layer: middle layer - thickness vector + displacement
        deformedcoefs3D.block(0, 0, greville.cols(), 3) = shell_coefs + MidPointDisp2.transpose() - deformed_thickness;
        
    
        // Middle layer: middle layer + displacement
        deformedcoefs3D.block(greville.cols(), 0, greville.cols(), 3) = shell_coefs + MidPointDisp2.transpose();
        

        // Top layer: middle layer + thickness vector + displacement
        deformedcoefs3D.block(2 * greville.cols(), 0, greville.cols(), 3) = shell_coefs + MidPointDisp2.transpose() + deformed_thickness;


        // Create and visualize the deformed volume
        gsGeometry<>::uPtr deformed_volume = vbasis.makeGeometry(deformedcoefs3D);
        gsWriteParaview(*deformed_volume, "deformed_volume", 1000, true);
        

        // Add solution field to the deformed volume for visualization
        gsMultiPatch<> mp_volume;
        mp_volume.addPatch(deformed_volume->clone());
        

        // Create volume displacement field
        gsMatrix<> volumeDisplacement = deformedcoefs3D;

        // Get original volume control points
        gsMatrix<> originalVolCoefs(3*greville.cols(), 3);
        originalVolCoefs.block(0, 0, greville.cols(), 3) = coefs-dcoefs;
        originalVolCoefs.block(greville.cols(), 0, greville.cols(), 3) = coefs;
        originalVolCoefs.block(2*greville.cols(), 0, greville.cols(), 3) = coefs+dcoefs;
        
        // Calculate displacement (deformed volume - original volume)
        volumeDisplacement -= originalVolCoefs;

        // Create geometry representing displacement
        gsGeometry<>::uPtr displacementGeo = vbasis.makeGeometry(volumeDisplacement);
        gsMultiPatch<> mp_displacement;
        mp_displacement.addPatch(displacementGeo->clone());
        
        
        // Create and visualize displacement field
        // gsField<> displacementField(mp_volume, mp_displacement);
        
        // Prepare displacement data for coupling
        gsMatrix<> displacementData(3, quad_xy.cols());
        displacementData << MidPointDisp, MidPointDisp;

        
        // Write displacement data to interface for FSI coupling
        participant.writeData(SolidMesh, DisplacementData, quadPointIDs, displacementData);

        if (get_writeTime)
            t_write += participant.writeTime();

        // do the coupling
        precice_dt = participant.advance(dt);

        if (participant.requiresReadingCheckpoint())
        {
            U = U_checkpoint;
            V = V_checkpoint;
            A = A_checkpoint;
            timestep = timestep_checkpoint;
        }
        else
        {
            time += dt;
            timestep++;

            gsField<> solField(mid_surface_geom, solution);
            if (timestep % plotmod == 0 && plot)
            {
                std::string fileName = dirname + "/solution" + util::to_string(timestep);
                gsWriteParaview<>(solField, fileName, 500);
                std::string fileNameSol = "solution" + util::to_string(timestep) + "0";
                collection.addTimestep(fileNameSol, time, ".vts");

                std::string fileName2 = dirname + "/deformed_volume" + util::to_string(timestep);
                gsWriteParaview(*deformed_volume, fileName2, 1000, true);
                std::string fileNameVol = "deformed_volume" + util::to_string(timestep);
                collection_volume.addTimestep(fileNameVol, time, ".vts");
            }

            otherDataMatrix << time;
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
        collection_volume.save();
    }

    delete timeIntegrator;
    return EXIT_SUCCESS;

}
