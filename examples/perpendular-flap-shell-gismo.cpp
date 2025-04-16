/**
 * @file perpendicular-flap-shell-gismo.cpp
 * @brief Calculate mechanical response of perpendicular flap shell element
 *
 * This program calculates the mechanical response of a perpendicular flap shell element,
 * including deformation, stress and strain. It uses the gsPreCICE library for structural 
 * analysis and gsThinShellAssembler for shell element analysis.
 *
 * For Kirchhoff–Love shells, a mapping is established between the mid-surface (P) 
 * and the top/bottom surfaces (Q). In this version, the mapping is updated every time step 
 * based on the current mid-surface geometry: The current unit normal is recalculated,
 * and the corresponding thickness offset is computed for updated 3D control points.
 *
 * This ensures that even for large deformations and rotations the generated volume 
 * (i.e. the top/bottom surfaces) remain perpendicular to the current mid-surface.
 *
 * (Author i
 */

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsLookupFunction.h>

#include <gsAssembler/gsExprEvaluator.h>

// #include <gsElasticity/gsMassAssembler.h>
// #include <gsElasticity/gsElasticityAssembler.h>

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
    ////////////////////////////////////////////////////////////////////////////
    // Parse command line options
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
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6: RK4
    std::string output("");
    std::string dirname = "./output";

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt("e", "degreeElevation",
               "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate);
    cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement loops", numRefine);
    cmd.addString("c", "config", "PreCICE config file", precice_config);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m", "method", "1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson", method);
    cmd.addSwitch("write", "Create a file with point data", write);
    cmd.addSwitch("readTime", "Get the read time", get_readTime);
    cmd.addSwitch("writeTime", "Get the write time", get_writeTime);
    cmd.addInt("n", "dof1", "Number of basis function in one direction", n);
    cmd.addInt("d", "degree", "Degree of a surface", degree);
    cmd.addString("o", "output", "Name of the output file.", output);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    GISMO_ASSERT(gsFileManager::fileExists(precice_config), "No precice config file has been defined");

    // Adjust values if necessary
    degree = math::max((index_t)(0), degree);
    n = math::max(n, degree + 1);
    m = math::max(m, degree + 1);

    gsInfo << "----------------------\n\n"
           << "n: " << n << "\n"
           << "degree: " << degree << "\n"
           << "output: " << output << "\n"
           << "----------------------\n\n";

    ////////////////////////////////////////////////////////////////////////////
    // Construct the 2D mid-surface geometry
    gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
    gsKnotVector<> kv2(0, 1, m - degree - 1, degree + 1);
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
    gsMatrix<> greville = basis.anchors();

    // Generate control points for the vertical flap.
    // Flap is positioned in the YZ plane with x = 0.
    gsMatrix<> coefs(greville.cols(), 3);
    for (index_t col = 0; col != greville.cols(); col++)
    {
        real_t x = greville(0, col);
        real_t y = greville(1, col);
        coefs(col, 0) = 0;    // x = 0 (fixed)
        coefs(col, 1) = y;    // y-coordinate
        coefs(col, 2) = x;    // Map parameter x to physical z for vertical orientation
    }
    gsTensorBSpline<2, real_t> surface(basis, coefs);
    gsExprEvaluator<> ev;
    ev.setIntegrationElements(basis);
    auto G = ev.getMap(surface);

    ////////////////////////////////////////////////////////////////////////////
    // Construct the 3D (volumetric) shell geometry.
    // Calculate the initial normal offset (dcoefs) and establish the initial mapping
    gsKnotVector<> kv3(0, 1, 0, 3);
    gsTensorBSplineBasis<3, real_t> vbasis(kv1, kv2, kv3);
    gsMatrix<> coefs3D(3 * greville.cols(), 3);

    // Define constant thickness, here 0.1 (can be modified as needed).
    gsFunctionExpr<> thickness("0.1", 2);
    gsMatrix<> T = thickness.eval(greville);
    gsMatrix<> dcoefs(greville.cols(), 3);  // Initial normal offset (fixed mapping method)
    dcoefs.setZero();
    gsMatrix<> N;
    for (index_t k = 0; k != greville.cols(); k++)
    {
        N = ev.eval(sn(G).normalized(), greville.col(k));
        dcoefs.row(k).transpose() = N * T(0, k);
    }

    // Construct 3D control points (Bottom, Middle, Top layers) using mid-surface control points and initial dcoefs
    coefs3D.block(0,                 0, greville.cols(), 3) = coefs - dcoefs; // Bottom layer
    coefs3D.block(greville.cols(),   0, greville.cols(), 3) = coefs;        // Middle layer
    coefs3D.block(2 * greville.cols(), 0, greville.cols(), 3) = coefs + dcoefs; // Top layer

    gsGeometry<>::uPtr volume = vbasis.makeGeometry(coefs3D);
    gsWriteParaview(*volume, "volume", 1000, true);

    ////////////////////////////////////////////////////////////////////////////
    // Construct mid-surface geometry for shell analysis.
    index_t componentSize = volume->coefs().rows() / 3;
    gsMatrix<> coefsSurf(componentSize, 3);
    gsMatrix<> dcoefsSurf(componentSize, 3);
    for (index_t k = 0; k != componentSize; k++)
    {
        coefsSurf.row(k) = (coefs3D.row(k) + coefs3D.row(k + 2 * componentSize)) / 2;  // Mid-surface points
        dcoefsSurf.row(k) = (coefs3D.row(k + 2 * componentSize) - coefs3D.row(k)) / 2; // Thickness direction vector
    }
    coefsSurf.transposeInPlace();
    dcoefsSurf.transposeInPlace();
    gsWriteParaviewPoints(coefsSurf, "coefsSurf");
    gsWriteParaviewPoints(dcoefsSurf, "dcoefsSurf");

    gsGeometry<>::uPtr surf_mid = basis.makeGeometry(coefsSurf.transpose());
    gsMultiPatch<> mid_surface_geom;
    mid_surface_geom.addPatch(std::move(surf_mid));
    gsWriteParaview(mid_surface_geom, "surf_mid", 1000, true);

    gsOptionList quadOptions = gsExprAssembler<>::defaultOptions();
    gsMatrix<> quad_shell_uv = gsQuadrature::getAllNodes(basis.basis(0), quadOptions);
    gsMatrix<> quad_shell_xy = mid_surface_geom.patch(0).eval(quad_shell_uv);

    // Embed 2D geometry into 3D for coupling.
    gsMultiPatch<> solutions;
    gsConstantFunction<> g(0., 0., 3);
    gsConstantFunction<> gravity(0., 0., 3);

    ////////////////////////////////////////////////////////////////////////////
    // PreCICE coupling initialization.
    std::string participantName = "Solid";
    gsPreCICE<real_t> participant(participantName, precice_config);
    gsDebugVar("Got here 1");

    std::string SolidMesh = "Solid-Mesh";
    std::string StressData = "Stress";
    std::string DisplacementData = "Displacement";

    std::vector<patchSide> couplingInterfaces(2);
    couplingInterfaces[0] = patchSide(0, boundary::front);
    couplingInterfaces[1] = patchSide(0, boundary::back);

    gsDebugVar("Got here 2");
    gsMatrix<> quad_uv = gsQuadrature::getAllNodes(vbasis.basis(0), quadOptions, couplingInterfaces);
    gsMatrix<> quad_xy = volume->eval(quad_uv);
    gsDebugVar(quad_xy);
    gsWriteParaviewPoints(quad_xy, "quadPointsAll");

    gsVector<index_t> quadPointIDs;
    participant.addMesh(SolidMesh, quad_xy, quadPointIDs);

    real_t precice_dt = participant.initialize();

    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 1);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, 0, 0, false, 2);
    bcInfo.addCondition(0, boundary::south, condition_type::clamped, 0, 0, false, -1);
    bcInfo.setGeoMap(mid_surface_geom);

    ////////////////////////////////////////////////////////////////////////////
    // Material properties and assembler setup.
    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;
    gsFunctionExpr<> E_modulus(std::to_string(E), 3);
    gsFunctionExpr<> PoissonRatio(std::to_string(nu), 3);
    gsFunctionExpr<> Density(std::to_string(rho), 3);
    gsFunctionExpr<> thickness_3d("0.1", 3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E_modulus;
    parameters[1] = &PoissonRatio;

    gsOptionList options;
    gsMatrix<> quadPointsData(quad_shell_xy.rows(), quad_shell_xy.cols());
    quadPointsData.setZero();
    gsLookupFunction<real_t> surfForce(quad_shell_xy, quadPointsData);

    gsMaterialMatrixBase<real_t>::uPtr materialMatrix;
    options.addInt("Material", "Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden", 0);
    options.addInt("Implementation", "Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral", 1);
    materialMatrix = getMaterialMatrix<3, real_t>(mid_surface_geom, thickness_3d, parameters, Density, options);

    gsMultiBasis<> bases(basis);
    gsThinShellAssemblerBase<real_t>* assembler;
    assembler = new gsThinShellAssembler<3, real_t, true>(mid_surface_geom, bases, bcInfo, surfForce, materialMatrix);
    gsOptionList assemblerOptions = options.wrapIntoGroup("Assembler");
    assembler->assemble();
    assembler->setOptions(assemblerOptions);

    // Assemble constant mass and stiffness matrices.
    assembler->assembleMass();
    gsSparseMatrix<> M = assembler->massMatrix();
    assembler->assemble();
    gsSparseMatrix<> K = assembler->matrix();

    gsFileManager::mkdir(dirname);
    gsParaviewCollection collection(dirname + "/solution");
    gsParaviewCollection collection_volume(dirname + "/deformed_volume");

    ////////////////////////////////////////////////////////////////////////////
    // Time stepping and dynamic solver setup.
    real_t t_read = 0, t_write = 0;
    real_t dt = precice_dt;

    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian =
        [&assembler, &solutions](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> &m)
        {
            assembler->constructSolution(x, solutions);
            ThinShellAssemblerStatus status = assembler->assembleMatrix(solutions);
            m = assembler->matrix();
            return true;
        };

    gsStructuralAnalysisOps<real_t>::TResidual_t Residual =
        [&assembler, &solutions](gsMatrix<real_t> const &x, real_t t, gsVector<real_t> &result)
        {
            assembler->constructSolution(x, solutions);
            ThinShellAssemblerStatus status = assembler->assembleVector(solutions);
            result = assembler->rhs();
            return true;
        };

    gsSparseMatrix<> C(assembler->numDofs(), assembler->numDofs());
    gsStructuralAnalysisOps<real_t>::Damping_t Damping =
        [&C](const gsVector<real_t>&, gsSparseMatrix<real_t>& m)
        { m = C; return true; };
    gsStructuralAnalysisOps<real_t>::Mass_t Mass =
        [&M](gsSparseMatrix<real_t>& m)
        { m = M; return true; };

    gsDynamicBase<real_t>* timeIntegrator;
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

    // Initialize displacement, velocity and acceleration (all zero).
    gsMatrix<> solVector = gsMatrix<real_t>::Zero(assembler->numDofs(), 1);
    gsVector<> F = assembler->rhs();

    gsVector<> F_checkpoint, U_checkpoint, V_checkpoint, A_checkpoint, U, V, A;
    F_checkpoint = F;
    U_checkpoint = U = gsVector<real_t>::Zero(assembler->numDofs(), 1);
    V_checkpoint = V = gsVector<real_t>::Zero(assembler->numDofs(), 1);
    A_checkpoint = A = gsVector<real_t>::Zero(assembler->numDofs(), 1);

    index_t timestep = 0, timestep_checkpoint = 0;
    real_t time = 0;

    ////////////////////////////////////////////////////////////////////////////
    // Main time stepping loop.
    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
        {
            U_checkpoint = U;
            V_checkpoint = V;
            A_checkpoint = A;
            gsInfo << "Checkpoint written:\n"
                   << "\t ||U|| = " << U.norm() << "\n"
                   << "\t ||V|| = " << V.norm() << "\n"
                   << "\t ||A|| = " << A.norm() << "\n";
            timestep_checkpoint = timestep;
        }
        
        // Read stress data from preCICE
        gsMatrix<> stressDataRead; // Declare matrix to hold read data
        participant.readData(SolidMesh, StressData, quadPointIDs, stressDataRead);
        // Assign the read data to quadPointsData. 
        // NOTE: Ensure dimensions/structure of stressDataRead match expectations for quadPointsData.
        // This might require reshaping or selecting specific components depending on the StressData format.
        quadPointsData = stressDataRead; 
        
        if (get_readTime)
            t_read += participant.readTime();

        assembler->assemble();
        F = assembler->rhs();

        gsInfo << "Solving timestep " << time << "...\n";
        timeIntegrator->step(time, dt, U, V, A);
        solVector = U;
        gsInfo << "Finished\n";

        // Adjust time step (sync with PreCICE)
        dt = std::min(dt, precice_dt);

        // Construct the deformed solution of the mid-surface from current displacement U (updated mid-surface geometry)
        gsMultiPatch<> solution = assembler->constructDisplacement(U);
        // Calculate mid-surface displacement at integration points and control points (for output only, not used for new normal calculation)
        gsMatrix<> MidPointDisp = solution.patch(0).eval(quad_shell_uv);
        gsMatrix<> MidPointDisp2 = solution.patch(0).eval(greville);

        // --- Update current mid-surface normal and calculate the latest thickness offset ---
        // Construct the updated mid-surface geometry to calculate the current normal.
        gsGeometry<>::uPtr updated_mid_geom = basis.makeGeometry(solution.patch(0).coefs());
        // Rebuild evaluator to calculate the updated unit normal
        gsExprEvaluator<> ev_shell_updated;
        ev_shell_updated.setIntegrationElements(basis);
        auto G_shell_updated = ev_shell_updated.getMap(*updated_mid_geom);

        gsMatrix<> updated_dcoefs(greville.cols(), 3);
        for (index_t k = 0; k < greville.cols(); ++k)
        {
            // Calculate the current unit normal at the parameter point greville.col(k)
            gsMatrix<> curNormal = ev_shell_updated.eval(sn(G_shell_updated).normalized(), greville.col(k));
            updated_dcoefs.row(k).transpose() = curNormal * T(0, k);
        }

        // --- Reconstruct 3D volume geometry using the updated mid-surface and thickness offset ---
        gsMatrix<> deformedcoefs3D(3 * greville.cols(), 3);
        // Bottom layer: Mid-surface control points minus the latest thickness offset
        deformedcoefs3D.block(0, 0, greville.cols(), 3) = updated_mid_geom->coefs() - updated_dcoefs;
        // Middle layer: Directly the mid-surface control points
        deformedcoefs3D.block(greville.cols(), 0, greville.cols(), 3) = updated_mid_geom->coefs();
        // Top layer: Mid-surface control points plus the latest thickness offset
        deformedcoefs3D.block(2 * greville.cols(), 0, greville.cols(), 3) = updated_mid_geom->coefs() + updated_dcoefs;

        // Construct the new volume geometry using the updated 3D control points
        gsGeometry<>::uPtr deformed_volume = vbasis.makeGeometry(deformedcoefs3D);
        gsWriteParaview(*deformed_volume, "deformed_volume", 1000, true);

        // Optional: Generate deformation field for visualization (calculate displacement = current volume geometry - initial volume geometry)
        gsMatrix<> volumeDisplacement = deformedcoefs3D;
        gsMatrix<> originalVolCoefs(3 * greville.cols(), 3);
        originalVolCoefs.block(0, 0, greville.cols(), 3) = coefs - dcoefs;
        originalVolCoefs.block(greville.cols(), 0, greville.cols(), 3) = coefs;
        originalVolCoefs.block(2 * greville.cols(), 0, greville.cols(), 3) = coefs + dcoefs;
        volumeDisplacement -= originalVolCoefs;
        gsGeometry<>::uPtr displacementGeo = vbasis.makeGeometry(volumeDisplacement);
        gsMultiPatch<> mp_displacement;
        mp_displacement.addPatch(displacementGeo->clone());

        // Prepare displacement data for coupling (send mid-surface displacement at integration points to the coupling interface)
        gsMatrix<> displacementData(3, quad_xy.cols());
        displacementData << MidPointDisp, MidPointDisp;
        participant.writeData(SolidMesh, DisplacementData, quadPointIDs, displacementData);

        if (get_writeTime)
            t_write += participant.writeTime();

        // Perform PreCICE coupling
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
        }
    } // end while coupling loop

    if (get_readTime)
        gsInfo << "Read time: " << t_read << "\n";
    if (get_writeTime)
        gsInfo << "Write time: " << t_write << "\n";

    if (plot)
    {
        collection.save();
        collection_volume.save();
    }

    delete timeIntegrator;
    return EXIT_SUCCESS;
}