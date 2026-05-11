/**
  @file perpendicular-flap-IGA-vertex-fluid.cpp
  @brief Fluid participant for preCICE FSI coupling with vertex data exchange.
*/

#include <gismo.h>

#include <gsIncompressibleFlow/src/gsINSSolver.h>
#include <gsIncompressibleFlow/src/gsINSSolverALE.h>
#include <gsIncompressibleFlow/src/gsINSMeshDeformer.h>

#include <gsPreCICE/gsPreCICE.h>
#include <map>
#include <algorithm>
#include <fstream>
#include <vector>
#include <cmath>
#include <memory>

using namespace gismo;

namespace
{
/// Replace non-finite entries (NaN/Inf) with zero and log a short summary.
static bool sanitizeFinite(gsMatrix<real_t>& M, const std::string& label, index_t step, index_t maxPrint = 5)
{
    index_t bad = 0;
    index_t printed = 0;
    real_t maxAbs = 0.0;
    for (index_t i = 0; i < M.rows(); ++i)
    {
        for (index_t j = 0; j < M.cols(); ++j)
        {
            const real_t v = M(i,j);
            if (!std::isfinite(v))
            {
                ++bad;
                if (printed < maxPrint)
                {
                    gsWarn << "[FSI-Fluid|Step " << step << "] Non-finite in " << label
                           << " at (" << i << "," << j << ") = " << v << " -> 0\n";
                    ++printed;
                }
                M(i,j) = 0.0;
            }
            else
            {
                maxAbs = std::max(maxAbs, std::abs(v));
            }
        }
    }
    if (bad > 0)
        gsWarn << "[FSI-Fluid|Step " << step << "] Sanitized " << bad << " non-finite values in "
               << label << " (max|finite|=" << maxAbs << ").\n";
    return bad == 0;
}

static gsMatrix<real_t> buildALEMeshDisplacement(
    const gsMultiPatch<>& curGeo,
    const gsMultiPatch<>& refGeo,
    const std::vector<index_t>& deformablePatches,
    const gsINSAssemblerUnsteadyALE<real_t, RowMajor>& aleAssembler)
{
    const index_t dim = curGeo.geoDim();
    const gsDofMapper& mapper = aleAssembler.getMappers()[0];
    const gsMultiBasis<>& velBasis = aleAssembler.getBases()[0];
    const index_t udofs = aleAssembler.getUdofs();

    gsMatrix<real_t> meshDisp(dim * udofs + aleAssembler.getPdofs(), 1);
    meshDisp.setZero();

    for (index_t p : deformablePatches)
    {
        if (p >= curGeo.nPatches() || p >= refGeo.nPatches()
            || p >= static_cast<index_t>(velBasis.nPieces()))
            continue;

        const gsBasis<>& basis = velBasis.piece(p);
        gsMatrix<real_t> anchors = basis.anchors();
        gsMatrix<real_t> curVals, refVals;
        curGeo.patch(p).eval_into(anchors, curVals);
        refGeo.patch(p).eval_into(anchors, refVals);

        for (index_t i = 0; i < anchors.cols(); ++i)
        {
            if (!mapper.is_free(i, p))
                continue;

            const index_t base = mapper.index(i, p);
            for (index_t d = 0; d < dim; ++d)
                meshDisp(base + d * udofs, 0) = curVals(d, i) - refVals(d, i);
        }
    }

    return meshDisp;
}

/// @brief Construct perpendicular flap fluid domain (8-patch)
/// Channel: x \in [-3,3], y \in [0,4] assembled from 8 patches:
/// 0-2 central 3-patch flap region in [-1,1]×[0,2]; plus extensions:
/// 3: [-3,-1]×[0,2], 4: [1,3]×[0,2], 5: [-3,-1]×[2,4], 6: [-1,1]×[2,4], 7: [1,3]×[2,4]
/// Flap: x \in [-0.05,0.05], y \in [0,1]
static gsMultiPatch<> makePerpendicularFlapFluidDomain()
{
    const real_t xl = -1.0, xr = 1.0;
    const real_t yb =  0.0, yt = 2.0;
    const real_t fxL = -0.05, fxR = 0.05; // flap thickness 0.1
    const real_t fyT =  1.0;              // flap height 1.0

    gsMultiPatch<> mp;
    gsKnotVector<> uKnots(0, 1, 0, 2), vKnots(0, 1, 0, 2); // Bilinear patches

    // Patch 0: Left trapezoid
    gsMatrix<> coefs0(4, 2);
    coefs0 << xl,  yb, fxL, yb, xl,  yt, fxL, fyT;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs0));

    // Patch 1: Right trapezoid
    gsMatrix<> coefs1(4, 2);
    coefs1 << fxR, yb, xr,  yb, fxR, fyT, xr,  yt;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs1));

    // Patch 2: Top trapezoid
    gsMatrix<> coefs2(4, 2);
    coefs2 << fxL, fyT, fxR, fyT, xl,  yt, xr,  yt;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs2));

    // Patch 3: left extension bottom [-3,-1]x[0,2]
    gsMatrix<> coefs3(4, 2);
    coefs3 << -3.0, yb,
              -1.0, yb,
              -3.0, yt,
              -1.0, yt;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs3));

    // Patch 4: right extension bottom [1,3]x[0,2]
    gsMatrix<> coefs4(4, 2);
    coefs4 << 1.0, yb,
              3.0, yb,
              1.0, yt,
              3.0, yt;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs4));

    // Patch 5: left extension top [-3,-1]x[2,4]
    gsMatrix<> coefs5(4, 2);
    coefs5 << -3.0, yt,
              -1.0, yt,
              -3.0, 4.0,
              -1.0, 4.0;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs5));

    // Patch 6: center top [-1,1]x[2,4]
    gsMatrix<> coefs6(4, 2);
    coefs6 << xl, yt,
              xr, yt,
              xl, 4.0,
              xr, 4.0;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs6));

    // Patch 7: right extension top [1,3]x[2,4]
    gsMatrix<> coefs7(4, 2);
    coefs7 << 1.0, yt,
              3.0, yt,
              1.0, 4.0,
              3.0, 4.0;
    mp.addPatch(gsTensorBSpline<2>(uKnots, vKnots, coefs7));

    mp.computeTopology();
    return mp;
}

/// @brief Build quadrature points on the FSI interface for vertex-based coupling
static void buildInterfaceQuadraturePoints(
    const gsMultiPatch<>& patches,
    const gsMultiBasis<>& basis,
    const std::vector<patchSide>& fsiInterfaceSides,
    gsMatrix<real_t>& quadPoints_xy,      // physical coordinates (2 x N)
    gsMatrix<real_t>& quadPoints_uv,      // parametric coordinates (2 x N)
    gsVector<real_t>& quadWeights,         // integration weights: quad_weight * |boundary_Jacobian|
    std::vector<std::pair<index_t, boxSide>>& quadSideLocations,
    std::vector<index_t>& patchOffsets)
{
    gsOptionList quadOptions = gsExprAssembler<>::defaultOptions();

    std::vector<gsMatrix<real_t>> patchQuad_uv;
    std::vector<gsMatrix<real_t>> patchQuad_xy;
    std::vector<gsVector<real_t>> patchQuad_wt;

    patchOffsets.resize(fsiInterfaceSides.size() + 1);
    patchOffsets[0] = 0;

    for (size_t i = 0; i < fsiInterfaceSides.size(); ++i)
    {
        const patchSide& ps = fsiInterfaceSides[i];
        const short_t fixDir = ps.side().direction();

        // Get quadrature rule for this boundary
        typename gsQuadRule<real_t>::uPtr QuRule =
            gsQuadrature::getPtr(*basis.basis(ps.patch).domain(), quadOptions, fixDir);

        // Count points
        typename gsBasis<real_t>::domainIter domIt  = basis.basis(ps.patch).domain()->beginBdr(ps.side());
        typename gsBasis<real_t>::domainIter domEnd = basis.basis(ps.patch).domain()->endBdr(ps.side());
        index_t nPts = 0;
        for (; domIt < domEnd; ++domIt)
            nPts += QuRule->numNodes();

        gsMatrix<real_t> side_uv(basis.basis(ps.patch).dim(), nPts);
        gsVector<real_t> side_wt(nPts);

        // Collect nodes and weights with boundary Jacobian
        domIt = basis.basis(ps.patch).domain()->beginBdr(ps.side());
        index_t off = 0;
        gsMatrix<real_t> nodes;
        gsVector<real_t> weights;
        for (; domIt < domEnd; ++domIt)
        {
            QuRule = gsQuadrature::getPtr(*basis.basis(ps.patch).domain(), quadOptions, fixDir);
            QuRule->mapTo(domIt.lowerCorner(), domIt.upperCorner(), nodes, weights);
            const index_t nn = QuRule->numNodes();
            side_uv.block(0, off, nodes.rows(), nn) = nodes;

            for (index_t k = 0; k < nn; ++k)
            {
                gsMatrix<real_t> jac;
                patches.patch(ps.patch).jacobian_into(nodes.col(k), jac);
                const short_t tangDir = 1 - fixDir; // 2D only
                const real_t bndJac = jac.col(tangDir).norm();
                side_wt(off + k) = weights(k) * bndJac;
            }
            off += nn;
        }

        gsMatrix<real_t> side_xy = patches.patch(ps.patch).eval(side_uv);

        patchQuad_uv.push_back(side_uv);
        patchQuad_xy.push_back(side_xy);
        patchQuad_wt.push_back(side_wt);
        patchOffsets[i + 1] = patchOffsets[i] + nPts;
    }

    index_t totalPoints = patchOffsets.back();
    quadPoints_xy.resize(patches.geoDim(), totalPoints);
    quadPoints_uv.resize(patches.domainDim(), totalPoints);
    quadWeights.resize(totalPoints);
    quadSideLocations.resize(totalPoints);

    for (size_t i = 0; i < fsiInterfaceSides.size(); ++i)
    {
        const patchSide& ps = fsiInterfaceSides[i];
        index_t offset = patchOffsets[i];
        index_t nPts = patchQuad_uv[i].cols();

        quadPoints_uv.block(0, offset, patches.domainDim(), nPts) = patchQuad_uv[i];
        quadPoints_xy.block(0, offset, patches.geoDim(), nPts) = patchQuad_xy[i];
        quadWeights.segment(offset, nPts) = patchQuad_wt[i];

        for (index_t j = 0; j < nPts; ++j)
            quadSideLocations[offset + j] = std::make_pair(ps.patch, ps.side());
    }
}

} // namespace

int main(int argc, char* argv[])
{
    //==========================================================================
    // Command-line Options
    //==========================================================================

    std::string preciceConfig;
    std::string geomPath = "auto";
    real_t viscosity = 1.0;
    real_t density   = 1.0;
    real_t totalTime = 5.0;
    real_t dt = 0.01;

    int numRefine = 0;
    int degreeElevate = 0;

    bool plot = false;
    real_t plotEvery = 1.0;
    index_t plotEverySteps = 10;
    bool fixPressure = true;
    bool usePiolaTraction = true;
    bool useCauchyTraction = false;

    // ALE method: 0=HE, 1=IHE, 2=LE, 6=BHE
    int aleMethod = mesh_deform_method::IHE;
    real_t alePoisson = 0.3;
    bool aleCheck = true;

    gsCmdLine cmd("Fluid participant (ALE, vertex-based) for perpendicular-flap FSI via preCICE");
    cmd.addString("c", "config",  "preCICE configuration file", preciceConfig);
    cmd.addString("g", "geom",    "Fluid geometry XML", geomPath);
    cmd.addReal  ("v", "nu",      "Kinematic viscosity", viscosity);
    cmd.addReal  ("R", "density", "Fluid density", density);
    cmd.addReal  ("t", "time",    "Total simulation time", totalTime);
    cmd.addReal  ("s", "step",    "Time step", dt);
    cmd.addInt   ("r", "refine",  "Uniform h-refinements", numRefine);
    cmd.addInt   ("p", "degree",  "Degree elevation", degreeElevate);
    cmd.addSwitch("plot",       "Write Paraview outputs", plot);
    cmd.addReal  ("E", "plot-every", "Plot interval in seconds", plotEvery);
    cmd.addInt   ("S", "plot-steps", "Plot every N steps", plotEverySteps);
    cmd.addSwitch("fixPressure", "Fix pressure gauge by setting p=0 at one corner", fixPressure);
    cmd.addSwitch("piolaTraction", "Apply Piola pull-back traction", usePiolaTraction);
    cmd.addSwitch("cauchyTraction", "Disable Piola pull-back (small-deformation)", useCauchyTraction);
    cmd.addInt   ("a", "aleMethod", "ALE mesh method (0=HE, 1=IHE, 2=LE, 6=BHE)", aleMethod);
    cmd.addReal  ("M", "alePoisson", "Poisson's ratio for ALE linear elasticity", alePoisson);
    cmd.addSwitch("C", "aleCheck", "Check ALE mesh quality", aleCheck);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if (useCauchyTraction)
        usePiolaTraction = false;

    const auto pdeMethod = static_cast<mesh_deform_method::method>(aleMethod);
    if (!gsINSMeshDeformer<real_t>::isSupportedMethod(pdeMethod))
    {
        gsWarn << "Invalid ALE method " << aleMethod
               << ". Valid values are 0=HE, 1=IHE, 2=LE, and 6=BHE.\n";
        return -1;
    }

    // Try to locate a default preCICE config
    if (preciceConfig.empty())
    {
        auto fileExists = [](const std::string &p)->bool{
            std::ifstream ifs(p.c_str());
            return static_cast<bool>(ifs);
        };

        std::vector<std::string> candidates = {
            "../optional/gsPreCICE/examples/perpendicular-flap/precice-config-vertex.xml",
            "../../optional/gsPreCICE/examples/perpendicular-flap/precice-config-vertex.xml",
            "../../../optional/gsPreCICE/examples/perpendicular-flap/precice-config-vertex.xml",
            "optional/gsPreCICE/examples/perpendicular-flap/precice-config-vertex.xml",
            "precice-config-vertex.xml",
            // Fallback to existing vertex config
            "../optional/gsPreCICE/examples/perpendicular-flap/precice-config.xml",
            "../../optional/gsPreCICE/examples/perpendicular-flap/precice-config.xml",
            "precice-config.xml"
        };

        for (const auto& c : candidates)
        {
            if (fileExists(c)) { preciceConfig = c; break; }
        }
    }

    if (preciceConfig.empty())
    {
        gsWarn << "No preCICE config passed. Use -c <precice-config.xml>.\n";
        return -1;
    }

    gsInfo << "[FSI-Fluid-Vertex] Geometry: " << geomPath << "\n";
    gsInfo << "[FSI-Fluid-Vertex] preCICE config: " << preciceConfig << "\n";

    //==========================================================================
    // Geometry Setup
    //==========================================================================

    gsMultiPatch<> fluidDomain;
    if (geomPath == "auto")
    {
        gsInfo << "[FSI-Fluid-Vertex] Building perpendicular-flap fluid domain (8 patches)...\n";
        fluidDomain = makePerpendicularFlapFluidDomain();
    }
    else
    {
        GISMO_ASSERT(gsFileManager::fileExists(geomPath), "Geometry file not found: " << geomPath);
        gsReadFile<>(geomPath, fluidDomain);
    }

    for (int d = 0; d < degreeElevate; ++d)
        fluidDomain.degreeElevate(1);

    gsMultiBasis<> basisVelocity(fluidDomain);
    gsMultiBasis<> basisPressure(fluidDomain);

    // Taylor–Hood
    basisVelocity.setDegree(basisVelocity.maxCwiseDegree() + 1);
    for (int r = 0; r < numRefine; ++r)
    {
        basisVelocity.uniformRefine();
        basisPressure.uniformRefine();
    }

    // FSI interface sides
    std::vector<patchSide> fsiInterfaceSides;
    fsiInterfaceSides.emplace_back(0, boundary::east);
    fsiInterfaceSides.emplace_back(2, boundary::south);
    fsiInterfaceSides.emplace_back(1, boundary::west);

    // Deformable sub-domain
    const std::vector<index_t> deformablePatches = {0, 1, 2};
    gsMultiPatch<> deformableDomain;
    for (index_t pid : deformablePatches)
        deformableDomain.addPatch(fluidDomain.patch(pid).clone());
    deformableDomain.computeTopology();

    // Mesh deformer setup
    std::unique_ptr<gsINSMeshDeformer<real_t>> meshDeformerPtr;
    {
        meshDeformerPtr = std::make_unique<gsINSMeshDeformer<real_t>>(
            deformableDomain, pdeMethod);
        meshDeformerPtr->setMovingBoundaries(fsiInterfaceSides);
        meshDeformerPtr->options().setReal("PoissonsRatio", alePoisson);
        meshDeformerPtr->options().setSwitch("Check", aleCheck);
        meshDeformerPtr->initialize();
        gsInfo << "[FSI-Fluid-Vertex] Using gsINSMeshDeformer (method=" << pdeMethod << ").\n";
    }

    // Boundary conditions
    gsBoundaryConditions<> bcInfo;
    gsConstantFunction<> inletVec(10.0, 0.0, 2);
    bcInfo.addCondition(3, boundary::west, condition_type::dirichlet, &inletVec, 0);
    bcInfo.addCondition(5, boundary::west, condition_type::dirichlet, &inletVec, 0);

    gsConstantFunction<> zeroVel(0.0, 0.0, 2);
    bcInfo.addCondition(3, boundary::south, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(4, boundary::south, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(5, boundary::north, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(6, boundary::north, condition_type::dirichlet, &zeroVel, 0);
    bcInfo.addCondition(7, boundary::north, condition_type::dirichlet, &zeroVel, 0);

    if (fixPressure)
        bcInfo.addCornerValue(boundary::southwest, 0.0, 3, 1);

    // FSI Interface BC: u = w (mesh velocity)
    std::vector<std::shared_ptr<gsFSIMeshVelocityFunction<>>> ifaceMeshVelFuncs;
    ifaceMeshVelFuncs.reserve(fsiInterfaceSides.size());
    for (const auto& ps : fsiInterfaceSides)
    {
        auto f = std::make_shared<gsFSIMeshVelocityFunction<>>(
            &fluidDomain.patch(ps.patch), &fluidDomain.patch(ps.patch), dt);
        ifaceMeshVelFuncs.push_back(f);
        bcInfo.addCondition(ps.patch, ps.side(), condition_type::dirichlet, f.get(), 0);
    }

    // PDE & solver
    gsConstantFunction<> fZero(0.0, 0.0, 2);
    gsNavStokesPde<real_t> nsPde(fluidDomain, bcInfo, &fZero, viscosity);

    std::vector<gsMultiBasis<>> spaces;
    spaces.push_back(basisVelocity);
    spaces.push_back(basisPressure);

    gsFlowSolverParams<real_t> params(nsPde, spaces);
    params.options().setReal("timeStep", dt);
    params.options().setInt ("nonlin.maxIt", 50);
    params.options().setReal("nonlin.tol", 1e-6);
    params.options().setInt ("lin.maxIt", 250);
    params.options().setReal("lin.tol", 1e-8);

    gsINSSolverUnsteadyALE<> solver(memory::make_shared_not_owned(&params));
    solver.setFSIInterfaceSides(fsiInterfaceSides);
    solver.initialize();
    solver.setALEActive(true);

    gsMultiPatch<> refFluidGeometry = solver.getAssembler()->getPatches();
    gsMultiPatch<> prevPatches = solver.getAssembler()->getPatches();

    //==========================================================================
    // Build VERTEX-BASED coupling mesh (quadrature points)
    //==========================================================================
    gsPreCICE<real_t> participant("Fluid", preciceConfig);

    const std::string FluidMesh = "Fluid-Mesh";
    const std::string StressData = "Stress";
    const std::string DisplData = "Displacement";

    // Build quadrature points on FSI interface
    gsMatrix<real_t> quadPoints_xy, quadPoints_uv;
    gsVector<real_t> quadWeights;  // integration weights: quad_weight * |boundary_Jacobian|
    std::vector<std::pair<index_t, boxSide>> quadSideLocations;
    std::vector<index_t> patchOffsets;

    buildInterfaceQuadraturePoints(fluidDomain, basisVelocity, fsiInterfaceSides,
                                   quadPoints_xy, quadPoints_uv, quadWeights, quadSideLocations, patchOffsets);

    gsInfo << "[FSI-Fluid-Vertex] Built " << quadPoints_xy.cols() << " quadrature points on FSI interface\n";
    gsInfo << "[FSI-Fluid-Vertex] Total weight (interface length): " << quadWeights.sum() << "\n";
    for (size_t i = 0; i < fsiInterfaceSides.size(); ++i)
    {
        gsInfo << "  - Side " << i << ": " << (patchOffsets[i+1] - patchOffsets[i]) << " points\n";
    }

    // Add mesh to preCICE
    gsVector<index_t> fluidMeshIDs;
    participant.addMesh(FluidMesh, quadPoints_xy, fluidMeshIDs);
    gsInfo << "[FSI-Fluid-Vertex] Added Fluid-Mesh with " << fluidMeshIDs.size() << " vertices\n";

    // Check if initial data is required
    const bool needsInitialData = participant.requiresInitialData();

    // Initialize preCICE
    real_t precice_dt = participant.initialize();
    gsInfo << "[FSI-Fluid-Vertex] preCICE initialized, dt_max = " << precice_dt << "\n";

    // Write initial data
    if (needsInitialData)
    {
        gsInfo << "[FSI-Fluid-Vertex] Writing initial stress (zero)...\n";
        gsMatrix<real_t> stress0(2, fluidMeshIDs.size());
        stress0.setZero();
        participant.writeData(FluidMesh, StressData, fluidMeshIDs, stress0);
    }

    //==========================================================================
    // Time loop
    //==========================================================================
    gsInfo << "[FSI-Fluid-Vertex] Entering main coupling loop...\n";

    gsMatrix<> Uchk, Pchk;
    gsINSSolverUnsteadyALE<>::ALEState solverALEChk;
    // Checkpoint: store coefficients only (not entire MultiPatch objects)
    // to avoid invalidating meshDeformer's internal references
    std::vector<gsMatrix<real_t>> prevPatchesCoefsChk(fluidDomain.nPatches());
    std::vector<gsMatrix<real_t>> fluidDomainCoefsChk(fluidDomain.nPatches());
    std::vector<gsMatrix<real_t>> deformableDomainCoefsChk(deformablePatches.size());
    real_t time = 0.0;
    index_t step = 0;

    gsParaviewCollection pvU("fsi_flap_fluid_vertex_U");
    gsParaviewCollection pvP("fsi_flap_fluid_vertex_P");

    while (participant.isCouplingOngoing())
    {
        // Checkpoint (write) - save coefficients only to preserve meshDeformer references
        if (participant.requiresWritingCheckpoint())
        {
            Uchk = solver.solutionCoefs(0);
            Pchk = solver.solutionCoefs(1);
            solverALEChk = solver.saveALEState();
            for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
            {
                prevPatchesCoefsChk[p] = prevPatches.patch(p).coefs();
                fluidDomainCoefsChk[p] = fluidDomain.patch(p).coefs();
            }
            for (size_t i = 0; i < deformablePatches.size(); ++i)
                deformableDomainCoefsChk[i] = deformableDomain.patch(i).coefs();
            gsInfo << "[FSI-Fluid-Vertex|Step " << step << "] Checkpoint written\n";
        }

        // Read displacement at quadrature points (via RBF mapping from Solid)
        gsMatrix<> disp_quad;
        participant.readData(FluidMesh, DisplData, fluidMeshIDs, disp_quad);

        if (step % 20 == 0)
            gsInfo << "[FSI-Fluid-Vertex|Step " << step << "] Read displacement: "
                   << disp_quad.rows() << "x" << disp_quad.cols() << "\n";

        // Update quadrature point positions (reference + displacement)
        gsMatrix<real_t> quadPoints_xy_current = quadPoints_xy + disp_quad;

        // Build boundary displacement for mesh deformer
        gsMultiPatch<> boundaryDispDeformer;
        std::map<index_t,index_t> patchToLocal;
        for (size_t i = 0; i < deformablePatches.size(); ++i)
        {
            const index_t pid = deformablePatches[i];
            patchToLocal[pid] = static_cast<index_t>(i);
            boundaryDispDeformer.addPatch(deformableDomain.patch(i).clone());
            boundaryDispDeformer.patch(i).coefs().setZero();
        }
        boundaryDispDeformer.computeTopology();

        // Interpolate displacement from quadrature points to boundary basis coefficients
        for (size_t i = 0; i < fsiInterfaceSides.size(); ++i)
        {
            const patchSide& ps = fsiInterfaceSides[i];
            auto it = patchToLocal.find(ps.patch);
            if (it == patchToLocal.end())
                continue;
            const index_t lp = it->second;

            // Get quadrature point data for this side
            index_t startIdx = patchOffsets[i];
            index_t endIdx = patchOffsets[i + 1];
            index_t nPts = endIdx - startIdx;

            gsMatrix<real_t> sideQuadUV = quadPoints_uv.block(0, startIdx, 2, nPts);
            gsMatrix<real_t> sideDisp = disp_quad.block(0, startIdx, 2, nPts);

            // Interpolate to boundary basis using least squares
            typename gsBasis<>::uPtr bndBasis = boundaryDispDeformer.patch(lp).basis().boundaryBasis(ps.side());
            gsMatrix<> anchors = bndBasis->anchors();

            // Convert 2D parametric coords to 1D boundary parameter
            gsMatrix<real_t> quadParam1D(1, nPts);
            for (index_t j = 0; j < nPts; ++j)
            {
                if (ps.side() == boundary::west || ps.side() == boundary::east)
                    quadParam1D(0, j) = sideQuadUV(1, j);
                else
                    quadParam1D(0, j) = sideQuadUV(0, j);
            }

            // Evaluate basis at quadrature points and solve least squares
            gsMatrix<real_t> Nvals;
            bndBasis->eval_into(quadParam1D, Nvals); // nBasis x nPts

            // Solve N^T N c = N^T d for each component
            gsMatrix<real_t> NtN = Nvals * Nvals.transpose();
            gsMatrix<real_t> Ntd = Nvals * sideDisp.transpose(); // nBasis x 2

            // Regularize if needed
            NtN += 1e-10 * gsMatrix<real_t>::Identity(NtN.rows(), NtN.cols());
            gsMatrix<real_t> coeffs = NtN.ldlt().solve(Ntd); // nBasis x 2

            // Scatter into boundary coefficients
            gsMatrix<index_t> bIdx = boundaryDispDeformer.patch(lp).basis().boundary(ps.side());
            for (index_t j = 0; j < bIdx.rows() && j < coeffs.rows(); ++j)
            {
                const index_t idx = bIdx(j,0);
                boundaryDispDeformer.patch(lp).coefs()(idx, 0) = coeffs(j, 0);
                if (coeffs.cols() > 1)
                    boundaryDispDeformer.patch(lp).coefs()(idx, 1) = coeffs(j, 1);
            }
        }

        // Deform ALE mesh
        real_t maxDispMag = (disp_quad.size() > 0) ? disp_quad.colwise().norm().maxCoeff() : 0.0;
        const bool hasMotion = (disp_quad.size() > 0) && (maxDispMag > 1e-12);

        if (hasMotion)
        {
            // Reset to reference configuration before deformation
            for (size_t i = 0; i < deformablePatches.size(); ++i)
                deformableDomain.patch(i).coefs() = refFluidGeometry.patch(deformablePatches[i]).coefs();

            try
            {
                meshDeformerPtr->updateMesh(boundaryDispDeformer);
            }
            catch (const std::exception& e)
            {
                gsWarn << "[FSI-Fluid-Vertex] Mesh deformation failed: " << e.what() << "\n";
            }
        }
        else
        {
            for (size_t i = 0; i < deformablePatches.size(); ++i)
                deformableDomain.patch(i).coefs() = refFluidGeometry.patch(deformablePatches[i]).coefs();
        }

        // Sync to full domain
        for (size_t i = 0; i < deformablePatches.size(); ++i)
        {
            const index_t pid = deformablePatches[i];
            fluidDomain.patch(pid).coefs() = deformableDomain.patch(i).coefs();
        }

        gsMultiPatch<>& patches = const_cast<gsMultiPatch<>&>(solver.getAssembler()->getPatches());
        for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
            patches.patch(p).coefs() = fluidDomain.patch(p).coefs();
        const gsMultiPatch<>& curGeo = patches;

        // Update boundary velocity functions
        const real_t dt_step = std::min(dt, precice_dt);
        for (index_t i = 0; i < static_cast<index_t>(fsiInterfaceSides.size()); ++i)
        {
            const index_t p = fsiInterfaceSides[i].patch;
            ifaceMeshVelFuncs[i]->setCurrent(&curGeo.patch(p));
            ifaceMeshVelFuncs[i]->setPrevious(&prevPatches.patch(p));
            ifaceMeshVelFuncs[i]->setDt(dt_step);
        }

        // Compute mesh displacement for ALE convection
        gsMatrix<real_t> meshDisp;
        try
        {
            GISMO_ASSERT(solver.getALEAssembler(), "ALE assembler not found");
            GISMO_ASSERT(!solver.getALEAssembler()->getMappers().empty(), "Velocity mapper not found");
            GISMO_ASSERT(!solver.getALEAssembler()->getBases().empty(), "Velocity basis not found");
            meshDisp = buildALEMeshDisplacement(
                curGeo, refFluidGeometry, deformablePatches, *solver.getALEAssembler());
        }
        catch (const std::exception& e)
        {
            gsWarn << "[FSI-Fluid-Vertex] Error computing mesh displacement: " << e.what() << "\n";
            meshDisp.resize(0,0);
        }

        params.options().setReal("timeStep", dt_step);
        if (meshDisp.size() > 0)
            solver.updateMesh(meshDisp);

        // Solve fluid
        solver.nextIteration();

        // Compute traction at quadrature points
        gsField<> velField = solver.constructSolution(0);
        gsField<> prsField = solver.constructSolution(1);
        const gsMultiPatch<> curFluidGeometry = solver.getAssembler()->getPatches();

        gsMatrix<real_t> stress_quad(2, quadPoints_xy.cols());
        stress_quad.setZero();

        const bool includeViscous = true;
        const bool pullBackToReference = usePiolaTraction;
        const gsVector<real_t> solidCenter = (gsVector<real_t>(2) << 0.0, 0.5).finished();

        for (index_t k = 0; k < quadPoints_xy.cols(); ++k)
        {
            const index_t patchId = quadSideLocations[k].first;
            const boxSide side = quadSideLocations[k].second;
            const gsVector<real_t> uv = quadPoints_uv.col(k);

            gsVector<real_t> traction = computeFluidTractionAtPoint(
                velField, prsField, refFluidGeometry, curFluidGeometry,
                uv, patchId, side, density, viscosity, includeViscous, pullBackToReference,
                &solidCenter);

            stress_quad.col(k) = traction;
        }

        // Validate exchanged stresses before writing them to preCICE.
        sanitizeFinite(stress_quad, "stress_quad", step);

        if (step % 20 == 0)
        {
            real_t maxTrac = stress_quad.cwiseAbs().maxCoeff();
            real_t meanTrac = stress_quad.cwiseAbs().mean();
            gsInfo << "[FSI-Fluid-Vertex|Step " << step << "] Traction stats: max=" << maxTrac
                   << ", mean=" << meanTrac << "\n";
        }

        // Write stress to preCICE
        participant.writeData(FluidMesh, StressData, fluidMeshIDs, stress_quad);

        // Advance coupling
        precice_dt = participant.advance(dt_step);

        // Checkpoint (read) - restore coefficients only to preserve meshDeformer references
        if (participant.requiresReadingCheckpoint())
        {
            solver.setSolutionCoefs(Uchk, 0);
            solver.setSolutionCoefs(Pchk, 1);
            solver.restoreALEState(solverALEChk);

            // Restore coefficients (not entire objects) to keep meshDeformer valid
            for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
            {
                prevPatches.patch(p).coefs() = prevPatchesCoefsChk[p];
                fluidDomain.patch(p).coefs() = fluidDomainCoefsChk[p];
            }
            for (size_t i = 0; i < deformablePatches.size(); ++i)
                deformableDomain.patch(i).coefs() = deformableDomainCoefsChk[i];

            gsMultiPatch<>& patches = const_cast<gsMultiPatch<>&>(solver.getAssembler()->getPatches());
            for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
                patches.patch(p).coefs() = fluidDomain.patch(p).coefs();

            gsInfo << "[FSI-Fluid-Vertex|Step " << step << "] Checkpoint restored\n";
        }
        else
        {
            time += dt_step;
            ++step;

            prevPatches = solver.getAssembler()->getPatches();

            if (plot && step % plotEverySteps == 0)
            {
                gsField<> vel = velField;
                gsField<> pre = prsField;
                std::string suf = util::to_string(step);
                gsWriteParaview(vel, "fsi_flap_vertex_U_" + suf, 400, false);
                gsWriteParaview(pre, "fsi_flap_vertex_P_" + suf, 400, false);

                for (index_t p = 0; p < fluidDomain.nPatches(); ++p)
                {
                    pvU.addTimestep("fsi_flap_vertex_U_" + suf + std::to_string(p), time, ".vts");
                    pvP.addTimestep("fsi_flap_vertex_P_" + suf + std::to_string(p), time, ".vts");
                }
            }
        }
    }

    if (plot)
    {
        pvU.save();
        pvP.save();
    }

    participant.finalize();
    gsInfo << "[FSI-Fluid-Vertex] Coupling completed successfully.\n";
    return 0;
}
