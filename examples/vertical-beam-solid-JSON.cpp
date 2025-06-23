/** @file vertical-beam-solid-JSON.cpp
 *
 * @brief JSON-based solid participant for the vertical beam example
 *
 * This file is part of the G+Smo library.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * Author(s): J. Li
 */

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsLookupFunction.h>
#include <gsJSON/gsJSON.h>
#include <gsPreCICE/gsJSONValidator.h>

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

// Helper function to safely read custom matrices from JSON
gsMatrix<> safeGetMatrix(const json& j) 
{
    try {
        if (!j.contains("rows") || !j.contains("cols") || !j.contains("data")) {
            gsWarn << "Matrix JSON missing required fields (rows, cols, data)\n";
            return gsMatrix<>();
        }
        
        index_t rows = j["rows"].get<index_t>();
        index_t cols = j["cols"].get<index_t>();
        
        if (!j["data"].is_array() || j["data"].size() != rows * cols) {
            gsWarn << "Matrix data has incorrect size. Expected: " << (rows * cols) 
                   << ", actual: " << j["data"].size() << "\n";
            return gsMatrix<>();
        }
        
        gsMatrix<> result(rows, cols);
        const auto& data = j["data"];
        
        for (index_t i = 0; i < rows; ++i) {
            for (index_t j = 0; j < cols; ++j) {
                result(i, j) = data[i * cols + j].get<real_t>();
            }
        }
        
        return result;
    } catch (const std::exception& e) {
        gsWarn << "Error parsing matrix from JSON: " << e.what() << "\n";
        return gsMatrix<>();
    }
}

// Helper function to apply quadrature options from JSON configuration
void applyQuadratureOptions(gsOptionList& quadOptions, const json& quadConfig, const std::string& source) {
    if (quadConfig.contains("quA")) {
        real_t quA = quadConfig["quA"].get<real_t>();
        quadOptions.setReal("quA", quA);
    }
    
    if (quadConfig.contains("quB")) {
        index_t quB = quadConfig["quB"].get<index_t>();
        quadOptions.setInt("quB", quB);
    }
    
    if (quadConfig.contains("quRule")) {
        index_t quRule = quadConfig["quRule"].get<index_t>();
        quadOptions.setInt("quRule", quRule);
    }
    
    if (quadConfig.contains("overInt")) {
        bool overInt = quadConfig["overInt"].get<bool>();
        quadOptions.setSwitch("overInt", overInt);
    }
    
    if (quadConfig.contains("quAb")) {
        real_t quAb = quadConfig["quAb"].get<real_t>();
        quadOptions.setReal("quAb", quAb);
    }
    
    if (quadConfig.contains("quBb")) {
        index_t quBb = quadConfig["quBb"].get<index_t>();
        quadOptions.setInt("quBb", quBb);
    }
    
    gsInfo << "Quadrature options loaded from " << source << ":\n";
    gsInfo << "  quA: " << quadOptions.getReal("quA") << "\n";
    gsInfo << "  quB: " << quadOptions.getInt("quB") << "\n";
    gsInfo << "  quRule: " << quadOptions.getInt("quRule") << "\n";
    gsInfo << "  overInt: " << quadOptions.getSwitch("overInt") << "\n";
}

// Helper function to create quadrature options from JSON
gsOptionList createQuadratureOptions(const json& config, const std::string& quadratureFile = "") {
    gsOptionList quadOptions = gsExprAssembler<>::defaultOptions();
    
    // First try to load from separate quadrature file if specified
    if (!quadratureFile.empty()) {
        std::ifstream quadFile(quadratureFile);
        if (quadFile.good()) {
            quadFile.close();
            gsJSON quadConfig(quadratureFile);
            
            if (quadConfig.contains("quadrature")) {
                applyQuadratureOptions(quadOptions, quadConfig["quadrature"], "separate file: " + quadratureFile);
                return quadOptions;
            } else {
                gsWarn << "Quadrature file " << quadratureFile << " does not contain 'quadrature' section\n";
            }
        } else {
            gsWarn << "Quadrature file not found: " << quadratureFile << "\n";
        }
    }
    
    // Fall back to main configuration file
    if (config.contains("quadrature")) {
        applyQuadratureOptions(quadOptions, config["quadrature"], "main configuration");
    } else {
        gsInfo << "No quadrature options found, using defaults\n";
    }
    
    return quadOptions;
}

// Helper function to convert boundary string to enum
boundary::side getBoundaryFromString(const std::string& str) {
    if (str == "north") return boundary::north;
    if (str == "south") return boundary::south;
    if (str == "east") return boundary::east;
    if (str == "west") return boundary::west;
    if (str == "front") return boundary::front;
    if (str == "back") return boundary::back;
    gsWarn << "Unknown boundary: " << str << "\n";
    return boundary::north; // default
}

// Create geometry from JSON config
gsMultiPatch<> createGeometry(const nlohmann::json& geomConfig) {
    gsMultiPatch<> patches;
    
    if (geomConfig["type"] == "rectangle") {
        real_t x0 = geomConfig["bounds"]["x0"].get<real_t>();
        real_t y0 = geomConfig["bounds"]["y0"].get<real_t>();
        real_t x1 = geomConfig["bounds"]["x1"].get<real_t>();
        real_t y1 = geomConfig["bounds"]["y1"].get<real_t>();
        
        // Create 2D rectangle then embed to 3D (like original program)
        patches.addPatch(gsNurbsCreator<>::BSplineRectangle(x0, y0, x1, y1));
        patches.addAutoBoundaries();
        patches.embed(3);  // Embed 2D geometry into 3D space (z=0)
    }
    
    return patches;
}

int main(int argc, char* argv[])
{
    // Process command line options
    std::string jsonFile = "vertical-beam-solid-config.json";
    std::string quadratureFile = "";
    
    // Simple command line argument parsing
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if ((arg == "-j" || arg == "--json") && i + 1 < argc) {
            jsonFile = argv[++i];
        }
        else if ((arg == "-q" || arg == "--quadrature") && i + 1 < argc) {
            quadratureFile = argv[++i];
        }
        else if (arg == "--help" || arg == "-h") {
            gsInfo << "Usage: " << argv[0] << " [options]\n"
                   << "Options:\n"
                   << "  -j, --json FILE      JSON configuration file\n"
                   << "  -q, --quadrature FILE JSON quadrature options file (optional)\n"
                   << "  -h, --help           Show this help message\n";
            return EXIT_SUCCESS;
        }
    }
    
    gsInfo << "Using configuration file: " << jsonFile << "\n";
    if (!quadratureFile.empty()) {
        gsInfo << "Using quadrature options file: " << quadratureFile << "\n";
    }

    // Check if files exist
    std::ifstream configCheck(jsonFile);
    bool configExists = configCheck.good();
    configCheck.close();

    if (!configExists) {
        gsWarn << "Configuration file not found: " << jsonFile << "\n";
        return EXIT_FAILURE;
    }

    // Load the configuration
    gsInfo << "Loading configuration from " << jsonFile << "\n";
    gsJSON config(jsonFile);

    // Extract participant information
    std::string participantName = config["participant"]["name"].get<std::string>();
    std::string preciceConfig = config["participant"]["preciceConfigFile"].get<std::string>();
    
    gsInfo << "Setting up participant: " << participantName << "\n";

    // Create geometry from JSON configuration
    gsMultiPatch<> patches = createGeometry(config["geometry"]);

    // Apply refinement from JSON
    gsMultiBasis<> bases(patches);
    if (config.contains("refinement")) {
        index_t numElevate = config["refinement"].value("numElevate", 0);
        index_t numRefine = config["refinement"].value("numRefine", 0);
        
        if (numElevate != 0) {
            patches.degreeElevate(numElevate);
        }
        
        for (int r = 0; r < numRefine; ++r) {
            patches.uniformRefine();
        }
        
        bases = gsMultiBasis<>(patches);
    }

    gsInfo << "Patches: " << patches.nPatches() << ", degree: " << bases.minCwiseDegree() << "\n";

    // Extract material parameters from JSON
    real_t rho = config["material"]["density"].get<real_t>();
    real_t E = config["material"]["youngsModulus"].get<real_t>();
    real_t nu = config["material"]["poissonsRatio"].get<real_t>();
    real_t thickness = config["material"]["thickness"].get<real_t>();

    // Set up coupling interfaces from JSON
    std::vector<patchSide> couplingInterfaces;
    if (config.contains("couplingInterfaces")) {
        for (const auto& interfaceConfig : config["couplingInterfaces"]) {
            index_t patch = interfaceConfig["patch"].get<index_t>();
            boundary::side side = getBoundaryFromString(interfaceConfig["boundary"].get<std::string>());
            couplingInterfaces.push_back(patchSide(patch, side));
            gsInfo << "Added coupling interface: patch " << patch << ", boundary " << interfaceConfig["boundary"].get<std::string>() << "\n";
        }
    } else {
        // Default: assume all boundaries are coupling interfaces (for backward compatibility)
        gsWarn << "No coupling interfaces specified, using default behavior\n";
    }

    // Initialize preCICE
    gsPreCICE<real_t> participant(participantName, preciceConfig);

    // Set up mesh with JSON-configured quadrature options
    gsOptionList quadOptions = createQuadratureOptions(config, quadratureFile);
    gsMatrix<> quadPoints;
    
    if (!couplingInterfaces.empty()) {
        // Use specified coupling interfaces
        quadPoints = gsQuadrature::getAllNodes(bases.basis(0), quadOptions, couplingInterfaces);
        gsInfo << "Generated quadrature points on specified coupling interfaces\n";
    } else {
        // Default: use all boundaries as coupling interfaces
        quadPoints = gsQuadrature::getAllNodes(bases.basis(0), quadOptions);
    }
    
    gsVector<index_t> quadPointIDs;
    participant.addMesh("Solid-Mesh", quadPoints, quadPointIDs);
    gsMatrix<> quadPointsData(3, quadPoints.cols());
    quadPointsData.setZero();
    
    gsInfo << "Created mesh 'Solid-Mesh' with " << quadPoints.cols() << " quadrature points\n";

    // Initialize the coupling
    real_t precice_dt = participant.initialize();

    // Define boundary conditions from JSON
    gsBoundaryConditions<> bcInfo;
    
    for (const auto& bcConfig : config["boundaryConditions"]) {
        index_t patch = bcConfig["patch"].get<index_t>();
        std::string boundaryStr = bcConfig["boundary"].get<std::string>();
        std::string type = bcConfig["type"].get<std::string>();
        
        boundary::side side = getBoundaryFromString(boundaryStr);
        
        if (type == "dirichlet") {
            bcInfo.addCondition(patch, side, condition_type::dirichlet, nullptr, -1);
        } else if (type == "clamped") {
            index_t component = bcConfig.value("component", 2);
            bcInfo.addCondition(patch, side, condition_type::clamped, nullptr, component);
        }
    }
    
    bcInfo.setGeoMap(patches);

    // Set up surface force function
    gsLookupFunction<real_t> surfForce(quadPoints, quadPointsData);

    // Set up material matrices
    gsFunctionExpr<> E_modulus(std::to_string(E), 3);
    gsFunctionExpr<> PoissonRatio(std::to_string(nu), 3);
    gsFunctionExpr<> Density(std::to_string(rho), 3);
    gsFunctionExpr<> t(std::to_string(thickness), 3);

    std::vector<gsFunctionSet<>*> parameters(2);
    parameters[0] = &E_modulus;
    parameters[1] = &PoissonRatio;

    gsOptionList options;
    options.addInt("Material", "Material model: (0): SvK | (1): NH | (2): NH_ext | (3): MR | (4): Ogden", 0);
    options.addInt("Implementation", "Implementation: (0): Composites | (1): Analytical | (2): Generalized | (3): Spectral", 1);
    
    auto materialMatrix = getMaterialMatrix<3, real_t>(patches, t, parameters, Density, options);

    // Create assembler
    gsThinShellAssembler<3, real_t, true> assembler(patches, bases, bcInfo, surfForce, materialMatrix.get());
    gsOptionList assemblerOptions = quadOptions.wrapIntoGroup("ExprAssembler");
    assembler.setOptions(assemblerOptions);

    // Initialize time stepping variables
    index_t timestep = 0;
    index_t timestep_checkpoint = 0;

    // Compute mass and stiffness matrices
    assembler.assembleMass();
    gsSparseMatrix<> M = assembler.massMatrix();
    assembler.assemble();

    // Set up output
    std::string dirname = config["output"].value("directory", "./output");
    gsFileManager::mkdir(dirname);
    gsParaviewCollection collection(dirname + "/solution");

    real_t dt = precice_dt;

    // Define solver functions
    gsMultiPatch<> solutions;
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = 
        [&assembler, &solutions](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> & m) {
        ThinShellAssemblerStatus status;
        assembler.constructSolution(x, solutions);
        status = assembler.assembleMatrix(solutions);
        m = assembler.matrix();
        return true;
    };

    gsStructuralAnalysisOps<real_t>::TResidual_t Residual = 
        [&assembler, &solutions](gsMatrix<real_t> const &x, real_t t, gsVector<real_t> & result) {
        ThinShellAssemblerStatus status;
        assembler.constructSolution(x, solutions);
        status = assembler.assembleVector(solutions);
        result = assembler.rhs();
        return true;
    };

    gsSparseMatrix<> C = gsSparseMatrix<>(assembler.numDofs(), assembler.numDofs());
    gsStructuralAnalysisOps<real_t>::Damping_t Damping = 
        [&C](const gsVector<real_t> &, gsSparseMatrix<real_t> & m) { m = C; return true; };
    gsStructuralAnalysisOps<real_t>::Mass_t Mass = 
        [&M](gsSparseMatrix<real_t> & m) { m = M; return true; };

    // Create time integrator based on JSON configuration
    gsDynamicBase<real_t>* timeIntegrator = nullptr;
    std::string method = config["solver"].value("method", "newmark");
    
    if (method == "explicit_euler") {
        timeIntegrator = new gsDynamicExplicitEuler<real_t, true>(Mass, Damping, Jacobian, Residual);
    } else if (method == "implicit_euler") {
        timeIntegrator = new gsDynamicImplicitEuler<real_t, true>(Mass, Damping, Jacobian, Residual);
    } else if (method == "newmark") {
        timeIntegrator = new gsDynamicNewmark<real_t, true>(Mass, Damping, Jacobian, Residual);
    } else if (method == "bathe") {
        timeIntegrator = new gsDynamicBathe<real_t, true>(Mass, Damping, Jacobian, Residual);
    } else if (method == "wilson") {
        timeIntegrator = new gsDynamicWilson<real_t, true>(Mass, Damping, Jacobian, Residual);
        timeIntegrator->options().setReal("gamma", 1.4);
    } else if (method == "rk4") {
        timeIntegrator = new gsDynamicRK4<real_t, true>(Mass, Damping, Jacobian, Residual);
    } else {
        gsWarn << "Unknown solver method: " << method << ", using Newmark\n";
        timeIntegrator = new gsDynamicNewmark<real_t, true>(Mass, Damping, Jacobian, Residual);
    }

    real_t tolerance = config["solver"].value("tolerance", 1e-3);
    bool verbose = config["solver"].value("verbose", true);
    
    timeIntegrator->options().setReal("DT", dt);
    timeIntegrator->options().setReal("TolU", tolerance);
    timeIntegrator->options().setSwitch("Verbose", verbose);

    // Initialize solution vectors
    gsVector<> F = assembler.rhs();
    gsVector<> F_checkpoint, U_checkpoint, V_checkpoint, A_checkpoint, U, V, A;

    F_checkpoint = F;
    U_checkpoint = U = gsVector<real_t>::Zero(assembler.numDofs(), 1);
    V_checkpoint = V = gsVector<real_t>::Zero(assembler.numDofs(), 1);
    A_checkpoint = A = gsVector<real_t>::Zero(assembler.numDofs(), 1);

    real_t time = 0;
    bool plot = config["output"].value("plot", true);
    index_t plotmod = config["output"].value("plotmod", 1);

    // Define solution variable outside the loop
    gsMultiPatch<> solution;

    // Time integration loop
    while (participant.isCouplingOngoing()) {
        timestep++;
        
        if (participant.requiresWritingCheckpoint()) {
            U_checkpoint = U;
            V_checkpoint = V;
            A_checkpoint = A;
            timestep_checkpoint = timestep;
        }

        // Read stress data from coupling (using JSON configuration)
        if (config.contains("readData")) {
            for (const auto& readData : config["readData"]) {
                std::string meshName = readData["mesh"].get<std::string>();
                std::string dataName = readData["name"].get<std::string>();
                
                if (meshName == "Solid-Mesh" && dataName == "Stress") {
                    participant.readData(meshName, dataName, quadPointIDs, quadPointsData);
                    gsInfo << "Read data: " << dataName << " from mesh: " << meshName << "\n";
                }
            }
        } else {
            // Backward compatibility: use hardcoded values
            participant.readData("Solid-Mesh", "Stress", quadPointIDs, quadPointsData);
        }
        
        assembler.assemble();
        F = assembler.rhs();

        // solve gismo timestep
        timeIntegrator->step(time, dt, U, V, A);

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt, precice_dt);

        // Write displacement data to coupling (using JSON configuration)
        gsVector<> displacements = U;
        solution = assembler.constructDisplacement(displacements);
        gsMatrix<> quadPointsDisplacements = solution.patch(0).eval(quadPoints);
        
        if (config.contains("writeData")) {
            for (const auto& writeData : config["writeData"]) {
                std::string meshName = writeData["mesh"].get<std::string>();
                std::string dataName = writeData["name"].get<std::string>();
                
                if (meshName == "Solid-Mesh" && dataName == "Displacement") {
                    // preCICE expects 3D data, write full displacement vector
                    participant.writeData(meshName, dataName, quadPointIDs, quadPointsDisplacements);
                    gsInfo << "Wrote data: " << dataName << " to mesh: " << meshName << "\n";
                }
            }
        } else {
            // Backward compatibility: use hardcoded values
            participant.writeData("Solid-Mesh", "Displacement", quadPointIDs, quadPointsDisplacements);
        }

        // do the coupling
        precice_dt = participant.advance(dt);

        if (participant.requiresReadingCheckpoint()) {
            U = U_checkpoint;
            V = V_checkpoint;
            A = A_checkpoint;
            timestep = timestep_checkpoint;
            // Reconstruct solution based on restored displacement
            gsVector<> displacements = U;
            solution = assembler.constructDisplacement(displacements);
        } else {
            // advance variables
            time += dt;
        }

        // Always save visualization based on current solution state
        gsField<> solField(patches, solution);
        if (timestep % plotmod == 0 && plot) {
            std::string fileName = dirname + "/solution" + util::to_string(timestep);
            gsWriteParaview<>(solField, fileName, 500);
            fileName = "solution" + util::to_string(timestep) + "0";
            collection.addTimestep(fileName, time, ".vts");
        }
    }

    // Save final solution and collection
    if (plot) {
        // Ensure final solution is based on current displacement state
        gsVector<> finalDisplacements = U;
        gsMultiPatch<> finalSolution = assembler.constructDisplacement(finalDisplacements);
        
        // Save final solution frame
        gsField<> finalSolField(patches, finalSolution);
        std::string finalFileName = dirname + "/solution" + util::to_string(timestep);
        gsWriteParaview<>(finalSolField, finalFileName, 500);
        
        // Add final frame to collection if not already added
        if (timestep % plotmod != 0) {
            std::string fileName = "solution" + util::to_string(timestep) + "0";
            collection.addTimestep(fileName, time, ".vts");
        }
        
        // Save the animation collection
        collection.save();
        gsInfo << "Saved animation collection to " << dirname << "/solution.pvd\n";
    }

    delete timeIntegrator;
    return EXIT_SUCCESS;
} 