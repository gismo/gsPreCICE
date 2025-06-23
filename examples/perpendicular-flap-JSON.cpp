/** @file perpendicular-flap-JSON.cpp
 *
 * @brief JSON-based elasticity participant for the preCICE example "perpendicular-flap"
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
#include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsLookupFunction.h>
#include <gsJSON/gsJSON.h>
#include <gsPreCICE/gsJSONValidator.h>

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

// Helper function to safely read custom matrices from JSON
gsMatrix<> safeGetMatrix(const json& j) 
{
    try {
        if (!j.contains("rows") || !j.contains("cols") || !j.contains("data")) {
            gsWarn << "Matrix JSON missing required fields (rows, cols, data)\n";
            return gsMatrix<>(); // Return empty matrix
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
    gsOptionList quadOptions = gsAssembler<>::defaultOptions();
    
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

// Helper function to create B-spline geometry from JSON
gsMultiPatch<> createBSplineGeometry(const json& geometryConfig) {
    gsMultiPatch<> patches;
    
    if (geometryConfig["type"] != "bspline") {
        gsWarn << "Only B-spline geometry is supported\n";
        return patches;
    }
    
    index_t dim = geometryConfig["dimension"].get<index_t>();
    
    if (dim == 2) {
        // Create knot vectors from JSON
        std::vector<gsKnotVector<>> knotVectors;
        
        for (const auto& kvConfig : geometryConfig["knotVectors"]) {
            std::vector<real_t> knots = kvConfig["knots"].get<std::vector<real_t>>();
            std::vector<index_t> multiplicities = kvConfig["multiplicities"].get<std::vector<index_t>>();
            
            // Handle special case for direction 1 with numElements
            if (kvConfig.contains("numElements")) {
                index_t numElements = kvConfig["numElements"].get<index_t>();
                gsKnotVector<> kv(knots[0], knots[1], numElements, multiplicities[0]);
                knotVectors.push_back(kv);
            } else {
                // Create knot vector with given knots and multiplicities  
                std::vector<real_t> expandedKnots;
                for (size_t i = 0; i < knots.size(); ++i) {
                    for (index_t m = 0; m < multiplicities[i]; ++m) {
                        expandedKnots.push_back(knots[i]);
                    }
                }
                gsKnotVector<> kv(expandedKnots);
                knotVectors.push_back(kv);
            }
        }
        
        // Create tensor B-spline basis
        gsTensorBSplineBasis<2> basis(knotVectors[0], knotVectors[1]);
        
        // Get control points
        gsMatrix<> coefs;
        
        if (geometryConfig.contains("controlPoints")) {
            const auto& cpConfig = geometryConfig["controlPoints"];
            
            // Check if coordinates are directly specified
            if (cpConfig.contains("coordinates")) {
                coefs = safeGetMatrix(cpConfig["coordinates"]);
                if (coefs.size() == 0) {
                    gsWarn << "Failed to read control point coordinates from JSON, using default anchors\n";
                    coefs = basis.anchors();
                    coefs.transposeInPlace();
                }
                gsInfo << "Using directly specified control point coordinates\n";
            } else {
                // Use default anchors and apply transformations
                coefs = basis.anchors();
                coefs.transposeInPlace();
                
                if (cpConfig.contains("scaleX")) {
                    real_t scaleX = cpConfig["scaleX"].get<real_t>();
                    coefs.col(0).array() *= scaleX;
                }
                
                if (cpConfig.contains("offsetX")) {
                    real_t offsetX = cpConfig["offsetX"].get<real_t>();
                    coefs.col(0).array() += offsetX;
                }
                gsInfo << "Using default anchors with transformations\n";
            }
        } else {
            // No control point configuration, use default anchors
            coefs = basis.anchors();
            coefs.transposeInPlace();
            gsInfo << "Using default anchor points\n";
        }
        
        patches.addPatch(basis.makeGeometry(give(coefs)));
    }
    
    return patches;
}

int main(int argc, char* argv[])
{
    // Process command line options
    std::string jsonFile = "perpendicular-flap-config.json";
    std::string schemaFile = "perpendicular-flap-schema.json";
    std::string quadratureFile = "";
    
    // Simple command line argument parsing
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if ((arg == "-j" || arg == "--json") && i + 1 < argc) {
            jsonFile = argv[++i];
        }
        else if ((arg == "-s" || arg == "--schema") && i + 1 < argc) {
            schemaFile = argv[++i];
        }
        else if ((arg == "-q" || arg == "--quadrature") && i + 1 < argc) {
            quadratureFile = argv[++i];
        }
        else if (arg == "--help" || arg == "-h") {
            gsInfo << "Usage: " << argv[0] << " [options]\n"
                   << "Options:\n"
                   << "  -j, --json FILE      JSON configuration file (default: perpendicular-flap-config.json)\n"
                   << "  -s, --schema FILE    JSON schema file (default: perpendicular-flap-schema.json)\n"
                   << "  -q, --quadrature FILE JSON quadrature options file (optional)\n"
                   << "  -h, --help           Show this help message\n";
            return EXIT_SUCCESS;
        }
    }
    
    gsInfo << "Using configuration file: " << jsonFile << "\n";
    gsInfo << "Using schema file: " << schemaFile << "\n";
    if (!quadratureFile.empty()) {
        gsInfo << "Using quadrature options file: " << quadratureFile << "\n";
    }

    // Check if files exist
    std::ifstream configCheck(jsonFile);
    bool configExists = configCheck.good();
    configCheck.close();

    std::ifstream schemaCheck(schemaFile);
    bool schemaExists = schemaCheck.good();
    schemaCheck.close();

    // Handle missing configuration file
    if (!configExists) {
        gsWarn << "Configuration file not found: " << jsonFile << "\n";
        return EXIT_FAILURE;
    }

    // Load the configuration
    gsInfo << "Loading configuration from " << jsonFile << "\n";
    gsJSON config(jsonFile);

    // Validate against schema if it exists
    if (schemaExists) {
        gsInfo << "Validating configuration against schema\n";
        gsJSONValidator validator(schemaFile);
        std::string validationError;
        
        if (!validator.validate(config, validationError)) {
            gsWarn << "Configuration validation failed: " << validationError << "\n";
            return EXIT_FAILURE;
        }
        
        gsInfo << "Configuration validated successfully\n";
    } else {
        gsWarn << "Schema file not found: " << schemaFile << "\n";
        gsWarn << "Skipping schema validation\n"; // The validation is not mendatory
    }

    // Extract participant information
    std::string participantName = config["participant"]["name"].get<std::string>();
    // Now it is still depend on the preciceConfigFile, but it is not necessary.
    std::string preciceConfig = config["participant"]["preciceConfigFile"].get<std::string>();

    
    gsInfo << "Setting up participant: " << participantName << "\n";
    gsInfo << "Using preCICE config file: " << preciceConfig << "\n";

    // Create geometry from JSON configuration
    gsInfo << "Creating B-spline geometry\n";
    gsMultiPatch<> patches = createBSplineGeometry(config["geometry"]);
    
    if (patches.nPatches() == 0) {
        gsWarn << "Failed to create geometry from configuration\n";
        return EXIT_FAILURE;
    }

    // Create bases
    gsMultiBasis<> bases(patches);

    // Apply refinement from JSON
    if (config.contains("refinement")) {
        index_t numElevate = config["refinement"].value("numElevate", 0);
        index_t numRefine = config["refinement"].value("numRefine", 0);
        
        // Set degree
        bases.setDegree(bases.maxCwiseDegree() + numElevate);

        // h-refine each basis
        for (int r = 0; r < numRefine; ++r)
            bases.uniformRefine();
    }

    gsInfo << "Patches: " << patches.nPatches() << ", degree: " << bases.minCwiseDegree() << "\n";

    // Extract material parameters from JSON
    real_t rho = config["material"]["density"].get<real_t>();
    real_t E = config["material"]["youngsModulus"].get<real_t>();
    real_t nu = config["material"]["poissonsRatio"].get<real_t>();
    
    gsInfo << "Material properties: rho=" << rho << ", E=" << E << ", nu=" << nu << "\n";

    // Set up coupling interfaces from JSON
    std::vector<patchSide> couplingInterfaces;
    for (const auto& interfaceConfig : config["couplingInterfaces"]) {
        index_t patch = interfaceConfig["patch"].get<index_t>();
        boundary::side side = getBoundaryFromString(interfaceConfig["boundary"].get<std::string>());
        couplingInterfaces.push_back(patchSide(patch, side));
    }

    // Initialize preCICE
    gsPreCICE<real_t> participant(participantName, preciceConfig);

    // Set up meshes from JSON configuration
    // Create quadrature options from JSON configuration
    gsOptionList quadOptions = createQuadratureOptions(config, quadratureFile);

    gsInfo << "Setting up meshes\n";
    std::map<std::string, gsVector<index_t>> meshVertexIDs;
    std::map<std::string, gsMatrix<>> meshCoordinates;
    
    for (const auto& meshConfig : config["meshes"]) {
        std::string meshName = meshConfig["name"].get<std::string>();
        std::string meshType = meshConfig["type"].get<std::string>();
        
        if (meshType == "quadrature") {
            // Get quadrature nodes on coupling interface using JSON-configured options
            gsMatrix<> quad_uv = gsQuadrature::getAllNodes(bases.basis(0), quadOptions, couplingInterfaces);
            gsMatrix<> quad_xy = patches.patch(0).eval(quad_uv);
            
            gsVector<index_t> vertexIDs;
            participant.addMesh(meshName, quad_xy, vertexIDs);
            meshVertexIDs[meshName] = vertexIDs;
            meshCoordinates[meshName] = quad_xy;
            
            gsInfo << "Added quadrature mesh: " << meshName << " with " << quad_xy.cols() << " points\n";
        }
    }

    // Initialize the coupling
    real_t precice_dt = participant.initialize();

    // Define boundary conditions from JSON
    gsBoundaryConditions<> bcInfo;
    gsConstantFunction<> g_D(0, patches.geoDim());
    
    // Set up stress lookup function for coupling
    gsMatrix<> quad_stress(2, meshCoordinates.begin()->second.cols());
    quad_stress.setZero();
    gsLookupFunction<real_t> g_L(meshCoordinates.begin()->second, quad_stress);
    
    // Add boundary conditions from JSON
    for (const auto& bcConfig : config["boundaryConditions"]) {
        index_t patch = bcConfig["patch"].get<index_t>();
        boundary::side side = getBoundaryFromString(bcConfig["boundary"].get<std::string>());
        std::string type = bcConfig["type"].get<std::string>();
        
        if (type == "dirichlet") {
            index_t component = bcConfig["component"].get<index_t>();
            real_t value = bcConfig["value"].get<real_t>();
            gsConstantFunction<> g_val(value, patches.geoDim());
            bcInfo.addCondition(patch, side, condition_type::dirichlet, &g_val, component);
        }
    }
    
    // Add coupling boundary conditions (Neumann)
    for (const auto& interfaceConfig : config["couplingInterfaces"]) {
        index_t patch = interfaceConfig["patch"].get<index_t>();
        boundary::side side = getBoundaryFromString(interfaceConfig["boundary"].get<std::string>());
        bcInfo.addCondition(patch, side, condition_type::neumann, &g_L, -1, true);
    }
    
    bcInfo.setGeoMap(patches);

    // Create source functions
    gsConstantFunction<> g(0., 0., 2);
    gsConstantFunction<> gravity(0., 0., 2);

    // Create mass assembler
    gsMassAssembler<real_t> massAssembler(patches, bases, bcInfo, gravity);
    massAssembler.options().setReal("Density", rho);
    massAssembler.assemble();

    // Create stiffness assembler
    gsElasticityAssembler<real_t> assembler(patches, bases, bcInfo, g);
    assembler.options().setReal("YoungsModulus", E);
    assembler.options().setReal("PoissonsRatio", nu);
    assembler.options().setInt("MaterialLaw", material_law::saint_venant_kirchhoff);
    assembler.assemble();

    gsSparseMatrix<> M = massAssembler.matrix();
    gsSparseMatrix<> K = assembler.matrix();
    gsSparseMatrix<> C = gsSparseMatrix<>(assembler.numDofs(), assembler.numDofs());

    // Time step
    real_t dt = precice_dt;

    // Initialize solution vectors
    gsMatrix<> solVector;
    solVector.setZero(assembler.numDofs(), 1);

    std::vector<gsMatrix<>> fixedDofs = assembler.allFixedDofs();
    gsVector<> U_checkpoint, V_checkpoint, A_checkpoint, U, V, A;

    U_checkpoint = U = gsVector<real_t>::Zero(assembler.numDofs(), 1);
    V_checkpoint = V = gsVector<real_t>::Zero(assembler.numDofs(), 1);
    A_checkpoint = A = gsVector<real_t>::Zero(assembler.numDofs(), 1);

    // Set up output
    gsParaviewCollection collection("./output/solution");
    index_t timestep = 0;
    index_t timestep_checkpoint = 0;
    
    bool plot = config["output"].value("plot", false);
    index_t plotmod = config["output"].value("plotmod", 10);

    // Define solver functions
    gsStructuralAnalysisOps<real_t>::Jacobian_t Jacobian = 
        [&assembler, &fixedDofs](gsMatrix<real_t> const &x, gsSparseMatrix<real_t> & m) {
        assembler.assemble(x, fixedDofs);
        m = assembler.matrix();
        return true;
    };

    gsStructuralAnalysisOps<real_t>::TResidual_t Residual = 
        [&assembler, &fixedDofs](gsMatrix<real_t> const &x, real_t t, gsVector<real_t> & result) {
        assembler.assemble(x, fixedDofs);
        result = assembler.rhs();
        return true;
    };

    gsStructuralAnalysisOps<real_t>::TForce_t TForce = 
        [&assembler](real_t, gsVector<real_t> & result) {
        assembler.assemble();
        result = assembler.rhs();
        return true;
    };

    gsStructuralAnalysisOps<real_t>::Damping_t Damping = 
        [&C](const gsVector<real_t> &, gsSparseMatrix<real_t> & m) { m = C; return true; };
    gsStructuralAnalysisOps<real_t>::Mass_t Mass = 
        [&M](gsSparseMatrix<real_t> & m) { m = M; return true; };
    gsStructuralAnalysisOps<real_t>::Stiffness_t Stiffness = 
        [&K](gsSparseMatrix<real_t> & m) { m = K; return true; };

    // Create time integrator based on JSON configuration
    gsDynamicBase<real_t>* timeIntegrator = nullptr;
    std::string method = config["solver"].value("method", "newmark");
    bool nonlinear = config["solver"].value("nonlinear", false);
    
    if (method == "explicit_euler") {
        timeIntegrator = new gsDynamicExplicitEuler<real_t, true>(Mass, Damping, Jacobian, Residual);
    } else if (method == "implicit_euler") {
        timeIntegrator = new gsDynamicImplicitEuler<real_t, true>(Mass, Damping, Jacobian, Residual);
    } else if (method == "newmark") {
        if (nonlinear) {
            timeIntegrator = new gsDynamicNewmark<real_t, true>(Mass, Damping, Jacobian, Residual);
        } else {
            timeIntegrator = new gsDynamicNewmark<real_t, false>(Mass, Damping, Stiffness, TForce);
        }
    } else if (method == "bathe") {
        timeIntegrator = new gsDynamicBathe<real_t, true>(Mass, Damping, Jacobian, Residual);
    } else if (method == "wilson") {
        timeIntegrator = new gsDynamicWilson<real_t, true>(Mass, Damping, Jacobian, Residual);
        timeIntegrator->options().setReal("gamma", 1.4);
    } else if (method == "rk4") {
        timeIntegrator = new gsDynamicRK4<real_t, true>(Mass, Damping, Jacobian, Residual);
    } else {
        gsWarn << "Unknown solver method: " << method << ", using Newmark\n";
        timeIntegrator = new gsDynamicNewmark<real_t, false>(Mass, Damping, Stiffness, TForce);
    }

    real_t tolerance = config["solver"].value("tolerance", 1e-6);
    bool verbose = config["solver"].value("verbose", true);
    
    timeIntegrator->options().setReal("DT", dt);
    timeIntegrator->options().setReal("TolU", tolerance);
    timeIntegrator->options().setSwitch("Verbose", verbose);

    real_t time = 0;

    // Set up point data output if configured
    gsStructuralAnalysisOutput<real_t>* writer = nullptr;
    if (config["output"].contains("pointData") && 
        config["output"]["pointData"].value("enabled", false)) {
        
        std::string filename = config["output"]["pointData"]["filename"].get<std::string>();
        gsMatrix<> points = safeGetMatrix(config["output"]["pointData"]["points"]);
        
        if (points.size() > 0) {
            writer = new gsStructuralAnalysisOutput<real_t>(filename, points);
            writer->init({"x", "y"}, {"time"});
        }
    }

    // Plot initial solution
    if (plot) {
        gsMultiPatch<> solution;
        assembler.constructSolution(solVector, fixedDofs, solution);

        gsField<> solField(patches, solution);
        std::string fileName = "./output/solution" + util::to_string(timestep);
        gsWriteParaview<>(solField, fileName, 500);
        fileName = "solution" + util::to_string(timestep) + "0";
        collection.addTimestep(fileName, time, ".vts");
    }

    gsInfo << "Starting time loop with initial dt = " << dt << "\n";

    // Time integration loop
    while (participant.isCouplingOngoing()) {
        timestep++;
        gsInfo << "Time step: " << timestep << ", t = " << time << "\n";
        
        if (participant.requiresWritingCheckpoint()) {
            U_checkpoint = U;
            V_checkpoint = V;
            A_checkpoint = A;
            timestep_checkpoint = timestep;
            
            gsInfo << "Checkpoint written:\n";
            gsInfo << "\t ||U|| = " << U.norm() << "\n";
            gsInfo << "\t ||V|| = " << V.norm() << "\n";
            gsInfo << "\t ||A|| = " << A.norm() << "\n";
        }

        // Read stress data from coupling
        for (const auto& readData : config["readData"]) {
            std::string meshName = readData["mesh"].get<std::string>();
            std::string dataName = readData["name"].get<std::string>();
            
            if (meshVertexIDs.count(meshName)) {
                participant.readData(meshName, dataName, meshVertexIDs[meshName], quad_stress);
                gsInfo << "Read data: " << dataName << " from mesh: " << meshName << "\n";
            }
        }

        // Solve elasticity timestep
        gsInfo << "Solving timestep " << time << "...\n";
        timeIntegrator->step(time, dt, U, V, A);
        solVector = U;
        gsInfo << "Finished\n";

        // Potentially adjust non-matching timestep sizes
        dt = std::min(dt, precice_dt);

        // Write displacement data to coupling
        for (const auto& writeData : config["writeData"]) {
            std::string meshName = writeData["mesh"].get<std::string>();
            std::string dataName = writeData["name"].get<std::string>();
            
            if (meshVertexIDs.count(meshName)) {
                gsMultiPatch<> solution;
                assembler.constructSolution(solVector, fixedDofs, solution);
                
                // Get quadrature points in parametric domain using JSON-configured options
                gsMatrix<> quad_uv = gsQuadrature::getAllNodes(bases.basis(0), quadOptions, couplingInterfaces);
                
                // Evaluate solution at quadrature points
                gsMatrix<> result(patches.geoDim(), quad_uv.cols());
                solution.patch(0).eval_into(quad_uv, result);
                
                participant.writeData(meshName, dataName, meshVertexIDs[meshName], result);
                gsInfo << "Wrote data: " << dataName << " to mesh: " << meshName << "\n";
            }
        }

        // Advance coupling
        precice_dt = participant.advance(dt);

        if (participant.requiresReadingCheckpoint()) {
            U = U_checkpoint;
            V = V_checkpoint;
            A = A_checkpoint;
            timestep = timestep_checkpoint;
        } else {
            // Advance variables
            time += dt;

            gsMultiPatch<> solution;
            assembler.constructSolution(solVector, fixedDofs, solution);
            gsField<> solField(patches, solution);
            
            if (timestep % plotmod == 0) {
                if (plot) {
                    std::string fileName = "./output/solution" + util::to_string(timestep);
                    gsWriteParaview<>(solField, fileName, 500);
                    fileName = "solution" + util::to_string(timestep) + "0";
                    collection.addTimestep(fileName, time, ".vts");
                }
                
                // Write point data
                if (writer) {
                    gsMatrix<> pointDataMatrix, otherDataMatrix(1, 1);
                    gsMatrix<> points = safeGetMatrix(config["output"]["pointData"]["points"]);
                    solution.patch(0).eval_into(points, pointDataMatrix);
                    otherDataMatrix << time;
                    writer->add(pointDataMatrix, otherDataMatrix);
                }
            }
        }
    }

    // Finalize the coupling
    participant.finalize();
    gsInfo << "Coupling finalized\n";

    if (plot) {
        collection.save();
    }

    // Clean up
    delete timeIntegrator;
    if (writer) delete writer;

    // Print timing information
    gsInfo << "Timing statistics:\n";
    gsInfo << "- Read time: " << participant.readTime() << " seconds\n";
    gsInfo << "- Write time: " << participant.writeTime() << " seconds\n";
    gsInfo << "- Initialize time: " << participant.initializeTime() << " seconds\n";

    return EXIT_SUCCESS;
} 