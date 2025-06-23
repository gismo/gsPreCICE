/** @file vertical-beam-fluid-JSON.cpp
 *
 * @brief JSON-based fluid participant for the vertical beam example
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
#include <gsJSON/gsJSON.h>
#include <gsPreCICE/gsJSONValidator.h>

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

// Helper function to create geometry from JSON
gsMultiPatch<> createGeometry(const json& geometryConfig) {
    gsMultiPatch<> patches;
    
    if (geometryConfig["type"] == "rectangle") {
        real_t x0 = geometryConfig["bounds"]["x0"].get<real_t>();
        real_t y0 = geometryConfig["bounds"]["y0"].get<real_t>();
        real_t x1 = geometryConfig["bounds"]["x1"].get<real_t>();
        real_t y1 = geometryConfig["bounds"]["y1"].get<real_t>();
        
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
    std::string jsonFile = "vertical-beam-fluid-config.json";
    std::string schemaFile = "vertical-beam-fluid-schema.json";
    
    // Simple command line argument parsing
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if ((arg == "-j" || arg == "--json") && i + 1 < argc) {
            jsonFile = argv[++i];
        }
        else if ((arg == "-s" || arg == "--schema") && i + 1 < argc) {
            schemaFile = argv[++i];
        }
        else if (arg == "--help" || arg == "-h") {
            gsInfo << "Usage: " << argv[0] << " [options]\n"
                   << "Options:\n"
                   << "  -j, --json FILE      JSON configuration file\n"
                   << "  -s, --schema FILE    JSON schema file\n"
                   << "  -h, --help           Show this help message\n";
            return EXIT_SUCCESS;
        }
    }
    
    gsInfo << "Using configuration file: " << jsonFile << "\n";

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
    gsInfo << "Using preCICE config file: " << preciceConfig << "\n";

    // Create geometry from JSON configuration
    gsInfo << "Creating geometry\n";
    gsMultiPatch<> patches = createGeometry(config["geometry"]);
    
    if (patches.nPatches() == 0) {
        gsWarn << "Failed to create geometry from configuration\n";
        return EXIT_FAILURE;
    }

    // Create discrete fluid mesh
    gsMatrix<> bbox = patches.patch(0).support();
    
    index_t numPoints = config["fluidMesh"].value("numPoints", 100);
    gsMatrix<> pointGrid = gsPointGrid<>(bbox, numPoints);
    gsMatrix<> mesh = patches.patch(0).eval(pointGrid);

    gsInfo << "Created fluid mesh with " << mesh.cols() << " points\n";

    // Initialize preCICE
    gsPreCICE<real_t> participant(participantName, preciceConfig);

    // For fluid participant, we don't create the mesh - we receive it from Solid
    gsInfo << "Fluid participant will receive mesh from Solid participant\n";
    std::map<std::string, gsVector<index_t>> meshVertexIDs;
    std::map<std::string, gsMatrix<>> meshCoordinates;

    // Initialize the coupling
    real_t precice_dt = participant.initialize();
    gsInfo << "Initialized participant with dt = " << precice_dt << "\n";
    
    // After initialization, get the mesh information from preCICE
    gsVector<index_t> meshIDs;
    gsMatrix<> meshPoints;
    participant.getMeshVertexIDsAndCoordinates("Solid-Mesh", meshIDs, meshPoints);
    
    meshVertexIDs["Solid-Mesh"] = meshIDs;
    meshCoordinates["Solid-Mesh"] = meshPoints;
    
    gsInfo << "Received mesh 'Solid-Mesh' with " << meshPoints.cols() << " points\n";

    // Initialize stress data based on received mesh size (3D for preCICE)
    gsMatrix<> StressPointData(3, meshPoints.cols());

    // Initialize time stepping variables
    real_t t = 0, dt = precice_dt;
    index_t timestep = 0;

    // Set up output
    std::string dirname = config["output"].value("directory", "./output");
    gsFileManager::mkdir(dirname);
    gsParaviewCollection collection(dirname + "/fluid_solution");

    // Load configuration
    std::string loadType = config["loading"].value("type", "constant");
    real_t loadMagnitude = config["loading"].value("magnitude", -5e10);
    
    gsInfo << "Load configuration: type=" << loadType 
           << ", magnitude=" << loadMagnitude << "\n";

    // Set up simple CSV output if configured
    std::ofstream csvFile;
    bool writePointData = config["output"].contains("pointData") && 
                         config["output"]["pointData"].value("enabled", false);
    
    if (writePointData) {
        std::string filename = config["output"]["pointData"]["filename"].get<std::string>();
        csvFile.open(filename);
        csvFile << "time,force_z\n";
    }

    bool plot = config["output"].value("plot", false);
    index_t plotmod = config["output"].value("plotmod", 10);

    gsInfo << "Starting time loop\n";

    // Time integration loop
    while (participant.isCouplingOngoing()) {
        timestep++;
        gsInfo << "Time step: " << timestep << ", t = " << t << "\n";
        
        if (participant.requiresWritingCheckpoint()) {
            gsInfo << "Writing Checkpoint\n";
        }

        // Read displacement data from coupling
        gsMatrix<> meshPointDisplacements;
        participant.readData("Solid-Mesh", "Displacement", meshVertexIDs["Solid-Mesh"], meshPointDisplacements);
        gsInfo << "Read Displacement data from Solid-Mesh (" << meshPointDisplacements.cols() << " points)\n";

        // Calculate stress data based on load configuration - using same logic as original vertical-beam-vertex-vertex-fluid.cpp
        if (loadType == "constant") {
            // loadCase==0 in original: Write data at the quadrature points
            StressPointData.setZero();
            StressPointData.row(2).setConstant(loadMagnitude);  // Apply force in z-direction as in original
        } else if (loadType == "spring") {
            // loadCase==1 in original: Spring load case
            // Original code: StressPointData(2,k) = -meshPointDisplacements(2,k);
            StressPointData.setZero();
            for (index_t k = 0; k != meshPointDisplacements.cols(); k++) {
                StressPointData(2, k) = -meshPointDisplacements(2, k);  // Use z-displacement as in original
            }
            // Impact loading at t==0 as in original
            if (t == 0) {
                StressPointData.row(2).array() += loadMagnitude;
            }
        } else if (loadType == "custom") {
            // Custom load - could be time-dependent or position-dependent
            StressPointData.setZero();
            real_t timeLoad = config["loading"].value("timeFunction", 1.0);
            StressPointData.row(2).setConstant(loadMagnitude * timeLoad);  // Apply in z-direction
        }

        // Write stress data to coupling
        participant.writeData("Solid-Mesh", "Stress", meshVertexIDs["Solid-Mesh"], StressPointData);
        gsInfo << "Wrote Stress data to Solid-Mesh\n";

        // Advance coupling
        precice_dt = participant.advance(dt);
        dt = std::min(precice_dt, dt);

        if (participant.requiresReadingCheckpoint()) {
            gsInfo << "Reading Checkpoint\n";
        } else {
            t += dt;

            // Output visualization and point data
            if (timestep % plotmod == 0 && plot) {
                // Create point cloud visualization with stress data
                gsMatrix<> currentMeshPoints = meshCoordinates["Solid-Mesh"];
                
                // Create a point cloud with stress magnitude data
                gsMatrix<> stressMagnitude(1, StressPointData.cols());
                for (index_t i = 0; i < StressPointData.cols(); ++i) {
                    stressMagnitude(0, i) = StressPointData.col(i).norm(); // Magnitude of stress vector
                }
                
                std::string fileName = dirname + "/fluid_solution" + util::to_string(timestep);
                
                // Write point cloud with stress data
                gsWriteParaviewPoints(currentMeshPoints, stressMagnitude, fileName);
                
                fileName = "fluid_solution" + util::to_string(timestep) + "0";
                collection.addTimestep(fileName, t, ".vts");
            }

            // Write point data
            if (writePointData && csvFile.is_open()) {
                // Write total force in z-direction (row 2)
                real_t totalForce = StressPointData.row(2).sum();
                csvFile << t << "," << totalForce << "\n";
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
    if (csvFile.is_open()) {
        csvFile.close();
    }

    return EXIT_SUCCESS;
} 