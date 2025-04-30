/** @file dummy-participant-JSON.cpp
 *
 * @brief Example of a preCICE participant using JSON configuration
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
 #include <gsJSON/gsJSON.h>
 #include <gsJSON/gsJSONValidator.h>
 
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
         
         gsInfo << "Matrix JSON structure: rows=" << rows 
                << ", cols=" << cols 
                << ", data.size=" << j["data"].size() << "\n";
         
         if (!j["data"].is_array() || j["data"].size() != rows * cols) {
             gsWarn << "Matrix data has incorrect size. Expected: " << (rows * cols) 
                    << ", actual: " << j["data"].size() << "\n";
             return gsMatrix<>();
         }
         
         gsMatrix<> result(rows, cols);
         const auto& data = j["data"];
         
         for (index_t i = 0; i < rows; ++i) {
             for (index_t j = 0; j < cols; ++j) {  // <-- Fixed condition
                 result(i, j) = data[i * cols + j].get<real_t>();
             }
         }
         
         return result;
     } catch (const std::exception& e) {
         gsWarn << "Error parsing matrix from JSON: " << e.what() << "\n";
         return gsMatrix<>();
     }
 }

 int main(int argc, char* argv[])
 {
     // Process command line options - simple parsing approach
     std::string jsonFile = "participant-config.json";
     std::string schemaFile = "participant-schema.json";
     
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
                    << "  -j, --json FILE    JSON configuration file (default: participant-config.json)\n"
                    << "  -s, --schema FILE  JSON schema file (default: participant-schema.json)\n"
                    << "  -h, --help         Show this help message\n";
             return EXIT_SUCCESS;
         }
     }
     
     gsInfo << "Using configuration file: " << jsonFile << "\n";
     gsInfo << "Using schema file: " << schemaFile << "\n";

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
         gsWarn << "Skipping schema validation\n";
     }

     // Extract participant information
     std::string participantName = config["participant"]["name"].get<std::string>();
     std::string configFile = config["participant"]["preciceConfigFile"].get<std::string>();
     
     gsInfo << "Setting up participant: " << participantName << "\n";
     gsInfo << "Using preCICE config file: " << configFile << "\n";
     
     // Initialize preCICE
     gsPreCICE<real_t> interface(participantName, configFile);
     
     // Set up meshes
     gsInfo << "Setting up meshes\n";
     for (const auto& mesh : config["meshes"]) {
         std::string meshName = mesh["name"].get<std::string>();
         
         // Create mesh points based on configuration
         gsMatrix<> points;
         if (mesh["type"].get<std::string>() == "structured") {
             // Generate structured grid
             int nx = mesh["nx"].get<int>();
             int ny = mesh["ny"].get<int>();
             int dim = mesh["dim"].get<int>();
             
             points.resize(dim, nx * ny);
             
             gsVector<> xRange(2);
             gsVector<> yRange(2); 
             
             if (mesh["x-range"].is_array()) {
                 xRange(0) = mesh["x-range"][0].get<real_t>();
                 xRange(1) = mesh["x-range"][1].get<real_t>();
             } else {
                 xRange(0) = mesh["x-range"]["data"][0].get<real_t>();
                 xRange(1) = mesh["x-range"]["data"][1].get<real_t>();
             }

             if (mesh["y-range"].is_array()) {
                 yRange(0) = mesh["y-range"][0].get<real_t>();
                 yRange(1) = mesh["y-range"][1].get<real_t>();
             } else {
                 yRange(0) = mesh["y-range"]["data"][0].get<real_t>();
                 yRange(1) = mesh["y-range"]["data"][1].get<real_t>();
             }
             
             real_t dx = (xRange[1] - xRange[0]) / (nx - 1);
             real_t dy = (yRange[1] - yRange[0]) / (ny - 1);
             
             int idx = 0;
             for (int i = 0; i < nx; i++) {
                 for (int j = 0; j < ny; j++) {
                     points(0, idx) = xRange[0] + i * dx;
                     points(1, idx) = yRange[0] + j * dy;
                     if (dim > 2) points(2, idx) = 0.0;
                     idx++;
                 }
             }
         } 
         else if (mesh["type"].get<std::string>() == "custom") {
             // Debug output to help identify structure issues
             gsInfo << "Processing custom mesh: " << meshName << "\n";
             
             // Load custom points from config safely
             try {
                 points = safeGetMatrix(mesh["points"]);
                 if (points.size() == 0) {
                     gsWarn << "Invalid custom mesh points for mesh: " << meshName << "\n";
                     continue; // Skip this mesh
                 }
             } catch (const std::exception& e) {
                 gsWarn << "Error parsing custom mesh points: " << e.what() << "\n";
                 continue; // Skip this mesh
             }
         }
         
         // Check if this participant should provide this mesh
         if (participantName == "Solver1" && meshName == "Mesh2") {
             gsInfo << "Skipping setup for Mesh2 (received from Solver2)\n";
         } else {
             // Attempt to add the mesh and catch any exceptions
             try {
                 gsVector<index_t> vertexIDs;
                 interface.addMesh(meshName, points, vertexIDs);
                 gsInfo << "Added mesh: " << meshName << " with " << points.cols() << " points\n";
             } catch (const std::exception& e) {
                 gsWarn << "Could not add mesh " << meshName << ": " << e.what() << "\n";
                 gsWarn << "This mesh is likely received from another participant\n";
             }
         }
     }
 
     // Initialize the coupling
     real_t dt = interface.initialize();
     real_t t = 0;
     
     gsInfo << "Starting time loop with initial dt = " << dt << "\n";
     
     // Time loop
     int timeStep = 0;
     while (interface.isCouplingOngoing()) {
         timeStep++;
         gsInfo << "Time step: " << timeStep << ", t = " << t << "\n";
         
         // Process write data actions
         for (const auto& writeData : config["writeData"]) {
             std::string meshName = writeData["mesh"].get<std::string>();
             std::string dataName = writeData["name"].get<std::string>();
             
             // Get points for writing data
             gsMatrix<> coords;
             if (writeData["coords"]["type"].get<std::string>() == "mesh") {
                 gsVector<index_t> ids;
                 interface.getMeshVertexIDsAndCoordinates(meshName, ids, coords);
             } else {
                 try {
                     if (writeData["coords"].contains("points")) {
                         coords = safeGetMatrix(writeData["coords"]["points"]);
                         if (coords.size() == 0) {
                             gsWarn << "Empty coordinate points for " << dataName << ", skipping\n";
                             continue;
                         }
                     } else {
                         gsWarn << "Missing 'points' field in coordinates\n";
                         continue;
                     }
                 } catch (const std::exception& e) {
                     gsWarn << "Error parsing coordinate points: " << e.what() << "\n";
                     continue;
                 }
             }
             
             // Generate sample data (this would be your actual data in a real application)
             gsMatrix<> values;
             std::string dataType = writeData["dataType"].get<std::string>();
             
             if (dataType == "constant") {
                 real_t value = writeData["value"].get<real_t>();
                 values.resize(1, coords.cols());
                 values.setConstant(value);
             } else if (dataType == "timeDependent") {
                 values.resize(1, coords.cols());
                 for (index_t i = 0; i < coords.cols(); i++) {
                     values(0, i) = sin(t + coords(0, i)); // Example function of time and position
                 }
             } else if (dataType == "custom") {
                 if (!writeData.contains("values")) {
                     gsWarn << "Custom data type missing 'values' field\n";
                     values.resize(1, coords.cols());
                     values.setZero();
                 } else {
                     values = safeGetMatrix(writeData["values"]);
                     
                     // Ensure values has appropriate dimensions for writing
                     if (values.cols() != coords.cols()) {
                         gsWarn << "Custom values matrix has incorrect dimensions. "
                                << "Expected cols: " << coords.cols() 
                                << ", actual: " << values.cols() << "\n";
                                
                         // Fix the matrix size
                         gsMatrix<> corrected(1, coords.cols());
                         corrected.setZero();
                         
                         // Copy available values
                         index_t copySize = std::min(values.cols(), coords.cols());
                         for (index_t i = 0; i < copySize; ++i) {
                             corrected(0, i) = (values.rows() > 0 && values.cols() > i) ? values(0, i) : 0.0;
                         }
                         
                         values = corrected;
                     }
                 }
             }
             
             // Write data to preCICE
             interface.writeData(meshName, dataName, coords, values);
             gsInfo << "Wrote data: " << dataName << " to mesh: " << meshName << "\n";
         }
         
         // Handle checkpoints if needed
         if (interface.requiresWritingCheckpoint()) 
         {
             gsInfo << "Writing checkpoint\n";
             // In a real application, you would save your state here
         }
         
         // Advance coupling
         real_t precice_dt = interface.advance(dt);
         dt = std::min(precice_dt, dt); // Use the smaller timestep
         
         // Process read data actions
         for (const auto& readData : config["readData"]) {

            try {
                 std::string meshName = readData["mesh"].get<std::string>();
                 std::string dataName = readData["name"].get<std::string>();
                 
                 // Check if coords field is properly defined
                 if (!readData.contains("coords")) {
                     gsWarn << "Missing 'coords' field in readData for " << dataName << "\n";
                     continue;
                 }
                 
                 // Get points for reading data
                 gsMatrix<> coords;
                 if (readData["coords"]["type"].get<std::string>() == "mesh") {
                     gsVector<index_t> ids;
                     try {
                         interface.getMeshVertexIDsAndCoordinates(meshName, ids, coords);
                     } catch (const std::exception& e) {
                         gsWarn << "Error getting mesh coordinates: " << e.what() << "\n";
                         continue;
                     }
                 } else {
                     try {
                         if (readData["coords"].contains("points")) {
                             coords = safeGetMatrix(readData["coords"]["points"]);
                         } else {
                             gsWarn << "Missing 'points' field in coordinates\n";
                             continue;
                         }
                     } catch (const std::exception& e) {
                         gsWarn << "Error parsing coordinate points: " << e.what() << "\n";
                         continue;
                     }
                 }
                 
                 // Skip if empty coordinates
                 if (coords.cols() == 0) {
                     gsWarn << "Empty coordinates for " << dataName << ", skipping\n";
                     continue;
                 }
                 
                 // Read data safely
                 gsMatrix<> values;
                 try {
                     interface.readData(meshName, dataName, coords, values);
                     
                     gsInfo << "Read data: " << dataName << " from mesh: " << meshName << "\n";
                     if (values.size() > 0) {
                         gsInfo << "First value: " << values(0, 0) << "\n";
                     }
                 } catch (const std::exception& e) {
                     gsWarn << "Error reading data: " << e.what() << "\n";
                 }
             } catch (const std::exception& e) {
                 gsWarn << "Error processing readData entry: " << e.what() << "\n";
             }
         }
         
         if (interface.requiresReadingCheckpoint()) {
             gsInfo << "Reading checkpoint\n";
             // In a real application, you would restore your state here
         } else {
             // Only update time if we're not resetting to a checkpoint
             t += dt;
         }
     }
 
     // Finalize the coupling
     interface.finalize();
     gsInfo << "Coupling finalized\n";
     
     // Print timing information
     gsInfo << "Timing statistics:\n";
     gsInfo << "- Read time: " << interface.readTime() << " seconds\n";
     gsInfo << "- Write time: " << interface.writeTime() << " seconds\n";
     gsInfo << "- Initialize time: " << interface.initializeTime() << " seconds\n";
     
     return EXIT_SUCCESS;
 }