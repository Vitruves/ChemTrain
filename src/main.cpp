#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <map>
#include <variant>
#include <algorithm>
#include <stdexcept>
#include <atomic>
#include <omp.h>

// RDKit Includes
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/types.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/PartialCharges/GasteigerParams.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/ShapeHelpers/ShapeEncoder.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <RDGeneral/BoostStartInclude.h>
#include <RDGeneral/BoostEndInclude.h>

// Local Descriptor Includes
#include "common.hpp" // Defining desfact::DescriptorResult
#include "registry.hpp"

// Platform specific includes for terminal
#include <strings.h>
#include <sys/ioctl.h>
#include <unistd.h>

// CSV parser includes
#include "csv_parser.h"
#include "csv_wrapper.hpp"

// Terminal colors for pretty output
namespace Color {
    const std::string Reset   = "\033[0m";
    const std::string Green   = "\033[32m";
    const std::string Yellow  = "\033[33m";
    const std::string Blue    = "\033[34m";
    const std::string Magenta = "\033[35m";
    const std::string Cyan    = "\033[36m";
    const std::string Red     = "\033[31m";
    const std::string Bold    = "\033[1m";
}

// Get terminal width for progress bar
int getTerminalWidth() {
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_col > 0 ? w.ws_col : 80;
}

// Format time in human-readable format
std::string formatTime(double seconds) {
    if (seconds < 60) return std::to_string(static_cast<int>(seconds)) + "s";
    if (seconds < 3600) {
        int mins = static_cast<int>(seconds) / 60;
        int secs = static_cast<int>(seconds) % 60;
        return std::to_string(mins) + "m " + std::to_string(secs) + "s";
    }
    int hours = static_cast<int>(seconds) / 3600;
    int mins = (static_cast<int>(seconds) % 3600) / 60;
    return std::to_string(hours) + "h " + std::to_string(mins) + "m";
}

// Calculate visual length of a string (accounting for ANSI color codes)
size_t visibleLength(const std::string &str) {
    size_t len = 0;
    bool inEscapeSequence = false;
    for (char c : str) {
        if (c == '\033') {
            inEscapeSequence = true;
            continue;
        }
        if (inEscapeSequence) {
            if (c == 'm')
                inEscapeSequence = false;
            continue;
        }
        len++;
    }
    return len;
}

// Draw a nice progress bar
void drawProgressBar(long long current, long long total, 
                     const std::chrono::steady_clock::time_point &startTime, 
                     bool quietMode = false) {
    // Skip in quiet mode
    if (quietMode) return;
    
    static auto lastUpdate = std::chrono::steady_clock::now();
    auto now = std::chrono::steady_clock::now();
    
    // Limit refresh rate to avoid flicker but ensure smooth updates
    // Use 100ms as update frequency (10 updates per second) 
    if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastUpdate).count() < 100) {
        return;
    }
    lastUpdate = now;
    
    // Clear the current line completely 
    std::cout << "\r" << std::string(getTerminalWidth(), ' ') << "\r";
    
    // Calculate progress for the bar
    bool knownTotal = total > 0;
    double progress = knownTotal ? static_cast<double>(current) / total : 0;
    if (progress > 1.0) progress = 1.0;
    
    // Calculate metrics
    double elapsed = std::chrono::duration<double>(now - startTime).count();
    double rate = elapsed > 0 ? current / elapsed : 0;
    std::string timeInfo;
    
    if (knownTotal && progress > 0) {
        double estimated = elapsed / progress;
        double remaining = estimated - elapsed;
        timeInfo = "ETA: " + formatTime(remaining);
    } else {
        timeInfo = "Time: " + formatTime(elapsed);
    }
    
    // Progress bar components
    const int barWidth = 25;
    
    // Format percentage with 0.01% precision
    std::string percentDisplay;
    if (knownTotal) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(1) << (progress * 100.0) << "%";
        percentDisplay = ss.str();
    } else {
        percentDisplay = "??%";
    }
    
    // Color based on progress
    std::string percentColor = progress < 0.33 ? Color::Yellow : 
                               (progress < 0.66 ? Color::Cyan : Color::Green);
    
    std::string prefix = Color::Bold + "⚙ Process:" + Color::Reset;
    std::string stats = " " + Color::Bold + (knownTotal ? 
                 std::to_string(current) + "/" + std::to_string(total) : 
                 std::to_string(current)) + Color::Reset;
    
    std::stringstream rateStream;
    rateStream << std::fixed << std::setprecision(1) << rate;
    std::string rateInfo = rateStream.str() + " mol/s";
    
    // Create progress bar
    std::string progressBar;
    if (knownTotal) {
        // Unicode blocks for smooth progress bar
        static const char* blocks[] = {"", "▏", "▎", "▍", "▌", "▋", "▊", "▉"};
        double barBlocks = barWidth * 1.0;
        int completeBlocks = static_cast<int>(progress * barBlocks);
        int partialBlock = static_cast<int>((progress * barBlocks - completeBlocks) * 8);
        
        progressBar = "[" + Color::Green;
        for (int i = 0; i < completeBlocks; i++) progressBar += "█";
        if (completeBlocks < barWidth) {
            progressBar += blocks[partialBlock];
            progressBar += Color::Reset;
            for (int i = completeBlocks + 1; i < barWidth; i++) progressBar += " ";
        }
        progressBar += Color::Reset + "]";
    } else {
        // Spinner for indeterminate progress
        static int spinnerPos = 0;
        static const std::vector<std::string> spinnerFrames = {"|", "/", "-", "\\"};
        std::string spinner = spinnerFrames[spinnerPos++ % spinnerFrames.size()];
        
        progressBar = "[" + std::string(barWidth / 2 - 1, ' ') + 
                        Color::Bold + spinner + Color::Reset + 
                        std::string(barWidth / 2, ' ') + "]";
    }
    
    // Format into fixed-width sections
    std::cout << std::left 
               << std::setw(15) << prefix 
               << " " << percentColor << std::setw(7) << percentDisplay << Color::Reset
               << " " << progressBar 
               << " " << std::setw(15) << stats
               << " " << std::setw(12) << rateInfo
               << " " << std::setw(12) << timeInfo
               << std::flush;
}

// Convert descriptor result to string in a safe way
std::string safeVariantToString(const desfact::DescriptorResult &result, size_t maxLen = 500) {
    try {
        std::string resultStr;
        if (std::holds_alternative<double>(result)) {
            double value = std::get<double>(result);
            if (std::isnan(value)) return "NA";
            if (std::isinf(value)) return value > 0 ? "INF" : "-INF";
            std::stringstream ss;
            ss << std::fixed << std::setprecision(4) << value;
            resultStr = ss.str();
        } else if (std::holds_alternative<int>(result))
            resultStr = std::to_string(std::get<int>(result));
        else if (std::holds_alternative<std::string>(result))
            resultStr = std::get<std::string>(result);
        else if (std::holds_alternative<std::vector<double>>(result)) {
            const auto &vec = std::get<std::vector<double>>(result);
            std::stringstream ss;
            ss << "[";
            for (size_t i = 0; i < vec.size(); ++i) {
                if (i > 0) ss << ",";
                if (std::isnan(vec[i])) ss << "NA";
                else if (std::isinf(vec[i])) ss << (vec[i] > 0 ? "INF" : "-INF");
                else ss << std::fixed << std::setprecision(4) << vec[i];
            }
            ss << "]";
            resultStr = ss.str();
        } else
            resultStr = "NA";
         
        if (resultStr.length() > maxLen) {
            return resultStr.substr(0, maxLen) + "...";
        }
        return resultStr;
    } catch (...) {
        return "ERROR";
    }
}

// Help message function
void printUsage() {
    const int colWidth = 28;
    std::cout << Color::Bold << Color::Cyan << "DESCRIPTOR FACTORY" << Color::Reset << "\n";
    std::cout << "Simple molecular descriptor calculation utility\n\n";
    
    std::cout << Color::Bold << "USAGE:" << Color::Reset << "\n";
    std::cout << "  " << Color::Cyan << "descriptorfactory [OPTIONS]" << Color::Reset << "\n\n";
    
    std::cout << Color::Bold << "OPTIONS:" << Color::Reset << "\n";
    std::cout << "  " << Color::Yellow << std::left << std::setw(colWidth) << "-h, --help" << Color::Reset << "Display this help message\n";
    std::cout << "  " << Color::Yellow << std::left << std::setw(colWidth) << "-i, --input <file>" << Color::Reset << "Input CSV file\n";
    std::cout << "  " << Color::Yellow << std::left << std::setw(colWidth) << "-o, --output <file>" << Color::Reset << "Output CSV file\n";
    std::cout << "  " << Color::Yellow << std::left << std::setw(colWidth) << "-v, --verbose" << Color::Reset << "Show more detailed output\n";
    std::cout << "  " << Color::Yellow << std::left << std::setw(colWidth) << "-q, --quiet" << Color::Reset << "Minimal output (no progress bar)\n";
    std::cout << "  " << Color::Yellow << std::left << std::setw(colWidth) << "--list-descriptors" << Color::Reset << "List all available descriptors\n";
    std::cout << "\n";
    
    std::cout << Color::Bold << "ADVANCED OPTIONS:" << Color::Reset << "\n";
    std::cout << "  " << Color::Yellow << std::left << std::setw(colWidth) << "-d, --descriptors <list>" << Color::Reset << "Comma-separated list of descriptors\n";
    std::cout << "  " << Color::Yellow << std::left << std::setw(colWidth) << "--smiles-col <name/idx>" << Color::Reset << "Name or index of SMILES column\n";
    std::cout << "\n";
    
    std::cout << Color::Bold << "EXAMPLES:" << Color::Reset << "\n";
    std::cout << "  " << Color::Cyan << "descriptorfactory -i input.csv -o output.csv" << Color::Reset << "\n";
    std::cout << "  " << Color::Cyan << "descriptorfactory --list-descriptors" << Color::Reset << "\n";
    std::cout << "  " << Color::Cyan << "descriptorfactory -i input.csv -o output.csv -d MolWt,LogP" << Color::Reset << "\n";
}

// Parse a CSV line handling quotes properly
std::vector<std::string> parseCSVLine(const std::string& line, char delimiter = ',') {
    std::vector<std::string> fields;
    std::string field;
    bool inQuotes = false;
    
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        
        if (c == '"') {
            if (inQuotes && i + 1 < line.size() && line[i + 1] == '"') {
                // Handle escaped quotes (double quotes)
                field += '"';
                ++i; // Skip the next quote
            } else {
                // Toggle quote mode
                inQuotes = !inQuotes;
            }
        } else if (c == delimiter && !inQuotes) {
            // End of field
            fields.push_back(field);
            field.clear();
        } else {
            field += c;
        }
    }
    
    // Add the last field
    fields.push_back(field);
    return fields;
}

// Constants for parallel processing
const size_t BATCH_SIZE = 5000;  // Larger batch size for better performance
const size_t MIN_BATCH_SIZE = 500; // Minimum batch size 
const int MAX_NUM_THREADS = 32;    // Allow using more threads if available

int main(int argc, char **argv) {
    std::string inputFile, outputFile;
    std::vector<std::string> requestedDescriptorNames;
    std::string smilesColSpec = "SMILES";
    bool listDescriptors = false, verboseMode = false;
    bool quietMode = false;
    
    // Check for help flag first
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printUsage();
            return 0;
        }
    }
    
    // Parse all other arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-i" || arg == "--input") && i + 1 < argc)
            inputFile = argv[++i];
        else if ((arg == "-o" || arg == "--output") && i + 1 < argc)
            outputFile = argv[++i];
        else if ((arg == "-d" || arg == "--descriptors") && i + 1 < argc) {
            std::string descriptorArg = argv[++i];
            size_t pos = 0;
            while ((pos = descriptorArg.find(',')) != std::string::npos) {
                requestedDescriptorNames.push_back(descriptorArg.substr(0, pos));
                descriptorArg.erase(0, pos + 1);
            }
            if (!descriptorArg.empty())
                requestedDescriptorNames.push_back(descriptorArg);
        }
        else if (arg == "--smiles-col" && i + 1 < argc)
            smilesColSpec = argv[++i];
        else if (arg == "--list-descriptors")
            listDescriptors = true;
        else if (arg == "-v" || arg == "--verbose")
            verboseMode = true;
        else if (arg == "-q" || arg == "--quiet")
            quietMode = true;
        else
            std::cerr << "Warning: Ignoring unknown argument: " << arg << std::endl;
    }
    
    // Validate required arguments
    if ((inputFile.empty() || outputFile.empty()) && !listDescriptors) {
        std::cerr << "Error: Input and output files are required.\n";
        printUsage();
        return 1;
    }

    // Initialize descriptor registry
    desfact::initializeRegistry();
    auto &registry = desfact::DescriptorRegistry::getInstance();
    
    // List descriptors if requested
    if (listDescriptors) {
        std::cout << Color::Bold << "Available Descriptors:" << Color::Reset << "\n";
        std::map<std::string, std::vector<std::pair<std::string, std::string>>> categorized;
        
        // Group descriptors by category
        for (const auto &name : registry.getDescriptors()) {
            auto descriptor = registry.getDescriptor(name);
            if (!descriptor) continue;
            categorized[descriptor->getCategory()].push_back({name, descriptor->getDescription()});
        }
        
        // Print descriptors by category
        for (const auto &pair : categorized) {
            std::cout << "  [" << Color::Yellow << pair.first << Color::Reset << "]\n";
            for (const auto &desc : pair.second) {
                std::cout << "    " << Color::Green << std::left << std::setw(25) << desc.first 
                          << Color::Reset << " " << desc.second << "\n";
            }
            std::cout << "\n";
        }
        return 0;
    }

    // Determine which descriptors to calculate
    std::vector<std::string> finalDescriptorNames;
    if (requestedDescriptorNames.empty()) {
        // If no descriptors specified, use all available ones
        for (const auto &name : registry.getDescriptors())
            finalDescriptorNames.push_back(name);
        std::sort(finalDescriptorNames.begin(), finalDescriptorNames.end());
    } else {
        // Use only requested descriptors
        for (const auto &reqName : requestedDescriptorNames) {
            if (registry.getDescriptor(reqName))
                finalDescriptorNames.push_back(reqName);
            else
                std::cerr << "Warning: Descriptor '" << reqName << "' not found, skipping.\n";
        }
    }

    std::cout << Color::Bold << Color::Cyan << "Calculating " << finalDescriptorNames.size() 
               << " descriptors for input: " << inputFile << Color::Reset << std::endl;

    // Initialize the CSV parser with error handling
    csv::Parser csvParser(inputFile, ',', true, [](CSVErrorCode error, size_t line, const std::string& msg) {
        std::cerr << "CSV Error (line " << line << "): " << msg << std::endl;
    });

    // Open output file
    std::ofstream outputStream(outputFile);
    if (!outputStream.is_open()) {
        std::cerr << "Error: Cannot open output file: " << outputFile << std::endl;
        return 1;
    }

    // Create a larger buffer for file operations to improve I/O performance
    const size_t FILE_BUFFER_SIZE = 1048576; // 1MB buffer
    std::vector<char> outputBuffer(FILE_BUFFER_SIZE);
    outputStream.rdbuf()->pubsetbuf(outputBuffer.data(), FILE_BUFFER_SIZE);

    // Get SMILES column index
    int smilesColIndex = 0;
    if (!smilesColSpec.empty()) {
        try {
            // Check if it's a number
            smilesColIndex = std::stoi(smilesColSpec);
        } catch (const std::invalid_argument&) {
            // It's a column name, look it up in the headers
            smilesColIndex = csvParser.getHeaderIndex(smilesColSpec);
            if (smilesColIndex < 0) {
                std::cerr << "Error: SMILES column name '" << smilesColSpec << "' not found in CSV header." << std::endl;
                return 1;
            }
        }
    } else {
        // Try to find a column named "SMILES"
        smilesColIndex = csvParser.getHeaderIndex("SMILES");
        if (smilesColIndex < 0) smilesColIndex = 0; // Default to first column
    }

    // Write original header + new descriptor columns
    outputStream << csvParser.getHeaderLine();
    for (const auto &name : finalDescriptorNames)
        outputStream << "," << name;
    outputStream << "\n";

    // Set the number of threads for OpenMP
    int numThreads = omp_get_max_threads();
    // Allow more aggressive threading
    if (numThreads > MAX_NUM_THREADS) {
        numThreads = MAX_NUM_THREADS;
    }
    omp_set_num_threads(numThreads);
    
    if (verboseMode) {
        std::cout << "Using " << numThreads << " threads for processing" << std::endl;
    }

    // Processing statistics
    std::atomic<long long> totalProcessed(0);
    std::atomic<long long> successCount(0);
    std::atomic<long long> errorCount(0);
    
    // Estimate total lines for progress bar
    long long totalLinesEstimate = 0;
    if (!quietMode) {
        // Rough estimation of total lines without reading the whole file
        std::ifstream counterStream(inputFile);
        if (counterStream) {
            std::string line;
            // Skip header
            std::getline(counterStream, line);
            // Count lines with a limit to avoid slow startup
            const int MAX_COUNT_LINES = 10000;
            int countedLines = 0;
            while (countedLines < MAX_COUNT_LINES && std::getline(counterStream, line))
                countedLines++;
            
            // If we read the whole file, use the exact count
            if (countedLines < MAX_COUNT_LINES) {
                totalLinesEstimate = countedLines;
            } else {
                // Estimate based on file size
                size_t currentPos = counterStream.tellg();
                counterStream.seekg(0, std::ios::end);
                size_t fileSize = counterStream.tellg();
                
                // Avoid division by zero
                if (currentPos > 0) {
                    double bytesPerLine = static_cast<double>(currentPos) / countedLines;
                    totalLinesEstimate = static_cast<long long>(fileSize / bytesPerLine);
                }
            }
            counterStream.close();
        }
    }

    // Start timing
    auto startTime = std::chrono::steady_clock::now();
    auto lastProgressUpdateTime = startTime;
    
    // Process file in batches for better efficiency
    std::vector<std::stringstream> threadBuffers(numThreads);
    const size_t OUTPUT_FLUSH_THRESHOLD = 10000;
    
    // Create a thread-local pipeline for better cache performance
    std::vector<desfact::DescriptorPipeline> threadPipelines(numThreads);
    for (int t = 0; t < numThreads; t++) {
        for (const auto &name : finalDescriptorNames) {
            threadPipelines[t].addDescriptor(name);
        }
    }
    
    // Process the file in parallel batches
    bool processing = true;
    #pragma omp parallel
    {
        int threadId = omp_get_thread_num();
        auto& localBuffer = threadBuffers[threadId];
        auto& localPipeline = threadPipelines[threadId];
        size_t localBufferLineCount = 0;
        
        while (processing) {
            std::vector<csv::CSVLineWrapper> batch;
            size_t batchSize = 0;
            
            // Each thread gets its own batch
            #pragma omp critical(csv_read)
            {
                batch = csvParser.nextBatch(BATCH_SIZE / numThreads);
                batchSize = batch.size();
                if (batch.empty()) {
                    processing = false;
                }
            }
            
            // Process this thread's batch
            if (batchSize > 0) {
                long long batchStartIdx = 0;
                
                // Replace the atomic capture with a simpler approach
                #pragma omp critical(batch_index)
                {
                    batchStartIdx = totalProcessed;
                    totalProcessed += batchSize;
                }
                
                // Process all the records in this thread's batch
                for (size_t i = 0; i < batchSize; i++) {
                    auto& line = batch[i];
                    
                    // Get SMILES from the CSV line
                    std::string smiles = line.getField(smilesColIndex);
                    
                    // Start building output line with original content
                    std::string outputLine;
                    for (size_t j = 0; j < line.getFieldCount(); j++) {
                        if (j > 0) outputLine += ",";
                        outputLine += line.getField(j);
                    }
                    
                    if (!smiles.empty()) {
                        try {
                            // Parse SMILES to molecule
                            std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
                            
                            if (mol) {
                                // Calculate descriptors
                                auto descriptorResults = localPipeline.process(mol.get());
                                
                                // Output descriptor values
                                for (const auto &name : finalDescriptorNames) {
                                    outputLine += ",";
                                    auto it = descriptorResults.find(name);
                                    if (it != descriptorResults.end()) {
                                        outputLine += safeVariantToString(it->second);
                                    } else {
                                        outputLine += "NA";
                                    }
                                }
                                
                                successCount++;
                            } else {
                                // Failed to parse SMILES
                                for (size_t j = 0; j < finalDescriptorNames.size(); j++) {
                                    outputLine += ",NA_ParseError";
                                }
                                errorCount++;
                            }
                        } catch (const std::exception& e) {
                            // Error during processing
                            for (size_t j = 0; j < finalDescriptorNames.size(); j++) {
                                outputLine += ",NA_Error";
                            }
                            errorCount++;
                            if (verboseMode) {
                                #pragma omp critical(error_log)
                                std::cerr << "Error processing line " << (batchStartIdx + i) << ": " << e.what() << std::endl;
                            }
                        }
                    } else {
                        // Empty SMILES
                        for (size_t j = 0; j < finalDescriptorNames.size(); j++) {
                            outputLine += ",NA_EmptySmiles";
                        }
                        errorCount++;
                    }
                    
                    outputLine += "\n";
                    localBuffer << outputLine;
                    localBufferLineCount++;
                }
                
                // Periodically flush the thread-local buffer to the output file
                if (localBufferLineCount >= OUTPUT_FLUSH_THRESHOLD) {
                    #pragma omp critical(file_write)
                    {
                        outputStream << localBuffer.str();
                    }
                    localBuffer.str("");
                    localBuffer.clear();
                    localBufferLineCount = 0;
                }
                
                // Update progress bar periodically - using time-based updates
                auto now = std::chrono::steady_clock::now();
                if (!quietMode && 
                    std::chrono::duration_cast<std::chrono::milliseconds>(now - lastProgressUpdateTime).count() >= 100) {
                    #pragma omp critical(progress_update)
                    {
                        if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastProgressUpdateTime).count() >= 100) {
                            drawProgressBar(totalProcessed, totalLinesEstimate, startTime, quietMode);
                            lastProgressUpdateTime = now;
                        }
                    }
                }
            }
        }
        
        // Flush any remaining output from this thread
        if (localBufferLineCount > 0) {
            #pragma omp critical(file_write)
            {
                outputStream << localBuffer.str();
            }
        }
    }
    
    // Final progress bar update
    if (!quietMode) {
        drawProgressBar(totalProcessed, totalLinesEstimate, startTime, quietMode);
        std::cout << std::endl;
    }
    
    // Calculate performance metrics
    auto endTime = std::chrono::steady_clock::now();
    double totalTime = std::chrono::duration<double>(endTime - startTime).count();
    double speed = (totalTime > 0 && successCount > 0) ? successCount / totalTime : 0;

    // Print summary
    std::cout << "\r" << std::string(getTerminalWidth(), ' ') << "\r"; // Clear line
    std::cout << Color::Bold << Color::Green << "✓ Processing complete!" << Color::Reset << " ";
    std::cout << Color::Bold << totalProcessed << Color::Reset << " molecules processed ";
    std::cout << "(" << Color::Green << successCount << Color::Reset << " successful, ";
    if (errorCount > 0) {
        std::cout << Color::Red << errorCount << Color::Reset << " errors) in ";
    } else {
        std::cout << "0 errors) in ";
    }
    std::cout << Color::Magenta << formatTime(totalTime) << Color::Reset;
    std::cout << " (" << std::fixed << std::setprecision(1) << speed << " mol/s)\n";
    std::cout << Color::Bold << Color::Green << "Program completed successfully." << Color::Reset << std::endl;

    return 0;
}