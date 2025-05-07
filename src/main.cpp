#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include <map>
#include <variant>
#include <algorithm>
#include <stdexcept>
#include <atomic>
#include <omp.h>
#include <sys/stat.h>

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
#include <RDGeneral/RDLog.h>

// Local Descriptor Includes
#include "common.hpp" // Defining desfact::DescriptorResult
#include "registry.hpp"

// Platform specific includes for terminal
#include <strings.h>
#include <sys/ioctl.h>
#include <unistd.h>

// CSV parser includes
#include "csv_parser.h"
#include "csv_writer.h"

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

// Modify drawProgressBar function to show rows/total and use percentage for the bar
void drawProgressBar(long long current, long long total, 
                     const std::chrono::steady_clock::time_point &startTime, 
                     bool quietMode = false) {
    // Skip in quiet mode
    if (quietMode) return;
    
    static auto lastUpdate = std::chrono::steady_clock::now();
    auto now = std::chrono::steady_clock::now();
    
    // Limit refresh rate to avoid flicker but ensure smooth updates
    // Use 1ms as update frequency (1000 updates per second) for ultra-smooth display
    if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastUpdate).count() < 1) {
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

// Remove batch processing constants
const int MAX_NUM_THREADS = 128;    // Allow using more threads if available

// Add a function to count total lines in a CSV file
size_t countCSVRows(const std::string& inputFile) {
    FILE* file = fopen(inputFile.c_str(), "rb");
    if (!file) return 0;
    
    size_t lineCount = 0;
    char buffer[8192];
    bool inQuote = false;
    
    while (!feof(file)) {
        size_t bytesRead = fread(buffer, 1, sizeof(buffer), file);
        for (size_t i = 0; i < bytesRead; i++) {
            if (buffer[i] == '"') {
                inQuote = !inQuote;
            } else if (buffer[i] == '\n' && !inQuote) {
                lineCount++;
            }
        }
    }
    
    fclose(file);
    
    // If file doesn't end with newline, add one for the last line
    if (lineCount > 0 && !feof(file)) {
        lineCount++;
    }
    
    // Subtract 1 for the header if present (assume header is always present in this app)
    return lineCount > 0 ? lineCount - 1 : 0;
}

// Modified ProgressTracker class to track total rows
class ProgressTracker {
private:
    std::chrono::steady_clock::time_point startTime;
    std::chrono::steady_clock::time_point lastUpdateTime;
    size_t totalRows;  // Total number of rows (excluding header)
    std::atomic<size_t> processedLines;
    std::atomic<size_t> successCount;
    std::atomic<size_t> errorCount;
    bool quietMode;
    
public:
    ProgressTracker(const std::string& inputFile, bool quiet = false) 
        : startTime(std::chrono::steady_clock::now())
        , lastUpdateTime(startTime)
        , processedLines(0)
        , successCount(0)
        , errorCount(0)
        , quietMode(quiet) {
        
        // Count total rows in the CSV file
        totalRows = countCSVRows(inputFile);
        
        if (!quietMode && totalRows > 0) {
            std::cout << "Found " << totalRows << " rows to process." << std::endl;
        }
    }
    
    void updateProgress(bool success = true) {
        processedLines++;
        if (success) {
            successCount++;
        } else {
            errorCount++;
        }
        
        // Always update progress for each molecule processed for ultra-smooth display
        if (!quietMode) {
            drawProgress();
            lastUpdateTime = std::chrono::steady_clock::now();
        }
    }
    
    void drawProgress() const {
        auto now = std::chrono::steady_clock::now();
        // Pass processed lines and total rows directly
        drawProgressBar(static_cast<long long>(processedLines.load()), 
                        static_cast<long long>(totalRows), 
                        startTime, quietMode);
    }
    
    void printSummary() const {
        auto endTime = std::chrono::steady_clock::now();
        double totalTime = std::chrono::duration<double>(endTime - startTime).count();
        double speed = totalTime > 0 ? processedLines / totalTime : 0;
        
        std::cout << "\r" << std::string(getTerminalWidth(), ' ') << "\r";
        std::cout << Color::Bold << Color::Green << "✓ Processing complete!" << Color::Reset << " ";
        std::cout << Color::Bold << processedLines << Color::Reset << " molecules processed ";
        std::cout << "(" << Color::Green << successCount << Color::Reset << " successful, ";
        if (errorCount > 0) {
            std::cout << Color::Red << errorCount << Color::Reset << " errors) in ";
        } else {
            std::cout << "0 errors) in ";
        }
        std::cout << Color::Magenta << formatTime(totalTime) << Color::Reset;
        std::cout << " (" << std::fixed << std::setprecision(1) << speed << " mol/s)\n";
    }
    
    size_t getProcessedLines() const { return processedLines; }
    size_t getSuccessCount() const { return successCount; }
    size_t getErrorCount() const { return errorCount; }
};

// Main processing function
void processCSV(const std::string& inputFile, const std::string& outputFile,
                const std::vector<std::string>& descriptorNames,
                int smilesColIndex, bool verboseMode, bool quietMode) {
    
    // Initialize progress tracking first (this will count total rows)
    ProgressTracker progress(inputFile, quietMode);
    
    // Initialize CSV parser
    CSVParser* parser = csv_parser_init(inputFile.c_str(), ',', true, 
        [](CSVErrorCode error, size_t line, const char* msg, void* ctx) {
            std::cerr << "CSV Error (line " << line << "): " << msg << std::endl;
        }, nullptr, nullptr, nullptr);
    
    if (!parser) {
        throw std::runtime_error("Failed to initialize CSV parser");
    }
    
    CSVWriter* writer = csv_writer_init(outputFile.c_str(), ',',
        [](CSVWriteErrorCode error, const char* msg, void*) {
            std::cerr << "CSV Write Error: " << msg << std::endl;
        }, nullptr);
    
    if (!writer) {
        csv_parser_free(parser);
        throw std::runtime_error("Failed to initialize CSV writer");
    }
    
    // Write headers
    std::vector<const char*> headers;
    for (size_t i = 0; i < csv_parser_get_header_count(parser); i++) {
        headers.push_back(csv_parser_get_header_name(parser, i));
    }
    for (const auto& desc : descriptorNames) {
        headers.push_back(desc.c_str());
    }
    csv_writer_write_headers(writer, headers.data(), headers.size());
    
    // Initialize descriptor pipelines for each thread
    int numThreads = std::min(omp_get_max_threads(), MAX_NUM_THREADS);
    omp_set_num_threads(numThreads);
    
    if (verboseMode) {
        std::cout << "Using " << numThreads << " threads for processing" << std::endl;
    }

    std::vector<desfact::DescriptorPipeline> threadPipelines(numThreads);
    for (int t = 0; t < numThreads; t++) {
        for (const auto& name : descriptorNames) {
            threadPipelines[t].addDescriptor(name);
        }
    }
    
    // Process file line by line
    #pragma omp parallel
    {
        int threadId = omp_get_thread_num();
        auto& pipeline = threadPipelines[threadId];
        // Store full strings first to manage lifetime
        std::vector<std::string> outputFieldsStrings;
        outputFieldsStrings.reserve(headers.size());
        // Vector to hold the const char* pointers for the C API call
        std::vector<const char*> outputFieldsPtrs;
        outputFieldsPtrs.reserve(headers.size());

        while (true) {
            CSVLine* line = nullptr;
            size_t lineSize = 0;

            // Read next line
            #pragma omp critical(csv_read)
            {
                line = csv_parser_next_line(parser);
                if (line) {
                    // Estimate line size based on current position for progress (less accurate but avoids full read)
                    // Note: This is an approximation as line->line_number isn't directly bytes read
                    lineSize = csv_parser_ftell(parser);
                }
            }

            if (!line) break;

            // Process the line
            std::string smiles = line->field_count > smilesColIndex ? csv_line_get_field(line, smilesColIndex) : "";
            bool success = false;
            std::map<std::string, desfact::DescriptorResult> results; // Hold results here

            outputFieldsStrings.clear(); // Clear for the new row

            // Copy original fields
            for (size_t i = 0; i < line->field_count; i++) {
                 // Ensure field is not NULL before adding
                 const char* field_ptr = csv_line_get_field(line, i);
                 outputFieldsStrings.push_back(field_ptr ? field_ptr : "");
            }

            if (!smiles.empty()) {
                try {
                    // Use RDKit functionalities within try-catch
                    RDKit::ROMol* molPtr = RDKit::SmilesToMol(smiles);
                    std::unique_ptr<RDKit::ROMol> mol(molPtr); // Manage lifetime

                    if (mol) {
                        // Calculate descriptors
                        results = pipeline.process(mol.get());
                        success = true; // Assume success unless exception
                    }
                } catch (const std::exception& e) {
                    success = false; // Mark as failed on exception
                    if (verboseMode) {
                        #pragma omp critical(error_log)
                        std::cerr << "Error processing line " << line->line_number << " (SMILES: " << smiles << "): " << e.what() << std::endl;
                    }
                } catch (...) {
                     success = false; // Mark as failed on unknown exception
                    if (verboseMode) {
                        #pragma omp critical(error_log)
                        std::cerr << "Unknown error processing line " << line->line_number << " (SMILES: " << smiles << ")" << std::endl;
                    }
                }
            }

            // Add descriptor results (or placeholders if failed)
            for (const auto& name : descriptorNames) {
                if (success) {
                    auto it = results.find(name);
                    if (it != results.end()) {
                        // Store the actual string
                        outputFieldsStrings.push_back(safeVariantToString(it->second));
                    } else {
                        outputFieldsStrings.push_back("NA_Missing"); // Indicate descriptor wasn't found in results
                    }
                } else {
                    // Add placeholders for failed lines
                    outputFieldsStrings.push_back("NA_Error");
                }
            }

            // Prepare the C-style array of pointers just before writing
            outputFieldsPtrs.clear();
            for (const auto& str : outputFieldsStrings) {
                outputFieldsPtrs.push_back(str.c_str());
            }

            // Write output safely using the pointers derived from stable strings
            #pragma omp critical(csv_write)
            {
                if (!csv_writer_write_row(writer, outputFieldsPtrs.data(), outputFieldsPtrs.size())) {
                     // Optional: Log write error
                     if(verboseMode) {
                         std::cerr << "CSV Write Error occurred for line " << line->line_number << std::endl;
                     }
                }
                // Update progress: use actual line number and success status
                progress.updateProgress(success); // Pass success flag, file position is tracked internally
            }

            csv_line_free(line);
        }
    }
    
    // Cleanup
    csv_writer_flush(writer);
    csv_writer_free(writer);
    csv_parser_free(parser);
    
    // Print final summary
    if (!quietMode) {
        progress.printSummary();
    }
}

int main(int argc, char **argv) {
    // Move these variable declarations to the top of main
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
                std::cout << "" << Color::Green << std::left << std::setw(25) << desc.first 
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

    // Get SMILES column index
    int smilesColIndex = 0;
    if (!smilesColSpec.empty()) {
        try {
            // Check if it's a number
            smilesColIndex = std::stoi(smilesColSpec);
        } catch (const std::invalid_argument&) {
            // It's a column name, look it up in the headers
            CSVParser* parser = csv_parser_init(inputFile.c_str(), ',', true, 
                [](CSVErrorCode error, size_t line, const char* msg, void* ctx) {
                    std::cerr << "CSV Error (line " << line << "): " << msg << std::endl;
                }, nullptr, nullptr, nullptr);
            if (parser) {
                smilesColIndex = csv_parser_get_header_index(parser, smilesColSpec.c_str());
                csv_parser_free(parser);
                if (smilesColIndex < 0) {
                    std::cerr << "Error: SMILES column name '" << smilesColSpec << "' not found in CSV header." << std::endl;
                    return 1;
                }
            } else {
                std::cerr << "Error: Failed to initialize CSV parser for column lookup." << std::endl;
                return 1;
            }
        }
    } else {
        // Try to find a column named "SMILES"
        CSVParser* parser = csv_parser_init(inputFile.c_str(), ',', true, 
            [](CSVErrorCode error, size_t line, const char* msg, void* ctx) {
                std::cerr << "CSV Error (line " << line << "): " << msg << std::endl;
            }, nullptr, nullptr, nullptr);
        if (parser) {
            smilesColIndex = csv_parser_get_header_index(parser, "SMILES");
            csv_parser_free(parser);
            if (smilesColIndex < 0) smilesColIndex = 0; // Default to first column
        } else {
            std::cerr << "Error: Failed to initialize CSV parser for default SMILES column lookup." << std::endl;
            return 1;
        }
    }

    try {
        processCSV(inputFile, outputFile, finalDescriptorNames, 
                  smilesColIndex, verboseMode, quietMode);
        std::cout << Color::Bold << Color::Green 
                  << "Program completed successfully." << Color::Reset << std::endl;
    } catch (const std::exception& e) {
        std::cerr << Color::Bold << Color::Red 
                  << "Error: " << e.what() << Color::Reset << std::endl;
        return 1;
    }
    
    return 0;
}