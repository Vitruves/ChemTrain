/**
 * C++ wrapper for the high-performance C CSV parser
 */

#ifndef CSV_WRAPPER_HPP
#define CSV_WRAPPER_HPP

#include <string>
#include <vector>
#include <functional>
#include <stdexcept>
#include "csv_parser.h"

namespace csv {

// Error handler callback type
using ErrorCallback = std::function<void(CSVErrorCode, size_t, const std::string&)>;

// RAII wrapper for CSVLine
class CSVLineWrapper {
private:
    CSVLine* line_;
    
public:
    explicit CSVLineWrapper(CSVLine* line) : line_(line) {}
    
    ~CSVLineWrapper() {
        if (line_) {
            csv_line_free(line_);
        }
    }
    
    // Prevent copying
    CSVLineWrapper(const CSVLineWrapper&) = delete;
    CSVLineWrapper& operator=(const CSVLineWrapper&) = delete;
    
    // Allow moving
    CSVLineWrapper(CSVLineWrapper&& other) noexcept : line_(other.line_) {
        other.line_ = nullptr;
    }
    
    CSVLineWrapper& operator=(CSVLineWrapper&& other) noexcept {
        if (this != &other) {
            if (line_) {
                csv_line_free(line_);
            }
            line_ = other.line_;
            other.line_ = nullptr;
        }
        return *this;
    }
    
    // Access the underlying CSVLine
    CSVLine* get() const { return line_; }
    
    // Check if line is valid
    explicit operator bool() const { return line_ != nullptr; }
    
    // Get field by index
    std::string getField(size_t index) const {
        const char* field = csv_line_get_field(line_, index);
        return field ? field : "";
    }
    
    // Get number of fields
    size_t getFieldCount() const {
        return line_ ? line_->field_count : 0;
    }
    
    // Get line number
    size_t getLineNumber() const {
        return line_ ? line_->line_number : 0;
    }
};

// C error handler that forwards to C++ callback
extern "C" {
    static void csv_error_handler_bridge(CSVErrorCode error, size_t line_number, 
                                       const char* message, void* context) {
        auto* callback = static_cast<ErrorCallback*>(context);
        if (callback && *callback) {
            (*callback)(error, line_number, message ? message : "");
        }
    }
}

// CSV Parser wrapper class
class Parser {
private:
    ::CSVParser* parser_;
    ErrorCallback error_callback_;
    
public:
    Parser(const std::string& filename, char delimiter = ',', bool has_header = true,
            ErrorCallback error_callback = nullptr)
        : parser_(nullptr), error_callback_(error_callback) {
        
        parser_ = csv_parser_init(
            filename.c_str(), 
            delimiter,
            has_header,
            error_callback_ ? csv_error_handler_bridge : nullptr,
            error_callback_ ? &error_callback_ : nullptr
        );
        
        if (!parser_) {
            throw std::runtime_error("Failed to initialize CSV parser for file: " + filename);
        }
    }
    
    ~Parser() {
        if (parser_) {
            csv_parser_free(parser_);
        }
    }
    
    // Prevent copying
    Parser(const Parser&) = delete;
    Parser& operator=(const Parser&) = delete;
    
    // Read the next line
    CSVLineWrapper nextLine() {
        return CSVLineWrapper(csv_parser_next_line(parser_));
    }
    
    // Get header index by name
    int getHeaderIndex(const std::string& header_name) const {
        return csv_parser_get_header_index(parser_, header_name.c_str());
    }
    
    // Get field by header name
    std::string getField(const CSVLineWrapper& line, const std::string& header_name) const {
        if (!line) return "";
        
        const char* field = csv_line_get_field_by_name(line.get(), parser_, header_name.c_str());
        return field ? field : "";
    }
    
    // Read all lines at once (for small files)
    std::vector<CSVLineWrapper> readAllLines(size_t max_lines = 0) {
        std::vector<CSVLineWrapper> lines;
        
        while (true) {
            CSVLineWrapper line(csv_parser_next_line(parser_));
            if (!line) break;
            
            lines.push_back(std::move(line));
            
            if (max_lines > 0 && lines.size() >= max_lines) {
                break;
            }
        }
        
        return lines;
    }
    
    // Get the original header line
    std::string getHeaderLine() const {
        std::string result;
        
        // We can't directly access the header count from CSVParser struct since it's an incomplete type
        // Use getHeaderCount() which calls an appropriate C function
        const size_t header_count = getHeaderCount();
        
        for (size_t i = 0; i < header_count; i++) {
            if (i > 0) result += ",";
            result += getHeaderName(i);
        }
        return result;
    }
    
    // Get header count - use accessor function instead of direct struct access
    size_t getHeaderCount() const {
        // We'll need to add a C function to expose this
        if (!parser_) return 0;
        // Call the C function instead of direct access
        return csv_parser_get_header_count(parser_);
    }
    
    // Get header name by index - use accessor function instead of direct struct access
    std::string getHeaderName(size_t index) const {
        if (!parser_) return "";
        // Call the C function instead of direct array access
        const char* header = csv_parser_get_header_name(parser_, index);
        return header ? header : "";
    }
    
    // Read a batch of lines for parallel processing
    std::vector<CSVLineWrapper> nextBatch(size_t batch_size) {
        std::vector<CSVLineWrapper> batch;
        batch.reserve(batch_size);
        
        std::vector<CSVLine*> lines(batch_size, nullptr);
        size_t count = csv_parser_next_batch(parser_, batch_size, lines.data());
        
        for (size_t i = 0; i < count; i++) {
            batch.emplace_back(CSVLineWrapper(lines[i]));
        }
        
        return batch;
    }
};

} // namespace csv

#endif // CSV_WRAPPER_HPP 