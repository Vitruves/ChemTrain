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
    CSVLine* line;
    bool owns_line;
    
public:
    CSVLineWrapper(CSVLine* l, bool take_ownership = true) 
        : line(l), owns_line(take_ownership) {}
    
    ~CSVLineWrapper() {
        if (owns_line && line) {
            csv_line_free(line);
        }
    }
    
    // Prevent copying
    CSVLineWrapper(const CSVLineWrapper&) = delete;
    CSVLineWrapper& operator=(const CSVLineWrapper&) = delete;
    
    // Allow moving
    CSVLineWrapper(CSVLineWrapper&& other) noexcept 
        : line(other.line), owns_line(other.owns_line) {
        other.line = nullptr;
        other.owns_line = false;
    }
    
    CSVLineWrapper& operator=(CSVLineWrapper&& other) noexcept {
        if (this != &other) {
            if (owns_line && line) {
                csv_line_free(line);
            }
            line = other.line;
            owns_line = other.owns_line;
            other.line = nullptr;
            other.owns_line = false;
        }
        return *this;
    }
    
    // Access the underlying CSVLine
    CSVLine* get() const { return line; }
    
    // Check if line is valid
    explicit operator bool() const { return line != nullptr; }
    
    // Get field by index
    const char* getField(size_t index) const {
        return line ? csv_line_get_field(line, index) : nullptr;
    }
    
    // Get number of fields
    size_t getFieldCount() const {
        return line ? line->field_count : 0;
    }
    
    // Get line number
    size_t getLineNumber() const {
        return line ? line->line_number : 0;
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
    CSVParser* parser;
    ErrorCallback error_handler;
    std::function<void(size_t, double)> progress_handler;
    
    static void c_error_handler(CSVErrorCode error, size_t line, const char* msg, void* context) {
        auto* self = static_cast<Parser*>(context);
        if (self && self->error_handler) {
            self->error_handler(error, line, msg ? msg : "");
        }
    }

    static void c_progress_handler(size_t line, double progress, void* context) {
        auto* self = static_cast<Parser*>(context);
        if (self && self->progress_handler) {
            self->progress_handler(line, progress);
        }
    }

public:
    Parser(const std::string& filename, char delimiter = ',', bool has_header = true,
           std::function<void(CSVErrorCode, size_t, const std::string&)> eh = nullptr,
           std::function<void(size_t, double)> ph = nullptr)
        : parser(nullptr), error_handler(eh), progress_handler(ph) {
        
        parser = csv_parser_init(filename.c_str(), delimiter, has_header,
                               c_error_handler, this,
                               c_progress_handler, this);
        
        if (!parser) {
            throw std::runtime_error("Failed to initialize CSV parser");
        }
    }
    
    ~Parser() {
        if (parser) {
            csv_parser_free(parser);
        }
    }
    
    // Prevent copying
    Parser(const Parser&) = delete;
    Parser& operator=(const Parser&) = delete;
    
    // Allow moving
    Parser(Parser&& other) noexcept : parser(other.parser) {
        other.parser = nullptr;
    }
    
    Parser& operator=(Parser&& other) noexcept {
        if (this != &other) {
            if (parser) {
                csv_parser_free(parser);
            }
            parser = other.parser;
            other.parser = nullptr;
        }
        return *this;
    }
    
    // Read the next line
    CSVLineWrapper nextLine() {
        if (!parser) return CSVLineWrapper(nullptr);
        return CSVLineWrapper(csv_parser_next_line(parser));
    }
    
    // Get header index by name
    int getHeaderIndex(const std::string& header_name) const {
        return parser ? csv_parser_get_header_index(parser, header_name.c_str()) : -1;
    }
    
    // Get field by header name
    std::string getField(const CSVLineWrapper& line, const std::string& header_name) const {
        if (!line) return "";
        
        const char* field = csv_line_get_field_by_name(line.get(), parser, header_name.c_str());
        return field ? field : "";
    }
    
    // Read all lines at once (for small files)
    std::vector<CSVLineWrapper> readAllLines(size_t max_lines = 0) {
        std::vector<CSVLineWrapper> lines;
        
        while (true) {
            CSVLineWrapper line(csv_parser_next_line(parser));
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
        return parser ? csv_parser_get_header_count(parser) : 0;
    }
    
    // Get header name by index - use accessor function instead of direct struct access
    const char* getHeaderName(size_t index) const {
        return parser ? csv_parser_get_header_name(parser, index) : nullptr;
    }
};

} // namespace csv

#endif // CSV_WRAPPER_HPP 