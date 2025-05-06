/**
 * High-performance CSV parser in C
 * Optimized for speed and memory efficiency
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <sys/stat.h>

#include "csv_parser.h"

#define CSV_PARSER_BUFFER_SIZE (1024 * 1024)  // 1MB buffer for best performance
#define CSV_MAX_FIELDS 1024                  // Maximum fields per line

// CSV parser state
struct CSVParser {
    FILE* file;
    char* buffer;
    size_t buffer_size;
    size_t buffer_pos;
    size_t buffer_end;
    size_t line_number;
    size_t total_bytes;
    size_t bytes_read;
    char delimiter;
    bool has_header;
    char** headers;
    size_t header_count;
    CSVErrorHandler error_handler;
    void* error_context;
    CSVProgressCallback progress_callback;
    void* progress_context;
};

// Initialize CSV parser
CSVParser* csv_parser_init(const char* filename, char delimiter, bool has_header,
                          CSVErrorHandler error_handler, void* error_context,
                          CSVProgressCallback progress_callback, void* progress_context) {
    CSVParser* parser = (CSVParser*)malloc(sizeof(CSVParser));
    if (!parser) {
        if (error_handler) {
            error_handler(CSV_ERROR_MEMORY, 0, "Failed to allocate parser", error_context);
        }
        return NULL;
    }
    
    parser->file = fopen(filename, "rb");
    if (!parser->file) {
        free(parser);
        if (error_handler) {
            error_handler(CSV_ERROR_FILE_NOT_FOUND, 0, "Failed to open file", error_context);
        }
        return NULL;
    }
    
    // Get file size for progress tracking
    struct stat st;
    if (fstat(fileno(parser->file), &st) == 0) {
        parser->total_bytes = st.st_size;
    } else {
        parser->total_bytes = 0;
    }
    parser->bytes_read = 0;
    
    parser->buffer = (char*)malloc(CSV_PARSER_BUFFER_SIZE);
    if (!parser->buffer) {
        fclose(parser->file);
        free(parser);
        if (error_handler) {
            error_handler(CSV_ERROR_MEMORY, 0, "Failed to allocate buffer", error_context);
        }
        return NULL;
    }
    
    parser->buffer_size = CSV_PARSER_BUFFER_SIZE;
    parser->buffer_pos = 0;
    parser->buffer_end = 0;
    parser->line_number = 0;
    parser->delimiter = delimiter;
    parser->has_header = has_header;
    parser->headers = NULL;
    parser->header_count = 0;
    parser->error_handler = error_handler;
    parser->error_context = error_context;
    parser->progress_callback = progress_callback;
    parser->progress_context = progress_context;
    
    // Read header if needed
    if (has_header) {
        CSVLine* header_line = csv_parser_next_line(parser);
        if (header_line) {
            parser->header_count = header_line->field_count;
            parser->headers = (char**)malloc(sizeof(char*) * header_line->field_count);
            if (parser->headers) {
                for (size_t i = 0; i < header_line->field_count; i++) {
                    parser->headers[i] = strdup(header_line->fields[i]);
                }
            }
            csv_line_free(header_line);
        }
    }
    
    return parser;
}

// Fill the buffer with more data from the file
static bool fill_buffer(CSVParser* parser) {
    if (parser->buffer_pos >= parser->buffer_end) {
        // Buffer is empty, refill it
        parser->buffer_pos = 0;
        parser->buffer_end = fread(parser->buffer, 1, parser->buffer_size, parser->file);
        parser->bytes_read += parser->buffer_end;
        
        // Update progress
        if (parser->progress_callback && parser->total_bytes > 0) {
            double progress = (double)parser->bytes_read / parser->total_bytes * 100.0;
            parser->progress_callback(parser->line_number, progress, parser->progress_context);
        }
        
        return parser->buffer_end > 0;
    } else if (parser->buffer_pos > 0) {
        // Move remaining data to start of buffer
        size_t remaining = parser->buffer_end - parser->buffer_pos;
        memmove(parser->buffer, parser->buffer + parser->buffer_pos, remaining);
        parser->buffer_pos = 0;
        parser->buffer_end = remaining;
        
        // Read more data
        size_t bytes_read = fread(parser->buffer + remaining, 1, 
                                parser->buffer_size - remaining, parser->file);
        parser->buffer_end += bytes_read;
        parser->bytes_read += bytes_read;
        
        // Update progress
        if (parser->progress_callback && parser->total_bytes > 0) {
            double progress = (double)parser->bytes_read / parser->total_bytes * 100.0;
            parser->progress_callback(parser->line_number, progress, parser->progress_context);
        }
        
        return bytes_read > 0 || remaining > 0;
    }
    
    return true;
}

// Parse a CSV line
static CSVLine* parse_line(CSVParser* parser) {
    if (!parser || parser->buffer_pos >= parser->buffer_end) {
        return NULL;
    }
    
    CSVLine* line = (CSVLine*)malloc(sizeof(CSVLine));
    if (!line) {
        if (parser->error_handler) {
            parser->error_handler(CSV_ERROR_MEMORY, parser->line_number,
                                "Failed to allocate line", parser->error_context);
        }
        return NULL;
    }
    
    line->fields = (char**)malloc(sizeof(char*) * CSV_MAX_FIELDS);
    if (!line->fields) {
        free(line);
        if (parser->error_handler) {
            parser->error_handler(CSV_ERROR_MEMORY, parser->line_number,
                                "Failed to allocate fields", parser->error_context);
        }
        return NULL;
    }
    
    line->field_count = 0;
    line->line_number = parser->line_number + 1;
    
    char* field_start = parser->buffer + parser->buffer_pos;
    bool in_quotes = false;
    size_t field_len = 0;
    
    while (true) {
        // Ensure buffer has data
        if (parser->buffer_pos >= parser->buffer_end) {
            if (!fill_buffer(parser)) {
                break;
            }
            field_start = parser->buffer + parser->buffer_pos;
        }
        
        char c = parser->buffer[parser->buffer_pos++];
        
        if (c == '"') {
            in_quotes = !in_quotes;
        } else if ((c == parser->delimiter || c == '\n' || c == '\r') && !in_quotes) {
            // End of field
            char* field = (char*)malloc(field_len + 1);
            if (!field) {
                if (parser->error_handler) {
                    parser->error_handler(CSV_ERROR_MEMORY, parser->line_number,
                                        "Failed to allocate field", parser->error_context);
                }
                csv_line_free(line);
                return NULL;
            }
            
            // Copy field data
            size_t pos = 0;
            bool skip_next = false;
            for (size_t i = 0; i < field_len; i++) {
                if (field_start[i] == '"') {
                    if (skip_next) {
                        field[pos++] = '"';
                        skip_next = false;
                    } else {
                        skip_next = true;
                    }
                } else {
                    field[pos++] = field_start[i];
                    skip_next = false;
                }
            }
            field[pos] = '\0';
            
            line->fields[line->field_count++] = field;
            
            if (line->field_count >= CSV_MAX_FIELDS) {
                if (parser->error_handler) {
                    parser->error_handler(CSV_ERROR_TOO_MANY_FIELDS, parser->line_number,
                                        "Too many fields in line", parser->error_context);
                }
                csv_line_free(line);
                return NULL;
            }
            
            if (c == '\n' || c == '\r') {
                // End of line
                if (c == '\r' && parser->buffer_pos < parser->buffer_end &&
                    parser->buffer[parser->buffer_pos] == '\n') {
                    parser->buffer_pos++; // Skip \n after \r
                }
                break;
            }
            
            field_start = parser->buffer + parser->buffer_pos;
            field_len = 0;
        } else {
            field_len++;
        }
    }
    
    parser->line_number++;
    return line;
}

// Read the next line from the CSV file
CSVLine* csv_parser_next_line(CSVParser* parser) {
    if (!parser || !parser->file) return NULL;
    
    if (parser->buffer_pos >= parser->buffer_end) {
        if (!fill_buffer(parser)) {
            return NULL; // EOF
        }
    }
    
    return parse_line(parser);
}

// Free CSV parser resources
void csv_parser_free(CSVParser* parser) {
    if (!parser) return;
    
    if (parser->headers) {
        for (size_t i = 0; i < parser->header_count; i++) {
            free(parser->headers[i]);
        }
        free(parser->headers);
    }
    
    if (parser->buffer) {
        free(parser->buffer);
    }
    
    if (parser->file) {
        fclose(parser->file);
    }
    
    free(parser);
}

// Free CSV line resources
void csv_line_free(CSVLine* line) {
    if (!line) return;
    
    if (line->fields) {
        for (size_t i = 0; i < line->field_count; i++) {
            free(line->fields[i]);
        }
        free(line->fields);
    }
    
    free(line);
}

// Get field from CSV line by index
const char* csv_line_get_field(CSVLine* line, size_t index) {
    if (!line || index >= line->field_count) return NULL;
    return line->fields[index];
}

// Get header index by name (case-insensitive)
int csv_parser_get_header_index(CSVParser* parser, const char* header_name) {
    if (!parser || !parser->headers || !header_name) return -1;
    
    for (size_t i = 0; i < parser->header_count; i++) {
        if (strcasecmp(parser->headers[i], header_name) == 0) {
            return (int)i;
        }
    }
    
    return -1;
}

// Get field from CSV line by header name
const char* csv_line_get_field_by_name(CSVLine* line, CSVParser* parser, const char* header_name) {
    if (!line || !parser || !header_name) return NULL;
    
    int index = csv_parser_get_header_index(parser, header_name);
    if (index < 0) return NULL;
    
    return csv_line_get_field(line, (size_t)index);
}

// Get the number of headers in the CSV file
size_t csv_parser_get_header_count(CSVParser* parser) {
    return parser ? parser->header_count : 0;
}

// Get header name by index
const char* csv_parser_get_header_name(CSVParser* parser, size_t index) {
    if (!parser || index >= parser->header_count) return NULL;
    return parser->headers[index];
}

long csv_parser_ftell(CSVParser* parser) {
    if (!parser || !parser->file) return -1;
    return ftell(parser->file);
} 