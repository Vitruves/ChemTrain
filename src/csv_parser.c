/**
 * High-performance CSV parser in C
 * Optimized for speed and memory efficiency
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

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
    char delimiter;
    bool has_header;
    char** headers;
    size_t header_count;
    CSVErrorHandler error_handler;
    void* error_context;
};

// Initialize CSV parser
CSVParser* csv_parser_init(const char* filename, char delimiter, bool has_header, 
                          CSVErrorHandler error_handler, void* error_context) {
    CSVParser* parser = (CSVParser*)malloc(sizeof(CSVParser));
    if (!parser) {
        return NULL;
    }
    
    parser->file = fopen(filename, "rb");
    if (!parser->file) {
        free(parser);
        return NULL;
    }
    
    // Allocate buffer
    parser->buffer_size = CSV_PARSER_BUFFER_SIZE;
    parser->buffer = (char*)malloc(parser->buffer_size);
    if (!parser->buffer) {
        fclose(parser->file);
        free(parser);
        return NULL;
    }
    
    parser->buffer_pos = 0;
    parser->buffer_end = 0;
    parser->line_number = 0;
    parser->delimiter = delimiter;
    parser->has_header = has_header;
    parser->headers = NULL;
    parser->header_count = 0;
    parser->error_handler = error_handler;
    parser->error_context = error_context;
    
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

// Free CSV parser resources
void csv_parser_free(CSVParser* parser) {
    if (!parser) return;
    
    if (parser->file) {
        fclose(parser->file);
    }
    
    if (parser->buffer) {
        free(parser->buffer);
    }
    
    if (parser->headers) {
        for (size_t i = 0; i < parser->header_count; i++) {
            free(parser->headers[i]);
        }
        free(parser->headers);
    }
    
    free(parser);
}

// Fill the buffer with more data from the file
static bool csv_parser_fill_buffer(CSVParser* parser) {
    if (parser->buffer_pos >= parser->buffer_end) {
        // Buffer is empty, refill it
        parser->buffer_pos = 0;
        parser->buffer_end = fread(parser->buffer, 1, parser->buffer_size, parser->file);
        return parser->buffer_end > 0;
    } else if (parser->buffer_pos > 0) {
        // Move remaining data to the start of the buffer
        size_t remaining = parser->buffer_end - parser->buffer_pos;
        memmove(parser->buffer, parser->buffer + parser->buffer_pos, remaining);
        parser->buffer_pos = 0;
        parser->buffer_end = remaining;
        
        // Read more data
        size_t bytes_read = fread(parser->buffer + remaining, 1, 
                                 parser->buffer_size - remaining, parser->file);
        parser->buffer_end += bytes_read;
        return bytes_read > 0 || remaining > 0;
    }
    
    return true;
}

// Parse a CSV line
static CSVLine* csv_parser_parse_line(CSVParser* parser, const char* line, size_t length) {
    CSVLine* csv_line = (CSVLine*)malloc(sizeof(CSVLine));
    if (!csv_line) return NULL;
    
    csv_line->fields = (char**)malloc(sizeof(char*) * CSV_MAX_FIELDS);
    if (!csv_line->fields) {
        free(csv_line);
        return NULL;
    }
    
    csv_line->field_count = 0;
    
    // Parse the line into fields
    bool in_quotes = false;
    size_t field_start = 0;
    size_t i;
    
    for (i = 0; i < length; i++) {
        char c = line[i];
        
        if (c == '"') {
            in_quotes = !in_quotes;
        } else if ((c == parser->delimiter || c == '\n' || c == '\r') && !in_quotes) {
            // End of field
            size_t field_len = i - field_start;
            char* field = (char*)malloc(field_len + 1);
            
            if (field) {
                // Copy field data
                size_t dst_pos = 0;
                for (size_t src_pos = field_start; src_pos < i; src_pos++) {
                    if (line[src_pos] == '"' && src_pos + 1 < i && line[src_pos + 1] == '"') {
                        // Handle escaped quotes (double quote)
                        field[dst_pos++] = '"';
                        src_pos++; // Skip the second quote
                    } else if (line[src_pos] != '"') {
                        field[dst_pos++] = line[src_pos];
                    }
                }
                field[dst_pos] = '\0';
                
                csv_line->fields[csv_line->field_count++] = field;
                
                if (csv_line->field_count >= CSV_MAX_FIELDS) {
                    break; // Too many fields
                }
            }
            
            field_start = i + 1;
        }
    }
    
    // Handle the last field
    if (field_start < length) {
        size_t field_len = length - field_start;
        char* field = (char*)malloc(field_len + 1);
        
        if (field) {
            // Copy field data
            size_t dst_pos = 0;
            for (size_t src_pos = field_start; src_pos < length; src_pos++) {
                if (line[src_pos] == '"' && src_pos + 1 < length && line[src_pos + 1] == '"') {
                    // Handle escaped quotes (double quote)
                    field[dst_pos++] = '"';
                    src_pos++; // Skip the second quote
                } else if (line[src_pos] != '"') {
                    field[dst_pos++] = line[src_pos];
                }
            }
            field[dst_pos] = '\0';
            
            csv_line->fields[csv_line->field_count++] = field;
        }
    }
    
    return csv_line;
}

// Read the next line from the CSV file
CSVLine* csv_parser_next_line(CSVParser* parser) {
    if (!parser || !parser->file) return NULL;
    
    // Make sure buffer has data
    if (parser->buffer_pos >= parser->buffer_end) {
        if (!csv_parser_fill_buffer(parser)) {
            return NULL; // EOF
        }
    }
    
    // Find the end of the line
    char* line_start = parser->buffer + parser->buffer_pos;
    char* line_end = NULL;
    bool in_quotes = false;
    size_t i = parser->buffer_pos;
    
    while (i < parser->buffer_end) {
        char c = parser->buffer[i];
        
        if (c == '"') {
            in_quotes = !in_quotes;
        } else if ((c == '\n' || c == '\r') && !in_quotes) {
            line_end = parser->buffer + i;
            break;
        }
        
        i++;
    }
    
    if (!line_end) {
        // Line extends beyond current buffer, need more data
        if (i < parser->buffer_size) {
            // End of file reached
            if (i > parser->buffer_pos) {
                line_end = parser->buffer + i;
            } else {
                return NULL; // Empty line at EOF
            }
        } else {
            // Buffer is full but line is not complete
            // Increase buffer size
            size_t new_size = parser->buffer_size * 2;
            char* new_buffer = (char*)realloc(parser->buffer, new_size);
            
            if (!new_buffer) {
                if (parser->error_handler) {
                    parser->error_handler(CSV_ERROR_MEMORY, parser->line_number, 
                                         "Out of memory", parser->error_context);
                }
                return NULL;
            }
            
            parser->buffer = new_buffer;
            parser->buffer_size = new_size;
            
            // Read more data
            if (!csv_parser_fill_buffer(parser)) {
                // End of file reached, use what we have
                line_end = parser->buffer + parser->buffer_end;
            }
            
            // Retry with the larger buffer
            return csv_parser_next_line(parser);
        }
    }
    
    // Process line
    size_t line_length = line_end - line_start;
    CSVLine* csv_line = csv_parser_parse_line(parser, line_start, line_length);
    
    if (csv_line) {
        csv_line->line_number = ++parser->line_number;
    }
    
    // Skip any CR/LF sequence
    size_t eol_length = 0;
    if (line_end < parser->buffer + parser->buffer_end) {
        if (*line_end == '\r' && line_end + 1 < parser->buffer + parser->buffer_end && *(line_end + 1) == '\n') {
            eol_length = 2;
        } else {
            eol_length = 1;
        }
    }
    
    parser->buffer_pos = (line_end - parser->buffer) + eol_length;
    
    return csv_line;
}

// Get header index by name
int csv_parser_get_header_index(CSVParser* parser, const char* header_name) {
    if (!parser || !parser->headers || !header_name) return -1;
    
    for (size_t i = 0; i < parser->header_count; i++) {
        if (strcasecmp(parser->headers[i], header_name) == 0) {
            return (int)i;
        }
    }
    
    return -1; // Header not found
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

// Get field from CSV line
const char* csv_line_get_field(CSVLine* line, size_t index) {
    if (!line || index >= line->field_count) return NULL;
    return line->fields[index];
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
    if (!parser) return 0;
    return parser->header_count;
}

// Get header name by index
const char* csv_parser_get_header_name(CSVParser* parser, size_t index) {
    if (!parser || index >= parser->header_count) return NULL;
    return parser->headers[index];
}

// Get a batch of CSV lines for multiprocessing
size_t csv_parser_next_batch(CSVParser* parser, size_t batch_size, CSVLine** lines) {
    if (!parser || !lines || batch_size == 0) return 0;
    
    size_t count = 0;
    for (size_t i = 0; i < batch_size; i++) {
        CSVLine* line = csv_parser_next_line(parser);
        if (!line) break;
        
        lines[count++] = line;
    }
    
    return count;
}

// Free a batch of CSV lines
void csv_lines_free_batch(CSVLine** lines, size_t count) {
    if (!lines) return;
    
    for (size_t i = 0; i < count; i++) {
        if (lines[i]) {
            csv_line_free(lines[i]);
        }
    }
} 