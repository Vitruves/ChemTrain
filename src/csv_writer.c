#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "csv_writer.h"

#define CSV_WRITER_BUFFER_SIZE (1024 * 1024)  // 1MB buffer for best performance

// CSV writer state
struct CSVWriter {
    FILE* file;
    char* buffer;
    size_t buffer_size;
    size_t buffer_pos;
    char delimiter;
    size_t row_count;
    CSVWriteErrorHandler error_handler;
    void* error_context;
};

// Helper function to report errors
static void report_error(CSVWriter* writer, CSVWriteErrorCode code, const char* message) {
    if (writer && writer->error_handler) {
        writer->error_handler(code, message, writer->error_context);
    }
}

// Helper function to flush buffer to file
static bool flush_buffer(CSVWriter* writer) {
    if (writer->buffer_pos > 0) {
        size_t written = fwrite(writer->buffer, 1, writer->buffer_pos, writer->file);
        if (written != writer->buffer_pos) {
            report_error(writer, CSV_WRITE_ERROR_IO, "Failed to write to file");
            return false;
        }
        writer->buffer_pos = 0;
    }
    return true;
}

// Helper function to write to buffer
static bool write_to_buffer(CSVWriter* writer, const char* data, size_t length) {
    // If data is larger than remaining buffer space, flush first
    if (writer->buffer_pos + length > writer->buffer_size) {
        if (!flush_buffer(writer)) {
            return false;
        }
        
        // If data is larger than buffer size, write directly
        if (length > writer->buffer_size) {
            size_t written = fwrite(data, 1, length, writer->file);
            if (written != length) {
                report_error(writer, CSV_WRITE_ERROR_IO, "Failed to write to file");
                return false;
            }
            return true;
        }
    }
    
    // Copy data to buffer
    memcpy(writer->buffer + writer->buffer_pos, data, length);
    writer->buffer_pos += length;
    return true;
}

// Helper function to write a field with proper escaping
static bool write_field(CSVWriter* writer, const char* field) {
    if (!field) field = "";
    
    bool needs_quotes = false;
    const char* p;
    size_t quote_count = 0;
    
    // Check if field needs quoting
    for (p = field; *p; p++) {
        if (*p == '"') {
            needs_quotes = true;
            quote_count++;
        } else if (*p == writer->delimiter || *p == '\n' || *p == '\r') {
            needs_quotes = true;
        }
    }
    
    size_t field_len = p - field;
    
    // Write opening quote if needed
    if (needs_quotes && !write_to_buffer(writer, "\"", 1)) {
        return false;
    }
    
    // Write field content with quote escaping
    if (quote_count > 0) {
        const char* start = field;
        for (p = field; *p; p++) {
            if (*p == '"') {
                // Write up to quote
                size_t chunk_len = p - start;
                if (chunk_len > 0 && !write_to_buffer(writer, start, chunk_len)) {
                    return false;
                }
                // Write escaped quote
                if (!write_to_buffer(writer, "\"\"", 2)) {
                    return false;
                }
                start = p + 1;
            }
        }
        // Write remaining chunk
        size_t chunk_len = p - start;
        if (chunk_len > 0 && !write_to_buffer(writer, start, chunk_len)) {
            return false;
        }
    } else {
        // Write field as is
        if (!write_to_buffer(writer, field, field_len)) {
            return false;
        }
    }
    
    // Write closing quote if needed
    if (needs_quotes && !write_to_buffer(writer, "\"", 1)) {
        return false;
    }
    
    return true;
}

CSVWriter* csv_writer_init(const char* filename, char delimiter,
                          CSVWriteErrorHandler error_handler, void* error_context) {
    CSVWriter* writer = (CSVWriter*)malloc(sizeof(CSVWriter));
    if (!writer) {
        if (error_handler) {
            error_handler(CSV_WRITE_ERROR_MEMORY, "Failed to allocate writer", error_context);
        }
        return NULL;
    }
    
    writer->file = fopen(filename, "wb");
    if (!writer->file) {
        free(writer);
        if (error_handler) {
            error_handler(CSV_WRITE_ERROR_FILE_NOT_FOUND, "Failed to open file", error_context);
        }
        return NULL;
    }
    
    writer->buffer = (char*)malloc(CSV_WRITER_BUFFER_SIZE);
    if (!writer->buffer) {
        fclose(writer->file);
        free(writer);
        if (error_handler) {
            error_handler(CSV_WRITE_ERROR_MEMORY, "Failed to allocate buffer", error_context);
        }
        return NULL;
    }
    
    writer->buffer_size = CSV_WRITER_BUFFER_SIZE;
    writer->buffer_pos = 0;
    writer->delimiter = delimiter;
    writer->row_count = 0;
    writer->error_handler = error_handler;
    writer->error_context = error_context;
    
    return writer;
}

void csv_writer_free(CSVWriter* writer) {
    if (!writer) return;
    
    if (writer->buffer) {
        flush_buffer(writer);
        free(writer->buffer);
    }
    
    if (writer->file) {
        fclose(writer->file);
    }
    
    free(writer);
}

bool csv_writer_write_row(CSVWriter* writer, const char** fields, size_t field_count) {
    if (!writer || !fields) return false;
    
    for (size_t i = 0; i < field_count; i++) {
        // Write delimiter before all fields except the first
        if (i > 0 && !write_to_buffer(writer, &writer->delimiter, 1)) {
            return false;
        }
        
        if (!write_field(writer, fields[i])) {
            return false;
        }
    }
    
    // Write newline
    if (!write_to_buffer(writer, "\n", 1)) {
        return false;
    }
    
    writer->row_count++;
    return true;
}

bool csv_writer_write_headers(CSVWriter* writer, const char** headers, size_t header_count) {
    return csv_writer_write_row(writer, headers, header_count);
}

bool csv_writer_flush(CSVWriter* writer) {
    if (!writer) return false;
    return flush_buffer(writer);
}

size_t csv_writer_get_row_count(CSVWriter* writer) {
    return writer ? writer->row_count : 0;
} 