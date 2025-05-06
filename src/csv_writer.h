/**
 * High-performance CSV writer header
 * Optimized C implementation for maximum speed
 */

#ifndef CSV_WRITER_H
#define CSV_WRITER_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Error codes for CSV writing
typedef enum {
    CSV_WRITE_ERROR_NONE = 0,
    CSV_WRITE_ERROR_FILE_NOT_FOUND,
    CSV_WRITE_ERROR_MEMORY,
    CSV_WRITE_ERROR_IO
} CSVWriteErrorCode;

// Error handler function type
typedef void (*CSVWriteErrorHandler)(CSVWriteErrorCode error, 
                                   const char* message, void* context);

// Forward declarations
typedef struct CSVWriter CSVWriter;

/**
 * Initialize a CSV writer
 * 
 * @param filename Path to output CSV file
 * @param delimiter Field delimiter character (typically ',')
 * @param error_handler Function to call on error (can be NULL)
 * @param error_context User-defined context for error handler
 * @return A new CSVWriter or NULL on error
 */
CSVWriter* csv_writer_init(const char* filename, char delimiter,
                          CSVWriteErrorHandler error_handler, void* error_context);

/**
 * Free CSV writer resources
 * 
 * @param writer The writer to free
 */
void csv_writer_free(CSVWriter* writer);

/**
 * Write a row of fields to the CSV file
 * 
 * @param writer The CSV writer
 * @param fields Array of field strings
 * @param field_count Number of fields
 * @return true on success, false on error
 */
bool csv_writer_write_row(CSVWriter* writer, const char** fields, size_t field_count);

/**
 * Write header row to the CSV file
 * 
 * @param writer The CSV writer
 * @param headers Array of header strings
 * @param header_count Number of headers
 * @return true on success, false on error
 */
bool csv_writer_write_headers(CSVWriter* writer, const char** headers, size_t header_count);

/**
 * Flush any buffered data to disk
 * 
 * @param writer The CSV writer
 * @return true on success, false on error
 */
bool csv_writer_flush(CSVWriter* writer);

/**
 * Get the current number of rows written
 * 
 * @param writer The CSV writer
 * @return Number of rows written so far
 */
size_t csv_writer_get_row_count(CSVWriter* writer);

#ifdef __cplusplus
}
#endif

#endif /* CSV_WRITER_H */ 