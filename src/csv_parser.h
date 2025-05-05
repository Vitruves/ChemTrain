/**
 * High-performance CSV parser header
 * Optimized C implementation for maximum speed
 */

#ifndef CSV_PARSER_H
#define CSV_PARSER_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Error codes for CSV parsing
typedef enum {
    CSV_ERROR_NONE = 0,
    CSV_ERROR_FILE_NOT_FOUND,
    CSV_ERROR_MEMORY,
    CSV_ERROR_PARSE,
    CSV_ERROR_TOO_MANY_FIELDS
} CSVErrorCode;

// Error handler function type
typedef void (*CSVErrorHandler)(CSVErrorCode error, size_t line_number, 
                              const char* message, void* context);

// Forward declarations
typedef struct CSVParser CSVParser;

// CSV Line structure
typedef struct {
    char** fields;
    size_t field_count;
    size_t line_number;
} CSVLine;

/**
 * Initialize a CSV parser
 * 
 * @param filename Path to CSV file
 * @param delimiter Field delimiter character (typically ',')
 * @param has_header Whether the file has a header row
 * @param error_handler Function to call on error (can be NULL)
 * @param error_context User-defined context for error handler
 * @return A new CSVParser or NULL on error
 */
CSVParser* csv_parser_init(const char* filename, char delimiter, bool has_header,
                         CSVErrorHandler error_handler, void* error_context);

/**
 * Free CSV parser resources
 * 
 * @param parser The parser to free
 */
void csv_parser_free(CSVParser* parser);

/**
 * Read the next line from the CSV file
 * 
 * @param parser The CSV parser
 * @return A new CSVLine or NULL on EOF or error
 */
CSVLine* csv_parser_next_line(CSVParser* parser);

/**
 * Get header index by name (case-insensitive)
 * 
 * @param parser The CSV parser
 * @param header_name Name of the header to find
 * @return Index of header or -1 if not found
 */
int csv_parser_get_header_index(CSVParser* parser, const char* header_name);

/**
 * Free CSV line resources
 * 
 * @param line The CSV line to free
 */
void csv_line_free(CSVLine* line);

/**
 * Get field from CSV line by index
 * 
 * @param line The CSV line
 * @param index Index of the field to get
 * @return Field value or NULL if index is out of range
 */
const char* csv_line_get_field(CSVLine* line, size_t index);

/**
 * Get field from CSV line by header name
 * 
 * @param line The CSV line
 * @param parser The CSV parser (for header lookup)
 * @param header_name Name of the header
 * @return Field value or NULL if header is not found
 */
const char* csv_line_get_field_by_name(CSVLine* line, CSVParser* parser, const char* header_name);

/**
 * Get the number of headers in the CSV file
 * 
 * @param parser The CSV parser
 * @return Number of headers or 0 if no headers or parser is NULL
 */
size_t csv_parser_get_header_count(CSVParser* parser);

/**
 * Get header name by index
 * 
 * @param parser The CSV parser
 * @param index Index of the header to get
 * @return Header name or NULL if index is out of range or parser is NULL
 */
const char* csv_parser_get_header_name(CSVParser* parser, size_t index);

/**
 * Get a batch of CSV lines for multiprocessing
 * 
 * @param parser The CSV parser
 * @param batch_size Maximum number of lines to read
 * @param lines Pointer to array of CSVLine pointers to fill
 * @return Number of lines read, 0 on EOF or error
 */
size_t csv_parser_next_batch(CSVParser* parser, size_t batch_size, CSVLine** lines);

/**
 * Free a batch of CSV lines
 * 
 * @param lines Array of CSVLine pointers
 * @param count Number of lines in the array
 */
void csv_lines_free_batch(CSVLine** lines, size_t count);

#ifdef __cplusplus
}
#endif

#endif /* CSV_PARSER_H */ 