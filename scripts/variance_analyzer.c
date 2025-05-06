#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <ctype.h>
#include <omp.h>

#define MAX_LINE_LENGTH 65536  // Increased for larger lines
#define MAX_COLUMNS 1024
#define MAX_ROWS 10000000      // Increased max rows
#define PROGRESS_INTERVAL 0.01 // Progress reporting every 0.01%

// ANSI color codes for terminal output
#define RED_COLOR "\033[1;31m"
#define GREEN_COLOR "\033[1;32m"
#define YELLOW_COLOR "\033[1;33m"
#define BLUE_COLOR "\033[1;34m"
#define MAGENTA_COLOR "\033[1;35m"
#define RESET_COLOR "\033[0m"

typedef struct {
    double* data;
    double sum;
    double sum_squares;
    int count;
    char* name;
    double variance;
    double std_dev;
    bool has_variance;
    double correlation;    // Correlation with reference column
} Column;

// Function to trim whitespace from a string
void trim(char* str) {
    if (!str) return;
    
    char* end;
    // Trim leading space
    while(isspace((unsigned char)*str)) str++;
    if(*str == 0) return;  // All spaces
    
    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;
    
    // Write new null terminator
    *(end + 1) = 0;
}

// Function to safely parse a CSV line, handling quoted fields
int parse_csv_line(char* line, char** fields, int max_fields) {
    if (!line || !fields || max_fields <= 0) return 0;
    
    int field_count = 0;
    char* p = line;
    char* field_start = p;
    bool in_quotes = false;
    
    while (*p && field_count < max_fields) {
        if (*p == '"') {
            in_quotes = !in_quotes;
        } else if (*p == ',' && !in_quotes) {
            *p = '\0';  // Terminate the field
            fields[field_count] = field_start;
            trim(fields[field_count]);
            field_count++;
            field_start = p + 1;
        }
        p++;
    }
    
    // Handle the last field
    if (field_start && field_count < max_fields) {
        fields[field_count] = field_start;
        trim(fields[field_count]);
        field_count++;
    }
    
    return field_count;
}

// Check if a string can be parsed as a number
bool is_numeric(const char* str) {
    if (!str || *str == '\0') return false;
    
    char* endptr;
    strtod(str, &endptr);
    
    // Check if conversion was successful and consumed the entire string
    return *endptr == '\0';
}

// Comparison function for sorting columns by correlation (absolute value)
int compare_columns_by_correlation(const void* a, const void* b) {
    const Column* col_a = (const Column*)a;
    const Column* col_b = (const Column*)b;
    
    // Sort by absolute correlation value (descending)
    double abs_corr_a = fabs(col_a->correlation);
    double abs_corr_b = fabs(col_b->correlation);
    
    if (abs_corr_a > abs_corr_b) return -1;
    if (abs_corr_a < abs_corr_b) return 1;
    
    // If correlations are equal, sort by name
    return strcmp(col_a->name, col_b->name);
}

// Comparison function for sorting columns by variance
int compare_columns_by_variance(const void* a, const void* b) {
    const Column* col_a = (const Column*)a;
    const Column* col_b = (const Column*)b;
    
    // First, sort by whether they have variance or not (non-zero variance first)
    if (col_a->has_variance && !col_b->has_variance) return -1;
    if (!col_a->has_variance && col_b->has_variance) return 1;
    
    // If both have variance or both don't, sort by name
    return strcmp(col_a->name, col_b->name);
}

// Function to count total lines in file (for progress reporting)
long count_lines(FILE* file) {
    char buffer[8192];
    long count = 0;
    
    // Save current position
    fpos_t position;
    fgetpos(file, &position);
    
    // Reset to beginning of file
    rewind(file);
    
    while (fgets(buffer, sizeof(buffer), file) != NULL) {
        count++;
    }
    
    // Restore position
    fsetpos(file, &position);
    
    return count;
}

// Function to find a column index by name
int find_column_by_name(Column* columns, int num_columns, const char* name) {
    for (int i = 0; i < num_columns; i++) {
        if (strcmp(columns[i].name, name) == 0) {
            return i;
        }
    }
    return -1;
}

// Function to display usage information
void print_usage(const char* program_name) {
    printf("Usage: %s <csv_file> [--reference-col <column_name>]\n\n", program_name);
    printf("Options:\n");
    printf("  --reference-col <column_name>  Calculate correlation between this column and all others\n");
    printf("\nExample:\n");
    printf("  %s data.csv --reference-col Temperature\n", program_name);
}

int main(int argc, char* argv[]) {
    char* filename = NULL;
    char* reference_column = NULL;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--reference-col") == 0) {
            if (i + 1 < argc) {
                reference_column = argv[i + 1];
                i++; // Skip the next argument
            } else {
                fprintf(stderr, "Error: --reference-col requires a column name\n");
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        } else if (filename == NULL) {
            filename = argv[i];
        } else {
            fprintf(stderr, "Error: Unexpected argument '%s'\n", argv[i]);
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }
    
    if (filename == NULL) {
        fprintf(stderr, "Error: CSV file not specified\n");
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    
    // Set number of threads for OpenMP
    int num_threads = omp_get_max_threads();
    printf("%sUsing %d threads for parallel processing%s\n", 
           GREEN_COLOR, num_threads, RESET_COLOR);
    
    // Start timing
    double start_time = omp_get_wtime();
    
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }
    
    printf("%sAnalyzing file size...%s\n", YELLOW_COLOR, RESET_COLOR);
    long total_lines = count_lines(file);
    printf("%sFile contains %ld lines%s\n", YELLOW_COLOR, total_lines, RESET_COLOR);
    
    char* line = (char*)malloc(MAX_LINE_LENGTH);
    if (!line) {
        fprintf(stderr, "Memory allocation failed for line buffer\n");
        fclose(file);
        return EXIT_FAILURE;
    }
    
    char** header_tokens = (char**)malloc(MAX_COLUMNS * sizeof(char*));
    if (!header_tokens) {
        fprintf(stderr, "Memory allocation failed for header tokens\n");
        free(line);
        fclose(file);
        return EXIT_FAILURE;
    }
    
    Column* columns = (Column*)calloc(MAX_COLUMNS, sizeof(Column));
    if (!columns) {
        fprintf(stderr, "Memory allocation failed for columns\n");
        free(header_tokens);
        free(line);
        fclose(file);
        return EXIT_FAILURE;
    }
    
    // Initialize columns
    for (int i = 0; i < MAX_COLUMNS; i++) {
        columns[i].data = (double*)malloc(MAX_ROWS * sizeof(double));
        if (!columns[i].data) {
            fprintf(stderr, "Memory allocation failed for column data\n");
            // Clean up previously allocated memory
            for (int j = 0; j < i; j++) {
                free(columns[j].data);
            }
            free(columns);
            free(header_tokens);
            free(line);
            fclose(file);
            return EXIT_FAILURE;
        }
        columns[i].sum = 0.0;
        columns[i].sum_squares = 0.0;
        columns[i].count = 0;
        columns[i].name = NULL;
        columns[i].variance = 0.0;
        columns[i].std_dev = 0.0;
        columns[i].has_variance = false;
        columns[i].correlation = 0.0;
    }
    
    int num_columns = 0;
    
    // Read header
    if (fgets(line, MAX_LINE_LENGTH, file)) {
        // Remove newline character
        line[strcspn(line, "\n")] = 0;
        line[strcspn(line, "\r")] = 0;  // Also handle CR
        
        // Parse header line
        num_columns = parse_csv_line(line, header_tokens, MAX_COLUMNS);
        
        if (num_columns == 0) {
            fprintf(stderr, "Failed to parse header or no columns found\n");
            // Clean up
            for (int i = 0; i < MAX_COLUMNS; i++) {
                free(columns[i].data);
            }
            free(columns);
            free(header_tokens);
            free(line);
            fclose(file);
            return EXIT_FAILURE;
        }
        
        // Store column names
        for (int i = 0; i < num_columns; i++) {
            columns[i].name = strdup(header_tokens[i]);
            if (!columns[i].name) {
                fprintf(stderr, "Memory allocation failed for column name\n");
                // Clean up
                for (int j = 0; j < i; j++) {
                    free(columns[j].name);
                }
                for (int j = 0; j < MAX_COLUMNS; j++) {
                    free(columns[j].data);
                }
                free(columns);
                free(header_tokens);
                free(line);
                fclose(file);
                return EXIT_FAILURE;
            }
        }
    } else {
        fprintf(stderr, "Empty file or error reading header\n");
        // Clean up
        for (int i = 0; i < MAX_COLUMNS; i++) {
            free(columns[i].data);
        }
        free(columns);
        free(header_tokens);
        free(line);
        fclose(file);
        return EXIT_FAILURE;
    }
    
    printf("%sFound %d columns in header%s\n", GREEN_COLOR, num_columns, RESET_COLOR);
    
    // Check if reference column exists
    int ref_col_idx = -1;
    if (reference_column != NULL) {
        ref_col_idx = find_column_by_name(columns, num_columns, reference_column);
        if (ref_col_idx == -1) {
            fprintf(stderr, "Error: Reference column '%s' not found in CSV header\n", reference_column);
            print_usage(argv[0]);
            // Clean up
            for (int i = 0; i < num_columns; i++) {
                free(columns[i].name);
            }
            for (int i = 0; i < MAX_COLUMNS; i++) {
                free(columns[i].data);
            }
            free(columns);
            free(header_tokens);
            free(line);
            fclose(file);
            return EXIT_FAILURE;
        }
        printf("%sReference column: %s (index: %d)%s\n", 
               BLUE_COLOR, reference_column, ref_col_idx, RESET_COLOR);
    }
    
    // Calculate progress interval for 0.01% reporting
    long progress_step = (long)(total_lines * PROGRESS_INTERVAL);
    if (progress_step < 1) progress_step = 1;
    
    // Allocate tokens for data rows
    char*** thread_data_tokens = (char***)malloc(num_threads * sizeof(char**));
    for (int t = 0; t < num_threads; t++) {
        thread_data_tokens[t] = (char**)malloc(MAX_COLUMNS * sizeof(char*));
        if (!thread_data_tokens[t]) {
            fprintf(stderr, "Memory allocation failed for thread data tokens\n");
            // Clean up
            for (int j = 0; j < t; j++) {
                free(thread_data_tokens[j]);
            }
            free(thread_data_tokens);
            for (int i = 0; i < num_columns; i++) {
                free(columns[i].name);
            }
            for (int i = 0; i < MAX_COLUMNS; i++) {
                free(columns[i].data);
            }
            free(columns);
            free(header_tokens);
            free(line);
            fclose(file);
            return EXIT_FAILURE;
        }
    }
    
    // Create thread-local storage for sums
    double** thread_sums = (double**)malloc(num_threads * sizeof(double*));
    double** thread_sum_squares = (double**)malloc(num_threads * sizeof(double*));
    int** thread_counts = (int**)malloc(num_threads * sizeof(int*));
    
    // For correlation calculation
    double*** thread_products = NULL;
    if (ref_col_idx >= 0) {
        thread_products = (double***)malloc(num_threads * sizeof(double**));
        for (int t = 0; t < num_threads; t++) {
            thread_products[t] = (double**)malloc(num_columns * sizeof(double*));
            for (int i = 0; i < num_columns; i++) {
                thread_products[t][i] = (double*)calloc(1, sizeof(double));
            }
        }
    }
    
    for (int t = 0; t < num_threads; t++) {
        thread_sums[t] = (double*)calloc(num_columns, sizeof(double));
        thread_sum_squares[t] = (double*)calloc(num_columns, sizeof(double));
        thread_counts[t] = (int*)calloc(num_columns, sizeof(int));
        
        if (!thread_sums[t] || !thread_sum_squares[t] || !thread_counts[t]) {
            fprintf(stderr, "Memory allocation failed for thread computation data\n");
            // Clean up (skipping detailed cleanup for brevity)
            exit(EXIT_FAILURE);
        }
    }
    
    // Process data rows
    printf("%sProcessing data rows...%s\n", YELLOW_COLOR, RESET_COLOR);
    
    int max_columns_found = 0;
    char** line_buffer = (char**)malloc(num_threads * sizeof(char*));
    for (int t = 0; t < num_threads; t++) {
        line_buffer[t] = (char*)malloc(MAX_LINE_LENGTH);
    }
    
    long row = 0; // Skip header
    long last_progress_report = 0;
    
    // Read all lines into an array for parallel processing
    char** all_lines = (char**)malloc(total_lines * sizeof(char*));
    if (!all_lines) {
        fprintf(stderr, "Memory allocation failed for lines array\n");
        exit(EXIT_FAILURE);
    }
    
    // Read lines
    rewind(file); // Start from beginning
    // Skip header, properly checking the return value
    if (fgets(line, MAX_LINE_LENGTH, file) == NULL) {
        fprintf(stderr, "Error skipping header line\n");
        // Clean up code here
        return EXIT_FAILURE;
    }
    
    long lines_read = 0;
    printf("%sReading file into memory...%s\n", YELLOW_COLOR, RESET_COLOR);
    
    while (lines_read < total_lines - 1 && fgets(line, MAX_LINE_LENGTH, file)) {
        all_lines[lines_read] = strdup(line);
        lines_read++;
        
        if (lines_read % progress_step == 0) {
            double progress_pct = 100.0 * lines_read / (total_lines - 1);
            printf("%sReading: %.2f%% complete (%ld/%ld lines)%s\r", 
                   YELLOW_COLOR, progress_pct, lines_read, total_lines - 1, RESET_COLOR);
            fflush(stdout);
        }
    }
    printf("\n%sRead %ld lines into memory%s\n", GREEN_COLOR, lines_read, RESET_COLOR);
    
    // Process lines in parallel
    printf("%sAnalyzing data in parallel...%s\n", YELLOW_COLOR, RESET_COLOR);
    
    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        int local_max_columns = 0;
        
        #pragma omp for schedule(dynamic, 1000)
        for (long i = 0; i < lines_read; i++) {
            char* current_line = all_lines[i];
            strcpy(line_buffer[tid], current_line);
            
            // Remove newline character
            line_buffer[tid][strcspn(line_buffer[tid], "\n")] = 0;
            line_buffer[tid][strcspn(line_buffer[tid], "\r")] = 0;
            
            // Parse data line
            int count = parse_csv_line(line_buffer[tid], thread_data_tokens[tid], MAX_COLUMNS);
            
            if (count > local_max_columns) {
                local_max_columns = count;
            }
            
            // Check if reference column has a valid value for this row
            bool ref_col_valid = ref_col_idx >= 0 && 
                                ref_col_idx < count && 
                                thread_data_tokens[tid][ref_col_idx] && 
                                is_numeric(thread_data_tokens[tid][ref_col_idx]);
            
            double ref_value = 0.0;
            if (ref_col_valid) {
                ref_value = atof(thread_data_tokens[tid][ref_col_idx]);
            }
            
            // Process each column in this row
            int process_count = count < num_columns ? count : num_columns;
            for (int j = 0; j < process_count; j++) {
                if (thread_data_tokens[tid][j] && is_numeric(thread_data_tokens[tid][j])) {
                    double value = atof(thread_data_tokens[tid][j]);
                    thread_sums[tid][j] += value;
                    thread_sum_squares[tid][j] += value * value;
                    thread_counts[tid][j]++;
                    
                    // For correlation calculation with reference column
                    if (ref_col_valid && j != ref_col_idx) {
                        *(thread_products[tid][j]) += value * ref_value;
                    }
                }
            }
            
            #pragma omp critical
            {
                row++;
                if (row - last_progress_report >= progress_step) {
                    double progress_pct = 100.0 * row / lines_read;
                    printf("%sProcessing: %.2f%% complete (%ld/%ld lines)%s\r", 
                           YELLOW_COLOR, progress_pct, row, lines_read, RESET_COLOR);
                    fflush(stdout);
                    last_progress_report = row;
                }
            }
        }
        
        #pragma omp critical
        {
            if (local_max_columns > max_columns_found) {
                max_columns_found = local_max_columns;
            }
        }
    }
    printf("\n");
    
    // Create arrays to store products with reference column for correlation calculation
    double* sum_products = NULL;
    if (ref_col_idx >= 0) {
        sum_products = (double*)calloc(num_columns, sizeof(double));
    }
    
    // Consolidate thread-local results
    for (int t = 0; t < num_threads; t++) {
        for (int i = 0; i < num_columns; i++) {
            columns[i].sum += thread_sums[t][i];
            columns[i].sum_squares += thread_sum_squares[t][i];
            columns[i].count += thread_counts[t][i];
            
            // For correlation
            if (ref_col_idx >= 0 && i != ref_col_idx) {
                sum_products[i] += *(thread_products[t][i]);
            }
        }
    }
    
    // Free the processed lines
    for (long i = 0; i < lines_read; i++) {
        free(all_lines[i]);
    }
    free(all_lines);
    
    // Free thread-local data
    for (int t = 0; t < num_threads; t++) {
        free(thread_sums[t]);
        free(thread_sum_squares[t]);
        free(thread_counts[t]);
        free(thread_data_tokens[t]);
        free(line_buffer[t]);
        
        if (ref_col_idx >= 0) {
            for (int i = 0; i < num_columns; i++) {
                free(thread_products[t][i]);
            }
            free(thread_products[t]);
        }
    }
    free(thread_sums);
    free(thread_sum_squares);
    free(thread_counts);
    free(thread_data_tokens);
    free(line_buffer);
    if (ref_col_idx >= 0) {
        free(thread_products);
    }
    
    printf("%sProcessed a total of %ld rows%s\n", GREEN_COLOR, row, RESET_COLOR);
    printf("%sMaximum columns found in any row: %d%s\n", GREEN_COLOR, max_columns_found, RESET_COLOR);
    
    // Calculate variance for each column
    printf("%sCalculating variance...%s\n", YELLOW_COLOR, RESET_COLOR);
    int zero_variance_count = 0;
    int no_data_count = 0;
    
    for (int i = 0; i < num_columns; i++) {
        if (columns[i].count > 1) {
            double mean = columns[i].sum / columns[i].count;
            columns[i].variance = (columns[i].sum_squares / columns[i].count) - (mean * mean);
            
            // Handle potential floating point errors for small variances
            if (columns[i].variance < 0 && columns[i].variance > -1e-10) {
                columns[i].variance = 0;
            }
            
            columns[i].std_dev = sqrt(fabs(columns[i].variance));
            
            // Check if variance is zero or very close to zero (floating point precision)
            if (fabs(columns[i].variance) < 1e-10) {
                columns[i].has_variance = false;
                zero_variance_count++;
            } else {
                columns[i].has_variance = true;
            }
        } else {
            columns[i].has_variance = false;
            no_data_count++;
        }
    }
    
    // Calculate correlations with reference column if specified
    if (ref_col_idx >= 0 && columns[ref_col_idx].has_variance) {
        printf("%sCalculating correlations with reference column '%s'...%s\n", 
               YELLOW_COLOR, reference_column, RESET_COLOR);
        
        double ref_mean = columns[ref_col_idx].sum / columns[ref_col_idx].count;
        double ref_std_dev = columns[ref_col_idx].std_dev;
        
        for (int i = 0; i < num_columns; i++) {
            if (i != ref_col_idx && columns[i].has_variance) {
                double mean_i = columns[i].sum / columns[i].count;
                double n = (double)(columns[i].count < columns[ref_col_idx].count ? 
                  columns[i].count : columns[ref_col_idx].count);
                
                double covariance = sum_products[i] / n - (mean_i * ref_mean);
                double denominator = columns[i].std_dev * ref_std_dev;
                
                if (denominator != 0) {
                    columns[i].correlation = covariance / denominator;
                } else {
                    columns[i].correlation = 0;
                }
            } else if (i == ref_col_idx) {
                columns[i].correlation = 1.0; // Reference column correlates perfectly with itself
            }
        }
        
        // Sort columns by correlation magnitude for better reporting
        qsort(columns, num_columns, sizeof(Column), compare_columns_by_correlation);
    } else {
        // If no reference column or it has zero variance, sort by variance presence
        qsort(columns, num_columns, sizeof(Column), compare_columns_by_variance);
    }
    
    // Free correlation calculation data
    if (sum_products) {
        free(sum_products);
    }
    
    // Print results based on whether we're showing correlations or just variance
    if (ref_col_idx >= 0) {
        printf("\n%sColumn,Count,Variance,Correlation with %s%s\n", 
               GREEN_COLOR, reference_column, RESET_COLOR);
        
        for (int i = 0; i < num_columns; i++) {
            char variance_str[50];
            if (columns[i].count > 1) {
                snprintf(variance_str, sizeof(variance_str), "%f", columns[i].variance);
            } else {
                strcpy(variance_str, "N/A");
            }
            
            char correlation_str[50];
            if (columns[i].count > 1 && columns[i].has_variance) {
                snprintf(correlation_str, sizeof(correlation_str), "%f", columns[i].correlation);
            } else {
                strcpy(correlation_str, "N/A");
            }
            
            // Highlight different correlation ranges with colors
            if (i == ref_col_idx) {
                // Reference column itself is highlighted differently
                printf("%s%s,%d,%s,%s%s\n", 
                       BLUE_COLOR, columns[i].name, columns[i].count, 
                       variance_str, correlation_str, RESET_COLOR);
            } else if (!columns[i].has_variance || columns[i].count <= 1) {
                // Red for zero variance or no data
                printf("%s%s,%d,%s,%s%s\n", 
                       RED_COLOR, columns[i].name, columns[i].count, 
                       variance_str, correlation_str, RESET_COLOR);
            } else {
                // Color based on correlation strength
                double abs_corr = fabs(columns[i].correlation);
                
                if (abs_corr > 0.7) {
                    // Strong correlation (positive or negative)
                    printf("%s%s,%d,%s,%s%s\n", 
                           MAGENTA_COLOR, columns[i].name, columns[i].count, 
                           variance_str, correlation_str, RESET_COLOR);
                } else if (abs_corr > 0.3) {
                    // Moderate correlation
                    printf("%s%s,%d,%s,%s%s\n", 
                           YELLOW_COLOR, columns[i].name, columns[i].count, 
                           variance_str, correlation_str, RESET_COLOR);
                } else {
                    // Weak or no correlation
                    printf("%s,%d,%s,%s\n", 
                           columns[i].name, columns[i].count, 
                           variance_str, correlation_str);
                }
            }
        }
        
        // Print correlation summary statistics
        int strong_pos_corr = 0;
        int strong_neg_corr = 0;
        int moderate_corr = 0;
        int weak_corr = 0;
        
        for (int i = 0; i < num_columns; i++) {
            if (i != ref_col_idx && columns[i].has_variance && columns[i].count > 1) {
                double corr = columns[i].correlation;
                if (corr > 0.7) strong_pos_corr++;
                else if (corr < -0.7) strong_neg_corr++;
                else if (fabs(corr) > 0.3) moderate_corr++;
                else weak_corr++;
            }
        }
        
        printf("\n%sCorrelation Summary:%s\n", GREEN_COLOR, RESET_COLOR);
        printf("%sStrong positive correlations (r > 0.7): %d%s\n", 
               MAGENTA_COLOR, strong_pos_corr, RESET_COLOR);
        printf("%sStrong negative correlations (r < -0.7): %d%s\n", 
               MAGENTA_COLOR, strong_neg_corr, RESET_COLOR);
        printf("%sModerate correlations (0.3 < |r| < 0.7): %d%s\n", 
               YELLOW_COLOR, moderate_corr, RESET_COLOR);
        printf("Weak correlations (|r| < 0.3): %d\n", weak_corr);
    } else {
        // Standard variance output without correlations
        printf("\n%sColumn,Count,Variance%s\n", GREEN_COLOR, RESET_COLOR);
        for (int i = 0; i < num_columns; i++) {
            if (columns[i].count > 1) {
                if (!columns[i].has_variance) {
                    printf("%s%s,%d,%f%s\n", RED_COLOR, columns[i].name, columns[i].count, columns[i].variance, RESET_COLOR);
                } else {
                    printf("%s,%d,%f\n", columns[i].name, columns[i].count, columns[i].variance);
                }
            } else {
                printf("%s%s,%d,N/A%s\n", RED_COLOR, columns[i].name, columns[i].count, RESET_COLOR);
            }
        }
    }
    
    // Print variance summary statistics
    printf("\n%sVariance Summary:%s\n", GREEN_COLOR, RESET_COLOR);
    printf("Total columns: %d\n", num_columns);
    printf("Columns with zero variance: %s%d (%.2f%%)%s\n", 
           RED_COLOR, zero_variance_count, 
           (float)zero_variance_count / num_columns * 100.0, 
           RESET_COLOR);
printf("Columns with no data: %s%d (%.2f%%)%s\n", 
           RED_COLOR, no_data_count, 
           (float)no_data_count / num_columns * 100.0, 
           RESET_COLOR);
    printf("Columns with non-zero variance: %d (%.2f%%)\n", 
           num_columns - zero_variance_count - no_data_count,
           (float)(num_columns - zero_variance_count - no_data_count) / num_columns * 100.0);
    
    // Print execution time
    double end_time = omp_get_wtime();
    printf("\n%sExecution time: %.2f seconds%s\n", GREEN_COLOR, end_time - start_time, RESET_COLOR);
    
    // Clean up
    for (int i = 0; i < num_columns; i++) {
        free(columns[i].name);
    }
    for (int i = 0; i < MAX_COLUMNS; i++) {
        free(columns[i].data);
    }
    free(columns);
    free(header_tokens);
    free(line);
    fclose(file);
    
    return EXIT_SUCCESS;
}