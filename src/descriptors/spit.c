#include "../cregistry.h" // Required by ChemTrain system
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

// -------- Constants and Macros --------
#define SPIT_DIMS 4 // Number of features in the SPIT tensor: is_same_token, is_bond_pair, is_ring_related, normalized_index_distance

// Access macro for the flattened SPIT tensor
#define SPIT_TENSOR_IDX(i, j, k, num_tokens) (((i) * (num_tokens) + (j)) * SPIT_DIMS + (k))

// -------- Type Definitions --------
typedef enum {
    TOKEN_TYPE_UNKNOWN = 0,
    TOKEN_TYPE_ATOM,              // General atoms like C, N, O, S
    TOKEN_TYPE_BRACKET_ATOM,      // Atoms in brackets like [nH+], [CH2]
    TOKEN_TYPE_AROMATIC_LOWERCASE_ATOM, // c, n, o, s, p (lowercase aromatic)
    TOKEN_TYPE_BOND,              // -, =, #, $, :
    TOKEN_TYPE_RING_CLOSURE,      // 1, 2, %10
    TOKEN_TYPE_BRANCH_OPEN,       // (
    TOKEN_TYPE_BRANCH_CLOSE,      // )
    TOKEN_TYPE_DOT                // . (for disconnected structures)
} TokenType;

typedef struct {
    char* value;       // Dynamically allocated string for the token
    TokenType type;    // Type of the token
} SmilesToken;

typedef struct {
    SmilesToken* tokens; // Dynamically allocated array of SmilesToken
    int count;           // Number of tokens
    int capacity;        // Current capacity of the tokens array
} TokenizedSmiles;

typedef struct {
    TokenizedSmiles tokenized_smiles;
    double* spit_tensor_flat; // Flattened N x N x D tensor
    int N;                    // Number of tokens (dimension of the tensor)
    int tensor_valid;         // Flag indicating if tensor was successfully created
} SpitContext;


// -------- Forward Declarations for Static Helper Functions --------
static void init_tokenized_smiles(TokenizedSmiles* ts);
static int add_smiles_token(TokenizedSmiles* ts, const char* value_start, int len, TokenType type);
static void classify_and_add_token(TokenizedSmiles* ts, const char* current_char, int len, int is_bracket_atom, int is_aromatic_lowercase);

// -------- SMILES Tokenizer Implementation --------

static void init_tokenized_smiles(TokenizedSmiles* ts) {
    ts->tokens = NULL;
    ts->count = 0;
    ts->capacity = 0;
}

static int add_smiles_token(TokenizedSmiles* ts, const char* value_start, int len, TokenType type) {
    if (ts->count >= ts->capacity) {
        int new_capacity = (ts->capacity == 0) ? 16 : ts->capacity * 2;
        SmilesToken* new_tokens = (SmilesToken*)realloc(ts->tokens, new_capacity * sizeof(SmilesToken));
        if (!new_tokens) {
            fprintf(stderr, "SPIT: Failed to reallocate memory for tokens\n");
            return 0; // Failure
        }
        ts->tokens = new_tokens;
        ts->capacity = new_capacity;
    }

    char* token_value = (char*)malloc(len + 1);
    if (!token_value) {
        fprintf(stderr, "SPIT: Failed to allocate memory for token value\n");
        return 0; // Failure
    }
    strncpy(token_value, value_start, len);
    token_value[len] = '\0';

    ts->tokens[ts->count].value = token_value;
    ts->tokens[ts->count].type = type;
    ts->count++;
    return 1; // Success
}

static int is_two_letter_element(char c1, char c2) {
    if (!isalpha(c1) || !islower(c2)) return 0;
    // Common two-letter elements
    const char* two_letter_elements[] = {
        "Ac", "Ag", "Al", "Am", "Ar", "As", "At", "Au", "Ba", "Bh", "Bi", "Bk", "Br",
        "Cd", "Ce", "Cf", "Cl", "Cm", "Cn", "Co", "Cr", "Cs", "Cu", "Db", "Ds", "Dy",
        "Er", "Es", "Eu", "Fe", "Fl", "Fm", "Fr", "Ga", "Gd", "Ge", "He", "Hf", "Hg",
        "Ho", "Hs", "In", "Ir", "Kr", "La", "Li", "Lr", "Lu", "Mc", "Md", "Mg", "Mn",
        "Mo", "Mt", "Na", "Nb", "Nd", "Ne", "Nh", "Ni", "No", "Np", "Os", "Pa", "Pb",
        "Pd", "Pm", "Po", "Pr", "Pt", "Pu", "Ra", "Rb", "Re", "Rf", "Rg", "Rh", "Rn",
        "Ru", "Sb", "Sc", "Se", "Sg", "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te",
        "Th", "Ti", "Tl", "Tm", "Ts", "Og", "Yb", "Zn", "Zr"
        // Note: "U" is single letter, this list is for two-letter ones starting with U like Uue (if they exist and are common)
        // but generally we are looking for Xy pattern.
    };
    char symbol[3] = {c1, c2, '\0'};
    for (size_t i = 0; i < sizeof(two_letter_elements) / sizeof(two_letter_elements[0]); ++i) {
        if (strcmp(symbol, two_letter_elements[i]) == 0) {
            return 1;
        }
    }
    return 0;
}


TokenizedSmiles tokenize_smiles(const char* smiles) {
    TokenizedSmiles ts;
    init_tokenized_smiles(&ts);
    if (!smiles) return ts;

    const char* p = smiles;
    while (*p) {
        int token_len = 0;
        TokenType current_type = TOKEN_TYPE_UNKNOWN;

        if (*p == '[') { // Bracket atom
            const char* end_bracket = strchr(p, ']');
            if (end_bracket) {
                token_len = (end_bracket - p) + 1;
                current_type = TOKEN_TYPE_BRACKET_ATOM;
            } else { // Unterminated bracket, treat as error or single char
                token_len = 1;
                current_type = TOKEN_TYPE_UNKNOWN; // Or handle as error
            }
        } else if (isalpha(*p)) { // Atom
            if (isupper(*p) && isalpha(*(p + 1)) && islower(*(p + 1)) && is_two_letter_element(*p, *(p+1))) {
                token_len = 2; // Two-letter element like Cl, Br
                current_type = TOKEN_TYPE_ATOM;
            } else {
                token_len = 1; // Single-letter element or aromatic
                if (islower(*p)) { // c, n, o, s, p
                     // Check if it's a known aromatic e.g. c,n,o,s,p (sometimes se, as)
                    if (strchr("cnosp", *p)) {
                        current_type = TOKEN_TYPE_AROMATIC_LOWERCASE_ATOM;
                    } else { // other lowercase, less common or could be part of two-letter like 'se', 'as' if not caught by two-letter
                        current_type = TOKEN_TYPE_ATOM; // Or unknown if strict
                    }
                } else {
                    current_type = TOKEN_TYPE_ATOM; // C, N, O
                }
            }
        } else if (isdigit(*p)) { // Ring closure
            token_len = 1;
            current_type = TOKEN_TYPE_RING_CLOSURE;
        } else if (*p == '%') { // Extended ring closure (e.g., %10)
            if (isdigit(*(p + 1)) && isdigit(*(p + 2))) {
                token_len = 3;
            } else { // malformed, treat as single char
                token_len = 1;
            }
            current_type = TOKEN_TYPE_RING_CLOSURE;
        } else if (strchr("=~#$:", *p)) { // Bonds
            token_len = 1;
            current_type = TOKEN_TYPE_BOND;
        } else if (*p == '-') { // Explicit single bond
             token_len = 1;
             current_type = TOKEN_TYPE_BOND;
        } else if (*p == '(') {
            token_len = 1;
            current_type = TOKEN_TYPE_BRANCH_OPEN;
        } else if (*p == ')') {
            token_len = 1;
            current_type = TOKEN_TYPE_BRANCH_CLOSE;
        } else if (*p == '.') {
            token_len = 1;
            current_type = TOKEN_TYPE_DOT;
        } else if (*p == '*') { // Wildcard atom
            token_len = 1;
            current_type = TOKEN_TYPE_ATOM; // Treat '*' as an atom type
        }
        else { // Other characters (stereochemistry, etc.) - currently marked UNKNOWN
            token_len = 1;
            current_type = TOKEN_TYPE_UNKNOWN;
        }

        if (token_len == 0) { // Should not happen if logic is correct
            p++; continue;
        }
        
        if (!add_smiles_token(&ts, p, token_len, current_type)) {
            // Allocation failed, free what's been allocated and return empty
            for (int i = 0; i < ts.count; ++i) free(ts.tokens[i].value);
            free(ts.tokens);
            init_tokenized_smiles(&ts); // Reset to empty state
            return ts;
        }
        p += token_len;
    }
    return ts;
}

void free_tokenized_smiles(TokenizedSmiles* ts) {
    if (!ts) return;
    for (int i = 0; i < ts->count; ++i) {
        free(ts->tokens[i].value);
    }
    free(ts->tokens);
    ts->tokens = NULL;
    ts->count = 0;
    ts->capacity = 0;
}

// -------- SPIT Context Management --------
static int is_spit_feature_atom_type(TokenType type) {
    return type == TOKEN_TYPE_ATOM || type == TOKEN_TYPE_BRACKET_ATOM || type == TOKEN_TYPE_AROMATIC_LOWERCASE_ATOM;
}

static int is_spit_feature_bond_type(TokenType type) {
    return type == TOKEN_TYPE_BOND;
}

static int is_spit_feature_ring_type(TokenType type) {
    return type == TOKEN_TYPE_RING_CLOSURE;
}


SpitContext* create_spit_context(const char* smiles) {
    SpitContext* ctx = (SpitContext*)malloc(sizeof(SpitContext));
    if (!ctx) {
        fprintf(stderr, "SPIT: Failed to allocate memory for SpitContext\n");
        return NULL;
    }
    ctx->spit_tensor_flat = NULL;
    ctx->N = 0;
    ctx->tensor_valid = 0;
    init_tokenized_smiles(&ctx->tokenized_smiles);

    if (!smiles || *smiles == '\0') {
        return ctx; // Return context with N=0, tensor not allocated
    }

    ctx->tokenized_smiles = tokenize_smiles(smiles);
    ctx->N = ctx->tokenized_smiles.count;

    if (ctx->N == 0) {
        // Tokenization might have failed or SMILES was effectively empty
        // Free tokens if any were partially made before failure in tokenize_smiles
        free_tokenized_smiles(&ctx->tokenized_smiles);
        return ctx; 
    }

    size_t tensor_size = (size_t)ctx->N * ctx->N * SPIT_DIMS;
    ctx->spit_tensor_flat = (double*)malloc(tensor_size * sizeof(double));
    if (!ctx->spit_tensor_flat) {
        fprintf(stderr, "SPIT: Failed to allocate memory for SPIT tensor\n");
        free_tokenized_smiles(&ctx->tokenized_smiles);
        // ctx->N remains, but tensor_valid will be 0
        return ctx;
    }
    // Initialize tensor to zero to be safe
    memset(ctx->spit_tensor_flat, 0, tensor_size * sizeof(double));


    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            SmilesToken* token_i = &ctx->tokenized_smiles.tokens[i];
            SmilesToken* token_j = &ctx->tokenized_smiles.tokens[j];

            // Feature 0: is_same_token
            ctx->spit_tensor_flat[SPIT_TENSOR_IDX(i, j, 0, ctx->N)] = (strcmp(token_i->value, token_j->value) == 0) ? 1.0 : 0.0;

            // Feature 1: is_bond_pair (atom-bond or bond-atom type interaction)
            int ti_is_atom = is_spit_feature_atom_type(token_i->type);
            int tj_is_atom = is_spit_feature_atom_type(token_j->type);
            int ti_is_bond = is_spit_feature_bond_type(token_i->type);
            int tj_is_bond = is_spit_feature_bond_type(token_j->type);
            ctx->spit_tensor_flat[SPIT_TENSOR_IDX(i, j, 1, ctx->N)] = ((ti_is_atom && tj_is_bond) || (ti_is_bond && tj_is_atom)) ? 1.0 : 0.0;
            
            // Feature 2: is_ring_related (both are same ring closure tokens)
            int ti_is_ring = is_spit_feature_ring_type(token_i->type);
            int tj_is_ring = is_spit_feature_ring_type(token_j->type);
            ctx->spit_tensor_flat[SPIT_TENSOR_IDX(i, j, 2, ctx->N)] = (ti_is_ring && tj_is_ring && strcmp(token_i->value, token_j->value) == 0) ? 1.0 : 0.0;

            // Feature 3: normalized_index_distance
            if (ctx->N > 1) {
                ctx->spit_tensor_flat[SPIT_TENSOR_IDX(i, j, 3, ctx->N)] = (double)abs(i - j) / (ctx->N - 1.0);
            } else {
                ctx->spit_tensor_flat[SPIT_TENSOR_IDX(i, j, 3, ctx->N)] = 0.0;
            }
        }
    }
    ctx->tensor_valid = 1;
    return ctx;
}

void free_spit_context(SpitContext* ctx) {
    if (!ctx) return;
    free_tokenized_smiles(&ctx->tokenized_smiles);
    free(ctx->spit_tensor_flat);
    free(ctx);
}

// Helper to safely get SPIT value (not strictly necessary if direct access is careful)
static double get_spit_value(const SpitContext* ctx, int r, int c, int feature_idx) {
    if (!ctx || !ctx->tensor_valid || r < 0 || r >= ctx->N || c < 0 || c >= ctx->N || feature_idx < 0 || feature_idx >= SPIT_DIMS) {
        // This case should ideally not be hit if N is checked by caller
        // and indices are always valid.
        // For robustness, could return NaN or a specific error value.
        // For now, behavior is undefined or might crash if tensor_valid is false.
        // If tensor_valid is true, then access is as below.
        if (!ctx->tensor_valid) return 0.0; // Or NaN
    }
    return ctx->spit_tensor_flat[SPIT_TENSOR_IDX(r, c, feature_idx, ctx->N)];
}


// -------- General Math Helper functions --------
static double mean_array(const double* arr, int count) {
    if (count == 0) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < count; ++i) sum += arr[i];
    return sum / count;
}

static double variance_array(const double* arr, int count, double mean_val) {
    if (count < 2) return 0.0; // Variance is typically undefined or 0 for less than 2 samples
    double sum_sq_diff = 0.0;
    for (int i = 0; i < count; ++i) {
        double diff = arr[i] - mean_val;
        sum_sq_diff += diff * diff;
    }
    return sum_sq_diff / (count -1); // Sample variance
}

// -------- Descriptor Implementations --------

// Category 1: Symmetry & Positional Analysis

// 1. Token symmetry index
// Interpretation: Density of identical token pairs (i != j)
double calculateSpitTokenSymmetryIndex(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double sum_identical_pairs = 0.0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = i + 1; j < ctx->N; ++j) { // Iterate over unique pairs (i < j)
            sum_identical_pairs += get_spit_value(ctx, i, j, 0); // spit[i][j][0] is is_same_token
        }
    }
    
    double total_possible_pairs = (double)ctx->N * (ctx->N - 1) / 2.0;
    double result = (total_possible_pairs > 0) ? (sum_identical_pairs / total_possible_pairs) : 0.0;
    
    free_spit_context(ctx);
    return result;
}

// 4. Mean token-index distance
// Average of normalized_index_distance over all pairs (i, j)
double calculateSpitMeanTokenIndexDistance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N == 0) {
        free_spit_context(ctx);
        return 0.0;
    }

    if (ctx->N < 2) { // No distance if less than 2 tokens
         free_spit_context(ctx);
         return 0.0;
    }

    double sum_distances = 0.0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            sum_distances += get_spit_value(ctx, i, j, 3); // spit[i][j][3] is normalized_index_distance
        }
    }
    
    double result = sum_distances / ( (double)ctx->N * ctx->N );
    free_spit_context(ctx);
    return result;
}


// Category 2: Pattern Frequency Descriptors

// 11. Count of matching token pairs
// Sum of is_same_token for i < j
double calculateSpitMatchingTokenPairCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double count = 0.0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = i + 1; j < ctx->N; ++j) {
            if (get_spit_value(ctx, i, j, 0) == 1.0) { // is_same_token
                count += 1.0;
            }
        }
    }
    free_spit_context(ctx);
    return count;
}

// 20. Activation density (non-zero ratio)
// Ratio of non-zero cells in the entire SPIT tensor to total cells
double calculateSpitActivationDensity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
     if (!ctx || !ctx->tensor_valid || ctx->N == 0) {
        free_spit_context(ctx);
        return 0.0;
    }

    double non_zero_cells = 0.0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            for (int k = 0; k < SPIT_DIMS; ++k) {
                if (fabs(get_spit_value(ctx, i, j, k)) > 1e-9) { // Check for non-zero (floating point)
                    non_zero_cells += 1.0;
                }
            }
        }
    }
    
    double total_cells = (double)ctx->N * ctx->N * SPIT_DIMS;
    double result = (total_cells > 0) ? (non_zero_cells / total_cells) : 0.0;
    
    free_spit_context(ctx);
    return result;
}

// Category 3: Statistical Distributions

// 21. Mean of each feature plane (averaged)
// Average of the means of the 4 feature planes
double calculateSpitMeanOfFeaturePlanes(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N == 0) {
        free_spit_context(ctx);
        return 0.0;
    }

    double plane_means[SPIT_DIMS] = {0.0};
    double total_elements_per_plane = (double)ctx->N * ctx->N;

    if (total_elements_per_plane == 0) {
        free_spit_context(ctx);
        return 0.0;
    }

    for (int k = 0; k < SPIT_DIMS; ++k) {
        double sum_plane_k = 0.0;
        for (int i = 0; i < ctx->N; ++i) {
            for (int j = 0; j < ctx->N; ++j) {
                sum_plane_k += get_spit_value(ctx, i, j, k);
            }
        }
        plane_means[k] = sum_plane_k / total_elements_per_plane;
    }
    
    double sum_of_plane_means = 0.0;
    for(int k=0; k < SPIT_DIMS; ++k) {
        sum_of_plane_means += plane_means[k];
    }
    double result = sum_of_plane_means / SPIT_DIMS;

    free_spit_context(ctx);
    return result;
}


// Category 4: Token Role Differentiation

// 31. Most "connected" token (row sum) - average over features
// Max of (sum over j and k of spit[i][j][k]) for each i
double calculateSpitMostConnectedTokenScore(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N == 0) {
        free_spit_context(ctx);
        return 0.0;
    }

    double max_row_sum = -1.0; // Assuming activations are non-negative or this needs adjustment

    for (int i = 0; i < ctx->N; ++i) {
        double current_row_sum = 0.0;
        for (int j = 0; j < ctx->N; ++j) {
            for (int k = 0; k < SPIT_DIMS; ++k) {
                current_row_sum += get_spit_value(ctx, i, j, k);
            }
        }
        if (i == 0 || current_row_sum > max_row_sum) {
            max_row_sum = current_row_sum;
        }
    }
    
    free_spit_context(ctx);
    return (max_row_sum < 0 && ctx->N > 0) ? 0.0 : max_row_sum ; // if all sums were 0 or N=0
}

// Category 5: Motif & Texture Descriptors

// 46. Matrix trace score (sum i==j) - averaged over features
// Average of traces of the 4 feature planes
double calculateSpitMatrixTraceScore(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N == 0) {
        free_spit_context(ctx);
        return 0.0;
    }

    double feature_traces[SPIT_DIMS] = {0.0};
    for (int k = 0; k < SPIT_DIMS; ++k) {
        for (int i = 0; i < ctx->N; ++i) {
            feature_traces[k] += get_spit_value(ctx, i, i, k);
        }
    }

    double sum_of_traces = 0.0;
    for(int k=0; k < SPIT_DIMS; ++k) {
        sum_of_traces += feature_traces[k];
    }
    double result = (SPIT_DIMS > 0) ? (sum_of_traces / SPIT_DIMS) : 0.0;
    
    free_spit_context(ctx);
    return result;
}


// -------- Stubs/Placeholders for remaining descriptors --------
// Please implement these based on their descriptions, using the SpitContext.
// The interpretations might require further refinement.

// Category 1: Symmetry & Positional Analysis (Continued)
double calculateSpitRingSymmetryIndex(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double ring_pair_count = 0.0;
    double total_ring_tokens = 0.0;
    
    // Count total ring tokens
    for (int i = 0; i < ctx->N; ++i) {
        if (is_spit_feature_ring_type(ctx->tokenized_smiles.tokens[i].type)) {
            total_ring_tokens += 1.0;
        }
    }
    
    // Count matching ring pairs
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = i + 1; j < ctx->N; ++j) {
            if (get_spit_value(ctx, i, j, 2) > 0.0) { // is_ring_related
                ring_pair_count += 1.0;
            }
        }
    }
    
    double result = (total_ring_tokens > 0) ? (2.0 * ring_pair_count / total_ring_tokens) : 0.0;
    free_spit_context(ctx);
    return result;
}

double calculateSpitBondPairAsymmetry(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 3) {
        free_spit_context(ctx);
        return 0.0;
    }
    
    // For each atom, measure asymmetry of bond interactions with other atoms
    double asymmetry_sum = 0.0;
    int atom_count = 0;
    
    for (int i = 0; i < ctx->N; ++i) {
        if (is_spit_feature_atom_type(ctx->tokenized_smiles.tokens[i].type)) {
            double left_bonds = 0.0, right_bonds = 0.0;
            
            // Count bonds to the left
            for (int j = 0; j < i; ++j) {
                left_bonds += get_spit_value(ctx, i, j, 1);
            }
            
            // Count bonds to the right
            for (int j = i + 1; j < ctx->N; ++j) {
                right_bonds += get_spit_value(ctx, i, j, 1);
            }
            
            // Calculate asymmetry for this atom
            double total_bonds = left_bonds + right_bonds;
            double atom_asymmetry = (total_bonds > 0) ? 
                fabs(left_bonds - right_bonds) / total_bonds : 0.0;
            
            asymmetry_sum += atom_asymmetry;
            atom_count++;
        }
    }
    
    double result = (atom_count > 0) ? (asymmetry_sum / atom_count) : 0.0;
    free_spit_context(ctx);
    return result;
}

double calculateSpitStdDevTokenIndexDistances(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    // First calculate mean
    double sum = 0.0;
    double count = 0.0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            sum += get_spit_value(ctx, i, j, 3); // normalized_index_distance
            count += 1.0;
        }
    }
    
    double mean = (count > 0) ? (sum / count) : 0.0;
    
    // Now calculate std dev
    double variance_sum = 0.0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            double diff = get_spit_value(ctx, i, j, 3) - mean;
            variance_sum += diff * diff;
        }
    }
    
    double variance = (count > 1) ? (variance_sum / (count - 1.0)) : 0.0;
    double result = sqrt(variance);
    
    free_spit_context(ctx);
    return result;
}

double calculateSpitMaxDiagonalActivity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 1) {
        free_spit_context(ctx);
        return 0.0;
    }

    double max_diagonal_activity = 0.0;
    for (int i = 0; i < ctx->N; ++i) {
        double sum_activity = 0.0;
        for (int k = 0; k < SPIT_DIMS; ++k) {
            sum_activity += get_spit_value(ctx, i, i, k);
        }
        
        if (i == 0 || sum_activity > max_diagonal_activity) {
            max_diagonal_activity = sum_activity;
        }
    }
    
    free_spit_context(ctx);
    return max_diagonal_activity;
}

double calculateSpitMeanOffDiagonalSimilarity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double sum_similarity = 0.0;
    double count = 0.0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            if (i != j) {
                sum_similarity += get_spit_value(ctx, i, j, 0); // is_same_token
                count += 1.0;
            }
        }
    }
    
    double result = (count > 0) ? (sum_similarity / count) : 0.0;
    free_spit_context(ctx);
    return result;
}

double calculateSpitLowerTriangleActivationRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double lower_sum = 0.0;
    double upper_sum = 0.0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            if (i > j) { // Lower triangle
                for (int k = 0; k < SPIT_DIMS; ++k) {
                    lower_sum += get_spit_value(ctx, i, j, k);
                }
            } else if (i < j) { // Upper triangle
                for (int k = 0; k < SPIT_DIMS; ++k) {
                    upper_sum += get_spit_value(ctx, i, j, k);
                }
            }
            // Skip diagonal (i == j)
        }
    }
    
    double result = (upper_sum > 0) ? (lower_sum / upper_sum) : (lower_sum > 0 ? 1.0 : 0.0);
    free_spit_context(ctx);
    return result;
}

double calculateSpitDiagonalDominanceIndex(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 1) {
        free_spit_context(ctx);
        return 0.0;
    }

    double diagonal_sum = 0.0;
    double off_diagonal_sum = 0.0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            if (i == j) { // Diagonal
                for (int k = 0; k < SPIT_DIMS; ++k) {
                    diagonal_sum += get_spit_value(ctx, i, j, k);
                }
            } else { // Off-diagonal
                for (int k = 0; k < SPIT_DIMS; ++k) {
                    off_diagonal_sum += get_spit_value(ctx, i, j, k);
                }
            }
        }
    }
    
    double result = (off_diagonal_sum > 0) ? (diagonal_sum / off_diagonal_sum) : (diagonal_sum > 0 ? 1.0 : 0.0);
    free_spit_context(ctx);
    return result;
}

double calculateSpitBandwidthTop10Activation(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 3) {
        free_spit_context(ctx);
        return 0.0;
    }

    // Focus on the first feature plane (is_same_token) for simplicity
    int feature_idx = 0;
    
    // Create array of activation values for sorting
    int tensor_size = ctx->N * ctx->N;
    double* all_values = (double*)malloc(tensor_size * sizeof(double));
    if (!all_values) {
        free_spit_context(ctx);
        return 0.0;
    }
    
    // Collect values
    int idx = 0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            all_values[idx++] = get_spit_value(ctx, i, j, feature_idx);
        }
    }
    
    // Simple bubble sort (inefficient but OK for small tensors)
    for (int i = 0; i < tensor_size - 1; ++i) {
        for (int j = 0; j < tensor_size - i - 1; ++j) {
            if (all_values[j] < all_values[j + 1]) {
                double temp = all_values[j];
                all_values[j] = all_values[j + 1];
                all_values[j + 1] = temp;
            }
        }
    }
    
    // Determine threshold for top 10%
    int top_10_percent_count = (tensor_size + 9) / 10; // Ceiling of 10%
    if (top_10_percent_count <= 0) top_10_percent_count = 1;
    if (top_10_percent_count > tensor_size) top_10_percent_count = tensor_size;
    
    double threshold = all_values[top_10_percent_count - 1];
    free(all_values);
    
    // Find max distance between cells above threshold
    int max_distance = 0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            if (get_spit_value(ctx, i, j, feature_idx) >= threshold) {
                for (int k = 0; k < ctx->N; ++k) {
                    for (int l = 0; l < ctx->N; ++l) {
                        if (get_spit_value(ctx, k, l, feature_idx) >= threshold) {
                            int i_dist = abs(i - k);
                            int j_dist = abs(j - l);
                            int manhattan_dist = i_dist + j_dist;
                            if (manhattan_dist > max_distance) {
                                max_distance = manhattan_dist;
                            }
                        }
                    }
                }
            }
        }
    }
    
    double result = (double)max_distance;
    free_spit_context(ctx);
    return result;
}

// Category 2: Pattern Frequency Descriptors (Continued)
double calculateSpitBondAtomInteractionCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double total_interactions = 0.0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = i + 1; j < ctx->N; ++j) { // Only count each pair once
            total_interactions += get_spit_value(ctx, i, j, 1); // is_bond_pair
        }
    }
    
    free_spit_context(ctx);
    return total_interactions;
}

double calculateSpitRingRelatedTokenPairCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double total_pairs = 0.0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = i + 1; j < ctx->N; ++j) { // Only count each pair once
            total_pairs += get_spit_value(ctx, i, j, 2); // is_ring_related
        }
    }
    
    free_spit_context(ctx);
    return total_pairs;
}

double calculateSpitSelfPairCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 1) {
        free_spit_context(ctx);
        return 0.0;
    }

    // Count tokens on diagonal that are important atom types
    double self_pair_count = 0.0;
    
    for (int i = 0; i < ctx->N; ++i) {
        if (is_spit_feature_atom_type(ctx->tokenized_smiles.tokens[i].type) || 
            is_spit_feature_ring_type(ctx->tokenized_smiles.tokens[i].type)) {
            self_pair_count += 1.0;
        }
    }
    
    free_spit_context(ctx);
    return self_pair_count;
}

double calculateSpitRepetitiveTokenMotifCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 3) {
        free_spit_context(ctx);
        return 0.0;
    }

    double motif_count = 0.0;
    
    // Look for sequences of at least 3 identical tokens
    for (int i = 0; i < ctx->N - 2; ++i) {
        for (int j = i + 1; j < ctx->N - 1; ++j) {
            for (int k = j + 1; k < ctx->N; ++k) {
                // If all three tokens are the same
                if (strcmp(ctx->tokenized_smiles.tokens[i].value, ctx->tokenized_smiles.tokens[j].value) == 0 &&
                    strcmp(ctx->tokenized_smiles.tokens[j].value, ctx->tokenized_smiles.tokens[k].value) == 0) {
                    
                    // Check if they are consecutive
                    if (j == i + 1 && k == j + 1) {
                        motif_count += 1.0;
                    }
                }
            }
        }
    }
    
    free_spit_context(ctx);
    return motif_count;
}

double calculateSpitTokenRepeatSpan(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    int max_span = 0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = i + 1; j < ctx->N; ++j) {
            // If tokens are the same
            if (get_spit_value(ctx, i, j, 0) > 0.0) {
                int span = j - i;
                if (span > max_span) {
                    max_span = span;
                }
            }
        }
    }
    
    free_spit_context(ctx);
    return (double)max_span;
}

double calculateSpitMostFrequentInteractionPattern(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    // This is a simplification - we'll just look at token type combinations rather than
    // full tensor patterns, as that would require complex data structures
    int atom_atom = 0, atom_bond = 0, atom_ring = 0;
    int bond_bond = 0, bond_ring = 0, ring_ring = 0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = i + 1; j < ctx->N; ++j) {
            TokenType type_i = ctx->tokenized_smiles.tokens[i].type;
            TokenType type_j = ctx->tokenized_smiles.tokens[j].type;
            
            if (is_spit_feature_atom_type(type_i) && is_spit_feature_atom_type(type_j)) atom_atom++;
            else if ((is_spit_feature_atom_type(type_i) && is_spit_feature_bond_type(type_j)) ||
                     (is_spit_feature_bond_type(type_i) && is_spit_feature_atom_type(type_j))) atom_bond++;
            else if ((is_spit_feature_atom_type(type_i) && is_spit_feature_ring_type(type_j)) ||
                     (is_spit_feature_ring_type(type_i) && is_spit_feature_atom_type(type_j))) atom_ring++;
            else if (is_spit_feature_bond_type(type_i) && is_spit_feature_bond_type(type_j)) bond_bond++;
            else if ((is_spit_feature_bond_type(type_i) && is_spit_feature_ring_type(type_j)) ||
                     (is_spit_feature_ring_type(type_i) && is_spit_feature_bond_type(type_j))) bond_ring++;
            else if (is_spit_feature_ring_type(type_i) && is_spit_feature_ring_type(type_j)) ring_ring++;
        }
    }
    
    // Find the maximum interaction pattern
    int max_count = atom_atom;
    if (atom_bond > max_count) max_count = atom_bond;
    if (atom_ring > max_count) max_count = atom_ring;
    if (bond_bond > max_count) max_count = bond_bond;
    if (bond_ring > max_count) max_count = bond_ring;
    if (ring_ring > max_count) max_count = ring_ring;
    
    free_spit_context(ctx);
    return (double)max_count;
}

double calculateSpitTokenTransitionEntropy(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    // Define transition matrix for token types (8x8)
    const int NUM_TOKEN_TYPES = 8; // Number of token types in the enum
    int transitions[NUM_TOKEN_TYPES][NUM_TOKEN_TYPES] = {0};
    int type_counts[NUM_TOKEN_TYPES] = {0};
    
    // Count transitions
    for (int i = 0; i < ctx->N - 1; ++i) {
        TokenType from_type = ctx->tokenized_smiles.tokens[i].type;
        TokenType to_type = ctx->tokenized_smiles.tokens[i + 1].type;
        
        transitions[from_type][to_type]++;
        type_counts[from_type]++;
    }
    
    // Calculate entropy
    double entropy = 0.0;
    for (int i = 0; i < NUM_TOKEN_TYPES; ++i) {
        if (type_counts[i] > 0) {
            for (int j = 0; j < NUM_TOKEN_TYPES; ++j) {
                if (transitions[i][j] > 0) {
                    double p = (double)transitions[i][j] / type_counts[i];
                    entropy -= p * log2(p);
                }
            }
        }
    }
    
    free_spit_context(ctx);
    return entropy;
}

double calculateSpitNonZeroCellsPerFeature(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 1) {
        free_spit_context(ctx);
        return 0.0;
    }

    double non_zero_counts[SPIT_DIMS] = {0.0};
    
    for (int k = 0; k < SPIT_DIMS; ++k) {
        for (int i = 0; i < ctx->N; ++i) {
            for (int j = 0; j < ctx->N; ++j) {
                if (fabs(get_spit_value(ctx, i, j, k)) > 1e-9) {
                    non_zero_counts[k] += 1.0;
                }
            }
        }
    }
    
    // Average non-zero count across features
    double avg_non_zero = 0.0;
    for (int k = 0; k < SPIT_DIMS; ++k) {
        avg_non_zero += non_zero_counts[k];
    }
    avg_non_zero /= SPIT_DIMS;
    
    free_spit_context(ctx);
    return avg_non_zero;
}

// Category 3: Statistical Distributions (Continued)
double calculateSpitVarianceOfFeaturePlanes(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double feature_variances[SPIT_DIMS] = {0.0};
    
    for (int k = 0; k < SPIT_DIMS; ++k) {
        // Calculate mean for this feature plane
        double sum = 0.0;
        for (int i = 0; i < ctx->N; ++i) {
            for (int j = 0; j < ctx->N; ++j) {
                sum += get_spit_value(ctx, i, j, k);
            }
        }
        double mean = sum / (ctx->N * ctx->N);
        
        // Calculate variance
        double sq_diff_sum = 0.0;
        for (int i = 0; i < ctx->N; ++i) {
            for (int j = 0; j < ctx->N; ++j) {
                double diff = get_spit_value(ctx, i, j, k) - mean;
                sq_diff_sum += diff * diff;
            }
        }
        feature_variances[k] = sq_diff_sum / (ctx->N * ctx->N - 1);
    }
    
    // Average variance across features
    double avg_variance = 0.0;
    for (int k = 0; k < SPIT_DIMS; ++k) {
        avg_variance += feature_variances[k];
    }
    avg_variance /= SPIT_DIMS;
    
    free_spit_context(ctx);
    return avg_variance;
}

double calculateSpitFeatureWiseEntropy(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double feature_entropies[SPIT_DIMS] = {0.0};
    
    for (int k = 0; k < SPIT_DIMS; ++k) {
        // Simple approach: count occurrences of distinct values
        // For binary features (0/1), we count how many 0s and 1s
        int zeros = 0, ones = 0;
        
        for (int i = 0; i < ctx->N; ++i) {
            for (int j = 0; j < ctx->N; ++j) {
                double val = get_spit_value(ctx, i, j, k);
                if (fabs(val) < 1e-9) zeros++;
                else if (fabs(val - 1.0) < 1e-9) ones++;
                // Other values ignored for simplicity
            }
        }
        
        double total = zeros + ones;
        double p_zero = zeros / total;
        double p_one = ones / total;
        
        // Calculate entropy using Shannon formula
        if (p_zero > 0) feature_entropies[k] -= p_zero * log2(p_zero);
        if (p_one > 0) feature_entropies[k] -= p_one * log2(p_one);
    }
    
    // Average entropy across features
    double avg_entropy = 0.0;
    for (int k = 0; k < SPIT_DIMS; ++k) {
        avg_entropy += feature_entropies[k];
    }
    avg_entropy /= SPIT_DIMS;
    
    free_spit_context(ctx);
    return avg_entropy;
}

double calculateSpitFeatureCooccurrenceCorrelation(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double correlation_sum = 0.0;
    int pair_count = 0;
    
    // For simplicity, we'll just use the correlation between features 0 and 1
    // Full implementation would compute correlation between all pairs
    int f1 = 0, f2 = 1;
    
    // Calculate means
    double mean1 = 0.0, mean2 = 0.0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            mean1 += get_spit_value(ctx, i, j, f1);
            mean2 += get_spit_value(ctx, i, j, f2);
        }
    }
    mean1 /= (ctx->N * ctx->N);
    mean2 /= (ctx->N * ctx->N);
    
    // Calculate covariance and variances
    double covariance = 0.0, var1 = 0.0, var2 = 0.0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            double v1 = get_spit_value(ctx, i, j, f1);
            double v2 = get_spit_value(ctx, i, j, f2);
            double d1 = v1 - mean1;
            double d2 = v2 - mean2;
            
            covariance += d1 * d2;
            var1 += d1 * d1;
            var2 += d2 * d2;
        }
    }
    
    // Compute correlation coefficient
    double correlation = 0.0;
    if (var1 > 0 && var2 > 0) {
        correlation = covariance / sqrt(var1 * var2);
    }
    
    free_spit_context(ctx);
    return fabs(correlation); // Return absolute correlation
}

double calculateSpitMeanActivationAcrossRows(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 1) {
        free_spit_context(ctx);
        return 0.0;
    }

    double* row_sums = (double*)calloc(ctx->N, sizeof(double));
    if (!row_sums) {
        free_spit_context(ctx);
        return 0.0;
    }
    
    // Calculate sum for each row
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            for (int k = 0; k < SPIT_DIMS; ++k) {
                row_sums[i] += get_spit_value(ctx, i, j, k);
            }
        }
    }
    
    // Calculate mean of row sums
    double sum_of_row_sums = 0.0;
    for (int i = 0; i < ctx->N; ++i) {
        sum_of_row_sums += row_sums[i];
    }
    double mean_row_sum = sum_of_row_sums / ctx->N;
    
    free(row_sums);
    free_spit_context(ctx);
    return mean_row_sum;
}

double calculateSpitSkewnessInteractionMatrix(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 3) { // Need at least 3 samples for meaningful skewness
        free_spit_context(ctx);
        return 0.0;
    }

    // Calculate mean across all tensor values
    double sum = 0.0;
    int count = ctx->N * ctx->N * SPIT_DIMS;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            for (int k = 0; k < SPIT_DIMS; ++k) {
                sum += get_spit_value(ctx, i, j, k);
            }
        }
    }
    double mean = sum / count;
    
    // Calculate variance and skewness
    double sum_sq_diff = 0.0;
    double sum_cube_diff = 0.0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            for (int k = 0; k < SPIT_DIMS; ++k) {
                double diff = get_spit_value(ctx, i, j, k) - mean;
                sum_sq_diff += diff * diff;
                sum_cube_diff += diff * diff * diff;
            }
        }
    }
    
    double variance = sum_sq_diff / count;
    double std_dev = sqrt(variance);
    double skewness = 0.0;
    
    if (std_dev > 0) {
        skewness = (sum_cube_diff / count) / (std_dev * std_dev * std_dev);
    }
    
    free_spit_context(ctx);
    return skewness;
}

double calculateSpitKurtosisBondPairScores(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 4) { // Need at least 4 samples for meaningful kurtosis
        free_spit_context(ctx);
        return 0.0;
    }

    // Focusing just on bond pair scores (feature 1)
    int feature_idx = 1;
    
    // Calculate mean
    double sum = 0.0;
    int count = ctx->N * ctx->N;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            sum += get_spit_value(ctx, i, j, feature_idx);
        }
    }
    double mean = sum / count;
    
    // Calculate variance and kurtosis
    double sum_sq_diff = 0.0;
    double sum_fourth_diff = 0.0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            double diff = get_spit_value(ctx, i, j, feature_idx) - mean;
            double sq_diff = diff * diff;
            sum_sq_diff += sq_diff;
            sum_fourth_diff += sq_diff * sq_diff;
        }
    }
    
    double variance = sum_sq_diff / count;
    double kurtosis = 0.0;
    
    if (variance > 0) {
        kurtosis = (sum_fourth_diff / count) / (variance * variance) - 3.0; // Excess kurtosis
    }
    
    free_spit_context(ctx);
    return kurtosis;
}

double calculateSpitAverageDistanceWeightedActivation(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 1) {
        free_spit_context(ctx);
        return 0.0;
    }

    double weighted_sum = 0.0;
    double count = ctx->N * ctx->N;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            double distance_weight = 1.0 - get_spit_value(ctx, i, j, 3); // 1 - normalized_index_distance
            double activation_sum = 0.0;
            
            for (int k = 0; k < SPIT_DIMS - 1; ++k) { // Exclude distance feature itself
                activation_sum += get_spit_value(ctx, i, j, k);
            }
            
            weighted_sum += distance_weight * activation_sum;
        }
    }
    
    double result = (count > 0) ? (weighted_sum / count) : 0.0;
    free_spit_context(ctx);
    return result;
}

double calculateSpitFeatureSpecificSparsity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 1) {
        free_spit_context(ctx);
        return 0.0;
    }

    double non_zero_counts[SPIT_DIMS] = {0.0};
    
    for (int k = 0; k < SPIT_DIMS; ++k) {
        for (int i = 0; i < ctx->N; ++i) {
            for (int j = 0; j < ctx->N; ++j) {
                if (fabs(get_spit_value(ctx, i, j, k)) > 1e-9) {
                    non_zero_counts[k] += 1.0;
                }
            }
        }
    }
    
    // Average non-zero count across features
    double avg_non_zero = 0.0;
    for (int k = 0; k < SPIT_DIMS; ++k) {
        avg_non_zero += non_zero_counts[k];
    }
    avg_non_zero /= SPIT_DIMS;
    
    free_spit_context(ctx);
    return avg_non_zero;
}

double calculateSpitFeatureImportanceIndexPCA(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    // This is a placeholder implementation. PCA implementation or library is required.
    // For now, return a default value.
    free_spit_context(ctx);
    return 0.0; // Placeholder return
}

// Category 4: Token Role Differentiation (Continued)
double calculateSpitTokenHighestBondInteractions(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 1) {
        free_spit_context(ctx);
        return 0.0;
    }

    double max_bond_interactions = 0.0;
    int bond_feature_idx = 1; // is_bond_pair feature
    
    // For each token, calculate total bond interactions
    for (int i = 0; i < ctx->N; ++i) {
        double sum_bonds = 0.0;
        for (int j = 0; j < ctx->N; ++j) {
            sum_bonds += get_spit_value(ctx, i, j, bond_feature_idx);
        }
        
        if (i == 0 || sum_bonds > max_bond_interactions) {
            max_bond_interactions = sum_bonds;
        }
    }
    
    free_spit_context(ctx);
    return max_bond_interactions;
}

double calculateSpitRingParticipationIndexPerToken(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 1) {
        free_spit_context(ctx);
        return 0.0;
    }

    double max_ring_interactions = 0.0;
    int ring_feature_idx = 2; // is_ring_related feature
    
    // For each token, calculate total ring-related interactions
    for (int i = 0; i < ctx->N; ++i) {
        double sum_rings = 0.0;
        for (int j = 0; j < ctx->N; ++j) {
            sum_rings += get_spit_value(ctx, i, j, ring_feature_idx);
        }
        
        if (i == 0 || sum_rings > max_ring_interactions) {
            max_ring_interactions = sum_rings;
        }
    }
    
    free_spit_context(ctx);
    return max_ring_interactions;
}

double calculateSpitTokenEntropy(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double sum_entropies = 0.0;
    double* token_feature_values = (double*)malloc(ctx->N * sizeof(double));
    if (!token_feature_values) {
        free_spit_context(ctx);
        return 0.0;
    }
    
    // For each token and feature, compute entropy of its interactions
    for (int i = 0; i < ctx->N; ++i) {
        double token_entropy = 0.0;
        
        for (int k = 0; k < SPIT_DIMS; ++k) {
            // Collect all values for this token-feature combination
            for (int j = 0; j < ctx->N; ++j) {
                token_feature_values[j] = get_spit_value(ctx, i, j, k);
            }
            
            // Count frequencies of 0 and 1 (binary features)
            int zeros = 0, ones = 0;
            for (int j = 0; j < ctx->N; ++j) {
                if (fabs(token_feature_values[j]) < 1e-9) zeros++;
                else if (fabs(token_feature_values[j] - 1.0) < 1e-9) ones++;
            }
            
            double entropy = 0.0;
            if (zeros + ones > 0) {
                double p_zero = (double)zeros / (zeros + ones);
                double p_one = (double)ones / (zeros + ones);
                
                if (p_zero > 0) entropy -= p_zero * log2(p_zero);
                if (p_one > 0) entropy -= p_one * log2(p_one);
            }
            
            token_entropy += entropy;
        }
        
        // Average entropy across features
        token_entropy /= SPIT_DIMS;
        sum_entropies += token_entropy;
    }
    
    free(token_feature_values);
    double avg_token_entropy = sum_entropies / ctx->N;
    
    free_spit_context(ctx);
    return avg_token_entropy;
}

double calculateSpitPositionallyStableTokens(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double* variances = (double*)malloc(ctx->N * sizeof(double));
    if (!variances) {
        free_spit_context(ctx);
        return 0.0;
    }
    
    // Calculate variance of feature values for each token
    for (int i = 0; i < ctx->N; ++i) {
        double sum_variance = 0.0;
        
        for (int k = 0; k < SPIT_DIMS; ++k) {
            // Calculate mean
            double sum = 0.0;
            for (int j = 0; j < ctx->N; ++j) {
                sum += get_spit_value(ctx, i, j, k);
            }
            double mean = sum / ctx->N;
            
            // Calculate variance
            double sq_diff_sum = 0.0;
            for (int j = 0; j < ctx->N; ++j) {
                double diff = get_spit_value(ctx, i, j, k) - mean;
                sq_diff_sum += diff * diff;
            }
            double variance = sq_diff_sum / (ctx->N - 1);
            sum_variance += variance;
        }
        
        // Average variance across features
        variances[i] = sum_variance / SPIT_DIMS;
    }
    
    // Count tokens with variance below threshold (stable tokens)
    double threshold = 0.1; // Arbitrary threshold for "stable"
    int stable_token_count = 0;
    for (int i = 0; i < ctx->N; ++i) {
        if (variances[i] < threshold) {
            stable_token_count++;
        }
    }
    
    free(variances);
    double result = (double)stable_token_count;
    
    free_spit_context(ctx);
    return result;
}

double calculateSpitTokenClusterEntropy(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    // Simple clustering approach: group tokens by their type
    int type_counts[8] = {0}; // 8 token types in our enum
    
    for (int i = 0; i < ctx->N; ++i) {
        TokenType type = ctx->tokenized_smiles.tokens[i].type;
        if (type >= 0 && type < 8) {
            type_counts[type]++;
        }
    }
    
    // Calculate entropy of cluster sizes
    double entropy = 0.0;
    for (int t = 0; t < 8; ++t) {
        if (type_counts[t] > 0) {
            double p = (double)type_counts[t] / ctx->N;
            entropy -= p * log2(p);
        }
    }
    
    free_spit_context(ctx);
    return entropy;
}

double calculateSpitHighestInteractiveToken(const void* context, GetSmilesFunc getSmilesFunc) {
    // This is the same as calculateSpitMostConnectedTokenScore
    return calculateSpitMostConnectedTokenScore(context, getSmilesFunc);
}

double calculateSpitIntraClassInteractionRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double intra_class_sum = 0.0;
    double total_sum = 0.0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            if (i != j) { // Skip self-interactions
                double cell_sum = 0.0;
                for (int k = 0; k < SPIT_DIMS; ++k) {
                    cell_sum += get_spit_value(ctx, i, j, k);
                }
                
                total_sum += cell_sum;
                
                // Check if tokens are of the same type (intra-class)
                if (ctx->tokenized_smiles.tokens[i].type == ctx->tokenized_smiles.tokens[j].type) {
                    intra_class_sum += cell_sum;
                }
            }
        }
    }
    
    double result = (total_sum > 0) ? (intra_class_sum / total_sum) : 0.0;
    free_spit_context(ctx);
    return result;
}

double calculateSpitInterClassInteractionRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double inter_class_sum = 0.0;
    double total_sum = 0.0;
    
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            if (i != j) { // Skip self-interactions
                double cell_sum = 0.0;
                for (int k = 0; k < SPIT_DIMS; ++k) {
                    cell_sum += get_spit_value(ctx, i, j, k);
                }
                
                total_sum += cell_sum;
                
                // Check if tokens are of different types (inter-class)
                if (ctx->tokenized_smiles.tokens[i].type != ctx->tokenized_smiles.tokens[j].type) {
                    inter_class_sum += cell_sum;
                }
            }
        }
    }
    
    double result = (total_sum > 0) ? (inter_class_sum / total_sum) : 0.0;
    free_spit_context(ctx);
    return result;
}

double calculateSpitBondVsRingSeparationScore(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double bond_sum = 0.0;
    double ring_sum = 0.0;
    
    // Sum activations in bond feature plane (1) and ring feature plane (2)
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            bond_sum += get_spit_value(ctx, i, j, 1);
            ring_sum += get_spit_value(ctx, i, j, 2);
        }
    }
    
    // Calculate ratio of bond to ring activations
    double result = 0.0;
    if (ring_sum > 0) {
        result = bond_sum / ring_sum;
    } else if (bond_sum > 0) {
        result = 9999.0; // Large value to indicate high bond/ring ratio
    }
    
    free_spit_context(ctx);
    return result;
}

double calculateSpitCheckerboardPatternStrength(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double checkerboard_sum = 0.0;
    double total_cells = ctx->N * ctx->N;
    
    // Calculate checkerboard pattern: sum (-1)^(i+j) * value
    // For feature plane 0 (is_same_token)
    int feature_idx = 0;
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            double value = get_spit_value(ctx, i, j, feature_idx);
            double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
            checkerboard_sum += sign * value;
        }
    }
    
    // Normalize by total number of cells
    double result = fabs(checkerboard_sum) / total_cells;
    
    free_spit_context(ctx);
    return result;
}

double calculateSpitGridUniformity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 4) { // Need at least 4 tokens for 2x2 blocks
        free_spit_context(ctx);
        return 0.0;
    }

    // Use a 2x2 block size for simplicity
    int block_size = 2;
    int num_blocks_row = ctx->N / block_size;
    int num_blocks_col = ctx->N / block_size;
    int num_blocks = num_blocks_row * num_blocks_col;
    
    if (num_blocks == 0) {
        free_spit_context(ctx);
        return 0.0;
    }
    
    // Calculate sum for each block in feature plane 0
    int feature_idx = 0;
    double* block_sums = (double*)calloc(num_blocks, sizeof(double));
    if (!block_sums) {
        free_spit_context(ctx);
        return 0.0;
    }
    
    int block_idx = 0;
    for (int bi = 0; bi < num_blocks_row; ++bi) {
        for (int bj = 0; bj < num_blocks_col; ++bj) {
            double sum = 0.0;
            for (int i = bi * block_size; i < (bi + 1) * block_size && i < ctx->N; ++i) {
                for (int j = bj * block_size; j < (bj + 1) * block_size && j < ctx->N; ++j) {
                    sum += get_spit_value(ctx, i, j, feature_idx);
                }
            }
            block_sums[block_idx++] = sum;
        }
    }
    
    // Calculate standard deviation of block sums
    double mean_block_sum = 0.0;
    for (int i = 0; i < num_blocks; ++i) {
        mean_block_sum += block_sums[i];
    }
    mean_block_sum /= num_blocks;
    
    double sum_sq_diff = 0.0;
    for (int i = 0; i < num_blocks; ++i) {
        double diff = block_sums[i] - mean_block_sum;
        sum_sq_diff += diff * diff;
    }
    
    double std_dev = sqrt(sum_sq_diff / num_blocks);
    
    free(block_sums);
    free_spit_context(ctx);
    return std_dev;
}

double calculateSpitMax3x3InteractionPatchScore(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 3) { // Need at least 3x3
        free_spit_context(ctx);
        return 0.0;
    }

    int patch_size = 3;
    double max_patch_sum = 0.0;
    int feature_idx = 0; // Use first feature plane (is_same_token)
    
    // Slide a 3x3 window over the tensor and find max sum
    for (int i = 0; i <= ctx->N - patch_size; ++i) {
        for (int j = 0; j <= ctx->N - patch_size; ++j) {
            double patch_sum = 0.0;
            
            // Sum values in the 3x3 patch
            for (int pi = 0; pi < patch_size; ++pi) {
                for (int pj = 0; pj < patch_size; ++pj) {
                    patch_sum += get_spit_value(ctx, i + pi, j + pj, feature_idx);
                }
            }
            
            if ((i == 0 && j == 0) || patch_sum > max_patch_sum) {
                max_patch_sum = patch_sum;
            }
        }
    }
    
    free_spit_context(ctx);
    return max_patch_sum;
}

double calculateSpitLPatternMotifCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 3) { // Need at least 3x3 for L-patterns
        free_spit_context(ctx);
        return 0.0;
    }

    int feature_idx = 0; // Use first feature plane (is_same_token)
    int l_pattern_count = 0;
    
    // Look for L-patterns: three adjacent 1's forming an L shape
    for (int i = 0; i < ctx->N - 1; ++i) {
        for (int j = 0; j < ctx->N - 1; ++j) {
            // Check vertical L pattern
            if (get_spit_value(ctx, i, j, feature_idx) > 0.5 &&
                get_spit_value(ctx, i + 1, j, feature_idx) > 0.5 &&
                get_spit_value(ctx, i + 1, j + 1, feature_idx) > 0.5) {
                l_pattern_count++;
            }
            
            // Check horizontal L pattern
            if (get_spit_value(ctx, i, j, feature_idx) > 0.5 &&
                get_spit_value(ctx, i, j + 1, feature_idx) > 0.5 &&
                get_spit_value(ctx, i + 1, j + 1, feature_idx) > 0.5) {
                l_pattern_count++;
            }
        }
    }
    
    free_spit_context(ctx);
    return (double)l_pattern_count;
}

double calculateSpitSymmetryUnder90Rotation(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double sum_diff = 0.0;
    double total_comparisons = 0.0;
    int feature_idx = 0; // Use first feature plane (is_same_token)
    
    // Compare original with 90-degree rotation: spit[i][j][k] vs spit[N-1-j][i][k]
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            int rot_i = ctx->N - 1 - j;
            int rot_j = i;
            
            if (rot_i >= 0 && rot_i < ctx->N && rot_j >= 0 && rot_j < ctx->N) {
                double diff = fabs(get_spit_value(ctx, i, j, feature_idx) - 
                                  get_spit_value(ctx, rot_i, rot_j, feature_idx));
                sum_diff += diff;
                total_comparisons += 1.0;
            }
        }
    }
    
    // Normalize: 0 = perfect symmetry, 1 = maximum asymmetry
    double result = (total_comparisons > 0) ? (1.0 - sum_diff / total_comparisons) : 0.0;
    
    free_spit_context(ctx);
    return result;
}

double calculateSpitAntiDiagonalAsymmetryScore(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double sum_diff = 0.0;
    double total_comparisons = 0.0;
    int feature_idx = 0; // Use first feature plane
    
    // Compare pairs across the anti-diagonal: spit[i][j][k] vs spit[N-1-j][N-1-i][k]
    // but only when i+j != N-1 (not on the anti-diagonal)
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            if (i + j != ctx->N - 1) { // Skip anti-diagonal elements
                int anti_i = ctx->N - 1 - j;
                int anti_j = ctx->N - 1 - i;
                
                double diff = fabs(get_spit_value(ctx, i, j, feature_idx) - 
                                  get_spit_value(ctx, anti_i, anti_j, feature_idx));
                sum_diff += diff;
                total_comparisons += 1.0;
            }
        }
    }
    
    double result = (total_comparisons > 0) ? (sum_diff / total_comparisons) : 0.0;
    
    free_spit_context(ctx);
    return result;
}

double calculateSpitTopLeftToBottomRightImbalance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double top_left_sum = 0.0;
    double bottom_right_sum = 0.0;
    int midpoint = ctx->N / 2;
    
    // Calculate sum for top-left and bottom-right quadrants
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            double sum_features = 0.0;
            for (int k = 0; k < SPIT_DIMS; ++k) {
                sum_features += get_spit_value(ctx, i, j, k);
            }
            
            if (i < midpoint && j < midpoint) { // Top-left quadrant
                top_left_sum += sum_features;
            } else if (i >= midpoint && j >= midpoint) { // Bottom-right quadrant
                bottom_right_sum += sum_features;
            }
        }
    }
    
    // Calculate imbalance ratio
    double result = 0.0;
    if (top_left_sum + bottom_right_sum > 0) {
        result = (top_left_sum - bottom_right_sum) / (top_left_sum + bottom_right_sum);
    }
    
    free_spit_context(ctx);
    return result;
}

double calculateSpitEdgeTokenInteractionBias(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 3) { // Need at least 3x3 to have edge vs core
        free_spit_context(ctx);
        return 0.0;
    }

    double edge_sum = 0.0;
    double core_sum = 0.0;
    int edge_count = 0;
    int core_count = 0;
    
    // Calculate sum for edge vs core
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            double sum_features = 0.0;
            for (int k = 0; k < SPIT_DIMS; ++k) {
                sum_features += get_spit_value(ctx, i, j, k);
            }
            
            // Edge elements are on the perimeter of the matrix
            if (i == 0 || i == ctx->N - 1 || j == 0 || j == ctx->N - 1) {
                edge_sum += sum_features;
                edge_count++;
            } else { // Core elements
                core_sum += sum_features;
                core_count++;
            }
        }
    }
    
    double edge_avg = (edge_count > 0) ? (edge_sum / edge_count) : 0.0;
    double core_avg = (core_count > 0) ? (core_sum / core_count) : 0.0;
    
    // Calculate bias: positive means edge > core, negative means core > edge
    double result = 0.0;
    if (edge_avg + core_avg > 0) {
        result = (edge_avg - core_avg) / (edge_avg + core_avg);
    }
    
    free_spit_context(ctx);
    return result;
}

double calculateSpitCenterOfMatrixActivationShift(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    SpitContext* ctx = create_spit_context(smiles);
    if (!ctx || !ctx->tensor_valid || ctx->N < 2) {
        free_spit_context(ctx);
        return 0.0;
    }

    double weighted_i_sum = 0.0;
    double weighted_j_sum = 0.0;
    double total_weight = 0.0;
    int feature_idx = 0; // Use first feature plane (is_same_token)
    
    // Calculate center of mass
    for (int i = 0; i < ctx->N; ++i) {
        for (int j = 0; j < ctx->N; ++j) {
            double value = get_spit_value(ctx, i, j, feature_idx);
            weighted_i_sum += i * value;
            weighted_j_sum += j * value;
            total_weight += value;
        }
    }
    
    if (total_weight < 1e-9) {
        free_spit_context(ctx);
        return 0.0;
    }
    
    double center_i = weighted_i_sum / total_weight;
    double center_j = weighted_j_sum / total_weight;
    
    // Expected center is (N-1)/2, (N-1)/2
    double expected_center = (ctx->N - 1) / 2.0;
    
    // Calculate distance from expected center
    double distance = sqrt(pow(center_i - expected_center, 2) + 
                          pow(center_j - expected_center, 2));
    
    // Normalize by matrix size
    double result = distance / ctx->N;
    
    free_spit_context(ctx);
    return result;
}
