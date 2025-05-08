// src/descriptors/regexcounts.c
#include "../cregistry.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <regex.h>
#include <math.h>
#include <pthread.h>
#include <limits.h>

// Increase the regex cache size
#define MAX_CACHED_REGEX 500

// Add timestamp to track recency of use
typedef struct {
    const char* pattern;
    regex_t regex;
    int initialized; // 0 = not initialized, 1 = initialized_ok, -1 = compilation_failed
    unsigned long last_used;  // Timestamp to implement LRU replacement
} CachedRegex;

// Global regex cache
static CachedRegex regex_cache[MAX_CACHED_REGEX];
static int regex_cache_count = 0;
static unsigned long access_counter = 0;  // Simple counter to track usage

// Add a mutex for thread-safe access to the regex cache
static pthread_mutex_t regex_cache_mutex = PTHREAD_MUTEX_INITIALIZER;

// Use thread-local caches to reduce lock contention
#define LOCAL_CACHE_SIZE 64 // Increased size
__thread regex_t local_regex_cache[LOCAL_CACHE_SIZE];
__thread const char* local_patterns[LOCAL_CACHE_SIZE];
__thread int local_initialized[LOCAL_CACHE_SIZE]; // Ensure this uses the new LOCAL_CACHE_SIZE if it was separate
__thread int local_cache_count = 0;

// Helper function to find least recently used cache entry
static int findLRUCacheEntry() {
    int lru_index = 0;
    unsigned long oldest_time = ULONG_MAX;
    
    for (int i = 0; i < regex_cache_count; i++) {
        if (regex_cache[i].last_used < oldest_time) {
            oldest_time = regex_cache[i].last_used;
            lru_index = i;
        }
    }
    
    return lru_index;
}

// Helper function to get compiled regex - with thread-local cache
regex_t* getCompiledRegex(const char* pattern) {
    if (!pattern) return NULL;
    
    // 1. First check thread-local cache (no lock needed)
    for (int i = 0; i < local_cache_count; i++) {
        if (local_patterns[i] && strcmp(local_patterns[i], pattern) == 0) {
            if (local_initialized[i] == 1) { // Compiled successfully
                return &local_regex_cache[i];
            } else if (local_initialized[i] == -1) { // Known compilation failure
                return NULL; 
            }
        }
    }
    
    // 2. If we have room in the local cache, try to compile
    if (local_cache_count < LOCAL_CACHE_SIZE) {
        int idx = local_cache_count++;
        
        char* pattern_copy = strdup(pattern);
        if (!pattern_copy) {
            // fprintf(stderr, "Failed to allocate memory for regex pattern in local cache\n"); // Less verbose
            return NULL;
        }
        
        local_patterns[idx] = pattern_copy;
        
        int ret = regcomp(&local_regex_cache[idx], pattern, REG_EXTENDED);
        if (ret != 0) {
            local_initialized[idx] = -1; // Mark as compilation failed
            return NULL;
        }
        
        local_initialized[idx] = 1; // Mark as compiled successfully
        return &local_regex_cache[idx];
    }
    
    // 3. Only if local cache is full, check global cache (with lock)
    regex_t* result = NULL;
    
    pthread_mutex_lock(&regex_cache_mutex);
    access_counter++;
    
    // Check global cache
    for (int i = 0; i < regex_cache_count; i++) {
        if (regex_cache[i].pattern && strcmp(regex_cache[i].pattern, pattern) == 0) {
            if (regex_cache[i].initialized == 1) { // Compiled successfully
                regex_cache[i].last_used = access_counter;
                result = &regex_cache[i].regex;
                pthread_mutex_unlock(&regex_cache_mutex);
                return result;
            } else if (regex_cache[i].initialized == -1) { // Known compilation failure
                pthread_mutex_unlock(&regex_cache_mutex);
                return NULL;
            }
        }
    }
    
    // If not in cache (or present but not successfully compiled), attempt to add/replace
    int idx;
    if (regex_cache_count < MAX_CACHED_REGEX) {
        idx = regex_cache_count++;
    } else {
        idx = findLRUCacheEntry();
        
        if (regex_cache[idx].initialized == 1) { // if it was a successfully compiled regex
            regfree(&regex_cache[idx].regex);
        }
        if (regex_cache[idx].pattern) {
            free((void*)regex_cache[idx].pattern);
            regex_cache[idx].pattern = NULL;
        }
        regex_cache[idx].initialized = 0; // Reset state
    }
    
    char* pattern_copy_global = strdup(pattern);
    if (!pattern_copy_global) {
        // fprintf(stderr, "Failed to allocate memory for regex pattern in global cache\n"); // Less verbose
        pthread_mutex_unlock(&regex_cache_mutex);
        return NULL;
    }
    
    regex_cache[idx].pattern = pattern_copy_global;
    regex_cache[idx].last_used = access_counter;
    
    int ret_global = regcomp(&regex_cache[idx].regex, pattern, REG_EXTENDED);
    if (ret_global != 0) {
        // Make error printing conditional or less frequent to reduce overhead.
        // For example, only print the error for a given pattern once.
        // This simple version just removes it for now to test impact.
        // A more sophisticated approach might involve a flag per pattern to track if error was printed.
        // char error_msg[100];
        // regerror(ret_global, &regex_cache[idx].regex, error_msg, sizeof(error_msg));
        // fprintf(stderr, "Failed to compile regex '%s' for global cache: %s\n", pattern, error_msg);
        regex_cache[idx].initialized = -1; 
        pthread_mutex_unlock(&regex_cache_mutex);
        return NULL;
    }
    
    regex_cache[idx].initialized = 1; 
    result = &regex_cache[idx].regex;
    
    pthread_mutex_unlock(&regex_cache_mutex);
    return result;
}

// Free the thread-local cache
void freeThreadLocalRegexCache() {
    for (int i = 0; i < local_cache_count; i++) {
        if (local_initialized[i]) {
            regfree(&local_regex_cache[i]);
            local_initialized[i] = 0;
        }
        if (local_patterns[i]) {
            free((void*)local_patterns[i]);
            local_patterns[i] = NULL;
        }
    }
    local_cache_count = 0;
}

// Update the main free function to also free thread-local cache
void freeRegexCache() {
    freeThreadLocalRegexCache();
    
    pthread_mutex_lock(&regex_cache_mutex);
    
    for (int i = 0; i < regex_cache_count; i++) {
        if (regex_cache[i].initialized) {
            regfree(&regex_cache[i].regex);
            regex_cache[i].initialized = 0;
        }
        if (regex_cache[i].pattern) {
            free((void*)regex_cache[i].pattern);
            regex_cache[i].pattern = NULL;
        }
    }
    regex_cache_count = 0;
    
    pthread_mutex_unlock(&regex_cache_mutex);
}

// Helper functions for all descriptor implementations that need to be added

// Helper function to count regex matches in SMILES
int countSubstr(const char* str, const char* sub);  // Forward declaration

int countMatches(const char* s, const char* pattern) {
    if (!s || !pattern) return 0;

    // Direct substring checks for common simple patterns
    if (strcmp(pattern, "CCC") == 0) return countSubstr(s, pattern);
    if (strcmp(pattern, "C#C") == 0 || strcmp(pattern, "C#N") == 0 || strcmp(pattern, "N#N") == 0) return countSubstr(s, pattern);
    if (strcmp(pattern, "NC(=O)N") == 0) return countSubstr(s, pattern);
    if (strcmp(pattern, "NC(=S)N") == 0) return countSubstr(s, pattern);
    if (strcmp(pattern, "CF3") == 0 || strcmp(pattern, "C(F)(F)F") == 0 || strcmp(pattern, "[CF3]") == 0) return countSubstr(s, pattern);
    
    regmatch_t pmatch[1];
    int count = 0;
    size_t offset = 0;
    size_t len = strlen(s);
    int status;
    
    regex_t* regex_ptr = getCompiledRegex(pattern); // Rely solely on this for compiled regex
    
    if (!regex_ptr) {
        // Pattern either failed to compile and was cached as such, or some other error.
        // fprintf(stderr, "Skipping regex matching for pattern (compilation failed or not found): %s\n", pattern); // Optional: for debugging
        return 0; 
    }
    
    // Find all matches using the cached regex_ptr
    while (offset < len) {
        status = regexec(regex_ptr, s + offset, 1, pmatch, 0);
        if (status != 0) break; // REG_NOMATCH or error
        
        if (pmatch[0].rm_eo == 0 && pmatch[0].rm_so == 0) { // Empty match, can cause infinite loop if not handled
             // Advance offset by 1 if the match is empty and at the current offset,
             // unless we are at the end of the string.
            if (offset < len) {
                offset++;
            } else {
                break; // End of string
            }
            // If the regex can match an empty string, it might be better to adjust the regex itself.
            // For now, just ensure progress.
            if (pmatch[0].rm_eo == pmatch[0].rm_so) continue; // Re-evaluate loop condition
        } else if (pmatch[0].rm_eo == pmatch[0].rm_so) { // Also an empty match.
             // This case can happen if the regex matches an empty string not at the current offset.
             // To avoid infinite loops on empty matches, we must advance.
             // A common strategy is to advance by 1 if an empty match occurs.
            offset += (pmatch[0].rm_so + 1 > offset) ? pmatch[0].rm_so + 1 - offset : 1;
            // If it was a zero-length match, we should still count it.
             count++; 
             continue;
        }


        count++;
        offset += pmatch[0].rm_eo; // Move to the end of the last match
    }
    
    // Do NOT free regex_ptr here, it's managed by the cache.
    
    return count;
}

// Helper function to count characters in a string
int countChar(const char* s, char c) {
    int count = 0;
    while (*s) {
        if (*s == c) count++;
        s++;
    }
    return count;
}

// Helper function to count digits in a string
int countDigits(const char* s) {
    int count = 0;
    while (*s) {
        if (isdigit((unsigned char)*s)) count++;
        s++;
    }
    return count;
}

// Helper function to count non-overlapping substrings
int countSubstr(const char* text, const char* pattern) {
    if (!text || !pattern || !pattern[0]) return 0;
    int count = 0;
    const char* pos = text;
    size_t pattern_len = strlen(pattern);
    while ((pos = strstr(pos, pattern)) != NULL) {
        ++count;
        pos += pattern_len;
    }
    return count;
}

// Count N not preceded by C
int countNNotPrecededByC(const char* s) {
    int count = 0;
    for (size_t i = 0; s[i] != '\0'; ++i) {
        if (s[i] == 'N' && (i == 0 || s[i-1] != 'C')) {
            ++count;
        }
    }
    return count;
}

// Simple check for positive and negative charges in a string
// Returns 1 if both are present (zwitterion), 0 otherwise
int hasZwitterion(const char* s) {
    int hasPos = 0;
    int hasNeg = 0;
    
    while (*s) {
        if (*s == '+') hasPos = 1;
        else if (*s == '-') hasNeg = 1;
        s++;
    }
    
    return (hasPos && hasNeg) ? 1 : 0;
}

// Carbon count
double calculateSmilesCarbonCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'C') + countChar(smiles, 'c');
}

// Nitrogen count
double calculateSmilesNitrogenCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'N') + countChar(smiles, 'n');
}

// Oxygen count
double calculateSmilesOxygenCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'O') + countChar(smiles, 'o');
}

// Sulfur count
double calculateSmilesSulfurCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'S') + countChar(smiles, 's');
}

// Phosphorus count
double calculateSmilesPhosphorusCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'P') + countChar(smiles, 'p');
}

// Halogen count
double calculateSmilesHalogenCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    count += countChar(smiles, 'F');
    count += countChar(smiles, 'I');
    count += countChar(smiles, 'f');
    count += countChar(smiles, 'i');
    count += countSubstr(smiles, "Cl");
    count += countSubstr(smiles, "Br");
    count += countSubstr(smiles, "cl");
    count += countSubstr(smiles, "br");
    return count;
}

// Complex Ring Closure Count
double calculateSmilesComplexRingClosureCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    for (size_t i = 0; i < strlen(smiles); i++) {
        // Check for ring closure with value greater than 9
        if (smiles[i] == '%' && i + 2 < strlen(smiles) && isdigit(smiles[i+1]) && isdigit(smiles[i+2])) {
            count++;
        }
    }
    return count;
}

// Conjugated Systems Count
double calculateSmilesConjugatedSystemCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    char* pos = (char*)smiles;
    // Look for alternating single and double bonds in aromatic systems
    while ((pos = strstr(pos, "c=c-c=c")) != NULL) {
        count++;
        pos++;
    }
    
    pos = (char*)smiles;
    while ((pos = strstr(pos, "C=C-C=C")) != NULL) {
        count++;
        pos++;
    }
    
    // Count other conjugated patterns
    pos = (char*)smiles;
    while ((pos = strstr(pos, "C=C-C=O")) != NULL) {
        count++;
        pos++;
    }
    
    return count;
}

// Cyclopropyl Count
double calculateSmilesCyclopropylCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    char* pos = (char*)smiles;
    // Look for cyclopropyl patterns
    while ((pos = strstr(pos, "C1CC1")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// Cyclobutyl Count
double calculateSmilesCyclobutylCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    char* pos = (char*)smiles;
    // Look for cyclobutyl patterns
    while ((pos = strstr(pos, "C1CCC1")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// Cyclopentyl Count
double calculateSmilesCyclopentylCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    char* pos = (char*)smiles;
    // Look for cyclopentyl patterns
    while ((pos = strstr(pos, "C1CCCC1")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// Cyclohexyl Count
double calculateSmilesCyclohexylCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    char* pos = (char*)smiles;
    // Look for cyclohexyl patterns
    while ((pos = strstr(pos, "C1CCCCC1")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// Dimethyl Count
double calculateSmilesDimethylCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    char* pos = (char*)smiles;
    // Look for dimethyl patterns (C(C)(C))
    while ((pos = strstr(pos, "C(C)(C)")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// Double Bond Ring Junction Count
double calculateSmilesDoubleBondRingJunctionCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // This is a simplification - would need more complex parsing for accurate detection
    char* pos = (char*)smiles;
    int count = 0;
    while ((pos = strstr(pos, "C=C1")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// Ether Count
double calculateSmilesEtherCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    // Count COC patterns but avoid counting alcohols (OH)
    for (size_t i = 1; i < strlen(smiles) - 1; i++) {
        if (smiles[i] == 'O' && 
            (isalpha(smiles[i-1]) || smiles[i-1] == ']') && 
            (isalpha(smiles[i+1]) || smiles[i+1] == '[') &&
            smiles[i+1] != 'H') {
            count++;
        }
    }
    return count;
}

// Thiol Count
double calculateSmilesThiolCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for thiol groups
    const char* pattern = "S[H]|\\[SH\\]|\\[S;H1\\]|\\[SH-\\]|\\[S-\\]|S\\(H\\)";
    return countMatches(smiles, pattern);
}

// Thioether Count
double calculateSmilesThioetherCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for thioether groups
    const char* pattern = "CSC|\\[CH2\\]S\\[CH3\\]|\\[CH2\\]S\\[CH2\\]|S\\(C\\)C";
    return countMatches(smiles, pattern);
}

// Carbonyl Heteroatom Count
double calculateSmilesCarbonylHeteroatomCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for carbonyl adjacent to heteroatom
    const char* pattern = "C\\(=O\\)[NOS]";
    return countMatches(smiles, pattern);
}

// Aromatic Carbon Heteroatom Count
double calculateSmilesAromaticCarbonHeteroatomCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for aromatic carbon adjacent to heteroatom
    const char* pattern = "c[nospb]|[nospb]c";
    return countMatches(smiles, pattern);
}

// Carbon Phosphorus Bond Count
double calculateSmilesCarbonPhosphorusBondCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for carbon-phosphorus bonds
    const char* pattern = "CP|PC|\\[C\\]P|C\\[P\\]";
    return countMatches(smiles, pattern);
}

// Carboxylate Count
double calculateSmilesCarboxylateCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for carboxylate anions
    const char* pattern = "C\\(=O\\)\\[O-\\]";
    return countMatches(smiles, pattern);
}

// Branched Amide Count
double calculateSmilesBranchedAmideCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for branched amides
    const char* pattern = "C\\(=O\\)N\\([^)]*C[^)]*\\)";
    return countMatches(smiles, pattern);
}

// Heteroatom Sequence Count
double calculateSmilesHeteroatomSequenceCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for heteroatom sequences
    const char* pattern = "[CNO][CNO]|[CNSO][CNSO]|[CNP][CNP]|[CNOFPS][CNOFPS]";
    return countMatches(smiles, pattern);
}

// Trifluoromethyl Count
double calculateSmilesTrifluoromethylCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for trifluoromethyl groups
    const char* pattern = "CF3|C\\(F\\)\\(F\\)F|\\[CF3\\]";
    return countMatches(smiles, pattern);
}

// Guanidine Basic Count
double calculateSmilesGuanidineBasicCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for guanidine groups
    const char* pattern = "N=C\\(N\\)N|N=C\\(\\[NH2\\]\\)\\[NH2\\]";
    return countMatches(smiles, pattern);
}

// Imidazolium Count
double calculateSmilesImidazoliumCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for imidazolium cations
    const char* pattern = "\\[n[H]?\\+\\]";
    return countMatches(smiles, pattern);
}

// Pyridinium Count
double calculateSmilesPyridiniumCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for pyridinium cations
    const char* pattern = "\\[n\\+\\]";
    return countMatches(smiles, pattern);
}

// Quaternary Ammonium Basic Count
double calculateSmilesQuaternaryAmmoniumBasicCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for quaternary ammonium ions
    const char* pattern = "\\[N\\+\\]\\([^H)][^)]*\\)\\([^)]*\\)\\([^)]*\\)\\([^)]*\\)|\\[N;H0;\\+;X4\\]";
    return countMatches(smiles, pattern);
}

// Ionizable Nitrogen Count
double calculateSmilesIonizableNitrogenCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for ionizable nitrogens
    const char* pattern = "\\[[nN].*\\+\\]|\\[NH2\\]|\\[NH\\]|\\[N-\\]|\\[NH-\\]";
    return countMatches(smiles, pattern);
}

// Ionizable Oxygen Count
double calculateSmilesIonizableOxygenCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for ionizable oxygens
    const char* pattern = "\\[O-\\]|\\[OH\\]";
    return countMatches(smiles, pattern);
}

// Triple Bond Count
double calculateSmilesTripleBondCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    int count = 0;
    count += countChar(smiles, '#');
    count += countSubstr(smiles, "C#C");
    count += countSubstr(smiles, "C#N");
    count += countSubstr(smiles, "N#N");
    return count;
}

// Phenol Count
double calculateSmilesPhenolCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for phenol groups
    const char* pattern = "c1[c]cccc1O|Oc1[c]cccc1";
    return countMatches(smiles, pattern);
}

// Tertiary Carbon Count
double calculateSmilesTertiaryCarbonCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    
    // Custom counter for tertiary carbon pattern
    int count = 0;
    for (size_t i = 0; i < strlen(smiles); i++) {
        if (smiles[i] == 'C' && 
            i+1 < strlen(smiles) && smiles[i+1] == '(' &&
            strchr(smiles+i+1, ')') != NULL) {
            
            const char* str = smiles + i + 1;
            int parenthesis_depth = 0;
            int branch_count = 0;
            
            while (*str && branch_count < 3) {
                if (*str == '(') {
                    if (parenthesis_depth == 0) branch_count++;
                    parenthesis_depth++;
                }
                else if (*str == ')') {
                    parenthesis_depth--;
                }
                str++;
            }
            
            if (branch_count >= 3) count++;
        }
    }
    return count;
}

// N-Alkyl Count
double calculateSmilesNAlkylCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for N-alkyl patterns
    const char* pattern = "NC|\\[NH2\\]C|\\[NH\\]C";
    return countMatches(smiles, pattern);
}

// Tertiary Amine Basic Count
double calculateSmilesTertiaryAmineBasicCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for tertiary amines
    const char* pattern = "N\\([^)]+\\)\\([^)]+\\)\\([^)]+\\)";
    return countMatches(smiles, pattern);
}

// Imidazole Basic Count
double calculateSmilesImidazoleBasicCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for imidazole rings
    const char* pattern = "c1nccn1|c1ncnc1";
    return countMatches(smiles, pattern);
}

// Ionizable Phenol Count
double calculateSmilesIonizablePhenolCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for ionizable phenols
    const char* pattern = "c1[c]cccc1\\[O-\\]|c1[c]cccc1\\[OH\\]";
    return countMatches(smiles, pattern);
}

// Ionizable Phosphorus Count
double calculateSmilesIonizablePhosphorusCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for ionizable phosphorus
    const char* pattern = "\\[P\\+\\]|\\[P-\\]|P\\(=O\\)\\(O[H]?\\)";
    return countMatches(smiles, pattern);
}

// Ionizable Pyridine Count
double calculateSmilesIonizablePyridineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for ionizable pyridines
    const char* pattern = "\\[n\\+\\]|n[H]?\\+";
    return countMatches(smiles, pattern);
}

// Ionizable Amine Count
double calculateSmilesIonizableAmineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for ionizable amines
    const char* pattern = "\\[NH3\\+\\]|\\[NH2\\+\\]|\\[NH\\+\\]|\\[N\\+\\]";
    return countMatches(smiles, pattern);
}

// Ionizable Imidazole Count
double calculateSmilesIonizableImidazoleCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for ionizable imidazoles
    const char* pattern = "c1ncc\\[nH\\+\\]1|c1ncc\\[n\\+\\]1";
    return countMatches(smiles, pattern);
}

// Isocyanide Count
double calculateSmilesIsocyanideCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for isocyanides
    const char* pattern = "N#C|C#N";
    return countMatches(smiles, pattern);
}

// Sulfonamide Count
double calculateSmilesIonizableSulfonamideCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for ionizable sulfonamides
    const char* pattern = "S\\(=O\\)\\(=O\\)\\(\\[N-\\]|N\\[H\\]\\)";
    return countMatches(smiles, pattern);
}

// Sulfonate Count
double calculateSmilesSulfonateCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for sulfonates
    const char* pattern = "S\\(=O\\)\\(=O\\)\\[O-\\]|OS\\(=O\\)\\(=O\\)O";
    return countMatches(smiles, pattern);
}

// Sulfonic Acid Count
double calculateSmilesSulfonicAcidCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for sulfonic acids
    const char* pattern = "S\\(=O\\)\\(=O\\)O|S\\(=O\\)\\(=O\\)OH";
    return countMatches(smiles, pattern);
}

// Phosphate Anion Count
double calculateSmilesPhosphateAnionCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for phosphate anions
    const char* pattern = "P\\(=O\\)\\(\\[O-\\]\\)\\(O\\)O|P\\(=O\\)\\(O\\)\\(O\\)\\[O-\\]";
    return countMatches(smiles, pattern);
}

// Fluorinated Fragment Count
double calculateSmilesFluorinatedFragmentCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for fluorinated fragments
    const char* pattern = "[Cc](F|\\[F\\])+";
    return countMatches(smiles, pattern);
}

// Zwitterion Count
double calculateSmilesZwitterionCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Check for presence of both positive and negative charges
    return hasZwitterion(smiles) ? 1.0 : 0.0;
}

// Ionizable Sulfur Count
double calculateSmilesIonizableSulfurCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for ionizable sulfur
    const char* pattern = "\\[S\\+\\]|\\[S-\\]|\\[SH\\]|S\\[H\\]";
    return countMatches(smiles, pattern);
}

// General Ammonium Count
double calculateSmilesGeneralAmmoniumCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Pattern for general ammonium ions
    const char* pattern = "\\[NH4\\+\\]|\\[N\\+\\]|\\[N\\+\\]\\([^)]*\\)";
    return countMatches(smiles, pattern);
}

// Additional necessary functions from the C++ implementation
double calculateSmilesBoronCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'B');
}

double calculateSmilesFluorineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'F');
}

double calculateSmilesChlorineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countSubstr(smiles, "Cl");
}

double calculateSmilesBromineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countSubstr(smiles, "Br");
}

double calculateSmilesIodineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'I');
}

double calculateSmilesAromaticCarbonCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'c');
}

double calculateSmilesAromaticNitrogenCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'n');
}

double calculateSmilesAromaticOxygenCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 'o');
}

double calculateSmilesAromaticSulfurCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, 's');
}

double calculateSmilesAromaticHalogenCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    count += countChar(smiles, 'f');
    count += countChar(smiles, 'i');
    count += countSubstr(smiles, "cl");
    count += countSubstr(smiles, "br");
    return count;
}

double calculateSmilesSingleBondCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, '-');
}

double calculateSmilesDoubleBondCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, '=');
}

double calculateSmilesAromaticBondCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "([cnospfib])";
    int count = countMatches(smiles, pattern);
    count += countChar(smiles, ':');
    return count;
}

double calculateSmilesRingClosureCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countDigits(smiles);
}

double calculateSmilesBranchOpenCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, '(');
}

double calculateSmilesBranchCloseCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, ')');
}

double calculateSmilesBracketAtomCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, '[');
}

double calculateSmilesChiralCenterCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, '@');
}

double calculateSmilesPositiveChargeCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, '+');
}

double calculateSmilesNegativeChargeCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, '-');
}

double calculateSmilesStereoForwardCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, '/');
}

double calculateSmilesStereoBackwardCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, '\\');
}

double calculateSmilesDotCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return countChar(smiles, '.');
}

// SmilesComplexRingClosureCount
double smiles_complex_ring_closure_count(const char* smiles, GetSmilesFunc get_smiles) {
    int count = 0;
    for (size_t i = 0; i < strlen(smiles); i++) {
        // Check for ring closure with value greater than 9
        if (smiles[i] == '%' && i + 2 < strlen(smiles) && isdigit(smiles[i+1]) && isdigit(smiles[i+2])) {
            count++;
        }
    }
    return count;
}

// SmilesConjugatedSystemCount
double smiles_conjugated_system_count(const char* smiles, GetSmilesFunc get_smiles) {
    int count = 0;
    char* pos = (char*)smiles;
    // Look for alternating single and double bonds in aromatic systems
    while ((pos = strstr(pos, "c=c-c=c")) != NULL) {
        count++;
        pos++;
    }
    
    pos = (char*)smiles;
    while ((pos = strstr(pos, "C=C-C=C")) != NULL) {
        count++;
        pos++;
    }
    
    // Count other conjugated patterns
    pos = (char*)smiles;
    while ((pos = strstr(pos, "C=C-C=O")) != NULL) {
        count++;
        pos++;
    }
    
    return count;
}

// SmilesCyclobutylCount
double smiles_cyclobutyl_count(const char* smiles, GetSmilesFunc get_smiles) {
    int count = 0;
    char* pos = (char*)smiles;
    // Look for cyclobutyl patterns
    while ((pos = strstr(pos, "C1CCC1")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// SmilesCyclohexylCount
double smiles_cyclohexyl_count(const char* smiles, GetSmilesFunc get_smiles) {
    int count = 0;
    char* pos = (char*)smiles;
    // Look for cyclohexyl patterns
    while ((pos = strstr(pos, "C1CCCCC1")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// SmilesCyclopentylCount
double smiles_cyclopentyl_count(const char* smiles, GetSmilesFunc get_smiles) {
    int count = 0;
    char* pos = (char*)smiles;
    // Look for cyclopentyl patterns
    while ((pos = strstr(pos, "C1CCCC1")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// SmilesCyclopropylCount
double smiles_cyclopropyl_count(const char* smiles, GetSmilesFunc get_smiles) {
    int count = 0;
    char* pos = (char*)smiles;
    // Look for cyclopropyl patterns
    while ((pos = strstr(pos, "C1CC1")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// SmilesDimethylCount
double smiles_dimethyl_count(const char* smiles, GetSmilesFunc get_smiles) {
    int count = 0;
    char* pos = (char*)smiles;
    // Look for dimethyl patterns (C(C)(C))
    while ((pos = strstr(pos, "C(C)(C)")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// SmilesDoubleBondRingJunctionCount
double smiles_double_bond_ring_junction_count(const char* smiles, GetSmilesFunc get_smiles) {
    int count = 0;
    // This is a simplification - would need more complex parsing for accurate detection
    char* pos = (char*)smiles;
    while ((pos = strstr(pos, "C=C1")) != NULL) {
        count++;
        pos++;
    }
    return count;
}

// SmilesEsterCount
double calculateSmilesEsterCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int count = 0;
    char* pos = (char*)smiles;
    
    // Look for ester patterns C(=O)O
    while ((pos = strstr(pos, "C(=O)O")) != NULL) {
        count++;
        pos++;
    }
    
    pos = (char*)smiles;
    while ((pos = strstr(pos, "OC(=O)")) != NULL) {
        count++;
        pos++;
    }
    
    return count;
}

// SmilesEtherCount
double smiles_ether_count(const char* smiles, GetSmilesFunc get_smiles) {
    int count = 0;
    // Count COC patterns but avoid counting alcohols (OH)
    for (size_t i = 1; i < strlen(smiles) - 1; i++) {
        if (smiles[i] == 'O' && 
            (isalpha(smiles[i-1]) || smiles[i-1] == ']') && 
            (isalpha(smiles[i+1]) || smiles[i+1] == '[') &&
            smiles[i+1] != 'H') {
            count++;
        }
    }
    return count;
}

// Ketone Count
double calculateSmilesKetoneCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "C\\(=O\\)C";
    return countMatches(smiles, pattern);
}

// Primary Amine Count
double calculateSmilesPrimaryAmineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "\\[NH2\\]|N\\([H]\\)\\([H]\\)";
    return countMatches(smiles, pattern);
}

// Secondary Amine Count
double calculateSmilesSecondaryAmineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "\\[NH\\]|N\\([H]\\)\\([^H][^)]*\\)";
    return countMatches(smiles, pattern);
}

// Amide Count
double calculateSmilesAmideCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "C\\(=O\\)N|NC\\(=O\\)";
    return countMatches(smiles, pattern);
}

// Nitro Group Count
double calculateSmilesNitroCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "\\[N\\+\\]\\(=O\\)\\[O-\\]|N\\(=O\\)=O";
    return countMatches(smiles, pattern);
}

// Alcohol Count
double calculateSmilesAlcoholCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "CO|\\[OH\\]|O\\[H\\]";
    return countMatches(smiles, pattern);
}

// Alkyne Count
double calculateSmilesAlkyneCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "C#C";
    return countMatches(smiles, pattern);
}

// Nitrile Count
double calculateSmilesNitrileCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "C#N";
    return countMatches(smiles, pattern);
}

// Isocyanate Count
double calculateSmilesIsocyanateCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "N=C=O";
    return countMatches(smiles, pattern);
}

// Thiocyanate Count
double calculateSmilesThiocyanateCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "S=C=N|N=C=S";
    return countMatches(smiles, pattern);
}

// Aldehyde Count
double calculateSmilesAldehydeCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "C\\(=O\\)\\[H\\]|C=O\\[H\\]"
    // Augmented: Covers [CH]=O, [cH]=O, C(=O)[H], c(=O)[H], and explicit H1 variants.
    const char* pattern = "[CH1]=O|[cH1]=O|C\\(=O\\)\\[H1?\\]|c\\(=O\\)\\[H1?\\]|C\\[H1?\\]=O|c\\[H1?\\]=O";
    return countMatches(smiles, pattern);
}

// Anhydride Count
double calculateSmilesAnhydrideCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "C\\(=O\\)OC\\(=O\\)"
    // Augmented: Include aromatic carbons
    const char* pattern = "[Cc]\\(=O\\)O[Cc]\\(=O\\)";
    return countMatches(smiles, pattern);
}


// Phenyl Count
double calculateSmilesPhenylCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "c1ccccc1";
    return countMatches(smiles, pattern);
}

// Benzyl Count
double calculateSmilesBenzylCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "c1ccccc1C";
    return countMatches(smiles, pattern);
}

// Pyridine Count
double calculateSmilesPyridineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "c1ccncc1|c1cnccc1|c1cncc1";
    return countMatches(smiles, pattern);
}

// Pyrimidine Count
double calculateSmilesPyrimidineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "c1cncnc1|c1ncncn1"
    // Augmented: Add common permutations of N within the 6-membered aromatic ring.
    const char* pattern = "c1cncnc1|c1ncncc1|c1nccnc1|c1ncncn1|c1ncnc[nH]1|c1nc[nH]cn1"; // Basic pyrimidines
    return countMatches(smiles, pattern);
}

// Pyrrole Count
double calculateSmilesPyrroleCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "c1ccn[cH]1|c1cc[nH]c1"
    // Augmented: Include more explicit hydrogen and attachment point variations.
    const char* pattern = "c1cc[nH]c1|c1c[nH]cc1|c1[nH]ccc1|\\[nH\\]1cccc1|c1ccc[nH]1";
    return countMatches(smiles, pattern);
}

// Furan Count
double calculateSmilesFuranCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "c1ccoc1|c1cocc1";
    return countMatches(smiles, pattern);
}

// Thiophene Count
double calculateSmilesThiopheneCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "c1ccsc1|c1sccc1";
    return countMatches(smiles, pattern);
}

// Piperidine Count
double calculateSmilesPiperidineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "C1CCNCC1|N1CCCCC1";
    return countMatches(smiles, pattern);
}

// Piperazine Count
double calculateSmilesPiperazineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "C1CNCCN1|N1CCNCC1";
    return countMatches(smiles, pattern);
}

// Morpholine Count
double calculateSmilesMorpholineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "C1COCCN1|N1CCOCC1|O1CCNCC1";
    return countMatches(smiles, pattern);
}

// Thiomorpholine Count
double calculateSmilesThiomorpholineCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "C1CSCN1|N1CSCC1|S1CCNCC1"
    // Augmented: Add variants with explicit Hydrogens if needed, ensure all heteroatom positions.
    // Saturated 6-membered ring with one S and one N.
    const char* pattern = "C1CSCN1|C1NCSC1|C1SCNC1|C1CNSC1|S1CCNC1|N1CCSC1|S1CNCC1|N1CSCC1"; // Covers more permutations
    return countMatches(smiles, pattern);
}

// Carbamate Count
double calculateSmilesCarbamatCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "NC\\(=O\\)O|OC\\(=O\\)N";
    return countMatches(smiles, pattern);
}

// Urea Count
double calculateSmilesUreaCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "NC\\(=O\\)N";
    return countMatches(smiles, pattern);
}

// Thiourea Count
double calculateSmilesThioureaCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "NC\\(=S\\)N";
    return countMatches(smiles, pattern);
}

// Oxime Count
double calculateSmilesOximeCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "C=NO|N=CO";
    return countMatches(smiles, pattern);
}

// Hydrazone Count
double calculateSmilesHydrazoneCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "C=NN|N=CN";
    return countMatches(smiles, pattern);
}

// Hydrazide Count
double calculateSmilesHydrazideCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "C\\(=O\\)NN|NNC\\(=O\\)";
    return countMatches(smiles, pattern);
}

// Azide Count
double calculateSmilesAzideCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "\\[N-\\]=[N+]=N|N=[N+]=[N-]"
    // Augmented: Add resonance forms like [N-]-[N+]#N
    const char* pattern = "\\[N-\\]=\\[N\\+\\]=N|N=\\[N\\+\\]=\\[N-\\]|\\[N-\\]\\[N\\+\\]#N|N#\\[N\\+\\]\\[N-\\]";
    return countMatches(smiles, pattern);
}

// Diazo Count
double calculateSmilesDiazoCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "C=N=N|N=N=C"
    // Augmented: Include resonance forms and aliphatic/aromatic carbon
    const char* pattern = "[Cc]=N=N|N=N=[Cc]|[Cc]=\\[N\\+\\]=\\[N-\\]|\\[N-\\]=\\[N\\+\\]=[Cc]|[Cc]-\\[N\\+\\]#N|N#\\[N\\+\\]-[Cc]";
    return countMatches(smiles, pattern);
}

// Disulfide Count
double calculateSmilesDisulfideCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "CSS|SSC" (Counts S-S bond with at least one C)
    // Augmented: Count general S-S single bonds. This assumes 'SS' isn't part of a multi-letter element e.g. 'OsS'.
    // Common SMILES for S-S is just "SS" or "S-S".
    const char* pattern = "S-S|SS"; // Matches explicit S-S or implicit SS
    return countMatches(smiles, pattern);
}


// Phosphate Count
double calculateSmilesPhosphateCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "P\\(=O\\)\\(O\\)\\(O\\)O|OP\\(=O\\)";
    return countMatches(smiles, pattern);
}


// Sulfonamide Count
double calculateSmilesSulfonamideCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "S\\(=O\\)\\(=O\\)N|NS\\(=O\\)\\(=O\\)";
    return countMatches(smiles, pattern);
}

// Spiro Compound Count
double calculateSmilesSpiroCompoundCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Match patterns with shared ring junction atoms
    const char* pattern = "\\[.\\]1[^1]*\\[.\\]2[^2]*1[^1]*2";
    return countMatches(smiles, pattern);
}

// Fused Ring Count
double calculateSmilesFusedRingCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Current pattern is a heuristic. Creating a truly robust regex for all fused ring systems is exceptionally complex.
    // This specific heuristic is: "\\[.\\]1[^1]*\\[.\\]2[^2]*\\[.\\][^1]*1[^2]*2"
    // No simple "synonym" augmentation for this kind of structural query with regex.
    // Leaving as is due to complexity.
    const char* pattern = "\\[.\\]1[^1]*\\[.\\]2[^2]*\\[.\\][^1]*1[^2]*2";
    return countMatches(smiles, pattern);
}

// Biphenyl Count
double calculateSmilesBiphenylCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "c1ccccc1-c1ccccc1";
    return countMatches(smiles, pattern);
}

// Naphthalene Count
double calculateSmilesNaphthaleneCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    const char* pattern = "c1ccc2ccccc2c1";
    return countMatches(smiles, pattern);
}

// Indole Count
double calculateSmilesIndoleCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "c1ccc2c\\(c1\\)[nH]cc2|c1ccc2c\\(c1\\)nc[cH]2"
    // Augmented: Add more common SMILES representations for the pyrrole part fused to benzene.
    const char* pattern = "c1ccc2c\\(c1\\)[nH]cc2|c1ccc2c\\(c1\\)nc[cH]2|c1ccc2c\\(c1\\)c[nH]c2|c1ccc2c\\(c1\\)cc[nH]c2|c1ccc2[nH]ccc2c1|c1ccc2c[nH]cc1";
    return countMatches(smiles, pattern);
}

// Benzimidazole Count
double calculateSmilesBenzimidazoleCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "c1ccc2c\\(c1\\)ncn2|c1ccc2nc[nH]c2c1"
    // Augmented: Add variations for [nH] placement and common SMILES forms.
    const char* pattern = "c1ccc2c\\(c1\\)nc[nH]2|c1ccc2c\\(c1\\)[nH]cn2|c1ccc2ncn[cH]2|c1nc2ccccc2[nH]1|c1[nH]c2ccccc2nc1";
    return countMatches(smiles, pattern);
}

// Benzothiazole Count
double calculateSmilesBenzothiazoleCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "c1ccc2c\\(c1\\)ncs2|c1ccc2scnc2c1"
    // Augmented: These two cover the main N-C-S and S-C-N arrangements in the thiazole part.
    // Adding minor variants if any for specific SMILES representations.
    const char* pattern = "c1ccc2c\\(c1\\)ncs2|c1ccc2c\\(c1\\)scn2|c1nc2ccccc2s1|c1sc2ccccc2n1";
    return countMatches(smiles, pattern);
}

// Benzoxazole Count
double calculateSmilesBenzoxazoleCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Original: "c1ccc2c\\(c1\\)nco2|c1ccc2ocnc2c1"
    // Augmented: Similar to benzothiazole, covering N-C-O and O-C-N.
    const char* pattern = "c1ccc2c\\(c1\\)nco2|c1ccc2c\\(c1\\)ocn2|c1nc2ccccc2o1|c1oc2ccccc2n1";
    return countMatches(smiles, pattern);
}

// Branching Index - Counts ratio of branch points to total atom count
double calculateSmilesBranchingIndex(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int branches = countChar(smiles, '(');
    int atoms = 0;
    for (const char* p = smiles; *p; p++) {
        if (isalpha(*p) && (p == smiles || !isalpha(*(p-1)))) atoms++;
    }
    return atoms > 0 ? (double)branches / atoms : 0.0;
}

// Heteroatom Clustering - Measures how clustered heteroatoms are vs. spread out
double calculateSmilesHeteroatomClustering(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int clusters = 0;
    int heteroatoms = 0;
    int in_cluster = 0;
    
    for (const char* p = smiles; *p; p++) {
        if (*p == 'N' || *p == 'O' || *p == 'S' || *p == 'P') {
            heteroatoms++;
            if (!in_cluster) {
                clusters++;
                in_cluster = 1;
            }
        } else if (isalpha(*p)) {
            in_cluster = 0;
        }
    }
    
    return heteroatoms > 0 ? (double)clusters / heteroatoms : 0.0;
}

// Alternating Element Pattern - Measures C-X-C-X patterns
double calculateSmilesAlternatingElementPattern(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int alternating = 0;
    size_t len = strlen(smiles);

    if (len < 4) return 0;

    for (size_t i = 0; i <= len - 4; i++) {
        char a1 = smiles[i];
        char a2 = smiles[i+1];
        char a3 = smiles[i+2];
        char a4 = smiles[i+3];

        // Check if a1 and a3 are carbon (aliphatic or aromatic)
        int is_a1_carbon = (a1 == 'C' || a1 == 'c');
        int is_a3_carbon = (a3 == 'C' || a3 == 'c');

        // Check if a2 and a4 are heteroatoms (N, O, S, P, n, o, s, p)
        int is_a2_hetero = (strchr("NOSPDnspd", a2) != NULL); // Added Boron 'B', 'b' could be added too if desired
        int is_a4_hetero = (strchr("NOSPDnspd", a4) != NULL);


        // Pattern: Carbon-Hetero-Carbon-Hetero
        if (is_a1_carbon && is_a2_hetero && is_a3_carbon && is_a4_hetero) {
            alternating++;
        }
        // Could add other patterns like Hetero-Carbon-Hetero-Carbon if the definition expands
        // e.g., if (is_a1_hetero && is_a2_carbon && is_a3_hetero && is_a4_carbon) { alternating++; }
    }
    return alternating;
}

// Aromatic Ring Edge Count - Counts aromatic atoms with non-aromatic neighbors
double calculateSmilesAromaticRingEdge(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int edges = 0;
    
    for (size_t i = 1; i < strlen(smiles) - 1; i++) {
        if (islower(smiles[i]) && (isupper(smiles[i-1]) || isupper(smiles[i+1]))) {
            edges++;
        }
    }
    
    return edges;
}

// Chiral Center Density - Ratio of chiral centers to carbon atoms
double calculateSmilesChiralCenterDensity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int chiral = countChar(smiles, '@');
    int carbons = countChar(smiles, 'C') + countChar(smiles, 'c');
    
    return carbons > 0 ? (double)chiral / carbons : 0.0;
}

// Structural Symmetry Estimate - Use repeating patterns as symmetry indicator
double calculateSmilesStructuralSymmetry(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int repeating_patterns = 0;
    char patterns[20][10] = {"CC", "CCC", "CCCC", "C1CC1", "c1ccc", "CNC", "COC"};
    
    for (int i = 0; i < 7; i++) {
        char* pattern = patterns[i];
        char* pos = (char*)smiles;
        int count = 0;
        while ((pos = strstr(pos, pattern)) != NULL) {
            count++;
            pos++;
        }
        if (count > 1) repeating_patterns += count;
    }
    
    return repeating_patterns;
}

// Ring Bridge Count - Counts bridges between ring systems
double calculateSmilesRingBridgeCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int bridges = 0;
    size_t len = strlen(smiles);
    
    // Need at least 6 characters for the pattern
    if (len < 6) {
        return 0;
    }
    
    // Look for patterns like ]1)-([
    for (size_t i = 0; i < len - 5; i++) {
        if (smiles[i] == ']' && 
            i + 1 < len && isdigit((unsigned char)smiles[i+1]) && 
            i + 2 < len && smiles[i+2] == ')' && 
            i + 3 < len && smiles[i+3] == '-' && 
            i + 4 < len && smiles[i+4] == '(' && 
            i + 5 < len && smiles[i+5] == '[') {
            bridges++;
        }
    }
    
    return bridges;
}

// Fingerprint Fragment Entropy - Measures the Shannon entropy of atom type distribution
double calculateSmilesAtomicEntropy(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int counts[26] = {0}; // For A-Z
    int total = 0;
    
    for (const char* p = smiles; *p; p++) {
        if (isalpha(*p) && isupper(*p)) {
            counts[*p - 'A']++;
            total++;
        }
    }
    
    double entropy = 0.0;
    for (int i = 0; i < 26; i++) {
        if (counts[i] > 0) {
            double p = (double)counts[i] / total;
            entropy -= p * log(p) / log(2.0);
        }
    }
    
    return entropy;
}

// Nitrogen Position Variance - Measures how consistently N appears in the SMILES
double calculateSmilesNitrogenPositionVariance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int positions[100] = {0};
    int count = 0;
    
    for (size_t i = 0; i < strlen(smiles) && count < 100; i++) {
        if (smiles[i] == 'N' || smiles[i] == 'n') {
            positions[count++] = i;
        }
    }
    
    if (count < 2) return 0;
    
    // Calculate mean position
    double mean = 0;
    for (int i = 0; i < count; i++) {
        mean += positions[i];
    }
    mean /= count;
    
    // Calculate variance
    double variance = 0;
    for (int i = 0; i < count; i++) {
        double diff = positions[i] - mean;
        variance += diff * diff;
    }
    variance /= count;
    
    return variance;
}

// Bond Type Diversity - Ratio of unique bond types to total bonds
double calculateSmilesBondTypeDiversity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int single_bonds = countChar(smiles, '-');
    int double_bonds = countChar(smiles, '=');
    int triple_bonds = countChar(smiles, '#');
    int aromatic_bonds = 0;
    
    // Approximate aromatic bonds by counting lowercase atoms
    for (const char* p = smiles; *p; p++) {
        if (islower(*p) && isalpha(*p)) {
            aromatic_bonds++;
        }
    }
    
    int unique_types = 0;
    if (single_bonds > 0) unique_types++;
    if (double_bonds > 0) unique_types++;
    if (triple_bonds > 0) unique_types++;
    if (aromatic_bonds > 0) unique_types++;
    
    int total_bonds = single_bonds + double_bonds + triple_bonds + aromatic_bonds;
    return total_bonds > 0 ? (double)unique_types / 4.0 : 0.0;
}

// Adjacent Heteroatom Ratio - Fraction of heteroatoms with adjacent heteroatoms
double calculateSmilesAdjacentHeteroatomRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int total_heteroatoms = 0;
    int adjacent_heteroatoms = 0;
    
    for (size_t i = 0; i < strlen(smiles); i++) {
        if (smiles[i] == 'N' || smiles[i] == 'O' || smiles[i] == 'S' || smiles[i] == 'P') {
            total_heteroatoms++;
            if ((i > 0 && (smiles[i-1] == 'N' || smiles[i-1] == 'O' || smiles[i-1] == 'S' || smiles[i-1] == 'P')) ||
                (i < strlen(smiles) - 1 && (smiles[i+1] == 'N' || smiles[i+1] == 'O' || smiles[i+1] == 'S' || smiles[i+1] == 'P'))) {
                adjacent_heteroatoms++;
            }
        }
    }
    
    return total_heteroatoms > 0 ? (double)adjacent_heteroatoms / total_heteroatoms : 0.0;
}

// Heteroatom to Carbon Ratio with Position Weighting
double calculateSmilesPositionWeightedHeteroRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    double weighted_hetero = 0;
    double weighted_carbon = 0;
    size_t len = strlen(smiles);
    
    for (size_t i = 0; i < len; i++) {
        double position_weight = 1.0 - ((double)i / len) * 0.5;  // Weight decreases with position
        if (smiles[i] == 'C') {
            weighted_carbon += position_weight;
        } else if (smiles[i] == 'N' || smiles[i] == 'O' || smiles[i] == 'S' || smiles[i] == 'P') {
            weighted_hetero += position_weight;
        }
    }
    
    return weighted_carbon > 0 ? weighted_hetero / weighted_carbon : 0.0;
}

// SMILES Fragment Complexity - Counts brackets, parentheses, and digits as complexity indicators
double calculateSmilesFragmentComplexity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int complexity = 0;
    
    complexity += countChar(smiles, '[') * 2;  // Brackets are weighted more
    complexity += countChar(smiles, '(') * 1.5;
    complexity += countDigits(smiles);
    complexity += countChar(smiles, '@') * 3;  // Chirality adds complexity
    
    int len = strlen(smiles);
    return len > 0 ? (double)complexity / len : 0.0;
}

// Oxygen Clustering Factor - Measures how much oxygens cluster together
double calculateSmilesOxygenClusteringFactor(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int total_oxygens = countChar(smiles, 'O');
    int clustered_oxygens = 0;
    
    // Count oxygen atoms adjacent to other oxygen atoms
    for (size_t i = 1; i < strlen(smiles); i++) {
        if (smiles[i] == 'O' && smiles[i-1] == 'O') {
            clustered_oxygens++;
        }
    }
    
    return total_oxygens > 0 ? (double)clustered_oxygens / total_oxygens : 0.0;
}

// Primary Amine to Total Nitrogen Ratio
double calculateSmilesPrimaryAmineRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int primary_amines = countMatches(smiles, "N\\([H]\\)\\([H]\\)|N\\([H]\\)[H]|\\[NH2\\]");
    int total_nitrogens = countChar(smiles, 'N') + countChar(smiles, 'n');
    
    return total_nitrogens > 0 ? (double)primary_amines / total_nitrogens : 0.0;
}

// Methylation Degree - Ratio of methyl groups to carbon atoms
double calculateSmilesMethylationDegree(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int methyl_groups = countMatches(smiles, "\\(C\\)|C\\)");
    int total_carbons = countChar(smiles, 'C') + countChar(smiles, 'c');
    
    return total_carbons > 0 ? (double)methyl_groups / total_carbons : 0.0;
}

// Non-Carbon Backbone Ratio - Measures proportion of backbone that's not carbon
double calculateSmilesNonCarbonBackboneRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int backbone_length = 0;
    int noncarbon_backbone = 0;
    
    // Simple approximation - count atoms not in branches
    int branch_depth = 0;
    for (const char* p = smiles; *p; p++) {
        if (*p == '(') {
            branch_depth++;
        } else if (*p == ')') {
            branch_depth--;
        } else if (branch_depth == 0 && isupper(*p)) {
            backbone_length++;
            if (*p != 'C') {
                noncarbon_backbone++;
            }
        }
    }
    
    return backbone_length > 0 ? (double)noncarbon_backbone / backbone_length : 0.0;
}

// Carbohydrate-likeness Score
double calculateSmilesCarbohydrateLikeness(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Carbohydrates typically have C:H:O ratio of 1:2:1 and hydroxyl groups
    int carbons = countChar(smiles, 'C') + countChar(smiles, 'c');
    int oxygens = countChar(smiles, 'O') + countChar(smiles, 'o');
    int hydroxyls = countMatches(smiles, "OH|\\[OH\\]");
    
    double ratio_score = 0;
    if (carbons > 0 && oxygens > 0) {
        double co_ratio = (double)carbons / oxygens;
        ratio_score = 1.0 - fabs(co_ratio - 1.0) * 0.5;  // Score decreases as we move from ideal 1:1 ratio
        if (ratio_score < 0) ratio_score = 0;
    }
    
    double hydroxyl_score = carbons > 0 ? (double)hydroxyls / carbons : 0;
    return (ratio_score * 0.5 + hydroxyl_score * 0.5);
}

// Trifunctional Carbon Ratio
double calculateSmilesTrifunctionalCarbonRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int trifunctional_carbons = countMatches(smiles, "C\\([^)]*\\)\\([^)]*\\)\\([^)]*\\)");
    int total_carbons = countChar(smiles, 'C');
    
    return total_carbons > 0 ? (double)trifunctional_carbons / total_carbons : 0.0;
}

// Ring Junction Density
double calculateSmilesRingJunctionDensity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int junctions = 0;
    
    // Look for digits that appear more than once (shared ring junction)
    for (int i = 1; i <= 9; i++) {
        char digit = '0' + i;
        int count = 0;
        for (const char* p = smiles; *p; p++) {
            if (*p == digit) count++;
        }
        if (count >= 2) junctions++;
    }
    
    int ring_count = 0;
    for (int i = 1; i <= 9; i++) {
        if (strchr(smiles, '0' + i) != NULL) {
            ring_count++;
        }
    }
    
    return ring_count > 0 ? (double)junctions / ring_count : 0;
}

// Sulfur Oxidation State Distribution
double calculateSmilesSulfurOxidationStateDistribution(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int sulfur_count = countChar(smiles, 'S');
    int sulfoxide_count = countMatches(smiles, "S\\(=O\\)");
    int sulfone_count = countMatches(smiles, "S\\(=O\\)\\(=O\\)");
    
    // Calculate weighted average oxidation state
    double weighted_sum = 0;
    if (sulfur_count > 0) {
        weighted_sum = (sulfur_count - sulfoxide_count - sulfone_count) * 0 + 
                       sulfoxide_count * 2 + 
                       sulfone_count * 4;
        return weighted_sum / sulfur_count; // Average oxidation state
    }
    return 0.0;
}

// Phosphorus Oxidation State Distribution
double calculateSmilesPhosphorusOxidationStateDistribution(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int phosphorus_count = countChar(smiles, 'P');
    int phosphine_oxide_count = countMatches(smiles, "P\\(=O\\)");
    int phosphate_count = countMatches(smiles, "P\\(=O\\)\\(O\\)\\(O\\)O");
    
    // Calculate weighted average oxidation state
    double weighted_sum = 0;
    if (phosphorus_count > 0) {
        weighted_sum = (phosphorus_count - phosphine_oxide_count - phosphate_count) * 0 + 
                      phosphine_oxide_count * 2 + 
                      phosphate_count * 5;
        return weighted_sum / phosphorus_count; // Average oxidation state
    }
    return 0.0;
}

// Hydrophilic Group Density
double calculateSmilesHydrophilicGroupDensity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Count hydrophilic groups
    int hydrophilic = 0;
    hydrophilic += countMatches(smiles, "OH|\\[OH\\]");  // Hydroxyls
    hydrophilic += countMatches(smiles, "NH2|\\[NH2\\]"); // Primary amines
    hydrophilic += countMatches(smiles, "C\\(=O\\)O"); // Carboxylic acids
    
    // Get total atom count (approximate)
    int atom_count = 0;
    for (const char* p = smiles; *p; p++) {
        if (isalpha(*p) && isupper(*p)) {
            atom_count++;
        }
    }
    
    return atom_count > 0 ? (double)hydrophilic / atom_count : 0.0;
}

// Terminal Group Diversity
double calculateSmilesTerminalGroupDiversity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Count unique terminal groups
    int methyl_terminals = countMatches(smiles, "\\)C");
    int hydroxyl_terminals = countMatches(smiles, "\\)O");
    int amino_terminals = countMatches(smiles, "\\)N");
    int halide_terminals = countMatches(smiles, "\\)[FClBrI]");
    
    // Count how many different terminal types are present
    int unique_terminals = 0;
    if (methyl_terminals > 0) unique_terminals++;
    if (hydroxyl_terminals > 0) unique_terminals++;
    if (amino_terminals > 0) unique_terminals++;
    if (halide_terminals > 0) unique_terminals++;
    
    int total_terminals = methyl_terminals + hydroxyl_terminals + amino_terminals + halide_terminals;
    return total_terminals > 0 ? (double)unique_terminals / 4.0 : 0.0;
}

// Nitrogen Environment Diversity
double calculateSmilesNitrogenEnvironmentDiversity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Count different nitrogen environments
    int amino_n = countMatches(smiles, "N\\([H]\\)");
    int amide_n = countMatches(smiles, "NC\\(=O\\)");
    int aromatic_n = countChar(smiles, 'n');
    int nitro_n = countMatches(smiles, "N\\(=O\\)=O|\\[N\\+\\]\\(=O\\)\\[O-\\]");
    
    // Count how many different environments are present
    int unique_envs = 0;
    if (amino_n > 0) unique_envs++;
    if (amide_n > 0) unique_envs++;
    if (aromatic_n > 0) unique_envs++;
    if (nitro_n > 0) unique_envs++;
    
    int total_n = countChar(smiles, 'N') + countChar(smiles, 'n');
    return total_n > 0 ? (double)unique_envs / 4.0 : 0.0;
}

// Chain Length Distribution
double calculateSmilesChainLengthDistribution(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Find longest continuous carbon chain
    int longest_chain = 0;
    int current_chain = 0;
    
    for (const char* p = smiles; *p; p++) {
        if (*p == 'C') {
            current_chain++;
            if (current_chain > longest_chain) {
                longest_chain = current_chain;
            }
        } else if (!isalpha(*p) || *p == 'c') {
            // Don't reset on bond symbols or aromatic carbons
            continue;
        } else {
            current_chain = 0;
        }
    }
    
    int total_carbons = countChar(smiles, 'C');
    return total_carbons > 0 ? (double)longest_chain / total_carbons : 0.0;
}

// Carbon Aromaticity Fraction
double calculateSmilesCarbonAromaticityFraction(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int aromatic_c = countChar(smiles, 'c');
    int aliphatic_c = countChar(smiles, 'C');
    int total_c = aromatic_c + aliphatic_c;
    
    return total_c > 0 ? (double)aromatic_c / total_c : 0.0;
}

// Nitrogen Aromaticity Fraction
double calculateSmilesNitrogenAromaticityFraction(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int aromatic_n = countChar(smiles, 'n');
    int aliphatic_n = countChar(smiles, 'N');
    int total_n = aromatic_n + aliphatic_n;
    
    return total_n > 0 ? (double)aromatic_n / total_n : 0.0;
}

// Chirality Center Distribution
double calculateSmilesChiralityDistribution(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Count chiral centers and their spread
    int chiral_count = countChar(smiles, '@');
    if (chiral_count < 2) return 0.0;  // Need at least 2 centers for distribution
    
    int positions[100] = {0};
    int count = 0;
    
    for (size_t i = 0; i < strlen(smiles) && count < chiral_count && count < 100; i++) {
        if (smiles[i] == '@') {
            positions[count++] = i;
        }
    }
    
    // Calculate average distance between chiral centers
    double total_distance = 0;
    for (int i = 1; i < count; i++) {
        total_distance += positions[i] - positions[i-1];
    }
    double avg_distance = total_distance / (count - 1);
    
    // Normalize by SMILES length
    return avg_distance / strlen(smiles);
}

// Fluorine to Carbon Ratio
double calculateSmilesFluorineToCarbonRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int fluorine_count = countChar(smiles, 'F');
    int carbon_count = countChar(smiles, 'C') + countChar(smiles, 'c');
    
    return carbon_count > 0 ? (double)fluorine_count / carbon_count : 0.0;
}

// Amide Fraction
double calculateSmilesAmideFraction(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int amide_count = countMatches(smiles, "C\\(=O\\)N|NC\\(=O\\)");
    int total_nitrogen = countChar(smiles, 'N') + countChar(smiles, 'n');
    
    return total_nitrogen > 0 ? (double)amide_count / total_nitrogen : 0.0;
}

// Heterocycle Count Normalized
double calculateSmilesHeterocycleCountNormalized(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Count nitrogen-containing, oxygen-containing, and sulfur-containing heterocycles
    int n_heterocycles = countMatches(smiles, "[a-z]*n[a-z]*[0-9]");
    int o_heterocycles = countMatches(smiles, "[a-z]*o[a-z]*[0-9]");
    int s_heterocycles = countMatches(smiles, "[a-z]*s[a-z]*[0-9]");
    
    int total_cycles = 0;
    for (int i = 1; i <= 9; i++) {
        char digit = '0' + i;
        if (strchr(smiles, digit) != NULL) {
            total_cycles++;
        }
    }
    
    return total_cycles > 0 ? (double)(n_heterocycles + o_heterocycles + s_heterocycles) / total_cycles : 0.0;
}

// Carbonyl Distribution
double calculateSmilesCarbonylDistribution(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    // Count different carbonyl-containing groups
    int ketone_count = countMatches(smiles, "C\\(=O\\)C");
    int aldehyde_count = countMatches(smiles, "C\\(=O\\)[H]");
    int amide_count = countMatches(smiles, "C\\(=O\\)N");
    int ester_count = countMatches(smiles, "C\\(=O\\)O");
    
    // Calculate entropy of distribution
    int total = ketone_count + aldehyde_count + amide_count + ester_count;
    if (total == 0) return 0.0;
    
    double entropy = 0.0;
    if (ketone_count > 0) {
        double p = (double)ketone_count / total;
        entropy -= p * log(p) / log(2.0);
    }
    if (aldehyde_count > 0) {
        double p = (double)aldehyde_count / total;
        entropy -= p * log(p) / log(2.0);
    }
    if (amide_count > 0) {
        double p = (double)amide_count / total;
        entropy -= p * log(p) / log(2.0);
    }
    if (ester_count > 0) {
        double p = (double)ester_count / total;
        entropy -= p * log(p) / log(2.0);
    }
    
    return entropy / 2.0;  // Normalize to 0-1 range (log2(4) = 2)
}

// Nitrogen Oxidation State Distribution
double calculateSmilesNitrogenOxidationStateDistribution(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int n_minus3 = countMatches(smiles, "N\\([H]\\)\\([H]\\)\\([H]\\)|\\[NH3\\]");  // NH3 (-3)
    int n_minus2 = countMatches(smiles, "N\\([H]\\)\\([H]\\)C|\\[NH2\\]C");        // NH2-C (-2)
    int n_minus1 = countMatches(smiles, "N\\([H]\\)C|\\[NH\\]C");                  // NH-C (-1)
    int n_plus1 = countMatches(smiles, "N=O|\\[N\\+\\]=O");                       // N=O (+1)
    int n_plus3 = countMatches(smiles, "N\\(=O\\)=O|\\[N\\+\\]\\(=O\\)\\[O-\\]"); // NO2 (+3)
    int n_plus5 = countMatches(smiles, "N\\(=O\\)\\(=O\\)=O");                    // NO3 (+5)
    
    int total_n = countChar(smiles, 'N') + countChar(smiles, 'n');
    if (total_n == 0) return 0.0;
    
    // Calculate weighted average oxidation state
    double weighted_sum = n_minus3 * (-3) + n_minus2 * (-2) + n_minus1 * (-1) + 
                          n_plus1 * 1 + n_plus3 * 3 + n_plus5 * 5;
    
    double remaining_n = total_n - (n_minus3 + n_minus2 + n_minus1 + n_plus1 + n_plus3 + n_plus5);
    weighted_sum += remaining_n * 0;  // Assume neutral state for remaining N
    
    return weighted_sum / total_n;
}

// Change this function's name to avoid duplication
// Old name: calculateSmilesAmideFraction
double calculateSmilesAmideToNitrogenRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    int amide_count = countMatches(smiles, "C\\(=O\\)N|NC\\(=O\\)");
    int total_nitrogen = countChar(smiles, 'N') + countChar(smiles, 'n');
    
    return total_nitrogen > 0 ? (double)amide_count / total_nitrogen : 0.0;
}

// Optimize performance by using specialized functions instead of regex for common patterns

// 1. Replace simple regex patterns with direct string functions
double calculateSmilesAliphaticChainCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    
    // Optimize by using direct substring counting for better performance
    return countSubstr(smiles, "CCC");
}

// Add a function to initialize all optimizations when the module loads
__attribute__((constructor))
void initRegexOptimizations() {
    // This is now a no-op as we're directly optimizing the implementations
}