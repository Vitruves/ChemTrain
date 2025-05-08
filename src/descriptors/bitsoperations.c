// src/descriptors/bitsoperations.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>

// --- Lookup Tables ---

// Popcount LUT
static const uint8_t POPCOUNT_LUT[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

// Leading Zeros LUT
static const uint8_t LEADING_ZEROS_LUT[256] = {
    8, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

// Trailing Zeros LUT
static const uint8_t TRAILING_ZEROS_LUT[256] = {
    8, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0
};

// CRC8 Table
static const uint8_t CRC8_TABLE[256] = {
    0x00, 0x07, 0x0E, 0x09, 0x1C, 0x1B, 0x12, 0x15,
    0x38, 0x3F, 0x36, 0x31, 0x24, 0x23, 0x2A, 0x2D,
    0x70, 0x77, 0x7E, 0x79, 0x6C, 0x6B, 0x62, 0x65,
    0x48, 0x4F, 0x46, 0x41, 0x54, 0x53, 0x5A, 0x5D,
    0xE0, 0xE7, 0xEE, 0xE9, 0xFC, 0xFB, 0xF2, 0xF5,
    0xD8, 0xDF, 0xD6, 0xD1, 0xC4, 0xC3, 0xCA, 0xCD,
    0x90, 0x97, 0x9E, 0x99, 0x8C, 0x8B, 0x82, 0x85,
    0xA8, 0xAF, 0xA6, 0xA1, 0xB4, 0xB3, 0xBA, 0xBD,
    0xC7, 0xC0, 0xC9, 0xCE, 0xDB, 0xDC, 0xD5, 0xD2,
    0xFF, 0xF8, 0xF1, 0xF6, 0xE3, 0xE4, 0xED, 0xEA,
    0xB7, 0xB0, 0xB9, 0xBE, 0xAB, 0xAC, 0xA5, 0xA2,
    0x8F, 0x88, 0x81, 0x86, 0x93, 0x94, 0x9D, 0x9A,
    0x27, 0x20, 0x29, 0x2E, 0x3B, 0x3C, 0x35, 0x32,
    0x1F, 0x18, 0x11, 0x16, 0x03, 0x04, 0x0D, 0x0A,
    0x57, 0x50, 0x59, 0x5E, 0x4B, 0x4C, 0x45, 0x42,
    0x6F, 0x68, 0x61, 0x66, 0x73, 0x74, 0x7D, 0x7A,
    0x89, 0x8E, 0x87, 0x80, 0x95, 0x92, 0x9B, 0x9C,
    0xB1, 0xB6, 0xBF, 0xB8, 0xAD, 0xAA, 0xA3, 0xA4,
    0xF9, 0xFE, 0xF7, 0xF0, 0xE5, 0xE2, 0xEB, 0xEC,
    0xC1, 0xC6, 0xCF, 0xC8, 0xDD, 0xDA, 0xD3, 0xD4,
    0x69, 0x6E, 0x67, 0x60, 0x75, 0x72, 0x7B, 0x7C,
    0x51, 0x56, 0x5F, 0x58, 0x4D, 0x4A, 0x43, 0x44,
    0x19, 0x1E, 0x17, 0x10, 0x05, 0x02, 0x0B, 0x0C,
    0x21, 0x26, 0x2F, 0x28, 0x3D, 0x3A, 0x33, 0x34,
    0x4E, 0x49, 0x40, 0x47, 0x52, 0x55, 0x5C, 0x5B,
    0x76, 0x71, 0x78, 0x7F, 0x6A, 0x6D, 0x64, 0x63,
    0x3E, 0x39, 0x30, 0x37, 0x22, 0x25, 0x2C, 0x2B,
    0x06, 0x01, 0x08, 0x0F, 0x1A, 0x1D, 0x14, 0x13,
    0xAE, 0xA9, 0xA0, 0xA7, 0xB2, 0xB5, 0xBC, 0xBB,
    0x96, 0x91, 0x98, 0x9F, 0x8A, 0x8D, 0x84, 0x83,
    0xDE, 0xD9, 0xD0, 0xD7, 0xC2, 0xC5, 0xCC, 0xCB,
    0xE6, 0xE1, 0xE8, 0xEF, 0xFA, 0xFD, 0xF4, 0xF3
};

// --- C Helper Functions ---

// Convert SMILES to byte array
uint8_t* smiles_to_bytes(const char* smiles, size_t* num_bytes) {
    *num_bytes = strlen(smiles);
    uint8_t* bytes = (uint8_t*)malloc(*num_bytes * sizeof(uint8_t));
    if (!bytes) return NULL;
    memcpy(bytes, smiles, *num_bytes);
    return bytes;
}

// Convert SMILES to bit array
bool* smiles_to_bits(const char* smiles, size_t* num_bits) {
    size_t len = strlen(smiles);
    *num_bits = len * 8;
    bool* bits = (bool*)malloc(*num_bits * sizeof(bool));
    if (!bits) return NULL;
    for (size_t i = 0; i < len; ++i) {
        uint8_t byte = (uint8_t)smiles[i];
        for (int j = 0; j < 8; ++j) {
            bits[i * 8 + j] = (byte >> (7 - j)) & 1;
        }
    }
    return bits;
}

// Population count for byte array
uint32_t popcount_bytes_c(const uint8_t* bytes, size_t num_bytes) {
    uint32_t count = 0;
    for (size_t i = 0; i < num_bytes; ++i) {
        count += POPCOUNT_LUT[bytes[i]];
    }
    return count;
}

// Bit transition count
uint32_t bit_transition_count_c(const bool* bits, size_t num_bits) {
    if (num_bits < 2) return 0;
    uint32_t count = 0;
    for (size_t i = 1; i < num_bits; ++i) {
        if (bits[i] != bits[i-1]) {
            count++;
        }
    }
    return count;
}

// Longest consecutive bit run
uint32_t longest_consecutive_bit_run_c(const bool* bits, size_t num_bits) {
    if (num_bits == 0) return 0;
    uint32_t maxRun = 0;
    uint32_t currentRun = 1;
    for (size_t i = 1; i < num_bits; ++i) {
        if (bits[i] == bits[i-1]) {
            currentRun++;
        } else {
            if (currentRun > maxRun) maxRun = currentRun;
            currentRun = 1;
        }
    }
    return (currentRun > maxRun) ? currentRun : maxRun;
}

// Byte XOR checksum
uint8_t byte_xor_checksum_c(const uint8_t* bytes, size_t num_bytes) {
    uint8_t checksum = 0;
    for (size_t i = 0; i < num_bytes; ++i) {
        checksum ^= bytes[i];
    }
    return checksum;
}

// Calculate entropy of a byte histogram
double calculate_entropy_c(const uint32_t* histogram, size_t totalCount) {
    if (totalCount == 0) return 0.0;
    double entropy = 0.0;
    for (int i = 0; i < 256; ++i) {
        if (histogram[i] > 0) {
            double probability = (double)histogram[i] / totalCount;
            entropy -= probability * log2(probability);
        }
    }
    return entropy;
}

// FNV-1a hash
uint32_t fnv1a_hash_c(const uint8_t* bytes, size_t num_bytes) {
    uint32_t hash = 2166136261u;
    for (size_t i = 0; i < num_bytes; ++i) {
        hash ^= bytes[i];
        hash *= 16777619u;
    }
    return hash;
}

// CRC8 checksum
uint8_t crc8_c(const uint8_t* bytes, size_t num_bytes) {
    uint8_t crc = 0;
    for (size_t i = 0; i < num_bytes; ++i) {
        crc = CRC8_TABLE[crc ^ bytes[i]];
    }
    return crc;
}

// Check if byte array is a palindrome
bool is_char_palindrome_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes <= 1) return true;
    for (size_t i = 0; i < num_bytes / 2; ++i) {
        if (bytes[i] != bytes[num_bytes - 1 - i]) return false;
    }
    return true;
}

// Check if bit array is a palindrome
bool is_bit_palindrome_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return true;
    for (size_t i = 0; i < num_bits / 2; ++i) {
        if (bits[i] != bits[num_bits - 1 - i]) return false;
    }
    return true;
}

// Bit distribution variance
double bit_distribution_variance_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes == 0) return 0.0;
    double popcount_sum = 0.0;
    for (size_t i = 0; i < num_bytes; ++i) {
        popcount_sum += POPCOUNT_LUT[bytes[i]];
    }
    double mean = popcount_sum / num_bytes;
    double variance = 0.0;
    for (size_t i = 0; i < num_bytes; ++i) {
        double diff = POPCOUNT_LUT[bytes[i]] - mean;
        variance += diff * diff;
    }
    return variance / num_bytes;
}

// Left/right bit balance
int left_right_bit_balance_c(const bool* bits, size_t num_bits) {
    if (num_bits == 0) return 0;
    size_t mid = num_bits / 2;
    int leftCount = 0, rightCount = 0;
    for (size_t i = 0; i < mid; ++i) if (bits[i]) leftCount++;
    for (size_t i = mid; i < num_bits; ++i) if (bits[i]) rightCount++;
    return leftCount - rightCount;
}

// qsort comparison function for uint8_t
int compare_uint8(const void* a, const void* b) {
    return (*(uint8_t*)a - *(uint8_t*)b);
}

// --- C Descriptor Implementations ---

// Hamming weight (bit count) of a byte array
uint32_t hamming_weight_c(const uint8_t* bytes, size_t num_bytes) {
    return popcount_bytes_c(bytes, num_bytes);
}

// Hamming distance between two byte arrays
uint32_t hamming_distance_c(const uint8_t* bytes1, const uint8_t* bytes2, size_t num_bytes) {
    uint32_t distance = 0;
    for (size_t i = 0; i < num_bytes; ++i) {
        distance += POPCOUNT_LUT[bytes1[i] ^ bytes2[i]];
    }
    return distance;
}

// Jaccard similarity between two bit arrays
double jaccard_similarity_c(const uint8_t* bytes1, const uint8_t* bytes2, size_t num_bytes) {
    uint32_t intersection = 0;
    uint32_t union_count = 0;
    
    for (size_t i = 0; i < num_bytes; ++i) {
        uint8_t intersect_byte = bytes1[i] & bytes2[i];
        uint8_t union_byte = bytes1[i] | bytes2[i];
        
        intersection += POPCOUNT_LUT[intersect_byte];
        union_count += POPCOUNT_LUT[union_byte];
    }
    
    return (union_count == 0) ? 0.0 : (double)intersection / union_count;
}

// Dice coefficient between two bit arrays
double dice_coefficient_c(const uint8_t* bytes1, const uint8_t* bytes2, size_t num_bytes) {
    uint32_t a = popcount_bytes_c(bytes1, num_bytes);
    uint32_t b = popcount_bytes_c(bytes2, num_bytes);
    
    uint32_t intersection = 0;
    for (size_t i = 0; i < num_bytes; ++i) {
        intersection += POPCOUNT_LUT[bytes1[i] & bytes2[i]];
    }
    
    return (a + b == 0) ? 0.0 : (2.0 * intersection) / (a + b);
}

// Tanimoto coefficient (same as Jaccard for binary data)
double tanimoto_coefficient_c(const uint8_t* bytes1, const uint8_t* bytes2, size_t num_bytes) {
    return jaccard_similarity_c(bytes1, bytes2, num_bytes);
}

// Count leading zeros in bit array
uint32_t count_leading_zeros_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes == 0) return 0;
    
    uint32_t count = 0;
    for (size_t i = 0; i < num_bytes; ++i) {
        if (bytes[i] == 0) {
            count += 8;
        } else {
            count += LEADING_ZEROS_LUT[bytes[i]];
            break;
        }
    }
    return count;
}

// Count trailing zeros in bit array
uint32_t count_trailing_zeros_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes == 0) return 0;
    
    uint32_t count = 0;
    for (size_t i = num_bytes; i > 0; --i) {
        if (bytes[i-1] == 0) {
            count += 8;
        } else {
            count += TRAILING_ZEROS_LUT[bytes[i-1]];
            break;
        }
    }
    return count;
}

// Bit clustering factor - measures how bits tend to cluster
double bit_clustering_factor_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0.0;
    
    uint32_t transitions = bit_transition_count_c(bits, num_bits);
    uint32_t ones = 0;
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) ones++;
    }
    
    double expected_transitions = 2.0 * ones * (num_bits - ones) / num_bits;
    return (expected_transitions == 0) ? 0.0 : (transitions / expected_transitions);
}

// Average bit run length
double average_bit_run_length_c(const bool* bits, size_t num_bits) {
    if (num_bits == 0) return 0.0;
    
    uint32_t transitions = bit_transition_count_c(bits, num_bits);
    return (transitions == 0) ? (double)num_bits : (double)num_bits / (transitions + 1);
}

// Bit alternation ratio (ratio of transitions to maximum possible)
double bit_alternation_ratio_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0.0;
    
    uint32_t transitions = bit_transition_count_c(bits, num_bits);
    return (double)transitions / (num_bits - 1);
}

// Longest sequence of zeros
uint32_t longest_zero_sequence_c(const bool* bits, size_t num_bits) {
    uint32_t max_length = 0;
    uint32_t current_length = 0;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (!bits[i]) {
            current_length++;
            if (current_length > max_length) max_length = current_length;
        } else {
            current_length = 0;
        }
    }
    
    return max_length;
}

// Longest sequence of ones
uint32_t longest_one_sequence_c(const bool* bits, size_t num_bits) {
    uint32_t max_length = 0;
    uint32_t current_length = 0;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            current_length++;
            if (current_length > max_length) max_length = current_length;
        } else {
            current_length = 0;
        }
    }
    
    return max_length;
}

// Bit dispersion - how evenly bits are distributed
double bit_dispersion_c(const bool* bits, size_t num_bits) {
    if (num_bits < 2) return 0.0;
    
    uint32_t ones = 0;
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) ones++;
    }
    
    if (ones == 0 || ones == num_bits) return 0.0;
    
    size_t* positions = (size_t*)malloc(ones * sizeof(size_t));
    if (!positions) return 0.0;
    
    size_t pos_idx = 0;
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) positions[pos_idx++] = i;
    }
    
    double avg_distance = 0.0;
    for (size_t i = 1; i < ones; ++i) {
        avg_distance += (positions[i] - positions[i-1]);
    }
    avg_distance /= (ones - 1);
    
    double variance = 0.0;
    for (size_t i = 1; i < ones; ++i) {
        double diff = (positions[i] - positions[i-1]) - avg_distance;
        variance += diff * diff;
    }
    variance /= (ones - 1);
    
    free(positions);
    return variance;
}

// Cosine similarity between bit vectors
double cosine_similarity_bits_c(const uint8_t* bytes1, const uint8_t* bytes2, size_t num_bytes) {
    uint32_t dot_product = 0;
    uint32_t magnitude1 = 0;
    uint32_t magnitude2 = 0;
    
    for (size_t i = 0; i < num_bytes; ++i) {
        dot_product += POPCOUNT_LUT[bytes1[i] & bytes2[i]];
        magnitude1 += POPCOUNT_LUT[bytes1[i]];
        magnitude2 += POPCOUNT_LUT[bytes2[i]];
    }
    
    return (magnitude1 == 0 || magnitude2 == 0) ? 0.0 : 
           (double)dot_product / (sqrt((double)magnitude1) * sqrt((double)magnitude2));
}

// Run-length encoding compression ratio
double rle_compression_ratio_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 1.0;
    
    uint32_t transitions = bit_transition_count_c(bits, num_bits);
    uint32_t runs = transitions + 1;
    
    double theoretical_size = runs * (log2(num_bits) + 1);
    return theoretical_size / num_bits;
}

// Bit density - ratio of set bits to total bits
double bit_density_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes == 0) return 0.0;
    
    uint32_t bit_count = popcount_bytes_c(bytes, num_bytes);
    return (double)bit_count / (num_bytes * 8);
}

// Bit periodicity - correlation between bits at distance d
double bit_periodicity_c(const bool* bits, size_t num_bits, size_t distance) {
    if (num_bits <= distance) return 0.0;
    
    size_t valid_positions = num_bits - distance;
    double correlation = 0.0;
    
    for (size_t i = 0; i < valid_positions; ++i) {
        correlation += (bits[i] == bits[i + distance]) ? 1.0 : 0.0;
    }
    
    return correlation / valid_positions;
}

// Autocorrelation - average correlation across multiple distances
double bit_autocorrelation_c(const bool* bits, size_t num_bits, size_t max_distance) {
    if (num_bits <= 1 || max_distance == 0) return 0.0;
    
    size_t actual_max = (max_distance < num_bits) ? max_distance : num_bits - 1;
    double sum_correlation = 0.0;
    
    for (size_t d = 1; d <= actual_max; ++d) {
        sum_correlation += bit_periodicity_c(bits, num_bits, d);
    }
    
    return sum_correlation / actual_max;
}

// Bit balancing - ratio of 1's to 0's
double bit_balance_ratio_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes == 0) return 0.0;
    
    uint32_t ones = popcount_bytes_c(bytes, num_bytes);
    uint32_t zeros = num_bytes * 8 - ones;
    
    return (zeros == 0) ? INFINITY : (double)ones / zeros;
}

// Bit entropy - information content of bit distribution
double bit_entropy_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes == 0) return 0.0;
    
    uint32_t ones = popcount_bytes_c(bytes, num_bytes);
    uint32_t total_bits = num_bytes * 8;
    uint32_t zeros = total_bits - ones;
    
    double p_one = (double)ones / total_bits;
    double p_zero = (double)zeros / total_bits;
    
    double entropy = 0.0;
    if (p_one > 0) entropy -= p_one * log2(p_one);
    if (p_zero > 0) entropy -= p_zero * log2(p_zero);
    
    return entropy;
}

// Byte diversity - number of unique bytes / total bytes
double byte_diversity_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes == 0) return 0.0;
    
    bool seen[256] = {false};
    uint32_t unique_count = 0;
    
    for (size_t i = 0; i < num_bytes; ++i) {
        if (!seen[bytes[i]]) {
            seen[bytes[i]] = true;
            unique_count++;
        }
    }
    
    return (double)unique_count / 256.0;
}

// Position-weighted bit sum
uint32_t position_weighted_bit_sum_c(const bool* bits, size_t num_bits) {
    uint32_t sum = 0;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) sum += i + 1;
    }
    
    return sum;
}

// Count of adjacent bit pairs (00, 01, 10, 11)
void adjacent_bit_pairs_c(const bool* bits, size_t num_bits, uint32_t counts[4]) {
    for (int i = 0; i < 4; i++) counts[i] = 0;
    
    if (num_bits < 2) return;
    
    for (size_t i = 0; i < num_bits - 1; ++i) {
        uint8_t pattern = (bits[i] << 1) | bits[i+1];
        counts[pattern]++;
    }
}

// Run gap ratio - ratio of largest gap between runs to average gap
double run_gap_ratio_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0.0;
    
    size_t* run_positions = (size_t*)malloc(num_bits * sizeof(size_t));
    if (!run_positions) return 0.0;
    
    size_t run_count = 0;
    bool current_bit = bits[0];
    run_positions[run_count++] = 0;
    
    for (size_t i = 1; i < num_bits; ++i) {
        if (bits[i] != current_bit) {
            current_bit = bits[i];
            run_positions[run_count++] = i;
        }
    }
    
    if (run_count <= 1) {
        free(run_positions);
        return 0.0;
    }
    
    size_t max_gap = 0;
    double total_gap = 0.0;
    
    for (size_t i = 1; i < run_count; ++i) {
        size_t gap = run_positions[i] - run_positions[i-1];
        if (gap > max_gap) max_gap = gap;
        total_gap += gap;
    }
    
    double avg_gap = total_gap / (run_count - 1);
    free(run_positions);
    
    return (avg_gap == 0) ? 0.0 : (double)max_gap / avg_gap;
}

// Morreau-Broto autocorrelation
double morreau_broto_autocorr_c(const bool* bits, size_t num_bits, size_t lag) {
    if (num_bits <= lag) return 0.0;
    
    double sum = 0.0;
    size_t count = num_bits - lag;
    
    for (size_t i = 0; i < count; ++i) {
        sum += (bits[i] ? 1.0 : -1.0) * (bits[i + lag] ? 1.0 : -1.0);
    }
    
    return sum / count;
}

// Bit interleaving difference
uint32_t bit_interleaving_diff_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes <= 1) return 0;
    
    uint32_t odd_popcount = 0;
    uint32_t even_popcount = 0;
    
    for (size_t i = 0; i < num_bytes; ++i) {
        for (int j = 0; j < 8; ++j) {
            bool bit = (bytes[i] >> (7 - j)) & 1;
            if ((i*8 + j) % 2 == 0) {
                even_popcount += bit ? 1 : 0;
            } else {
                odd_popcount += bit ? 1 : 0;
            }
        }
    }
    
    return (even_popcount > odd_popcount) ? 
           (even_popcount - odd_popcount) : (odd_popcount - even_popcount);
}

// First moment of bit distribution
double bit_first_moment_c(const bool* bits, size_t num_bits) {
    if (num_bits == 0) return 0.0;
    
    double sum = 0.0;
    double total_bits = 0.0;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            sum += i;
            total_bits += 1.0;
        }
    }
    
    return (total_bits == 0) ? 0.0 : sum / total_bits;
}

// Second moment (variance) of bit distribution
double bit_second_moment_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0.0;
    
    double mean = bit_first_moment_c(bits, num_bits);
    double sum_sq_diff = 0.0;
    double total_bits = 0.0;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            double diff = i - mean;
            sum_sq_diff += diff * diff;
            total_bits += 1.0;
        }
    }
    
    return (total_bits <= 1) ? 0.0 : sum_sq_diff / total_bits;
}

// Bit skewness - third moment of bit distribution
double bit_skewness_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 2) return 0.0;
    
    double mean = bit_first_moment_c(bits, num_bits);
    double variance = bit_second_moment_c(bits, num_bits);
    if (variance == 0.0) return 0.0;
    
    double std_dev = sqrt(variance);
    double sum_cube_diff = 0.0;
    double total_bits = 0.0;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            double diff = (i - mean) / std_dev;
            sum_cube_diff += diff * diff * diff;
            total_bits += 1.0;
        }
    }
    
    return (total_bits <= 2) ? 0.0 : sum_cube_diff / total_bits;
}

// Bit kurtosis - fourth moment of bit distribution
double bit_kurtosis_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 3) return 0.0;
    
    double mean = bit_first_moment_c(bits, num_bits);
    double variance = bit_second_moment_c(bits, num_bits);
    if (variance == 0.0) return 0.0;
    
    double std_dev = sqrt(variance);
    double sum_fourth_diff = 0.0;
    double total_bits = 0.0;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            double diff = (i - mean) / std_dev;
            double diff_sq = diff * diff;
            sum_fourth_diff += diff_sq * diff_sq;
            total_bits += 1.0;
        }
    }
    
    return (total_bits <= 3) ? 0.0 : (sum_fourth_diff / total_bits) - 3.0; // Excess kurtosis
}

// Bit segregation index - measures bit clustering by position
double bit_segregation_index_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0.0;
    
    uint32_t ones = 0;
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) ones++;
    }
    
    if (ones == 0 || ones == num_bits) return 0.0;
    
    double p_one = (double)ones / num_bits;
    double p_zero = 1.0 - p_one;
    double expected_different = 2.0 * p_one * p_zero * (num_bits - 1);
    
    uint32_t different_neighbors = 0;
    for (size_t i = 1; i < num_bits; ++i) {
        if (bits[i] != bits[i-1]) different_neighbors++;
    }
    
    return 1.0 - ((double)different_neighbors / expected_different);
}

// Fuzzy bit entropy - entropy with transition probabilities
double fuzzy_bit_entropy_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0.0;
    
    uint32_t counts[4] = {0}; // 00, 01, 10, 11
    for (size_t i = 0; i < num_bits - 1; ++i) {
        uint8_t pattern = (bits[i] << 1) | bits[i+1];
        counts[pattern]++;
    }
    
    double total = num_bits - 1;
    double entropy = 0.0;
    
    for (int i = 0; i < 4; ++i) {
        if (counts[i] > 0) {
            double p = counts[i] / total;
            entropy -= p * log2(p);
        }
    }
    
    return entropy;
}

// Shannon equitability - normalized entropy
double shannon_equitability_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes == 0) return 0.0;
    
    uint32_t histogram[256] = {0};
    for (size_t i = 0; i < num_bytes; ++i) {
        histogram[bytes[i]]++;
    }
    
    double entropy = calculate_entropy_c(histogram, num_bytes);
    double max_entropy = log2(256.0); // Maximum possible entropy with 256 symbols
    
    return entropy / max_entropy;
}

// Minimum distance between set bits
uint32_t min_distance_between_set_bits_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return UINT32_MAX;
    
    uint32_t min_distance = UINT32_MAX;
    int64_t last_set_pos = -1;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            if (last_set_pos >= 0) {
                uint32_t distance = i - last_set_pos;
                if (distance < min_distance) min_distance = distance;
            }
            last_set_pos = i;
        }
    }
    
    return (min_distance == UINT32_MAX) ? 0 : min_distance;
}

// Maximum distance between set bits
uint32_t max_distance_between_set_bits_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0;
    
    uint32_t max_distance = 0;
    int64_t last_set_pos = -1;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            if (last_set_pos >= 0) {
                uint32_t distance = i - last_set_pos;
                if (distance > max_distance) max_distance = distance;
            }
            last_set_pos = i;
        }
    }
    
    return max_distance;
}

// Median position of set bits
double median_position_of_set_bits_c(const bool* bits, size_t num_bits) {
    if (num_bits == 0) return 0.0;
    
    uint32_t ones = 0;
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) ones++;
    }
    
    if (ones == 0) return 0.0;
    
    size_t* positions = (size_t*)malloc(ones * sizeof(size_t));
    if (!positions) return 0.0;
    
    size_t pos_idx = 0;
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) positions[pos_idx++] = i;
    }
    
    double median;
    if (ones % 2 == 0) {
        median = (positions[ones/2 - 1] + positions[ones/2]) / 2.0;
    } else {
        median = positions[ones/2];
    }
    
    free(positions);
    return median;
}

// Moving average of bit density
double moving_avg_bit_density_c(const bool* bits, size_t num_bits, size_t window_size) {
    if (num_bits == 0 || window_size == 0 || window_size > num_bits) return 0.0;
    
    double max_density = 0.0;
    
    for (size_t i = 0; i <= num_bits - window_size; ++i) {
        uint32_t count = 0;
        for (size_t j = 0; j < window_size; ++j) {
            if (bits[i + j]) count++;
        }
        
        double density = (double)count / window_size;
        if (density > max_density) max_density = density;
    }
    
    return max_density;
}

// Bit avalanche factor - measures bit independence
double bit_avalanche_factor_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes <= 1) return 0.0;
    
    uint32_t total_flips = 0;
    uint32_t bit_pairs = 0;
    
    for (size_t i = 0; i < num_bytes; ++i) {
        for (int j = 0; j < 8; ++j) {
            uint8_t modified[num_bytes];
            memcpy(modified, bytes, num_bytes);
            
            // Flip this bit
            modified[i] ^= (1 << (7 - j));
            
            // Count how many other bits change in a simple hash
            uint8_t orig_hash = crc8_c(bytes, num_bytes);
            uint8_t mod_hash = crc8_c(modified, num_bytes);
            
            total_flips += POPCOUNT_LUT[orig_hash ^ mod_hash];
            bit_pairs++;
        }
    }
    
    return (bit_pairs == 0) ? 0.0 : (double)total_flips / (bit_pairs * 8);
}

// Approximate entropy - predictability measure
double approx_entropy_c(const bool* bits, size_t num_bits, size_t m) {
    if (num_bits <= m + 1) return 0.0;
    
    size_t count1 = num_bits - m + 1;
    size_t count2 = num_bits - m;
    
    double phi1 = 0.0;
    double phi2 = 0.0;
    
    // For pattern length m
    for (size_t i = 0; i < count1; ++i) {
        size_t matches = 0;
        for (size_t j = 0; j < count1; ++j) {
            bool match = true;
            for (size_t k = 0; k < m; ++k) {
                if (bits[i + k] != bits[j + k]) {
                    match = false;
                    break;
                }
            }
            if (match) matches++;
        }
        
        phi1 += log((double)matches / count1);
    }
    phi1 /= count1;
    
    // For pattern length m+1
    for (size_t i = 0; i < count2; ++i) {
        size_t matches = 0;
        for (size_t j = 0; j < count2; ++j) {
            bool match = true;
            for (size_t k = 0; k < m + 1; ++k) {
                if (bits[i + k] != bits[j + k]) {
                    match = false;
                    break;
                }
            }
            if (match) matches++;
        }
        
        phi2 += log((double)matches / count2);
    }
    phi2 /= count2;
    
    return phi1 - phi2;
}

// Bit quartet distribution - patterns of 4 consecutive bits
void bit_quartet_distribution_c(const bool* bits, size_t num_bits, uint32_t counts[16]) {
    for (int i = 0; i < 16; i++) counts[i] = 0;
    
    if (num_bits < 4) return;
    
    for (size_t i = 0; i <= num_bits - 4; ++i) {
        uint8_t pattern = 0;
        for (size_t j = 0; j < 4; ++j) {
            pattern = (pattern << 1) | (bits[i + j] ? 1 : 0);
        }
        counts[pattern]++;
    }
}

// Bit run count - number of runs of consecutive similar bits
uint32_t bit_run_count_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return num_bits;
    
    uint32_t runs = 1;
    for (size_t i = 1; i < num_bits; ++i) {
        if (bits[i] != bits[i-1]) runs++;
    }
    
    return runs;
}

// Bit run entropy - entropy of run lengths distribution
double bit_run_entropy_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0.0;
    
    uint32_t* run_lengths = (uint32_t*)malloc(num_bits * sizeof(uint32_t));
    if (!run_lengths) return 0.0;
    
    uint32_t run_count = 0;
    uint32_t current_length = 1;
    bool current_bit = bits[0];
    
    for (size_t i = 1; i < num_bits; ++i) {
        if (bits[i] == current_bit) {
            current_length++;
        } else {
            run_lengths[run_count++] = current_length;
            current_bit = bits[i];
            current_length = 1;
        }
    }
    run_lengths[run_count++] = current_length;
    
    // Calculate frequency of each run length
    uint32_t max_length = 0;
    for (size_t i = 0; i < run_count; ++i) {
        if (run_lengths[i] > max_length) max_length = run_lengths[i];
    }
    
    uint32_t* length_hist = (uint32_t*)calloc(max_length + 1, sizeof(uint32_t));
    if (!length_hist) {
        free(run_lengths);
        return 0.0;
    }
    
    for (size_t i = 0; i < run_count; ++i) {
        length_hist[run_lengths[i]]++;
    }
    
    double entropy = calculate_entropy_c(length_hist, run_count);
    
    free(run_lengths);
    free(length_hist);
    return entropy;
}

// Bit burst factor - ratio of longest run to average run
double bit_burst_factor_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 1.0;
    
    uint32_t longest_run = longest_consecutive_bit_run_c(bits, num_bits);
    double avg_run = average_bit_run_length_c(bits, num_bits);
    
    return (avg_run == 0.0) ? 0.0 : (double)longest_run / avg_run;
}

// Spectral flatness - Wiener entropy measure
double spectral_flatness_c(const bool* bits, size_t num_bits) {
    if (num_bits < 2) return 1.0; // By definition, a single bit has perfect flatness
    
    // Compute simple frequency spectrum by counting transitions at different lags
    size_t max_lag = num_bits / 2;
    double* spectrum = (double*)malloc(max_lag * sizeof(double));
    if (!spectrum) return 0.0;
    
    for (size_t lag = 1; lag <= max_lag; ++lag) {
        size_t matches = 0;
        size_t comparisons = num_bits - lag;
        
        for (size_t i = 0; i < comparisons; ++i) {
            if (bits[i] == bits[i + lag]) matches++;
        }
        
        spectrum[lag-1] = (double)matches / comparisons;
    }
    
    // Calculate geometric mean and arithmetic mean
    double log_sum = 0.0;
    double sum = 0.0;
    
    for (size_t i = 0; i < max_lag; ++i) {
        if (spectrum[i] > 0.0) {
            log_sum += log(spectrum[i]);
        } else {
            log_sum -= 10.0; // Substitute small value for log calculation
        }
        sum += spectrum[i];
    }
    
    double geo_mean = exp(log_sum / max_lag);
    double arith_mean = sum / max_lag;
    
    free(spectrum);
    return (arith_mean == 0.0) ? 0.0 : geo_mean / arith_mean;
}

// Central position-weighted bits sum
int32_t central_weighted_bit_sum_c(const bool* bits, size_t num_bits) {
    if (num_bits == 0) return 0;
    
    int32_t sum = 0;
    int64_t mid = num_bits / 2;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            // Weight is inverse of distance from center
            int64_t dist = labs((int64_t)i - mid);
            sum += (dist == 0) ? num_bits : (int32_t)((double)num_bits / (dist + 1));
        }
    }
    
    return sum;
}

// Bit asymmetry - difference between bit patterns in first/second half
double bit_asymmetry_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0.0;
    
    size_t mid = num_bits / 2;
    uint32_t different_bits = 0;
    
    for (size_t i = 0; i < mid; ++i) {
        if (bits[i] != bits[num_bits - 1 - i]) different_bits++;
    }
    
    return (double)different_bits / mid;
}

// Bit pattern persistence - how long patterns repeat
double bit_pattern_persistence_c(const bool* bits, size_t num_bits, size_t pattern_size) {
    if (num_bits < 2 * pattern_size) return 0.0;
    
    uint32_t max_repeat = 0;
    
    for (size_t i = 0; i <= num_bits - pattern_size; ++i) {
        uint32_t repeat_count = 0;
        
        for (size_t j = i + pattern_size; j <= num_bits - pattern_size; j += pattern_size) {
            bool match = true;
            for (size_t k = 0; k < pattern_size; ++k) {
                if (bits[i + k] != bits[j + k]) {
                    match = false;
                    break;
                }
            }
            
            if (match) {
                repeat_count++;
            } else {
                break;
            }
        }
        
        if (repeat_count > max_repeat) max_repeat = repeat_count;
    }
    
    return max_repeat;
}

// Average bit positional spread
double avg_bit_positional_spread_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0.0;
    
    uint32_t ones = 0;
    double sum_positions = 0.0;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            ones++;
            sum_positions += i;
        }
    }
    
    if (ones <= 1) return 0.0;
    
    double mean_pos = sum_positions / ones;
    double sum_squared_dev = 0.0;
    
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            double dev = i - mean_pos;
            sum_squared_dev += dev * dev;
        }
    }
    
    return sqrt(sum_squared_dev / ones);
}

// Bit pattern Lempel-Ziv complexity
uint32_t lempel_ziv_complexity_c(const bool* bits, size_t num_bits) {
    if (num_bits == 0) return 0;
    
    uint32_t complexity = 1;
    size_t i = 0, j, k, l;
    
    for (i = 1; i < num_bits; i++) {
        j = 0;
        k = i;
        l = 1;
        
        while (k < num_bits && j < i) {
            if (bits[k] == bits[j]) {
                k++;
                j++;
                l++;
            } else {
                j = j - l + 1;
                l = 1;
                if (j == i) {
                    complexity++;
                    break;
                }
            }
            
            if (j == i) {
                if (k == num_bits) break;
                j = 0;
                i = k;
                complexity++;
                break;
            }
        }
    }
    
    return complexity;
}

// Kolmogorov complexity estimation (using compress-and-count)
double kolmogorov_complexity_estimate_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return (double)num_bits;
    
    // Use Lempel-Ziv as a proxy for Kolmogorov complexity
    uint32_t lz_complexity = lempel_ziv_complexity_c(bits, num_bits);
    
    // Normalize by theoretical maximum complexity (log base 2)
    return (double)lz_complexity / ((double)num_bits / log2((double)num_bits));
}

// Modified Shannon entropy - uses overlapping bit sequences
double modified_shannon_entropy_c(const bool* bits, size_t num_bits, size_t block_size) {
    if (num_bits < block_size) return 0.0;
    
    size_t num_patterns = 1 << block_size;
    uint32_t* pattern_count = (uint32_t*)calloc(num_patterns, sizeof(uint32_t));
    if (!pattern_count) return 0.0;
    
    size_t window_count = num_bits - block_size + 1;
    
    for (size_t i = 0; i < window_count; ++i) {
        uint32_t pattern = 0;
        for (size_t j = 0; j < block_size; ++j) {
            pattern = (pattern << 1) | (bits[i + j] ? 1 : 0);
        }
        pattern_count[pattern]++;
    }
    
    double entropy = 0.0;
    for (size_t i = 0; i < num_patterns; ++i) {
        if (pattern_count[i] > 0) {
            double prob = (double)pattern_count[i] / window_count;
            entropy -= prob * log2(prob);
        }
    }
    
    free(pattern_count);
    return entropy / block_size; // Normalize by block size
}

// Bit sequence predictability
double bit_sequence_predictability_c(const bool* bits, size_t num_bits, size_t context_len) {
    if (num_bits <= context_len) return 0.5; // Default to unpredictable
    
    // Build context-based prediction model
    size_t num_contexts = 1 << context_len;
    uint32_t* zeros_after = (uint32_t*)calloc(num_contexts, sizeof(uint32_t));
    uint32_t* ones_after = (uint32_t*)calloc(num_contexts, sizeof(uint32_t));
    
    if (!zeros_after || !ones_after) {
        free(zeros_after);
        free(ones_after);
        return 0.5;
    }
    
    // Count occurrences for each context
    for (size_t i = 0; i < num_bits - context_len; ++i) {
        uint32_t context = 0;
        for (size_t j = 0; j < context_len; ++j) {
            context = (context << 1) | (bits[i + j] ? 1 : 0);
        }
        
        if (bits[i + context_len]) {
            ones_after[context]++;
        } else {
            zeros_after[context]++;
        }
    }
    
    // Calculate predictability as average certainty across contexts
    double total_certainty = 0.0;
    uint32_t active_contexts = 0;
    
    for (size_t i = 0; i < num_contexts; ++i) {
        uint32_t total = zeros_after[i] + ones_after[i];
        if (total > 0) {
            double p_zero = (double)zeros_after[i] / total;
            double p_one = (double)ones_after[i] / total;
            double certainty = (p_zero > p_one) ? p_zero : p_one;
            total_certainty += certainty;
            active_contexts++;
        }
    }
    
    free(zeros_after);
    free(ones_after);
    
    return (active_contexts == 0) ? 0.5 : total_certainty / active_contexts;
}

// Compression efficiency - ratio of compressed to original size
double compression_efficiency_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes == 0) return 0.0;
    
    // Simple RLE compression
    size_t compressed_size = 0;
    uint8_t current = bytes[0];
    uint8_t count = 1;
    
    for (size_t i = 1; i < num_bytes; ++i) {
        if (bytes[i] == current && count < 255) {
            count++;
        } else {
            compressed_size += 2; // Store (count, byte)
            current = bytes[i];
            count = 1;
        }
    }
    compressed_size += 2; // Last run
    
    return (double)compressed_size / num_bytes;
}

// NIST bit frequency test statistic
double nist_frequency_test_c(const bool* bits, size_t num_bits) {
    if (num_bits == 0) return 0.0;
    
    int32_t sum = 0;
    for (size_t i = 0; i < num_bits; ++i) {
        sum += bits[i] ? 1 : -1;
    }
    
    double s_obs = (double)abs(sum) / sqrt((double)num_bits);
    return s_obs;
}

// NIST runs test statistic
double nist_runs_test_c(const bool* bits, size_t num_bits) {
    if (num_bits <= 1) return 0.0;
    
    // Count ones and compute pi
    uint32_t ones = 0;
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) ones++;
    }
    
    double pi = (double)ones / num_bits;
    
    // Skip if pi is out of reasonable range
    if (fabs(pi - 0.5) >= (2.0 / sqrt(num_bits))) return 0.0;
    
    // Count runs
    uint32_t runs = 1;
    for (size_t i = 1; i < num_bits; ++i) {
        if (bits[i] != bits[i-1]) runs++;
    }
    
    // Calculate test statistic
    double runs_exp = 2.0 * num_bits * pi * (1.0 - pi);
    double std_dev = sqrt(2.0 * num_bits * pi * (1.0 - pi) * (2.0 * pi - 1.0));
    double v_obs = fabs(runs - runs_exp) / std_dev;
    
    return v_obs;
}

// Hamming weight difference from expected
double hamming_weight_deviation_c(const uint8_t* bytes, size_t num_bytes) {
    if (num_bytes == 0) return 0.0;
    
    uint32_t bit_count = popcount_bytes_c(bytes, num_bytes);
    uint32_t total_bits = num_bytes * 8;
    uint32_t expected = total_bits / 2; // Expect half the bits to be set
    
    return fabs((double)(bit_count - expected)) / total_bits;
}

// Count of specific byte values (e.g., count of zero bytes)
uint32_t specific_byte_count_c(const uint8_t* bytes, size_t num_bytes, uint8_t target) {
    uint32_t count = 0;
    for (size_t i = 0; i < num_bytes; ++i) {
        if (bytes[i] == target) count++;
    }
    return count;
}
