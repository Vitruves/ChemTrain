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
