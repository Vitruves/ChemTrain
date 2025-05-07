// src/descriptors/bitsoperations.c
#include "../cregistry.h"
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

double calculateStringLength(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return smiles ? (double)strlen(smiles) : 0.0;
}

double calculateBitCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes) return 0.0; // Memory allocation failed
    double result = (double)popcount_bytes_c(bytes, num_bytes);
    free(bytes);
    return result;
}

double calculateBitDensity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes || num_bytes == 0) {
        if(bytes) free(bytes);
        return 0.0;
    }
    double total_bits = (double)popcount_bytes_c(bytes, num_bytes);
    free(bytes);
    return total_bits / num_bytes;
}

double calculateAsciiSum(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    uint32_t sum = 0;
    for (size_t i = 0; smiles[i] != '\0'; ++i) {
        sum += (uint8_t)smiles[i];
    }
    return (double)sum;
}

double calculateLongestBitRun(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits) return 0.0;
    double result = (double)longest_consecutive_bit_run_c(bits, num_bits);
    free(bits);
    return result;
}

double calculateByteXorChecksum(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes) return 0.0;
    double result = (double)byte_xor_checksum_c(bytes, num_bytes);
    free(bytes);
    return result;
}

double calculateByteEntropy(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes || num_bytes == 0) {
         if(bytes) free(bytes);
         return 0.0;
    }
    uint32_t histogram[256] = {0};
    for (size_t i = 0; i < num_bytes; ++i) {
        histogram[bytes[i]]++;
    }
    double result = calculate_entropy_c(histogram, num_bytes);
    free(bytes);
    return result;
}

double calculateFirstByteValue(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    return (smiles && smiles[0] != '\0') ? (double)((uint8_t)smiles[0]) : 0.0;
}

double calculateLastByteValue(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t len = strlen(smiles);
    return (len > 0) ? (double)((uint8_t)smiles[len - 1]) : 0.0;
}

double calculateOddEvenBitRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits) return 0.0;
    uint32_t oddCount = 0, evenCount = 0;
    for (size_t i = 0; i < num_bits; ++i) {
        if (bits[i]) {
            if (i % 2 == 0) evenCount++;
            else oddCount++;
        }
    }
    free(bits);
    return (evenCount == 0) ? 0.0 : (double)oddCount / evenCount;
}

double calculateCharPalindromeCheck(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes) return 0.0;
    double result = is_char_palindrome_c(bytes, num_bytes) ? 1.0 : 0.0;
    free(bytes);
    return result;
}

double calculateBitTransitionCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits) return 0.0;
    double result = (double)bit_transition_count_c(bits, num_bits);
    free(bits);
    return result;
}

double calculateLeftRightBitBalance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits) return 0.0;
    double result = (double)left_right_bit_balance_c(bits, num_bits);
    free(bits);
    return result;
}

double calculateHammingWeightVariance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes) return 0.0;
    double result = bit_distribution_variance_c(bytes, num_bytes);
    free(bytes);
    return result;
}

double calculateFastHash(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes) return 0.0;
    double result = (double)fnv1a_hash_c(bytes, num_bytes);
    free(bytes);
    return result;
}

double calculateCrc8Checksum(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes) return 0.0;
    double result = (double)crc8_c(bytes, num_bytes);
    free(bytes);
    return result;
}

double calculateMostFrequentByte(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes || num_bytes == 0) {
        if(bytes) free(bytes);
        return 0.0;
    }
    uint32_t histogram[256] = {0};
    for (size_t i = 0; i < num_bytes; ++i) histogram[bytes[i]]++;
    uint8_t mostFrequent = 0;
    uint32_t maxCount = 0;
    for (int i = 0; i < 256; ++i) {
        if (histogram[i] > maxCount) {
            maxCount = histogram[i];
            mostFrequent = (uint8_t)i;
        }
    }
    free(bytes);
    return (double)mostFrequent;
}

double calculateLeastFrequentByte(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
     if (!bytes || num_bytes == 0) {
        if(bytes) free(bytes);
        return 0.0;
    }
    uint32_t histogram[256] = {0};
    for (size_t i = 0; i < num_bytes; ++i) histogram[bytes[i]]++;
    uint8_t leastFrequent = 0;
    uint32_t minCount = UINT32_MAX;
    bool found = false;
    for (int i = 0; i < 256; ++i) {
        if (histogram[i] > 0 && histogram[i] < minCount) {
            minCount = histogram[i];
            leastFrequent = (uint8_t)i;
            found = true;
        }
    }
    free(bytes);
    return found ? (double)leastFrequent : 0.0;
}

double calculateUniqueByteCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes) return 0.0;
    uint32_t histogram[256] = {0};
    for (size_t i = 0; i < num_bytes; ++i) histogram[bytes[i]]++;
    uint32_t uniqueCount = 0;
    for (int i = 0; i < 256; ++i) if (histogram[i] > 0) uniqueCount++;
    free(bytes);
    return (double)uniqueCount;
}

double calculateEvenIndexedBitCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits) return 0.0;
    uint32_t count = 0;
    for (size_t i = 0; i < num_bits; i += 2) if (bits[i]) count++;
    free(bits);
    return (double)count;
}

double calculateOddIndexedBitCount(const void* context, GetSmilesFunc getSmilesFunc) {
     const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits) return 0.0;
    uint32_t count = 0;
    for (size_t i = 1; i < num_bits; i += 2) if (bits[i]) count++;
    free(bits);
    return (double)count;
}

double calculateByteWiseMinimum(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles || smiles[0] == '\0') return 0.0;
    uint8_t minVal = (uint8_t)smiles[0];
    for (size_t i = 1; smiles[i] != '\0'; ++i) {
        if ((uint8_t)smiles[i] < minVal) minVal = (uint8_t)smiles[i];
    }
    return (double)minVal;
}

double calculateByteWiseMaximum(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles || smiles[0] == '\0') return 0.0;
    uint8_t maxVal = (uint8_t)smiles[0];
    for (size_t i = 1; smiles[i] != '\0'; ++i) {
        if ((uint8_t)smiles[i] > maxVal) maxVal = (uint8_t)smiles[i];
    }
    return (double)maxVal;
}

double calculateByteWiseMedian(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes || num_bytes == 0) {
        if(bytes) free(bytes);
        return 0.0;
    }
    qsort(bytes, num_bytes, sizeof(uint8_t), compare_uint8);
    double median;
    if (num_bytes % 2 == 0) {
        size_t mid1 = num_bytes / 2 - 1;
        size_t mid2 = num_bytes / 2;
        median = ((double)bytes[mid1] + (double)bytes[mid2]) / 2.0;
    } else {
        size_t mid = num_bytes / 2;
        median = (double)bytes[mid];
    }
    free(bytes);
    return median;
}

double calculateHighLowBitRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits || num_bits == 0) {
        if(bits) free(bits);
        return 0.0;
    }
    uint32_t highCount = 0;
    for (size_t i = 0; i < num_bits; ++i) if (bits[i]) highCount++;
    uint32_t lowCount = num_bits - highCount;
    free(bits);
    return (lowCount == 0) ? (double)num_bits : (double)highCount / lowCount;
}

double calculateCharSequenceComplexity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t len = strlen(smiles);
    if (len == 0) return 0.0;
    uint8_t seen[256] = {0};
    uint32_t uniqueCount = 0;
    for (size_t i = 0; i < len; ++i) {
        if (seen[(uint8_t)smiles[i]] == 0) {
            seen[(uint8_t)smiles[i]] = 1;
            uniqueCount++;
        }
    }
    return (double)uniqueCount / len;
}

double calculateAvgBitPopulationCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes || num_bytes == 0) {
         if(bytes) free(bytes);
         return 0.0;
    }
    double totalPopcount = (double)popcount_bytes_c(bytes, num_bytes);
    free(bytes);
    return totalPopcount / num_bytes;
}

double calculateAvgLeadingZeroCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
     if (!bytes || num_bytes == 0) {
        if(bytes) free(bytes);
        return 0.0;
    }
    uint32_t totalLeadingZeros = 0;
    for (size_t i = 0; i < num_bytes; ++i) totalLeadingZeros += LEADING_ZEROS_LUT[bytes[i]];
    free(bytes);
    return (double)totalLeadingZeros / num_bytes;
}

double calculateAvgTrailingZeroCount(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes || num_bytes == 0) {
        if(bytes) free(bytes);
        return 0.0;
    }
    uint32_t totalTrailingZeros = 0;
    for (size_t i = 0; i < num_bytes; ++i) totalTrailingZeros += TRAILING_ZEROS_LUT[bytes[i]];
    free(bytes);
    return (double)totalTrailingZeros / num_bytes;
}

double calculatePopcountVariance(const void* context, GetSmilesFunc getSmilesFunc) {
    // This is the same calculation as HammingWeightVariance
    return calculateHammingWeightVariance(context, getSmilesFunc);
}

double calculateByteAndReduction(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles || smiles[0] == '\0') return 0.0;
    uint8_t result = (uint8_t)smiles[0];
    for (size_t i = 1; smiles[i] != '\0'; ++i) result &= (uint8_t)smiles[i];
    return (double)result;
}

double calculateByteOrReduction(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles || smiles[0] == '\0') return 0.0;
    uint8_t result = (uint8_t)smiles[0];
    for (size_t i = 1; smiles[i] != '\0'; ++i) result |= (uint8_t)smiles[i];
    return (double)result;
}

double calculateByteXorReduction(const void* context, GetSmilesFunc getSmilesFunc) {
    // This is the same calculation as ByteXorChecksum
    return calculateByteXorChecksum(context, getSmilesFunc);
}

double calculateNonZeroByteDensity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t len = strlen(smiles);
    if (len == 0) return 0.0;
    uint32_t nonZeroCount = 0;
    for (size_t i = 0; i < len; ++i) if (smiles[i] != '\0') nonZeroCount++; // Should always be true unless embedded nulls
    return (double)nonZeroCount / len; // Simplified: Assumes no embedded nulls
}

double calculateAsciiProductMod256(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    uint32_t product = 1;
    for (size_t i = 0; smiles[i] != '\0'; ++i) {
        product = (product * (uint8_t)smiles[i]) % 256;
    }
    return (double)product;
}

double calculateBitAlternationFreq(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits || num_bits <= 1) {
        if(bits) free(bits);
        return 0.0;
    }
    uint32_t transitions = bit_transition_count_c(bits, num_bits);
    free(bits);
    return (double)transitions / (num_bits - 1);
}

double calculateCharDistributionVariance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes = strlen(smiles);
     if (num_bytes == 0) return 0.0;
    uint32_t histogram[256] = {0};
    for (size_t i = 0; i < num_bytes; ++i) histogram[(uint8_t)smiles[i]]++;
    double nonZeroFreqs[256];
    int nonZeroCount = 0;
    for(int i=0; i<256; ++i) {
        if (histogram[i] > 0) {
            nonZeroFreqs[nonZeroCount++] = (double)histogram[i] / num_bytes;
        }
    }
    if (nonZeroCount == 0) return 0.0;
    double mean = 1.0 / nonZeroCount;
    double variance = 0.0;
    for (int i=0; i<nonZeroCount; ++i) {
        double diff = nonZeroFreqs[i] - mean;
        variance += diff * diff;
    }
    return variance / nonZeroCount;
}

double calculateLowHighNibbleRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes = strlen(smiles);
    if (num_bytes == 0) return 0.0;
    uint32_t lowCount = 0, highCount = 0;
    for (size_t i = 0; i < num_bytes; ++i) {
        uint8_t byte = (uint8_t)smiles[i];
        uint8_t highNibble = (byte >> 4) & 0x0F;
        uint8_t lowNibble = byte & 0x0F;
        if (highNibble < 8) lowCount++; else highCount++;
        if (lowNibble < 8) lowCount++; else highCount++;
    }
    return (highCount == 0) ? (double)(num_bytes * 2) : (double)lowCount / highCount;
}

double calculateCharPosWeightedSum(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    double weightedSum = 0.0;
    for (size_t i = 0; smiles[i] != '\0'; ++i) {
        weightedSum += (double)((uint8_t)smiles[i]) * (i + 1);
    }
    return weightedSum;
}

double calculateNibbleEntropy(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes = strlen(smiles);
    if (num_bytes == 0) return 0.0;
    uint32_t nibbleHistogram[16] = {0};
    size_t totalNibbles = num_bytes * 2;
     for (size_t i = 0; i < num_bytes; ++i) {
        uint8_t byte = (uint8_t)smiles[i];
        nibbleHistogram[(byte >> 4) & 0x0F]++;
        nibbleHistogram[byte & 0x0F]++;
    }
    double entropy = 0.0;
    for (int i = 0; i < 16; ++i) {
        if (nibbleHistogram[i] > 0) {
            double probability = (double)nibbleHistogram[i] / totalNibbles;
            entropy -= probability * log2(probability);
        }
    }
    return entropy;
}

double calculateCharFrequencyEntropy(const void* context, GetSmilesFunc getSmilesFunc) {
    // Same as ByteEntropy
    return calculateByteEntropy(context, getSmilesFunc);
}

double calculateAsciiRangeCoverage(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles || smiles[0] == '\0') return 0.0;
    uint8_t minAscii = (uint8_t)smiles[0], maxAscii = (uint8_t)smiles[0];
    for (size_t i = 1; smiles[i] != '\0'; ++i) {
        uint8_t current = (uint8_t)smiles[i];
        if (current < minAscii) minAscii = current;
        if (current > maxAscii) maxAscii = current;
    }
    return (double)(maxAscii - minAscii + 1) / 256.0;
}

double calculateNibbleHistogramBalance(const void* context, GetSmilesFunc getSmilesFunc) {
     const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes = strlen(smiles);
    if (num_bytes == 0) return 0.0;
    uint32_t nibbleHistogram[16] = {0};
    size_t totalNibbles = num_bytes * 2;
     for (size_t i = 0; i < num_bytes; ++i) {
        uint8_t byte = (uint8_t)smiles[i];
        nibbleHistogram[(byte >> 4) & 0x0F]++;
        nibbleHistogram[byte & 0x0F]++;
    }
    double nonZeroFreqs[16];
    int nonZeroCount = 0;
    for(int i=0; i<16; ++i) {
        if (nibbleHistogram[i] > 0) {
            nonZeroFreqs[nonZeroCount++] = (double)nibbleHistogram[i] / totalNibbles;
        }
    }
    if (nonZeroCount <= 1) return 0.0;
    double mean = 1.0 / nonZeroCount;
    double variance = 0.0;
    for (int i = 0; i < nonZeroCount; ++i) {
        double diff = nonZeroFreqs[i] - mean;
        variance += diff * diff;
    }
    double stdDev = sqrt(variance / nonZeroCount);
    double cv = stdDev / mean;
    return 1.0 / (1.0 + cv);
}

double calculateBinaryCompressionRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t len = strlen(smiles);
    if (len < 2) return 1.0;
    double byteEntropy = calculateByteEntropy(context, getSmilesFunc);
    double normalizedEntropy = byteEntropy / 8.0;
    return 1.0 - (normalizedEntropy * 0.9);
}

double calculateBitReversalDistance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits || num_bits <= 1) {
        if(bits) free(bits);
        return 0.0;
    }
    uint32_t distance = 0;
    for (size_t i = 0; i < num_bits / 2; ++i) {
        if (bits[i] != bits[num_bits - 1 - i]) distance++;
    }
    free(bits);
    return (double)distance;
}

double calculateByteReversalDistance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes || num_bytes <= 1) {
         if(bytes) free(bytes);
         return 0.0;
    }
    uint32_t distance = 0;
    for (size_t i = 0; i < num_bytes / 2; ++i) {
        if (bytes[i] != bytes[num_bytes - 1 - i]) distance++;
    }
    free(bytes);
    return (double)distance;
}

double calculateMurmurHash3(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t len = strlen(smiles);
    if (len == 0) return 0.0;

    const uint32_t seed = 0x9747b28c;
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
    uint32_t h1 = seed;
    size_t nblocks = len / 4;

    for (size_t i = 0; i < nblocks; ++i) {
        uint32_t k1 = 0;
        memcpy(&k1, smiles + i * 4, 4); // Assumes little-endian, common case
        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> 17);
        k1 *= c2;
        h1 ^= k1;
        h1 = (h1 << 13) | (h1 >> 19);
        h1 = h1 * 5 + 0xe6546b64;
    }

    uint32_t k1 = 0;
    const uint8_t* tail = (const uint8_t*)(smiles + nblocks * 4);
    switch (len & 3) {
        case 3: k1 ^= tail[2] << 16; // Fallthrough
        case 2: k1 ^= tail[1] << 8;  // Fallthrough
        case 1: k1 ^= tail[0];
                k1 *= c1; k1 = (k1 << 15) | (k1 >> 17); k1 *= c2; h1 ^= k1;
    }

    h1 ^= len;
    h1 ^= (h1 >> 16);
    h1 *= 0x85ebca6b;
    h1 ^= (h1 >> 13);
    h1 *= 0xc2b2ae35;
    h1 ^= (h1 >> 16);
    return (double)h1;
}

double calculateCharRunLength(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles || smiles[0] == '\0') return 0.0;
    uint32_t runCount = 0;
    uint32_t totalRunLengthSum = 0;
    char prevChar = smiles[0];
    uint32_t currentRun = 1;
    for (size_t i = 1; smiles[i] != '\0'; ++i) {
        if (smiles[i] == prevChar) {
            currentRun++;
        } else {
            totalRunLengthSum += currentRun;
            runCount++;
            currentRun = 1;
            prevChar = smiles[i];
        }
    }
    totalRunLengthSum += currentRun;
    runCount++;
    return (double)totalRunLengthSum / runCount;
}

// Comparison function for qsort based on counts
typedef struct { uint32_t key; uint32_t value; } CountPair;
int compare_counts(const void* a, const void* b) {
    return ((CountPair*)b)->value - ((CountPair*)a)->value; // Descending order
}

double calculateByteRunEntropy(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
     size_t len = strlen(smiles);
     if (len <= 1) return 0.0;

    // Use a fixed-size array or hash map if keys can be large
    CountPair runLengthHist[256] = {0}; // Assume max run length won't exceed ~255 reasonably
    int hist_size = 0;

    char prevChar = smiles[0];
    uint32_t currentRun = 1;
    for (size_t i = 1; i < len; ++i) {
        if (smiles[i] == prevChar) {
            currentRun++;
        } else {
            bool found = false;
            for(int k=0; k<hist_size; ++k) {
                if (runLengthHist[k].key == currentRun) {
                    runLengthHist[k].value++;
                    found = true;
                    break;
                }
            }
            if (!found && hist_size < 256) {
                 runLengthHist[hist_size].key = currentRun;
                 runLengthHist[hist_size].value = 1;
                 hist_size++;
            }
            currentRun = 1;
            prevChar = smiles[i];
        }
    }
    // Add the last run
     bool found = false;
    for(int k=0; k<hist_size; ++k) {
        if (runLengthHist[k].key == currentRun) {
            runLengthHist[k].value++;
            found = true;
            break;
        }
    }
     if (!found && hist_size < 256) {
         runLengthHist[hist_size].key = currentRun;
         runLengthHist[hist_size].value = 1;
         hist_size++;
    }

    double entropy = 0.0;
    double totalRuns = 0.0;
    for (int i=0; i<hist_size; ++i) totalRuns += runLengthHist[i].value;

    if (totalRuns == 0.0) return 0.0;

    for (int i=0; i<hist_size; ++i) {
        double probability = (double)runLengthHist[i].value / totalRuns;
        entropy -= probability * log2(probability);
    }

    return entropy;
}

double calculateSingleBitPatternMatch(const void* context, GetSmilesFunc getSmilesFunc) {
     const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
     if (!bits || num_bits < 2) {
        if(bits) free(bits);
        return 0.0;
    }
    uint32_t patternCount = 0;
    for (size_t i = 0; i < num_bits - 1; ++i) {
        if (bits[i] && !bits[i+1]) patternCount++;
    }
    free(bits);
    return (double)patternCount / (num_bits - 1);
}

// Simplified BytePairTransitions - just counts unique pairs
double calculateBytePairTransitions(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t len = strlen(smiles);
    if (len < 2) return 0.0;

    // Estimate max possible pairs (256*256) - use a hash map or large array if memory allows
    // For simplicity, use a fixed large array and hope for no collisions or limited range
    #define PAIR_HASH_SIZE 65536
    uint32_t pairCounts[PAIR_HASH_SIZE] = {0};
    size_t totalPairs = len - 1;

    for (size_t i = 0; i < len - 1; ++i) {
        uint16_t pair = ((uint16_t)((uint8_t)smiles[i]) << 8) | (uint8_t)smiles[i+1];
        pairCounts[pair % PAIR_HASH_SIZE]++; // Simple hash
    }

    double entropy = 0.0;
    for(int i=0; i<PAIR_HASH_SIZE; ++i) {
        if(pairCounts[i] > 0) {
             double probability = (double)pairCounts[i] / totalPairs;
             entropy -= probability * log2(probability);
        }
    }
    return entropy;
}

double calculateBitDifferenceSum(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
     if (!bytes || num_bytes < 2) {
        if(bytes) free(bytes);
        return 0.0;
    }
    uint32_t diffSum = 0;
    for (size_t i = 0; i < num_bytes - 1; ++i) {
        diffSum += POPCOUNT_LUT[bytes[i] ^ bytes[i+1]];
    }
    free(bytes);
    return (double)diffSum;
}

double calculateFirstLastByteXor(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t len = strlen(smiles);
    if (len == 0) return 0.0;
    if (len == 1) return (double)((uint8_t)smiles[0]);
    return (double)(((uint8_t)smiles[0]) ^ ((uint8_t)smiles[len-1]));
}

// ByteMarkovProperty - placeholder, complex calculation omitted for C simplicity
double calculateByteMarkovProperty(const void* context, GetSmilesFunc getSmilesFunc) {
    return 0.0; // Placeholder
}

double calculateSwarPopcount(const void* context, GetSmilesFunc getSmilesFunc) {
    // SWAR Popcount implementation is complex and platform-dependent.
    // Reverting to the simpler LUT-based popcount.
    return calculateBitCount(context, getSmilesFunc);
}

// FastBitHistogram - uses nibble balance calculation as proxy
double calculateFastBitHistogram(const void* context, GetSmilesFunc getSmilesFunc) {
    return calculateNibbleHistogramBalance(context, getSmilesFunc);
}

double calculateBitDispersionIndex(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
     if (!bits || num_bits <= 1) {
        if(bits) free(bits);
        return 0.0;
    }
    uint32_t runLengths[num_bits]; // Max possible runs = num_bits
    int runCount = 0;
    bool currentBit = bits[0];
    uint32_t currentRun = 1;
    for (size_t i = 1; i < num_bits; ++i) {
        if (bits[i] == currentBit) {
            currentRun++;
        } else {
            if (runCount < num_bits) runLengths[runCount++] = currentRun;
            currentRun = 1;
            currentBit = bits[i];
        }
    }
     if (runCount < num_bits) runLengths[runCount++] = currentRun;
    free(bits);

    if (runCount == 0) return 0.0;
    double mean = (double)num_bits / runCount;
    double variance = 0.0;
    for (int i=0; i<runCount; ++i) {
        double diff = runLengths[i] - mean;
        variance += diff * diff;
    }
    variance /= runCount;
    return 1.0 / (1.0 + variance);
}

double calculateHalfByteRotationHash(const void* context, GetSmilesFunc getSmilesFunc) {
     const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    uint32_t hash = 0x8F1BBCDC;
    for (size_t i = 0; smiles[i] != '\0'; ++i) {
        uint8_t byte = (uint8_t)smiles[i];
        uint8_t high = (byte >> 4) & 0x0F;
        uint8_t low = byte & 0x0F;
        hash = ((hash << 5) | (hash >> 27)) ^ (high << 4 | low);
        hash = ((hash << 3) | (hash >> 29)) + ((low << 4) | high);
    }
    return (double)hash;
}

double calculateAvalanchePatternSensitivity(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
     if (!bytes || num_bytes == 0) {
        if(bytes) free(bytes);
        return 0.0;
    }
    uint32_t origHash = fnv1a_hash_c(bytes, num_bytes);
    uint8_t* modifiedBytes = (uint8_t*)malloc(num_bytes * sizeof(uint8_t));
     if (!modifiedBytes) { free(bytes); return 0.0; }
    memcpy(modifiedBytes, bytes, num_bytes);
    size_t midIndex = num_bytes / 2;
    modifiedBytes[midIndex] ^= 0x10; // Flip bit 4
    uint32_t modHash = fnv1a_hash_c(modifiedBytes, num_bytes);
    uint32_t hashDiff = origHash ^ modHash;
    uint32_t diffBits = 0;
    for(int i=0; i<4; ++i) diffBits += POPCOUNT_LUT[(hashDiff >> (i*8)) & 0xFF];
    free(bytes);
    free(modifiedBytes);
    return (double)diffBits / 32.0;
}

double calculateSequentialXorDiff(const void* context, GetSmilesFunc getSmilesFunc) {
     const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
     if (!bytes || num_bytes < 2) {
        if(bytes) free(bytes);
        return 0.0;
    }
    uint32_t xorSum = 0;
    for (size_t i = 0; i < num_bytes - 1; ++i) {
        xorSum += POPCOUNT_LUT[bytes[i] ^ bytes[i+1]];
    }
    free(bytes);
    return (double)xorSum / (num_bytes - 1);
}

double calculateByteDeltaEncodingSize(const void* context, GetSmilesFunc getSmilesFunc) {
     const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
     if (!bytes || num_bytes < 2) {
        if(bytes) free(bytes);
        return num_bytes;
    }
    uint32_t totalBits = 8;
    for (size_t i = 1; i < num_bytes; ++i) {
        int8_t delta = (int8_t)(bytes[i] - bytes[i-1]);
        if (delta == 0) totalBits += 1;
        else if (delta >= -4 && delta <= 3) totalBits += 3;
        else if (delta >= -16 && delta <= 15) totalBits += 5;
        else if (delta >= -64 && delta <= 63) totalBits += 7;
        else totalBits += 9;
    }
    free(bytes);
    return (double)totalBits / 8.0;
}

double calculateBitCorrelationCoef(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
     if (!bits || num_bits < 4) {
        if(bits) free(bits);
        return 0.0;
    }
    size_t n = num_bits / 2; // Use length of shorter sequence
    double sumX=0, sumY=0, sumXY=0, sumX2=0, sumY2=0;
    for (size_t i = 0; i < n; ++i) {
        double x = bits[2*i] ? 1.0 : 0.0;
        double y = bits[2*i+1] ? 1.0 : 0.0;
        sumX += x; sumY += y; sumXY += x*y; sumX2 += x*x; sumY2 += y*y;
    }
    free(bits);
    double numerator = n * sumXY - sumX * sumY;
    double denominator = sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));
    if (denominator <= 1e-9) return (numerator == 0) ? 1.0 : 0.0; // Handle zero denominator
    return numerator / denominator;
}

double calculateBitRunLengthSize(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits || num_bits == 0) {
        if(bits) free(bits);
        return 0.0;
    }
    uint32_t encodingSize = 0;
    uint32_t runCount = 0;
     if (num_bits > 0) {
        bool currentBit = bits[0];
        uint32_t currentRun = 1;
        for (size_t i = 1; i < num_bits; ++i) {
            if (bits[i] == currentBit) {
                currentRun++;
            } else {
                runCount++;
                encodingSize += 1 + (uint32_t)floor(log2((double)currentRun + 1e-9)); // Add epsilon for log2(1)
                currentRun = 1;
                currentBit = bits[i];
            }
        }
        runCount++;
        encodingSize += 1 + (uint32_t)floor(log2((double)currentRun + 1e-9));
        encodingSize += runCount; // Bit for each run's type (0/1)
    }
    free(bits);
    return (double)encodingSize / num_bits;
}

double calculateByteMirroringDistance(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
     if (!bytes || num_bytes == 0) {
        if(bytes) free(bytes);
        return 0.0;
    }
    static const uint8_t bitReversalTable[256] = { /* ... precomputed table ... */ }; // Need to define this table
    uint32_t totalDistance = 0;
    for (size_t i = 0; i < num_bytes; ++i) {
        uint8_t reversed = 0; // Simplified reversal
        uint8_t b = bytes[i];
        for(int j=0; j<8; ++j) reversed |= ((b >> j) & 1) << (7-j);
        totalDistance += POPCOUNT_LUT[bytes[i] ^ reversed];
    }
    free(bytes);
    return (double)totalDistance / num_bytes;
}

double calculateRollingByteChecksum(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    uint32_t a = 1, b = 0;
    for (size_t i = 0; smiles[i] != '\0'; ++i) {
        a = (a + (uint8_t)smiles[i]) % 65521;
        b = (b + a) % 65521;
    }
    return (double)((b << 16) | a);
}

double calculateZeroBitRatio(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
     if (!bits || num_bits == 0) {
        if(bits) free(bits);
        return 0.0;
    }
    uint32_t zeroBits = 0;
    for (size_t i = 0; i < num_bits; ++i) if (!bits[i]) zeroBits++;
    free(bits);
    return (double)zeroBits / num_bits;
}

double calculateBitAutocorrelation(const void* context, GetSmilesFunc getSmilesFunc) {
     const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
     if (!bits || num_bits < 2) {
        if(bits) free(bits);
        return 0.0;
    }
    double sum = 0.0; for (size_t i = 0; i < num_bits; ++i) sum += bits[i];
    double mean = sum / num_bits;
    double sumProduct = 0.0, sumSquares = 0.0;
    for (size_t i = 0; i < num_bits - 1; ++i) {
        double v1 = bits[i] - mean;
        double v2 = bits[i+1] - mean;
        sumProduct += v1 * v2;
    }
     for (size_t i = 0; i < num_bits; ++i) {
        double v = bits[i] - mean;
        sumSquares += v*v;
     }
    free(bits);
    return (sumSquares == 0.0) ? 1.0 : sumProduct / sumSquares;
}

double calculateDiagonalBitAccumulation(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bytes = strlen(smiles);
    if (num_bytes == 0) return 0.0;
    uint32_t sum_acc = 0;
     for (size_t i = 0; i < num_bytes; ++i) {
        uint8_t byte = (uint8_t)smiles[i];
        for (int j = 0; j < 8; ++j) {
            if ((byte & (1 << j)) && (((i * 8 + (7-j))) % 256 == byte)) { // Correct bit index
                sum_acc++;
            }
        }
    }
    return (double)sum_acc;
}

double calculateBitwiseMajorityFunction(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
    if (!smiles) return 0.0;
    size_t num_bits;
    bool* bits = smiles_to_bits(smiles, &num_bits);
    if (!bits || num_bits < 3) {
        if(bits) free(bits);
        return 0.0;
    }
    uint32_t majorityCount = 0;
    for (size_t i = 1; i < num_bits - 1; ++i) {
        int sum = bits[i-1] + bits[i] + bits[i+1];
        if (bits[i] == (sum >= 2)) majorityCount++;
    }
    free(bits);
    return (double)majorityCount / (num_bits - 2);
}

double calculateLzcntPattern(const void* context, GetSmilesFunc getSmilesFunc) {
     const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
     if (!bytes || num_bytes == 0) {
        if(bytes) free(bytes);
        return 0.0;
    }
    uint32_t lzcntHist[9] = {0};
    for (size_t i = 0; i < num_bytes; ++i) lzcntHist[LEADING_ZEROS_LUT[bytes[i]]]++;
    double entropy = 0.0;
    for (int i = 0; i < 9; ++i) {
        if (lzcntHist[i] > 0) {
            double prob = (double)lzcntHist[i] / num_bytes;
            entropy -= prob * log2(prob);
        }
    }
    free(bytes);
    return entropy;
}

double calculateTzcntPattern(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes || num_bytes == 0) {
        if(bytes) free(bytes);
        return 0.0;
    }
    uint32_t tzcntHist[9] = {0};
    for (size_t i = 0; i < num_bytes; ++i) tzcntHist[TRAILING_ZEROS_LUT[bytes[i]]]++;
    double entropy = 0.0;
    for (int i = 0; i < 9; ++i) {
        if (tzcntHist[i] > 0) {
            double prob = (double)tzcntHist[i] / num_bytes;
            entropy -= prob * log2(prob);
        }
    }
    free(bytes);
    return entropy;
}

double calculateByteGeometricMean(const void* context, GetSmilesFunc getSmilesFunc) {
    const char* smiles = getSmilesFunc(context);
     if (!smiles) return 0.0;
    size_t num_bytes;
    uint8_t* bytes = smiles_to_bytes(smiles, &num_bytes);
    if (!bytes || num_bytes == 0) {
         if(bytes) free(bytes);
         return 0.0;
    }
    double logSum = 0.0;
    uint32_t validCount = 0;
    for (size_t i = 0; i < num_bytes; ++i) {
        if (bytes[i] > 0) {
            logSum += log((double)bytes[i]);
            validCount++;
        }
    }
    free(bytes);
    if (validCount == 0) return 0.0;
    return exp(logSum / validCount);
}