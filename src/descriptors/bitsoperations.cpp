#include "../common.hpp"
#include <string>
#include <algorithm>
#include <memory>
#include <vector>
#include <numeric>
#include <map>
#include <cmath>
#include <set>
#include <unordered_map>
#include <array>

namespace desfact {

// --- Bit operation helper functions ---

// Thread-local cache to avoid repeated allocations
struct BitOperationsCache {
    std::vector<uint8_t> bytes;
    std::vector<uint8_t> nibbles;
    std::array<uint32_t, 256> byteHistogram;
    std::array<uint32_t, 16> nibbleHistogram;
    std::vector<bool> bits;
    
    void reset() {
        bytes.clear();
        nibbles.clear();
        byteHistogram.fill(0);
        nibbleHistogram.fill(0);
        bits.clear();
    }
    
    void computeFromString(const std::string& smiles) {
        reset();
        
        // Convert string to bytes
        bytes.reserve(smiles.size());
        for (char c : smiles) {
            bytes.push_back(static_cast<uint8_t>(c));
        }
        
        // Create nibbles (4-bit values)
        nibbles.reserve(smiles.size() * 2);
        for (uint8_t byte : bytes) {
            nibbles.push_back((byte >> 4) & 0x0F);  // High nibble
            nibbles.push_back(byte & 0x0F);         // Low nibble
        }
        
        // Compute byte histogram
        for (uint8_t byte : bytes) {
            byteHistogram[byte]++;
        }
        
        // Compute nibble histogram
        for (uint8_t nibble : nibbles) {
            nibbleHistogram[nibble]++;
        }
        
        // Convert string to bits
        bits.reserve(smiles.size() * 8);
        for (char c : smiles) {
            uint8_t byte = static_cast<uint8_t>(c);
            for (int i = 7; i >= 0; --i) {
                bits.push_back((byte >> i) & 1);
            }
        }
    }
};

static thread_local BitOperationsCache tlsBitCache;

// --- Base Class ---
class BitOperationsDescriptor : public Descriptor {
public:
    using Descriptor::Descriptor; // Inherit constructor
    std::string getCategory() const override { return "BitOperations"; }
};

// Helper functions

// Count of 1 bits in a byte using lookup table
static const std::array<uint8_t, 256> POPCOUNT_LUT = []() {
    std::array<uint8_t, 256> table{};
    for (int i = 0; i < 256; i++) {
        uint8_t count = 0;
        uint8_t val = i;
        while (val) {
            count += val & 1;
            val >>= 1;
        }
        table[i] = count;
    }
    return table;
}();

// Fast byte-wise population count
inline uint32_t popcountBytes(const std::vector<uint8_t>& bytes) {
    uint32_t count = 0;
    for (uint8_t byte : bytes) {
        count += POPCOUNT_LUT[byte];
    }
    return count;
}

// Fast bit transition count (0->1, 1->0)
inline uint32_t bitTransitionCount(const std::vector<bool>& bits) {
    if (bits.size() < 2) return 0;
    
    uint32_t count = 0;
    bool prev = bits[0];
    
    for (size_t i = 1; i < bits.size(); ++i) {
        if (bits[i] != prev) {
            count++;
            prev = bits[i];
        }
    }
    
    return count;
}

// Leading/trailing zero count lookups
static const std::array<uint8_t, 256> LEADING_ZEROS_LUT = []() {
    std::array<uint8_t, 256> table{};
    for (int i = 0; i < 256; i++) {
        if (i == 0) {
            table[i] = 8;
        } else {
            uint8_t count = 0;
            uint8_t val = i;
            while ((val & 0x80) == 0) {
                count++;
                val <<= 1;
            }
            table[i] = count;
        }
    }
    return table;
}();

static const std::array<uint8_t, 256> TRAILING_ZEROS_LUT = []() {
    std::array<uint8_t, 256> table{};
    for (int i = 0; i < 256; i++) {
        if (i == 0) {
            table[i] = 8;
        } else {
            uint8_t count = 0;
            uint8_t val = i;
            while ((val & 0x01) == 0) {
                count++;
                val >>= 1;
            }
            table[i] = count;
        }
    }
    return table;
}();

// Find longest consecutive bit run
inline uint32_t longestConsecutiveBitRun(const std::vector<bool>& bits) {
    if (bits.empty()) return 0;
    
    uint32_t maxRun = 0;
    uint32_t currentRun = 1;
    bool currentBit = bits[0];
    
    for (size_t i = 1; i < bits.size(); ++i) {
        if (bits[i] == currentBit) {
            currentRun++;
        } else {
            maxRun = std::max(maxRun, currentRun);
            currentRun = 1;
            currentBit = bits[i];
        }
    }
    
    return std::max(maxRun, currentRun);
}

// Compute byte-wise XOR checksum
inline uint8_t byteXorChecksum(const std::vector<uint8_t>& bytes) {
    uint8_t checksum = 0;
    for (uint8_t byte : bytes) {
        checksum ^= byte;
    }
    return checksum;
}

// Fast entropy calculation
inline double calculateEntropy(const std::array<uint32_t, 256>& histogram, size_t totalCount) {
    if (totalCount == 0) return 0.0;
    
    double entropy = 0.0;
    for (uint32_t count : histogram) {
        if (count > 0) {
            double probability = static_cast<double>(count) / totalCount;
            entropy -= probability * std::log2(probability);
        }
    }
    
    return entropy;
}

// FNV-1a hash function (fast hash)
inline uint32_t fnv1aHash(const std::vector<uint8_t>& bytes) {
    uint32_t hash = 2166136261u; // FNV offset basis
    for (uint8_t byte : bytes) {
        hash ^= byte;
        hash *= 16777619u; // FNV prime
    }
    return hash;
}

// CRC8 calculation
inline uint8_t crc8(const std::vector<uint8_t>& bytes) {
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
    
    uint8_t crc = 0;
    for (uint8_t byte : bytes) {
        crc = CRC8_TABLE[crc ^ byte];
    }
    return crc;
}

// Check if sequence is a palindrome
inline bool isCharPalindrome(const std::vector<uint8_t>& bytes) {
    size_t size = bytes.size();
    if (size <= 1) return true;
    
    for (size_t i = 0; i < size / 2; ++i) {
        if (bytes[i] != bytes[size - 1 - i]) {
            return false;
        }
    }
    return true;
}

// Check if bit sequence is a palindrome
inline bool isBitPalindrome(const std::vector<bool>& bits) {
    size_t size = bits.size();
    if (size <= 1) return true;
    
    for (size_t i = 0; i < size / 2; ++i) {
        if (bits[i] != bits[size - 1 - i]) {
            return false;
        }
    }
    return true;
}

// Calculate bit distribution variance
inline double bitDistributionVariance(const std::vector<uint8_t>& bytes) {
    if (bytes.empty()) return 0.0;
    
    std::vector<uint8_t> popcounts;
    popcounts.reserve(bytes.size());
    
    for (uint8_t byte : bytes) {
        popcounts.push_back(POPCOUNT_LUT[byte]);
    }
    
    double mean = std::accumulate(popcounts.begin(), popcounts.end(), 0.0) / popcounts.size();
    double variance = 0.0;
    
    for (uint8_t count : popcounts) {
        double diff = count - mean;
        variance += diff * diff;
    }
    
    return variance / popcounts.size();
}

// Compute left/right bit balance (positive if more bits in left half)
inline int leftRightBitBalance(const std::vector<bool>& bits) {
    size_t size = bits.size();
    if (size == 0) return 0;
    
    size_t mid = size / 2;
    int leftCount = 0, rightCount = 0;
    
    for (size_t i = 0; i < mid; ++i) {
        if (bits[i]) leftCount++;
    }
    
    for (size_t i = mid; i < size; ++i) {
        if (bits[i]) rightCount++;
    }
    
    return leftCount - rightCount;
}

// --- Descriptor Implementations ---

// 1. Raw String Length
DECLARE_DESCRIPTOR(StringLength, BitOperationsDescriptor, "Raw length of the SMILES string")
DESCRIPTOR_DEPENDENCIES(StringLength) { return {}; }
DescriptorResult StringLengthDescriptor::calculate(Context& context) const {
    return static_cast<double>(context.getSmiles().length());
}

// 2. Bit Count (total binary 1s)
DECLARE_DESCRIPTOR(BitCount, BitOperationsDescriptor, "Total count of binary 1 bits in the SMILES string")
DESCRIPTOR_DEPENDENCIES(BitCount) { return {}; }
DescriptorResult BitCountDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    return static_cast<double>(popcountBytes(tlsBitCache.bytes));
}

// 3. Bit Density (1s per byte)
DECLARE_DESCRIPTOR(BitDensity, BitOperationsDescriptor, "Density of 1 bits per byte")
DESCRIPTOR_DEPENDENCIES(BitDensity) { return {}; }
DescriptorResult BitDensityDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    double totalBits = popcountBytes(tlsBitCache.bytes);
    return totalBits / smiles.length();
}

// 4. Ascii Value Sum
DECLARE_DESCRIPTOR(AsciiSum, BitOperationsDescriptor, "Sum of ASCII values of all characters")
DESCRIPTOR_DEPENDENCIES(AsciiSum) { return {}; }
DescriptorResult AsciiSumDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    uint32_t sum = 0;
    for (char c : smiles) {
        sum += static_cast<uint8_t>(c);
    }
    return static_cast<double>(sum);
}

// 5. Longest Consecutive Bit Run
DECLARE_DESCRIPTOR(LongestBitRun, BitOperationsDescriptor, "Length of longest consecutive bit run (0s or 1s)")
DESCRIPTOR_DEPENDENCIES(LongestBitRun) { return {}; }
DescriptorResult LongestBitRunDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    return static_cast<double>(longestConsecutiveBitRun(tlsBitCache.bits));
}

// 6. Byte-wise XOR Checksum
DECLARE_DESCRIPTOR(ByteXorChecksum, BitOperationsDescriptor, "XOR of all bytes in the string")
DESCRIPTOR_DEPENDENCIES(ByteXorChecksum) { return {}; }
DescriptorResult ByteXorChecksumDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    return static_cast<double>(byteXorChecksum(tlsBitCache.bytes));
}

// 7. Byte Entropy
DECLARE_DESCRIPTOR(ByteEntropy, BitOperationsDescriptor, "Shannon entropy of byte distribution")
DESCRIPTOR_DEPENDENCIES(ByteEntropy) { return {}; }
DescriptorResult ByteEntropyDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    return calculateEntropy(tlsBitCache.byteHistogram, smiles.length());
}

// 8. First-byte Value
DECLARE_DESCRIPTOR(FirstByteValue, BitOperationsDescriptor, "Value of the first byte in the string")
DESCRIPTOR_DEPENDENCIES(FirstByteValue) { return {}; }
DescriptorResult FirstByteValueDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    return static_cast<double>(static_cast<uint8_t>(smiles[0]));
}

// 9. Last-byte Value
DECLARE_DESCRIPTOR(LastByteValue, BitOperationsDescriptor, "Value of the last byte in the string")
DESCRIPTOR_DEPENDENCIES(LastByteValue) { return {}; }
DescriptorResult LastByteValueDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    return static_cast<double>(static_cast<uint8_t>(smiles.back()));
}

// 10. Odd/Even Bit Ratio
DECLARE_DESCRIPTOR(OddEvenBitRatio, BitOperationsDescriptor, "Ratio of 1 bits at odd positions to even positions")
DESCRIPTOR_DEPENDENCIES(OddEvenBitRatio) { return {}; }
DescriptorResult OddEvenBitRatioDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    uint32_t oddCount = 0, evenCount = 0;
    
    for (size_t i = 0; i < bits.size(); ++i) {
        if (bits[i]) {
            if (i % 2 == 0) evenCount++;
            else oddCount++;
        }
    }
    
    return (evenCount == 0) ? 0.0 : static_cast<double>(oddCount) / evenCount;
}

// 11. Character Palindrome Check
DECLARE_DESCRIPTOR(CharPalindromeCheck, BitOperationsDescriptor, "Boolean indicating if string is a character palindrome")
DESCRIPTOR_DEPENDENCIES(CharPalindromeCheck) { return {}; }
DescriptorResult CharPalindromeCheckDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    return isCharPalindrome(tlsBitCache.bytes) ? 1.0 : 0.0;
}

// 13. Bit Transition Count
DECLARE_DESCRIPTOR(BitTransitionCount, BitOperationsDescriptor, "Count of bit transitions (0→1, 1→0)")
DESCRIPTOR_DEPENDENCIES(BitTransitionCount) { return {}; }
DescriptorResult BitTransitionCountDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    return static_cast<double>(bitTransitionCount(tlsBitCache.bits));
}

// 14. Left/Right Bit Balance
DECLARE_DESCRIPTOR(LeftRightBitBalance, BitOperationsDescriptor, "Difference between bit counts in left and right halves")
DESCRIPTOR_DEPENDENCIES(LeftRightBitBalance) { return {}; }
DescriptorResult LeftRightBitBalanceDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    return static_cast<double>(leftRightBitBalance(tlsBitCache.bits));
}

// 15. Hamming Weight Variance
DECLARE_DESCRIPTOR(HammingWeightVariance, BitOperationsDescriptor, "Variance of population counts across bytes")
DESCRIPTOR_DEPENDENCIES(HammingWeightVariance) { return {}; }
DescriptorResult HammingWeightVarianceDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    return bitDistributionVariance(tlsBitCache.bytes);
}

// 16. Fast Hash (FNV-1a)
DECLARE_DESCRIPTOR(FastHash, BitOperationsDescriptor, "FNV-1a hash value of the string")
DESCRIPTOR_DEPENDENCIES(FastHash) { return {}; }
DescriptorResult FastHashDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    return static_cast<double>(fnv1aHash(tlsBitCache.bytes));
}

// 17. CRC8 Checksum
DECLARE_DESCRIPTOR(Crc8Checksum, BitOperationsDescriptor, "CRC8 checksum of the string")
DESCRIPTOR_DEPENDENCIES(Crc8Checksum) { return {}; }
DescriptorResult Crc8ChecksumDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    return static_cast<double>(crc8(tlsBitCache.bytes));
}

// 18. Most Frequent Byte
DECLARE_DESCRIPTOR(MostFrequentByte, BitOperationsDescriptor, "Value of the most frequently occurring byte")
DESCRIPTOR_DEPENDENCIES(MostFrequentByte) { return {}; }
DescriptorResult MostFrequentByteDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    
    uint8_t mostFrequent = 0;
    uint32_t maxCount = 0;
    
    for (size_t i = 0; i < 256; ++i) {
        if (tlsBitCache.byteHistogram[i] > maxCount) {
            maxCount = tlsBitCache.byteHistogram[i];
            mostFrequent = static_cast<uint8_t>(i);
        }
    }
    
    return static_cast<double>(mostFrequent);
}

// 19. Least Frequent Byte
DECLARE_DESCRIPTOR(LeastFrequentByte, BitOperationsDescriptor, "Value of the least frequently occurring byte (that appears at least once)")
DESCRIPTOR_DEPENDENCIES(LeastFrequentByte) { return {}; }
DescriptorResult LeastFrequentByteDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    
    uint8_t leastFrequent = 0;
    uint32_t minCount = UINT32_MAX;
    
    for (size_t i = 0; i < 256; ++i) {
        if (tlsBitCache.byteHistogram[i] > 0 && tlsBitCache.byteHistogram[i] < minCount) {
            minCount = tlsBitCache.byteHistogram[i];
            leastFrequent = static_cast<uint8_t>(i);
        }
    }
    
    return static_cast<double>(leastFrequent);
}

// 20. Unique Byte Count
DECLARE_DESCRIPTOR(UniqueByteCount, BitOperationsDescriptor, "Count of unique byte values in the string")
DESCRIPTOR_DEPENDENCIES(UniqueByteCount) { return {}; }
DescriptorResult UniqueByteCountDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    
    uint32_t uniqueCount = 0;
    for (uint32_t count : tlsBitCache.byteHistogram) {
        if (count > 0) uniqueCount++;
    }
    
    return static_cast<double>(uniqueCount);
}

// 21. Even-indexed Bit Count
DECLARE_DESCRIPTOR(EvenIndexedBitCount, BitOperationsDescriptor, "Count of 1 bits at even-indexed positions")
DESCRIPTOR_DEPENDENCIES(EvenIndexedBitCount) { return {}; }
DescriptorResult EvenIndexedBitCountDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    uint32_t count = 0;
    for (size_t i = 0; i < bits.size(); i += 2) {
        if (bits[i]) count++;
    }
    
    return static_cast<double>(count);
}

// 22. Odd-indexed Bit Count
DECLARE_DESCRIPTOR(OddIndexedBitCount, BitOperationsDescriptor, "Count of 1 bits at odd-indexed positions")
DESCRIPTOR_DEPENDENCIES(OddIndexedBitCount) { return {}; }
DescriptorResult OddIndexedBitCountDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    uint32_t count = 0;
    for (size_t i = 1; i < bits.size(); i += 2) {
        if (bits[i]) count++;
    }
    
    return static_cast<double>(count);
}

// 23. Byte-wise Minimum
DECLARE_DESCRIPTOR(ByteWiseMinimum, BitOperationsDescriptor, "Minimum byte value in the string")
DESCRIPTOR_DEPENDENCIES(ByteWiseMinimum) { return {}; }
DescriptorResult ByteWiseMinimumDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    uint8_t minValue = 255;
    for (char c : smiles) {
        uint8_t value = static_cast<uint8_t>(c);
        if (value < minValue) {
            minValue = value;
        }
    }
    
    return static_cast<double>(minValue);
}

// 24. Byte-wise Maximum
DECLARE_DESCRIPTOR(ByteWiseMaximum, BitOperationsDescriptor, "Maximum byte value in the string")
DESCRIPTOR_DEPENDENCIES(ByteWiseMaximum) { return {}; }
DescriptorResult ByteWiseMaximumDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    uint8_t maxValue = 0;
    for (char c : smiles) {
        uint8_t value = static_cast<uint8_t>(c);
        if (value > maxValue) {
            maxValue = value;
        }
    }
    
    return static_cast<double>(maxValue);
}

// 25. Byte-wise Median
DECLARE_DESCRIPTOR(ByteWiseMedian, BitOperationsDescriptor, "Median byte value in the string")
DESCRIPTOR_DEPENDENCIES(ByteWiseMedian) { return {}; }
DescriptorResult ByteWiseMedianDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    std::vector<uint8_t> values;
    values.reserve(smiles.size());
    
    for (char c : smiles) {
        values.push_back(static_cast<uint8_t>(c));
    }
    
    std::sort(values.begin(), values.end());
    
    size_t mid = values.size() / 2;
    if (values.size() % 2 == 0) {
        return (values[mid - 1] + values[mid]) / 2.0;
    } else {
        return static_cast<double>(values[mid]);
    }
}

// 26. High/Low Bit Ratio
DECLARE_DESCRIPTOR(HighLowBitRatio, BitOperationsDescriptor, "Ratio of high bits (values 1) to low bits (values 0)")
DESCRIPTOR_DEPENDENCIES(HighLowBitRatio) { return {}; }
DescriptorResult HighLowBitRatioDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    if (bits.empty()) return 0.0;
    
    uint32_t highCount = 0;
    
    for (bool bit : bits) {
        if (bit) highCount++;
    }
    
    uint32_t lowCount = bits.size() - highCount;
    return (lowCount == 0) ? static_cast<double>(bits.size()) : static_cast<double>(highCount) / lowCount;
}

// 27. Character Sequence Complexity
DECLARE_DESCRIPTOR(CharSequenceComplexity, BitOperationsDescriptor, "Ratio of unique characters to total length")
DESCRIPTOR_DEPENDENCIES(CharSequenceComplexity) { return {}; }
DescriptorResult CharSequenceComplexityDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    std::set<char> uniqueChars(smiles.begin(), smiles.end());
    return static_cast<double>(uniqueChars.size()) / smiles.length();
}

// 28. Bit Population Count per Byte (Average)
DECLARE_DESCRIPTOR(AvgBitPopulationCount, BitOperationsDescriptor, "Average number of set bits per byte")
DESCRIPTOR_DEPENDENCIES(AvgBitPopulationCount) { return {}; }
DescriptorResult AvgBitPopulationCountDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    double totalPopcount = popcountBytes(tlsBitCache.bytes);
    
    return totalPopcount / smiles.length();
}

// 29. Leading Zero Count (Average)
DECLARE_DESCRIPTOR(AvgLeadingZeroCount, BitOperationsDescriptor, "Average number of leading zeros per byte")
DESCRIPTOR_DEPENDENCIES(AvgLeadingZeroCount) { return {}; }
DescriptorResult AvgLeadingZeroCountDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    
    uint32_t totalLeadingZeros = 0;
    for (uint8_t byte : tlsBitCache.bytes) {
        totalLeadingZeros += LEADING_ZEROS_LUT[byte];
    }
    
    return static_cast<double>(totalLeadingZeros) / smiles.length();
}

// 30. Trailing Zero Count (Average)
DECLARE_DESCRIPTOR(AvgTrailingZeroCount, BitOperationsDescriptor, "Average number of trailing zeros per byte")
DESCRIPTOR_DEPENDENCIES(AvgTrailingZeroCount) { return {}; }
DescriptorResult AvgTrailingZeroCountDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    
    uint32_t totalTrailingZeros = 0;
    for (uint8_t byte : tlsBitCache.bytes) {
        totalTrailingZeros += TRAILING_ZEROS_LUT[byte];
    }
    
    return static_cast<double>(totalTrailingZeros) / smiles.length();
}

// 31. Popcount Variance Across Bytes
DECLARE_DESCRIPTOR(PopcountVariance, BitOperationsDescriptor, "Variance of the bit count across bytes")
DESCRIPTOR_DEPENDENCIES(PopcountVariance) { return {}; }
DescriptorResult PopcountVarianceDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    
    std::vector<uint8_t> popcounts;
    popcounts.reserve(tlsBitCache.bytes.size());
    
    for (uint8_t byte : tlsBitCache.bytes) {
        popcounts.push_back(POPCOUNT_LUT[byte]);
    }
    
    double mean = std::accumulate(popcounts.begin(), popcounts.end(), 0.0) / popcounts.size();
    double variance = 0.0;
    
    for (uint8_t count : popcounts) {
        double diff = count - mean;
        variance += diff * diff;
    }
    
    return variance / popcounts.size();
}

// 32. Byte-wise AND Reduction
DECLARE_DESCRIPTOR(ByteAndReduction, BitOperationsDescriptor, "Result of ANDing all bytes together")
DESCRIPTOR_DEPENDENCIES(ByteAndReduction) { return {}; }
DescriptorResult ByteAndReductionDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    
    uint8_t result = 255;  // All bits set initially
    for (uint8_t byte : tlsBitCache.bytes) {
        result &= byte;
    }
    
    return static_cast<double>(result);
}

// 33. Byte-wise OR Reduction
DECLARE_DESCRIPTOR(ByteOrReduction, BitOperationsDescriptor, "Result of ORing all bytes together")
DESCRIPTOR_DEPENDENCIES(ByteOrReduction) { return {}; }
DescriptorResult ByteOrReductionDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    
    uint8_t result = 0;  // No bits set initially
    for (uint8_t byte : tlsBitCache.bytes) {
        result |= byte;
    }
    
    return static_cast<double>(result);
}

// 34. Byte-wise XOR Reduction
DECLARE_DESCRIPTOR(ByteXorReduction, BitOperationsDescriptor, "Result of XORing all bytes together")
DESCRIPTOR_DEPENDENCIES(ByteXorReduction) { return {}; }
DescriptorResult ByteXorReductionDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;

    tlsBitCache.computeFromString(smiles);

    uint8_t result = 0;  // No bits set initially
    for (uint8_t byte : tlsBitCache.bytes) {
        result ^= byte;
    }

    return static_cast<double>(result);
}

// 36. Non-zero Byte Density
DECLARE_DESCRIPTOR(NonZeroByteDensity, BitOperationsDescriptor, "Fraction of bytes that are non-zero")
DESCRIPTOR_DEPENDENCIES(NonZeroByteDensity) { return {}; }
DescriptorResult NonZeroByteDensityDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;

    tlsBitCache.computeFromString(smiles);

    uint32_t nonZeroCount = tlsBitCache.bytes.size() - tlsBitCache.byteHistogram[0];
    return static_cast<double>(nonZeroCount) / tlsBitCache.bytes.size();
}

// 37. ASCII Value Product (mod 256)
DECLARE_DESCRIPTOR(AsciiProductMod256, BitOperationsDescriptor, "Product of all ASCII values, modulo 256")
DESCRIPTOR_DEPENDENCIES(AsciiProductMod256) { return {}; }
DescriptorResult AsciiProductMod256Descriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    uint32_t product = 1;
    for (char c : smiles) {
        product = (product * static_cast<uint8_t>(c)) % 256;
    }
    
    return static_cast<double>(product);
}

// 38. Bit Alternation Frequency
DECLARE_DESCRIPTOR(BitAlternationFreq, BitOperationsDescriptor, "Frequency of bit flips (transitions) per bit")
DESCRIPTOR_DEPENDENCIES(BitAlternationFreq) { return {}; }
DescriptorResult BitAlternationFreqDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    const auto& bits = tlsBitCache.bits;
    
    if (bits.size() <= 1) return 0.0;
    
    uint32_t transitions = bitTransitionCount(bits);
    return static_cast<double>(transitions) / (bits.size() - 1);
}

// 39. Character Distribution Variance
DECLARE_DESCRIPTOR(CharDistributionVariance, BitOperationsDescriptor, "Variance of character frequencies")
DESCRIPTOR_DEPENDENCIES(CharDistributionVariance) { return {}; }
DescriptorResult CharDistributionVarianceDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    
    // Get non-zero frequencies
    std::vector<double> nonZeroFreqs;
    for (uint32_t count : tlsBitCache.byteHistogram) {
        if (count > 0) {
            nonZeroFreqs.push_back(static_cast<double>(count) / smiles.length());
        }
    }
    
    if (nonZeroFreqs.empty()) return 0.0;
    
    double mean = 1.0 / nonZeroFreqs.size();  // Theoretical mean for uniform distribution
    double variance = 0.0;
    
    for (double freq : nonZeroFreqs) {
        double diff = freq - mean;
        variance += diff * diff;
    }
    
    return variance / nonZeroFreqs.size();
}

// 40. Low/High Nibble Ratio
DECLARE_DESCRIPTOR(LowHighNibbleRatio, BitOperationsDescriptor, "Ratio of frequency of low nibbles (0-7) to high nibbles (8-15)")
DESCRIPTOR_DEPENDENCIES(LowHighNibbleRatio) { return {}; }
DescriptorResult LowHighNibbleRatioDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    
    uint32_t lowCount = 0, highCount = 0;
    
    for (uint8_t nibble : tlsBitCache.nibbles) {
        if (nibble < 8) lowCount++;
        else highCount++;
    }
    
    return (highCount == 0) ? static_cast<double>(tlsBitCache.nibbles.size()) : static_cast<double>(lowCount) / highCount;
}

// 41. Character Position Weighted Sum
DECLARE_DESCRIPTOR(CharPosWeightedSum, BitOperationsDescriptor, "Sum of byte values weighted by their position")
DESCRIPTOR_DEPENDENCIES(CharPosWeightedSum) { return {}; }
DescriptorResult CharPosWeightedSumDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    double weightedSum = 0.0;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        weightedSum += static_cast<uint8_t>(smiles[i]) * (i + 1);
    }
    
    return weightedSum;
}

// 43. Nibble Histogram Entropy
DECLARE_DESCRIPTOR(NibbleEntropy, BitOperationsDescriptor, "Shannon entropy of the nibble (4-bit) distribution")
DESCRIPTOR_DEPENDENCIES(NibbleEntropy) { return {}; }
DescriptorResult NibbleEntropyDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    
    // Calculate entropy using nibble histogram
    double entropy = 0.0;
    size_t totalNibbles = tlsBitCache.nibbles.size();
    
    for (uint32_t count : tlsBitCache.nibbleHistogram) {
        if (count > 0) {
            double probability = static_cast<double>(count) / totalNibbles;
            entropy -= probability * std::log2(probability);
        }
    }
    
    return entropy;
}

// 44. Character Frequency Histogram (Shannon entropy)
DECLARE_DESCRIPTOR(CharFrequencyEntropy, BitOperationsDescriptor, "Shannon entropy of character frequencies")
DESCRIPTOR_DEPENDENCIES(CharFrequencyEntropy) { return {}; }
DescriptorResult CharFrequencyEntropyDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    return calculateEntropy(tlsBitCache.byteHistogram, smiles.length());
}

// 45. ASCII Range Coverage
DECLARE_DESCRIPTOR(AsciiRangeCoverage, BitOperationsDescriptor, "Ratio of ASCII range covered by the string")
DESCRIPTOR_DEPENDENCIES(AsciiRangeCoverage) { return {}; }
DescriptorResult AsciiRangeCoverageDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    
    uint8_t minAscii = 255, maxAscii = 0;
    bool hasValues = false;
    
    for (size_t i = 0; i < 256; ++i) {
        if (tlsBitCache.byteHistogram[i] > 0) {
            minAscii = std::min(minAscii, static_cast<uint8_t>(i));
            maxAscii = std::max(maxAscii, static_cast<uint8_t>(i));
            hasValues = true;
        }
    }
    
    if (!hasValues) return 0.0;
    return static_cast<double>(maxAscii - minAscii + 1) / 256.0;
}

// 46. Nibble (4-bit) Histogram Analysis
DECLARE_DESCRIPTOR(NibbleHistogramBalance, BitOperationsDescriptor, "Balance of nibble frequencies (higher = more uniform)")
DESCRIPTOR_DEPENDENCIES(NibbleHistogramBalance) { return {}; }
DescriptorResult NibbleHistogramBalanceDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    
    // Get non-zero frequencies
    std::vector<double> nonZeroFreqs;
    for (uint32_t count : tlsBitCache.nibbleHistogram) {
        if (count > 0) {
            nonZeroFreqs.push_back(static_cast<double>(count) / tlsBitCache.nibbles.size());
        }
    }
    
    if (nonZeroFreqs.empty()) return 0.0;
    
    // Calculate coefficient of variation (lower = more uniform)
    double mean = 1.0 / nonZeroFreqs.size();
    double variance = 0.0;
    
    for (double freq : nonZeroFreqs) {
        double diff = freq - mean;
        variance += diff * diff;
    }
    
    double stdDev = std::sqrt(variance / nonZeroFreqs.size());
    double cv = stdDev / mean;
    
    // Transform to "balance" where higher = more uniform
    return 1.0 / (1.0 + cv);
}

// 47. Binary Compression Ratio Estimate
DECLARE_DESCRIPTOR(BinaryCompressionRatio, BitOperationsDescriptor, "Estimate of compressibility based on bit patterns")
DESCRIPTOR_DEPENDENCIES(BinaryCompressionRatio) { return {}; }
DescriptorResult BinaryCompressionRatioDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.length() < 2) return 1.0; // No compression possible for very short strings
    
    tlsBitCache.computeFromString(smiles);
    
    // Use entropy as a proxy for compressibility
    double byteEntropy = calculateEntropy(tlsBitCache.byteHistogram, smiles.length());
    
    // Normalize by maximum possible entropy (8 bits)
    double normalizedEntropy = byteEntropy / 8.0;
    
    // Estimate compression ratio (lower entropy = better compression)
    return 1.0 - (normalizedEntropy * 0.9); // Adjust factor for realistic compression ratio
}

// 48. Bit Reversal Distance
DECLARE_DESCRIPTOR(BitReversalDistance, BitOperationsDescriptor, "Hamming distance between bit pattern and its reversal")
DESCRIPTOR_DEPENDENCIES(BitReversalDistance) { return {}; }
DescriptorResult BitReversalDistanceDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    if (bits.size() <= 1) return 0.0;
    
    uint32_t distance = 0;
    size_t size = bits.size();
    
    for (size_t i = 0; i < size / 2; ++i) {
        if (bits[i] != bits[size - 1 - i]) {
            distance++;
        }
    }
    
    return static_cast<double>(distance);
}

// 49. Byte Reversal Distance
DECLARE_DESCRIPTOR(ByteReversalDistance, BitOperationsDescriptor, "Hamming distance between byte sequence and its reversal")
DESCRIPTOR_DEPENDENCIES(ByteReversalDistance) { return {}; }
DescriptorResult ByteReversalDistanceDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bytes = tlsBitCache.bytes;
    
    if (bytes.size() <= 1) return 0.0;
    
    uint32_t distance = 0;
    size_t size = bytes.size();
    
    for (size_t i = 0; i < size / 2; ++i) {
        if (bytes[i] != bytes[size - 1 - i]) {
            distance++;
        }
    }
    
    return static_cast<double>(distance);
}

// 50. MurmurHash3 Index
DECLARE_DESCRIPTOR(MurmurHash3, BitOperationsDescriptor, "MurmurHash3 hash value of the string")
DESCRIPTOR_DEPENDENCIES(MurmurHash3) { return {}; }
DescriptorResult MurmurHash3Descriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    // MurmurHash3 implementation (32-bit)
    const uint32_t seed = 0x9747b28c;
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
    
    uint32_t h1 = seed;
    
    // Process 4 bytes at a time
    size_t len = smiles.length();
    size_t nblocks = len / 4;
    
    for (size_t i = 0; i < nblocks; ++i) {
        uint32_t k1 = 0;
        for (int j = 0; j < 4; ++j) {
            k1 |= (static_cast<uint32_t>(static_cast<uint8_t>(smiles[i*4 + j])) << (j*8));
        }
        
        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> 17);
        k1 *= c2;
        
        h1 ^= k1;
        h1 = (h1 << 13) | (h1 >> 19);
        h1 = h1 * 5 + 0xe6546b64;
    }
    
    // Process remaining bytes
    uint32_t k1 = 0;
    for (size_t i = nblocks * 4; i < len; ++i) {
        k1 |= (static_cast<uint32_t>(static_cast<uint8_t>(smiles[i])) << ((i - nblocks * 4) * 8));
    }
    
    k1 *= c1;
    k1 = (k1 << 15) | (k1 >> 17);
    k1 *= c2;
    h1 ^= k1;
    
    // Finalization
    h1 ^= len;
    h1 ^= (h1 >> 16);
    h1 *= 0x85ebca6b;
    h1 ^= (h1 >> 13);
    h1 *= 0xc2b2ae35;
    h1 ^= (h1 >> 16);
    
    return static_cast<double>(h1);
}

// 51. Character Run Length
DECLARE_DESCRIPTOR(CharRunLength, BitOperationsDescriptor, "Average length of runs of identical characters")
DESCRIPTOR_DEPENDENCIES(CharRunLength) { return {}; }
DescriptorResult CharRunLengthDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    uint32_t runCount = 0;
    uint32_t totalRuns = 0;
    
    char prevChar = smiles[0];
    uint32_t currentRun = 1;
    
    for (size_t i = 1; i < smiles.length(); ++i) {
        if (smiles[i] == prevChar) {
            currentRun++;
        } else {
            totalRuns += currentRun;
            runCount++;
            currentRun = 1;
            prevChar = smiles[i];
        }
    }
    
    // Don't forget the last run
    totalRuns += currentRun;
    runCount++;
    
    return static_cast<double>(totalRuns) / runCount;
}

// 52. Byte-wise Run Entropy
DECLARE_DESCRIPTOR(ByteRunEntropy, BitOperationsDescriptor, "Entropy of run lengths of bytes")
DESCRIPTOR_DEPENDENCIES(ByteRunEntropy) { return {}; }
DescriptorResult ByteRunEntropyDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.length() <= 1) return 0.0;
    
    // Count runs
    std::unordered_map<uint32_t, uint32_t> runLengthHist;
    
    char prevChar = smiles[0];
    uint32_t currentRun = 1;
    
    for (size_t i = 1; i < smiles.length(); ++i) {
        if (smiles[i] == prevChar) {
            currentRun++;
        } else {
            runLengthHist[currentRun]++;
            currentRun = 1;
            prevChar = smiles[i];
        }
    }
    
    // Add the last run
    runLengthHist[currentRun]++;
    
    // Calculate entropy
    double entropy = 0.0;
    double totalRuns = 0.0;
    
    for (const auto& pair : runLengthHist) {
        totalRuns += pair.second;
    }
    
    for (const auto& pair : runLengthHist) {
        double probability = static_cast<double>(pair.second) / totalRuns;
        entropy -= probability * std::log2(probability);
    }
    
    return entropy;
}

// 53. Single-bit Pattern Match
DECLARE_DESCRIPTOR(SingleBitPatternMatch, BitOperationsDescriptor, "Frequency of single bit patterns ('1' followed by '0')")
DESCRIPTOR_DEPENDENCIES(SingleBitPatternMatch) { return {}; }
DescriptorResult SingleBitPatternMatchDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    if (bits.size() < 2) return 0.0;
    
    uint32_t patternCount = 0;
    
    for (size_t i = 0; i < bits.size() - 1; ++i) {
        if (bits[i] && !bits[i + 1]) {
            patternCount++;
        }
    }
    
    return static_cast<double>(patternCount) / (bits.size() - 1);
}

// 54. Byte Pair Transitions
DECLARE_DESCRIPTOR(BytePairTransitions, BitOperationsDescriptor, "Shannon entropy of byte pair transitions")
DESCRIPTOR_DEPENDENCIES(BytePairTransitions) { return {}; }
DescriptorResult BytePairTransitionsDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.length() < 2) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    const auto& bytes = tlsBitCache.bytes;
    
    // Count transitions between pairs
    std::unordered_map<uint16_t, uint32_t> pairCounts;
    
    for (size_t i = 0; i < bytes.size() - 1; ++i) {
        uint16_t pair = (static_cast<uint16_t>(bytes[i]) << 8) | bytes[i + 1];
        pairCounts[pair]++;
    }
    
    // Calculate entropy
    double entropy = 0.0;
    size_t totalPairs = bytes.size() - 1;
    
    for (const auto& pair : pairCounts) {
        double probability = static_cast<double>(pair.second) / totalPairs;
        entropy -= probability * std::log2(probability);
    }
    
    return entropy;
}

// 55. First-Order Bit Difference Sum
DECLARE_DESCRIPTOR(BitDifferenceSum, BitOperationsDescriptor, "Sum of bit differences between consecutive bytes")
DESCRIPTOR_DEPENDENCIES(BitDifferenceSum) { return {}; }
DescriptorResult BitDifferenceSumDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.length() < 2) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    const auto& bytes = tlsBitCache.bytes;
    
    uint32_t diffSum = 0;
    
    for (size_t i = 0; i < bytes.size() - 1; ++i) {
        // XOR gives 1 for different bits
        uint8_t bitDiff = bytes[i] ^ bytes[i + 1];
        diffSum += POPCOUNT_LUT[bitDiff];
    }
    
    return static_cast<double>(diffSum);
}

// 56. First/Last N-byte XOR
DECLARE_DESCRIPTOR(FirstLastByteXor, BitOperationsDescriptor, "XOR of first and last bytes")
DESCRIPTOR_DEPENDENCIES(FirstLastByteXor) { return {}; }
DescriptorResult FirstLastByteXorDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    const auto& bytes = tlsBitCache.bytes;
    
    if (bytes.size() < 2) return static_cast<double>(bytes[0]);
    
    return static_cast<double>(bytes[0] ^ bytes.back());
}

// 57. Byte-level Markov Property
DECLARE_DESCRIPTOR(ByteMarkovProperty, BitOperationsDescriptor, "Measure of how well bytes predict subsequent bytes")
DESCRIPTOR_DEPENDENCIES(ByteMarkovProperty) { return {}; }
DescriptorResult ByteMarkovPropertyDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.length() < 3) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    const auto& bytes = tlsBitCache.bytes;
    
    // Count byte frequencies and transitions
    std::array<uint32_t, 256> byteFreq{};
    std::map<std::pair<uint8_t, uint8_t>, uint32_t> transitions;
    
    for (size_t i = 0; i < bytes.size() - 1; ++i) {
        byteFreq[bytes[i]]++;
        transitions[{bytes[i], bytes[i+1]}]++;
    }
    byteFreq[bytes.back()]++; // Count the last byte
    
    // Calculate conditional entropy
    double conditionalEntropy = 0.0;
    size_t totalTransitions = bytes.size() - 1;
    
    for (const auto& tr : transitions) {
        double jointProb = static_cast<double>(tr.second) / totalTransitions;
        double marginalProb = static_cast<double>(byteFreq[tr.first.first]) / bytes.size();
        double conditionalProb = jointProb / marginalProb;
        
        conditionalEntropy -= jointProb * std::log2(conditionalProb);
    }
    
    // Normalize to [0,1] where 1 indicates high predictability
    double maxEntropy = 8.0; // Max possible is 8 bits
    return 1.0 - (conditionalEntropy / maxEntropy);
}

// 58. SWAR (SIMD within a register) popcount
DECLARE_DESCRIPTOR(SwarPopcount, BitOperationsDescriptor, "SIMD-within-a-register population count")
DESCRIPTOR_DEPENDENCIES(SwarPopcount) { return {}; }
DescriptorResult SwarPopcountDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bytes = tlsBitCache.bytes;
    
    if (bytes.empty()) return 0.0;
    
    uint32_t totalBits = 0;
    
    // Process in 32-bit chunks for efficiency
    for (size_t i = 0; i < bytes.size(); i += 4) {
        uint32_t chunk = 0;
        
        for (size_t j = 0; j < 4 && i + j < bytes.size(); ++j) {
            chunk |= (static_cast<uint32_t>(bytes[i + j]) << (j * 8));
        }
        
        // SWAR algorithm for popcount
        chunk = chunk - ((chunk >> 1) & 0x55555555);
        chunk = (chunk & 0x33333333) + ((chunk >> 2) & 0x33333333);
        chunk = (chunk + (chunk >> 4)) & 0x0F0F0F0F;
        chunk = chunk + (chunk >> 8);
        chunk = chunk + (chunk >> 16);
        
        totalBits += chunk & 0x3F;
    }
    
    return static_cast<double>(totalBits);
}

// 59. Fast Bit Histogram Calculation
DECLARE_DESCRIPTOR(FastBitHistogram, BitOperationsDescriptor, "Histogram uniformity of bit patterns")
DESCRIPTOR_DEPENDENCIES(FastBitHistogram) { return {}; }
DescriptorResult FastBitHistogramDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    
    // Count 4-bit patterns (nibbles)
    const auto& nibbles = tlsBitCache.nibbles;
    std::array<uint32_t, 16> nibbleHist{};
    
    for (uint8_t nibble : nibbles) {
        nibbleHist[nibble]++;
    }
    
    // Calculate uniformity (0 = all in one bin, 1 = perfectly uniform)
    double maxCount = 0;
    double nonZeroBins = 0;
    
    for (uint32_t count : nibbleHist) {
        if (count > 0) {
            nonZeroBins++;
            maxCount = std::max(maxCount, static_cast<double>(count));
        }
    }
    
    if (nonZeroBins <= 1) return 0.0;
    
    double expectedCount = nibbles.size() / nonZeroBins;
    double uniformity = 1.0 - ((maxCount - expectedCount) / (nibbles.size() - expectedCount));
    
    return uniformity;
}

// 60. Bit Dispersion Index
DECLARE_DESCRIPTOR(BitDispersionIndex, BitOperationsDescriptor, "Measure of how evenly distributed the bits are")
DESCRIPTOR_DEPENDENCIES(BitDispersionIndex) { return {}; }
DescriptorResult BitDispersionIndexDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    if (bits.size() <= 1) return 0.0;
    
    // Count runs of 0s and 1s
    std::vector<uint32_t> runLengths;
    uint32_t currentRun = 1;
    bool currentBit = bits[0];
    
    for (size_t i = 1; i < bits.size(); ++i) {
        if (bits[i] == currentBit) {
            currentRun++;
        } else {
            runLengths.push_back(currentRun);
            currentRun = 1;
            currentBit = bits[i];
        }
    }
    runLengths.push_back(currentRun);
    
    // Calculate variance of run lengths
    double mean = static_cast<double>(bits.size()) / runLengths.size();
    double variance = 0.0;
    
    for (uint32_t length : runLengths) {
        double diff = length - mean;
        variance += diff * diff;
    }
    
    variance /= runLengths.size();
    
    // Normalize: lower variance = more dispersed
    return 1.0 / (1.0 + variance);
}

// 61. Half-byte Rotation Hash
DECLARE_DESCRIPTOR(HalfByteRotationHash, BitOperationsDescriptor, "Hash value based on rotating half-bytes")
DESCRIPTOR_DEPENDENCIES(HalfByteRotationHash) { return {}; }
DescriptorResult HalfByteRotationHashDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    uint32_t hash = 0x8F1BBCDC;
    
    for (char c : smiles) {
        uint8_t byte = static_cast<uint8_t>(c);
        uint8_t high = (byte >> 4) & 0x0F;
        uint8_t low = byte & 0x0F;
        
        // Rotate and mix
        hash = ((hash << 5) | (hash >> 27)) ^ (high << 4 | low);
        hash = ((hash << 3) | (hash >> 29)) + ((low << 4) | high);
    }
    
    return static_cast<double>(hash);
}

// 62. Avalanche Pattern Sensitivity
DECLARE_DESCRIPTOR(AvalanchePatternSensitivity, BitOperationsDescriptor, "Measure of bit pattern change propagation")
DESCRIPTOR_DEPENDENCIES(AvalanchePatternSensitivity) { return {}; }
DescriptorResult AvalanchePatternSensitivityDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    const auto& bytes = tlsBitCache.bytes;
    
    // Simple avalanche test: compare hash of original vs hash with one bit changed
    uint32_t origHash = fnv1aHash(bytes);
    
    if (bytes.empty()) return 0.0;
    
    // Change middle byte, middle bit
    std::vector<uint8_t> modifiedBytes = bytes;
    size_t midIndex = bytes.size() / 2;
    modifiedBytes[midIndex] ^= 0x10; // Flip bit 4
    
    uint32_t modHash = fnv1aHash(modifiedBytes);
    
    // Count differing bits in hash
    uint32_t hashDiff = origHash ^ modHash;
    uint32_t diffBits = POPCOUNT_LUT[hashDiff & 0xFF] + 
                       POPCOUNT_LUT[(hashDiff >> 8) & 0xFF] + 
                       POPCOUNT_LUT[(hashDiff >> 16) & 0xFF] + 
                       POPCOUNT_LUT[(hashDiff >> 24) & 0xFF];
    
    return static_cast<double>(diffBits) / 32.0; // Normalize to [0,1]
}

// 63. Sequential XOR Difference
DECLARE_DESCRIPTOR(SequentialXorDiff, BitOperationsDescriptor, "Sum of sequential XOR differences between bytes")
DESCRIPTOR_DEPENDENCIES(SequentialXorDiff) { return {}; }
DescriptorResult SequentialXorDiffDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.length() < 2) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    const auto& bytes = tlsBitCache.bytes;
    
    uint32_t xorSum = 0;
    
    for (size_t i = 0; i < bytes.size() - 1; ++i) {
        xorSum += POPCOUNT_LUT[bytes[i] ^ bytes[i + 1]];
    }
    
    return static_cast<double>(xorSum) / (bytes.size() - 1);
}

// 64. Byte Delta Encoding Size
DECLARE_DESCRIPTOR(ByteDeltaEncodingSize, BitOperationsDescriptor, "Estimate of delta encoding efficiency")
DESCRIPTOR_DEPENDENCIES(ByteDeltaEncodingSize) { return {}; }
DescriptorResult ByteDeltaEncodingSizeDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.length() < 2) return static_cast<double>(smiles.length());
    
    tlsBitCache.computeFromString(smiles);
    const auto& bytes = tlsBitCache.bytes;
    
    // Calculate delta values
    std::vector<int8_t> deltas;
    deltas.reserve(bytes.size() - 1);
    
    for (size_t i = 1; i < bytes.size(); ++i) {
        deltas.push_back(static_cast<int8_t>(bytes[i] - bytes[i-1]));
    }
    
    // Estimate encoding size based on delta value range
    uint32_t totalBits = 8; // First byte stored as-is
    
    for (int8_t delta : deltas) {
        if (delta == 0) {
            totalBits += 1; // Just need a flag bit
        } else if (delta >= -4 && delta <= 3) {
            totalBits += 3; // 3 bits for small range
        } else if (delta >= -16 && delta <= 15) {
            totalBits += 5; // 5 bits for medium range
        } else if (delta >= -64 && delta <= 63) {
            totalBits += 7; // 7 bits for larger range
        } else {
            totalBits += 9; // Full delta plus flag
        }
    }
    
    return static_cast<double>(totalBits) / 8.0; // Return in bytes
}

// 65. Bit Correlation Coefficient
DECLARE_DESCRIPTOR(BitCorrelationCoef, BitOperationsDescriptor, "Correlation coefficient between even and odd bits")
DESCRIPTOR_DEPENDENCIES(BitCorrelationCoef) { return {}; }
DescriptorResult BitCorrelationCoefDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    if (bits.size() < 4) return 0.0;
    
    // Separate even and odd bits
    std::vector<bool> evenBits, oddBits;
    
    for (size_t i = 0; i < bits.size(); i += 2) {
        evenBits.push_back(bits[i]);
        if (i + 1 < bits.size()) {
            oddBits.push_back(bits[i + 1]);
        }
    }
    
    if (oddBits.empty()) return 0.0;
    
    // Calculate means
    double evenMean = 0.0, oddMean = 0.0;
    
    for (bool bit : evenBits) evenMean += bit ? 1.0 : 0.0;
    for (bool bit : oddBits) oddMean += bit ? 1.0 : 0.0;
    
    evenMean /= evenBits.size();
    oddMean /= oddBits.size();
    
    // Calculate correlation
    double numerator = 0.0, denom1 = 0.0, denom2 = 0.0;
    
    size_t minSize = std::min(evenBits.size(), oddBits.size());
    
    for (size_t i = 0; i < minSize; ++i) {
        double evenVal = evenBits[i] ? 1.0 : 0.0;
        double oddVal = oddBits[i] ? 1.0 : 0.0;
        
        double evenDiff = evenVal - evenMean;
        double oddDiff = oddVal - oddMean;
        
        numerator += evenDiff * oddDiff;
        denom1 += evenDiff * evenDiff;
        denom2 += oddDiff * oddDiff;
    }
    
    if (denom1 == 0.0 || denom2 == 0.0) return 0.0;
    
    return numerator / std::sqrt(denom1 * denom2);
}

// 66. Bit-level Run Length Encoding Size
DECLARE_DESCRIPTOR(BitRunLengthSize, BitOperationsDescriptor, "Estimate of bit-level run length encoding efficiency")
DESCRIPTOR_DEPENDENCIES(BitRunLengthSize) { return {}; }
DescriptorResult BitRunLengthSizeDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    if (bits.empty()) return 0.0;
    
    // Count runs
    std::vector<uint32_t> runLengths;
    
    bool currentBit = bits[0];
    uint32_t currentRun = 1;
    
    for (size_t i = 1; i < bits.size(); ++i) {
        if (bits[i] == currentBit) {
            currentRun++;
        } else {
            runLengths.push_back(currentRun);
            currentBit = bits[i];
            currentRun = 1;
        }
    }
    
    // Add the last run
    runLengths.push_back(currentRun);
    
    // Calculate encoding size
    uint32_t encodingSize = runLengths.size(); // One bit for each run to indicate 0/1
    
    for (uint32_t length : runLengths) {
        // Calculate bits needed to encode run length
        encodingSize += 1 + static_cast<uint32_t>(std::log2(length));
    }
    
    // Effectiveness ratio (smaller is better)
    return static_cast<double>(encodingSize) / bits.size();
}

// 67. Byte Mirroring Distance
DECLARE_DESCRIPTOR(ByteMirroringDistance, BitOperationsDescriptor, "Hamming distance of bytes from their bit-reversed versions")
DESCRIPTOR_DEPENDENCIES(ByteMirroringDistance) { return {}; }
DescriptorResult ByteMirroringDistanceDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bytes = tlsBitCache.bytes;
    
    if (bytes.empty()) return 0.0;
    
    // Pre-compute bit reversal lookup table
    static const std::array<uint8_t, 256> bitReversalTable = []() {
        std::array<uint8_t, 256> table{};
        for (int i = 0; i < 256; i++) {
            uint8_t reversed = 0;
            for (int j = 0; j < 8; j++) {
                reversed |= ((i >> j) & 1) << (7 - j);
            }
            table[i] = reversed;
        }
        return table;
    }();
    
    uint32_t totalDistance = 0;
    
    for (uint8_t byte : bytes) {
        uint8_t reversed = bitReversalTable[byte];
        uint8_t diff = byte ^ reversed;
        totalDistance += POPCOUNT_LUT[diff];
    }
    
    return static_cast<double>(totalDistance) / bytes.size();
}

// 68. Rolling Byte Checksum
DECLARE_DESCRIPTOR(RollingByteChecksum, BitOperationsDescriptor, "Adler-32 type rolling checksum")
DESCRIPTOR_DEPENDENCIES(RollingByteChecksum) { return {}; }
DescriptorResult RollingByteChecksumDescriptor::calculate(Context& context) const {
    const auto& smiles = context.getSmiles();
    if (smiles.empty()) return 0.0;
    
    tlsBitCache.computeFromString(smiles);
    const auto& bytes = tlsBitCache.bytes;
    
    uint32_t a = 1, b = 0; // Adler-32 algorithm components
    
    for (uint8_t byte : bytes) {
        a = (a + byte) % 65521; // Prime number for Adler-32
        b = (b + a) % 65521;
    }
    
    return static_cast<double>((b << 16) | a);
}

// 69. Zero Bit Ratio
DECLARE_DESCRIPTOR(ZeroBitRatio, BitOperationsDescriptor, "Ratio of 0 bits to total bits")
DESCRIPTOR_DEPENDENCIES(ZeroBitRatio) { return {}; }
DescriptorResult ZeroBitRatioDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    if (bits.empty()) return 0.0;
    
    uint32_t zeroBits = 0;
    for (bool bit : bits) {
        if (!bit) zeroBits++;
    }
    
    return static_cast<double>(zeroBits) / bits.size();
}

// 70. Bit-level Autocorrelation
DECLARE_DESCRIPTOR(BitAutocorrelation, BitOperationsDescriptor, "Autocorrelation of bit sequence with lag 1")
DESCRIPTOR_DEPENDENCIES(BitAutocorrelation) { return {}; }
DescriptorResult BitAutocorrelationDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    if (bits.size() < 2) return 0.0;
    
    double mean = 0.0;
    for (bool bit : bits) {
        mean += bit ? 1.0 : 0.0;
    }
    mean /= bits.size();
    
    double sumProduct = 0.0, sumSquares = 0.0;
    
    for (size_t i = 0; i < bits.size() - 1; ++i) {
        double val1 = bits[i] ? 1.0 : 0.0;
        double val2 = bits[i + 1] ? 1.0 : 0.0;
        
        sumProduct += (val1 - mean) * (val2 - mean);
    }
    
    for (size_t i = 0; i < bits.size(); ++i) {
        double val = bits[i] ? 1.0 : 0.0;
        sumSquares += (val - mean) * (val - mean);
    }
    
    return (sumSquares == 0.0) ? 0.0 : sumProduct / sumSquares;
}

// 71. Diagonal Bit Accumulation
DECLARE_DESCRIPTOR(DiagonalBitAccumulation, BitOperationsDescriptor, "Sum of bits where index equals value")
DESCRIPTOR_DEPENDENCIES(DiagonalBitAccumulation) { return {}; }
DescriptorResult DiagonalBitAccumulationDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bytes = tlsBitCache.bytes;
    
    if (bytes.empty()) return 0.0;
    
    uint32_t sum = 0;
    for (size_t i = 0; i < bytes.size(); ++i) {
        uint8_t byte = bytes[i];
        for (int j = 0; j < 8; ++j) {
            if ((byte & (1 << j)) && ((i * 8 + j) % 256 == byte)) {
                sum++;
            }
        }
    }
    
    return static_cast<double>(sum);
}

// 72. Bit-wise Majority Function
DECLARE_DESCRIPTOR(BitwiseMajorityFunction, BitOperationsDescriptor, "Majority vote from adjacent bits")
DESCRIPTOR_DEPENDENCIES(BitwiseMajorityFunction) { return {}; }
DescriptorResult BitwiseMajorityFunctionDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bits = tlsBitCache.bits;
    
    if (bits.size() < 3) return 0.0;
    
    uint32_t majorityCount = 0;
    
    for (size_t i = 1; i < bits.size() - 1; ++i) {
        // Check if the middle bit equals the majority of the 3-bit window
        bool majority = (bits[i-1] && bits[i]) || (bits[i] && bits[i+1]) || (bits[i-1] && bits[i+1]);
        if (bits[i] == majority) {
            majorityCount++;
        }
    }
    
    return static_cast<double>(majorityCount) / (bits.size() - 2);
}

// 73. LZCNT Pattern
DECLARE_DESCRIPTOR(LzcntPattern, BitOperationsDescriptor, "Pattern of leading zero counts")
DESCRIPTOR_DEPENDENCIES(LzcntPattern) { return {}; }
DescriptorResult LzcntPatternDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bytes = tlsBitCache.bytes;
    
    if (bytes.empty()) return 0.0;
    
    // Calculate entropy of leading zero counts
    std::array<uint32_t, 9> lzcntHist{}; // 0-8 leading zeros possible
    
    for (uint8_t byte : bytes) {
        uint8_t lzcnt = LEADING_ZEROS_LUT[byte];
        lzcntHist[lzcnt]++;
    }
    
    double entropy = 0.0;
    for (uint32_t count : lzcntHist) {
        if (count > 0) {
            double prob = static_cast<double>(count) / bytes.size();
            entropy -= prob * std::log2(prob);
        }
    }
    
    return entropy;
}

// 74. TZCNT Pattern
DECLARE_DESCRIPTOR(TzcntPattern, BitOperationsDescriptor, "Pattern of trailing zero counts")
DESCRIPTOR_DEPENDENCIES(TzcntPattern) { return {}; }
DescriptorResult TzcntPatternDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bytes = tlsBitCache.bytes;
    
    if (bytes.empty()) return 0.0;
    
    // Calculate entropy of trailing zero counts
    std::array<uint32_t, 9> tzcntHist{}; // 0-8 trailing zeros possible
    
    for (uint8_t byte : bytes) {
        uint8_t tzcnt = TRAILING_ZEROS_LUT[byte];
        tzcntHist[tzcnt]++;
    }
    
    double entropy = 0.0;
    for (uint32_t count : tzcntHist) {
        if (count > 0) {
            double prob = static_cast<double>(count) / bytes.size();
            entropy -= prob * std::log2(prob);
        }
    }
    
    return entropy;
}

// 75. Byte Value Geometric Mean
DECLARE_DESCRIPTOR(ByteGeometricMean, BitOperationsDescriptor, "Geometric mean of byte values")
DESCRIPTOR_DEPENDENCIES(ByteGeometricMean) { return {}; }
DescriptorResult ByteGeometricMeanDescriptor::calculate(Context& context) const {
    tlsBitCache.computeFromString(context.getSmiles());
    const auto& bytes = tlsBitCache.bytes;
    
    if (bytes.empty()) return 0.0;
    
    double logSum = 0.0;
    uint32_t validCount = 0;
    
    for (uint8_t byte : bytes) {
        if (byte > 0) {
            logSum += std::log(byte);
            validCount++;
        }
    }
    
    if (validCount == 0) return 0.0;
    return std::exp(logSum / validCount);
}

// Register the descriptors
void register_BitOperationsDescriptors() {
    auto& registry = DescriptorRegistry::getInstance();
    
    registry.registerDescriptor(std::make_shared<StringLengthDescriptor>());
    registry.registerDescriptor(std::make_shared<BitCountDescriptor>());
    registry.registerDescriptor(std::make_shared<BitDensityDescriptor>());
    registry.registerDescriptor(std::make_shared<AsciiSumDescriptor>());
    registry.registerDescriptor(std::make_shared<LongestBitRunDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteXorChecksumDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteEntropyDescriptor>());
    registry.registerDescriptor(std::make_shared<FirstByteValueDescriptor>());
    registry.registerDescriptor(std::make_shared<LastByteValueDescriptor>());
    registry.registerDescriptor(std::make_shared<OddEvenBitRatioDescriptor>());
    registry.registerDescriptor(std::make_shared<CharPalindromeCheckDescriptor>());
    registry.registerDescriptor(std::make_shared<BitTransitionCountDescriptor>());
    registry.registerDescriptor(std::make_shared<LeftRightBitBalanceDescriptor>());
    registry.registerDescriptor(std::make_shared<HammingWeightVarianceDescriptor>());
    registry.registerDescriptor(std::make_shared<FastHashDescriptor>());
    registry.registerDescriptor(std::make_shared<Crc8ChecksumDescriptor>());
    registry.registerDescriptor(std::make_shared<MostFrequentByteDescriptor>());
    registry.registerDescriptor(std::make_shared<LeastFrequentByteDescriptor>());
    registry.registerDescriptor(std::make_shared<UniqueByteCountDescriptor>());
    registry.registerDescriptor(std::make_shared<EvenIndexedBitCountDescriptor>());
    registry.registerDescriptor(std::make_shared<OddIndexedBitCountDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteWiseMinimumDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteWiseMaximumDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteWiseMedianDescriptor>());
    registry.registerDescriptor(std::make_shared<HighLowBitRatioDescriptor>());
    registry.registerDescriptor(std::make_shared<CharSequenceComplexityDescriptor>());
    registry.registerDescriptor(std::make_shared<AvgBitPopulationCountDescriptor>());
    registry.registerDescriptor(std::make_shared<AvgLeadingZeroCountDescriptor>());
    registry.registerDescriptor(std::make_shared<AvgTrailingZeroCountDescriptor>());
    registry.registerDescriptor(std::make_shared<PopcountVarianceDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteAndReductionDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteOrReductionDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteXorReductionDescriptor>());
    registry.registerDescriptor(std::make_shared<AsciiProductMod256Descriptor>());
    registry.registerDescriptor(std::make_shared<BitAlternationFreqDescriptor>());
    registry.registerDescriptor(std::make_shared<CharDistributionVarianceDescriptor>());
    registry.registerDescriptor(std::make_shared<LowHighNibbleRatioDescriptor>());
    registry.registerDescriptor(std::make_shared<CharPosWeightedSumDescriptor>());
    registry.registerDescriptor(std::make_shared<NibbleEntropyDescriptor>());
    registry.registerDescriptor(std::make_shared<CharFrequencyEntropyDescriptor>());
    registry.registerDescriptor(std::make_shared<AsciiRangeCoverageDescriptor>());
    registry.registerDescriptor(std::make_shared<NibbleHistogramBalanceDescriptor>());
    registry.registerDescriptor(std::make_shared<BinaryCompressionRatioDescriptor>());
    registry.registerDescriptor(std::make_shared<BitReversalDistanceDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteReversalDistanceDescriptor>());
    registry.registerDescriptor(std::make_shared<MurmurHash3Descriptor>());
    registry.registerDescriptor(std::make_shared<CharRunLengthDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteRunEntropyDescriptor>());
    registry.registerDescriptor(std::make_shared<SingleBitPatternMatchDescriptor>());
    registry.registerDescriptor(std::make_shared<BytePairTransitionsDescriptor>());
    registry.registerDescriptor(std::make_shared<BitDifferenceSumDescriptor>());
    registry.registerDescriptor(std::make_shared<FirstLastByteXorDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteMarkovPropertyDescriptor>());
    registry.registerDescriptor(std::make_shared<SwarPopcountDescriptor>());
    registry.registerDescriptor(std::make_shared<FastBitHistogramDescriptor>());
    registry.registerDescriptor(std::make_shared<BitDispersionIndexDescriptor>());
    registry.registerDescriptor(std::make_shared<HalfByteRotationHashDescriptor>());
    registry.registerDescriptor(std::make_shared<AvalanchePatternSensitivityDescriptor>());
    registry.registerDescriptor(std::make_shared<SequentialXorDiffDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteDeltaEncodingSizeDescriptor>());
    registry.registerDescriptor(std::make_shared<BitCorrelationCoefDescriptor>());
    registry.registerDescriptor(std::make_shared<BitRunLengthSizeDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteMirroringDistanceDescriptor>());
    registry.registerDescriptor(std::make_shared<RollingByteChecksumDescriptor>());
    registry.registerDescriptor(std::make_shared<ZeroBitRatioDescriptor>());
    registry.registerDescriptor(std::make_shared<BitAutocorrelationDescriptor>());
    registry.registerDescriptor(std::make_shared<DiagonalBitAccumulationDescriptor>());
    registry.registerDescriptor(std::make_shared<BitwiseMajorityFunctionDescriptor>());
    registry.registerDescriptor(std::make_shared<LzcntPatternDescriptor>());
    registry.registerDescriptor(std::make_shared<TzcntPatternDescriptor>());
    registry.registerDescriptor(std::make_shared<ByteGeometricMeanDescriptor>());
}


void register_StringLengthDescriptor() {
    auto descriptor = std::make_shared<StringLengthDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitCountDescriptor() {
    auto descriptor = std::make_shared<BitCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitDensityDescriptor() {
    auto descriptor = std::make_shared<BitDensityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_AsciiSumDescriptor() {
    auto descriptor = std::make_shared<AsciiSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_LongestBitRunDescriptor() {
    auto descriptor = std::make_shared<LongestBitRunDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteXorChecksumDescriptor() {
    auto descriptor = std::make_shared<ByteXorChecksumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteEntropyDescriptor() {
    auto descriptor = std::make_shared<ByteEntropyDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FirstByteValueDescriptor() {
    auto descriptor = std::make_shared<FirstByteValueDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_LastByteValueDescriptor() {
    auto descriptor = std::make_shared<LastByteValueDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_OddEvenBitRatioDescriptor() {
    auto descriptor = std::make_shared<OddEvenBitRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_CharPalindromeCheckDescriptor() {
    auto descriptor = std::make_shared<CharPalindromeCheckDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitTransitionCountDescriptor() {
    auto descriptor = std::make_shared<BitTransitionCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_LeftRightBitBalanceDescriptor() {
    auto descriptor = std::make_shared<LeftRightBitBalanceDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_HammingWeightVarianceDescriptor() {
    auto descriptor = std::make_shared<HammingWeightVarianceDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FastHashDescriptor() {
    auto descriptor = std::make_shared<FastHashDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_Crc8ChecksumDescriptor() {
    auto descriptor = std::make_shared<Crc8ChecksumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MostFrequentByteDescriptor() {
    auto descriptor = std::make_shared<MostFrequentByteDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_LeastFrequentByteDescriptor() {
    auto descriptor = std::make_shared<LeastFrequentByteDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_UniqueByteCountDescriptor() {
    auto descriptor = std::make_shared<UniqueByteCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_EvenIndexedBitCountDescriptor() {
    auto descriptor = std::make_shared<EvenIndexedBitCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_OddIndexedBitCountDescriptor() {
    auto descriptor = std::make_shared<OddIndexedBitCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteWiseMinimumDescriptor() {
    auto descriptor = std::make_shared<ByteWiseMinimumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteWiseMaximumDescriptor() {
    auto descriptor = std::make_shared<ByteWiseMaximumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteWiseMedianDescriptor() {
    auto descriptor = std::make_shared<ByteWiseMedianDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_HighLowBitRatioDescriptor() {
    auto descriptor = std::make_shared<HighLowBitRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_CharSequenceComplexityDescriptor() {
    auto descriptor = std::make_shared<CharSequenceComplexityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_AvgBitPopulationCountDescriptor() {
    auto descriptor = std::make_shared<AvgBitPopulationCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_AvgLeadingZeroCountDescriptor() {
    auto descriptor = std::make_shared<AvgLeadingZeroCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_AvgTrailingZeroCountDescriptor() {
    auto descriptor = std::make_shared<AvgTrailingZeroCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_PopcountVarianceDescriptor() {
    auto descriptor = std::make_shared<PopcountVarianceDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteAndReductionDescriptor() {
    auto descriptor = std::make_shared<ByteAndReductionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteOrReductionDescriptor() {
    auto descriptor = std::make_shared<ByteOrReductionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteXorReductionDescriptor() {
    auto descriptor = std::make_shared<ByteXorReductionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_AsciiProductMod256Descriptor() {
    auto descriptor = std::make_shared<AsciiProductMod256Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitAlternationFreqDescriptor() {
    auto descriptor = std::make_shared<BitAlternationFreqDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_CharDistributionVarianceDescriptor() {
    auto descriptor = std::make_shared<CharDistributionVarianceDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_LowHighNibbleRatioDescriptor() {
    auto descriptor = std::make_shared<LowHighNibbleRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_CharPosWeightedSumDescriptor() {
    auto descriptor = std::make_shared<CharPosWeightedSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_NibbleEntropyDescriptor() {
    auto descriptor = std::make_shared<NibbleEntropyDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_CharFrequencyEntropyDescriptor() {
    auto descriptor = std::make_shared<CharFrequencyEntropyDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_AsciiRangeCoverageDescriptor() {
    auto descriptor = std::make_shared<AsciiRangeCoverageDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NibbleHistogramBalanceDescriptor() {
    auto descriptor = std::make_shared<NibbleHistogramBalanceDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BinaryCompressionRatioDescriptor() {
    auto descriptor = std::make_shared<BinaryCompressionRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitReversalDistanceDescriptor() {
    auto descriptor = std::make_shared<BitReversalDistanceDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteReversalDistanceDescriptor() {
    auto descriptor = std::make_shared<ByteReversalDistanceDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_MurmurHash3Descriptor() {
    auto descriptor = std::make_shared<MurmurHash3Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_CharRunLengthDescriptor() {
    auto descriptor = std::make_shared<CharRunLengthDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteRunEntropyDescriptor() {
    auto descriptor = std::make_shared<ByteRunEntropyDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SingleBitPatternMatchDescriptor() {
    auto descriptor = std::make_shared<SingleBitPatternMatchDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BytePairTransitionsDescriptor() {
    auto descriptor = std::make_shared<BytePairTransitionsDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitDifferenceSumDescriptor() {
    auto descriptor = std::make_shared<BitDifferenceSumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FirstLastByteXorDescriptor() {
    auto descriptor = std::make_shared<FirstLastByteXorDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteMarkovPropertyDescriptor() {
    auto descriptor = std::make_shared<ByteMarkovPropertyDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SwarPopcountDescriptor() {
    auto descriptor = std::make_shared<SwarPopcountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_FastBitHistogramDescriptor() {
    auto descriptor = std::make_shared<FastBitHistogramDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitDispersionIndexDescriptor() {
    auto descriptor = std::make_shared<BitDispersionIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_HalfByteRotationHashDescriptor() {
    auto descriptor = std::make_shared<HalfByteRotationHashDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_AvalanchePatternSensitivityDescriptor() {
    auto descriptor = std::make_shared<AvalanchePatternSensitivityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SequentialXorDiffDescriptor() {
    auto descriptor = std::make_shared<SequentialXorDiffDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteDeltaEncodingSizeDescriptor() {
    auto descriptor = std::make_shared<ByteDeltaEncodingSizeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitCorrelationCoefDescriptor() {
    auto descriptor = std::make_shared<BitCorrelationCoefDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitRunLengthSizeDescriptor() {
    auto descriptor = std::make_shared<BitRunLengthSizeDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteMirroringDistanceDescriptor() {
    auto descriptor = std::make_shared<ByteMirroringDistanceDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_RollingByteChecksumDescriptor() {
    auto descriptor = std::make_shared<RollingByteChecksumDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ZeroBitRatioDescriptor() {
    auto descriptor = std::make_shared<ZeroBitRatioDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitAutocorrelationDescriptor() {
    auto descriptor = std::make_shared<BitAutocorrelationDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_DiagonalBitAccumulationDescriptor() {
    auto descriptor = std::make_shared<DiagonalBitAccumulationDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_BitwiseMajorityFunctionDescriptor() {
    auto descriptor = std::make_shared<BitwiseMajorityFunctionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_LzcntPatternDescriptor() {
    auto descriptor = std::make_shared<LzcntPatternDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_TzcntPatternDescriptor() {
    auto descriptor = std::make_shared<TzcntPatternDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_ByteGeometricMeanDescriptor() {
    auto descriptor = std::make_shared<ByteGeometricMeanDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_NonZeroByteDensityDescriptor() {
    auto descriptor = std::make_shared<NonZeroByteDensityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

} // namespace desfact