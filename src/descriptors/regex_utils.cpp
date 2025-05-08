#include "../common.hpp"
#include <re2/re2.h>
#include <string>
#include <algorithm>

// This file bridges C functions to C++ code that depends on them

namespace desfact {

// Forward declarations for C functions we need to wrap
extern "C" {
    int countMatches(const char* s, const char* pattern);
    int countChar(const char* s, char c);
    int countDigits(const char* s);
    int countSubstr(const char* text, const char* pattern);
    int countNNotPrecededByC(const char* s);
}

// C++ wrapper for countMatches that works with RE2
int countMatches(const std::string& smiles, const RE2& pattern) {
    re2::StringPiece input(smiles);
    int count = 0;
    std::string match;
    while (RE2::FindAndConsume(&input, pattern, &match)) {
        count++;
    }
    return count;
}

// C++ wrapper for countChar
int countChar(const std::string& s, char c) {
    return std::count(s.begin(), s.end(), c);
}

// C++ wrapper for countDigits
int countDigits(const std::string& s) {
    return std::count_if(s.begin(), s.end(), ::isdigit);
}

// C++ wrapper for countSubstr
int countSubstr(const std::string& text, const std::string& pattern) {
    if (pattern.empty()) return 0;
    int count = 0;
    size_t pos = 0;
    while ((pos = text.find(pattern, pos)) != std::string::npos) {
        ++count;
        pos += pattern.length();
    }
    return count;
}

// C++ wrapper for countNNotPrecededByC
int countNNotPrecededByC(const std::string& s) {
    int count = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == 'N' && (i == 0 || s[i-1] != 'C')) {
            ++count;
        }
    }
    return count;
}

} // namespace desfact 