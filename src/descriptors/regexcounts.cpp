#include "../common.hpp"
// #include <GraphMol/SmilesParse/SmilesWrite.h> // No longer needed here
#include <re2/re2.h>
#include <string>
#include <algorithm>
#include <memory>
#include <vector>

namespace desfact {

// Helper function to count regex matches in SMILES
int countMatches(const std::string& smiles, const RE2& pattern) {
    re2::StringPiece input(smiles);
    int count = 0;
    // Change back to std::string for compatibility with RE2::FindAndConsume
    std::string match; // Use std::string instead of std::string_view
    while (RE2::FindAndConsume(&input, pattern, &match)) { // Pass address of std::string
        count++;
    }
    return count;
}

// Helper function to count characters in a string
int countChar(const std::string& s, char c) {
    return std::count(s.begin(), s.end(), c);
}

// Helper function to count digits in a string
int countDigits(const std::string& s) {
    return std::count_if(s.begin(), s.end(), ::isdigit);
}

// Helper function to count non-overlapping substrings
int countSubstr(const std::string& text, const std::string& pattern) {
    if (pattern.empty()) return 0; // Avoid infinite loop
    int count = 0;
    size_t pos = 0;
    while ((pos = text.find(pattern, pos)) != std::string::npos) {
        ++count;
        pos += pattern.length(); // Move past the found pattern
    }
    return count;
}

int countNNotPrecededByC(const std::string& s) {
    int count = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == 'N' && (i == 0 || s[i-1] != 'C')) {
            ++count;
        }
    }
    return count;
}

// Carbon count
DECLARE_DESCRIPTOR(SmilesCarbonCount, RegexCounts, "Count of carbon atoms (C, c) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesCarbonCount) { return {}; }
DescriptorResult SmilesCarbonCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Count both C and c for total carbon
    return countChar(smiles, 'C') + countChar(smiles, 'c');
}

// Nitrogen count
DECLARE_DESCRIPTOR(SmilesNitrogenCount, RegexCounts, "Count of nitrogen atoms (N, n) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesNitrogenCount) { return {}; }
DescriptorResult SmilesNitrogenCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Count both N and n for total nitrogen
    return countChar(smiles, 'N') + countChar(smiles, 'n');
}

// Oxygen count
DECLARE_DESCRIPTOR(SmilesOxygenCount, RegexCounts, "Count of oxygen atoms (O, o) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesOxygenCount) { return {}; }
DescriptorResult SmilesOxygenCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Count both O and o for total oxygen
    return countChar(smiles, 'O') + countChar(smiles, 'o');
}

// Sulfur count
DECLARE_DESCRIPTOR(SmilesSulfurCount, RegexCounts, "Count of sulfur atoms (S, s) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSulfurCount) { return {}; }
DescriptorResult SmilesSulfurCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Count both S and s for total sulfur
    return countChar(smiles, 'S') + countChar(smiles, 's');
}

// Phosphorus count
DECLARE_DESCRIPTOR(SmilesPhosphorusCount, RegexCounts, "Count of phosphorus atoms (P, p) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPhosphorusCount) { return {}; }
DescriptorResult SmilesPhosphorusCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Count both P and p for total phosphorus
    return countChar(smiles, 'P') + countChar(smiles, 'p');
}

// Halogen count
DECLARE_DESCRIPTOR(SmilesHalogenCount, RegexCounts, "Count of halogen atoms (F, Cl, Br, I, f, cl, br, i) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesHalogenCount) { return {}; }
DescriptorResult SmilesHalogenCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    int count = 0;
    count += countChar(smiles, 'F');
    count += countChar(smiles, 'I');
    count += countChar(smiles, 'f'); // Add lowercase aromatic
    count += countChar(smiles, 'i'); // Add lowercase aromatic
    count += countSubstr(smiles, "Cl");
    count += countSubstr(smiles, "Br");
    count += countSubstr(smiles, "cl"); // Add lowercase aromatic
    count += countSubstr(smiles, "br"); // Add lowercase aromatic
    return count;
}

// Boron count
DECLARE_DESCRIPTOR(SmilesBoronCount, RegexCounts, "Count of boron atoms in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBoronCount) { return {}; }
DescriptorResult SmilesBoronCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, 'B');
}

// Silicon count
DECLARE_DESCRIPTOR(SmilesSiliconCount, RegexCounts, "Count of silicon atoms in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSiliconCount) { return {}; }
DescriptorResult SmilesSiliconCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countSubstr(smiles, "Si");
}

// Fluorine count
DECLARE_DESCRIPTOR(SmilesFluorineCount, RegexCounts, "Count of fluorine atoms in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesFluorineCount) { return {}; }
DescriptorResult SmilesFluorineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, 'F');
}

// Chlorine count
DECLARE_DESCRIPTOR(SmilesChlorineCount, RegexCounts, "Count of chlorine atoms in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesChlorineCount) { return {}; }
DescriptorResult SmilesChlorineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countSubstr(smiles, "Cl");
}

// Bromine count
DECLARE_DESCRIPTOR(SmilesBromineCount, RegexCounts, "Count of bromine atoms in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBromineCount) { return {}; }
DescriptorResult SmilesBromineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countSubstr(smiles, "Br");
}

// Iodine count
DECLARE_DESCRIPTOR(SmilesIodineCount, RegexCounts, "Count of iodine atoms in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIodineCount) { return {}; }
DescriptorResult SmilesIodineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, 'I');
}

// Aromatic carbon count
DECLARE_DESCRIPTOR(SmilesAromaticCarbonCount, RegexCounts, "Count of aromatic carbon atoms 'c' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticCarbonCount) { return {}; }
DescriptorResult SmilesAromaticCarbonCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, 'c');
}

// Aromatic nitrogen count
DECLARE_DESCRIPTOR(SmilesAromaticNitrogenCount, RegexCounts, "Count of aromatic nitrogen atoms 'n' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticNitrogenCount) { return {}; }
DescriptorResult SmilesAromaticNitrogenCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, 'n');
}

// Aromatic oxygen count
DECLARE_DESCRIPTOR(SmilesAromaticOxygenCount, RegexCounts, "Count of aromatic oxygen atoms 'o' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticOxygenCount) { return {}; }
DescriptorResult SmilesAromaticOxygenCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, 'o');
}

// Aromatic sulfur count
DECLARE_DESCRIPTOR(SmilesAromaticSulfurCount, RegexCounts, "Count of aromatic sulfur atoms 's' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticSulfurCount) { return {}; }
DescriptorResult SmilesAromaticSulfurCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, 's');
}


// Aromatic halogen count
DECLARE_DESCRIPTOR(SmilesAromaticHalogenCount, RegexCounts, "Count of aromatic halogen atoms (f, cl, br, i) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticHalogenCount) { return {}; }
DescriptorResult SmilesAromaticHalogenCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    int count = 0;
    count += countChar(smiles, 'f');
    count += countChar(smiles, 'i');
    // Use RE2 for "cl" and "br" as substrings to avoid counting 'c' or 'b' in other contexts
    static const RE2 pattern("(cl|br)");
    re2::StringPiece input(smiles);
    std::string match;
    while (RE2::FindAndConsume(&input, pattern, &match)) {
        count++;
    }
    // Also count bracketed aromatic halogens explicitly marked as aromatic (less common)
    static const RE2 bracket_pattern("\\[(F|Cl|Br|I);a\\]"); // Match [X;a]
    re2::StringPiece bracket_input(smiles);
    while(RE2::FindAndConsume(&bracket_input, bracket_pattern, &match)){
        count++;
    }
    return count;
}

// Single bond count (implicit bonds are not counted here)
DECLARE_DESCRIPTOR(SmilesSingleBondCount, RegexCounts, "Count of explicit single bonds '-' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSingleBondCount) { return {}; }
DescriptorResult SmilesSingleBondCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, '-');
}

// Double bond count (No change needed)
DECLARE_DESCRIPTOR(SmilesDoubleBondCount, RegexCounts, "Count of double bonds '=' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesDoubleBondCount) { return {}; }
DescriptorResult SmilesDoubleBondCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, '=');
}

// Triple bond count (No change needed)
DECLARE_DESCRIPTOR(SmilesTripleBondCount, RegexCounts, "Count of triple bonds '#' in SMILES")
DESCRIPTOR_DEPENDENCIES(SmilesTripleBondCount) { return {}; }
DescriptorResult SmilesTripleBondCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("#");
    return countMatches(smiles, pattern);
}

// Aromatic bond count (Approximation: counts lowercase letters often involved)
DECLARE_DESCRIPTOR(SmilesAromaticBondCount, RegexCounts, "Approx. count of aromatic bonds based on lowercase atoms in SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticBondCount) { return {}; }
DescriptorResult SmilesAromaticBondCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Counts lowercase letters typically involved in aromatic systems + explicit ':'
    static const RE2 pattern("([cnospfib])"); // More inclusive set
    int count = countMatches(smiles, pattern);
    count += countChar(smiles, ':'); // Add explicit aromatic bonds
    // This is a rough count, as implicit aromatic bonds are not directly represented
    return count;
}

// Ring closure count
DECLARE_DESCRIPTOR(SmilesRingClosureCount, RegexCounts, "Count of ring closures in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesRingClosureCount) { return {}; }
DescriptorResult SmilesRingClosureCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countDigits(smiles);
}

// Branch open count
DECLARE_DESCRIPTOR(SmilesBranchOpenCount, RegexCounts, "Count of branch openings '(' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBranchOpenCount) { return {}; }
DescriptorResult SmilesBranchOpenCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, '(');
}

// Branch close count
DECLARE_DESCRIPTOR(SmilesBranchCloseCount, RegexCounts, "Count of branch closings ')' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBranchCloseCount) { return {}; }
DescriptorResult SmilesBranchCloseCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, ')');
}

// Square bracket atom count
DECLARE_DESCRIPTOR(SmilesBracketAtomCount, RegexCounts, "Count of bracketed atoms in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBracketAtomCount) { return {}; }
DescriptorResult SmilesBracketAtomCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, '[');
}

// Chiral center count
DECLARE_DESCRIPTOR(SmilesChiralCenterCount, RegexCounts, "Count of chiral centers '@' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesChiralCenterCount) { return {}; }
DescriptorResult SmilesChiralCenterCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, '@');
}

// Positive charge count
DECLARE_DESCRIPTOR(SmilesPositiveChargeCount, RegexCounts, "Count of positive charges '+' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPositiveChargeCount) { return {}; }
DescriptorResult SmilesPositiveChargeCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, '+');
}

// Negative charge count
DECLARE_DESCRIPTOR(SmilesNegativeChargeCount, RegexCounts, "Count of negative charges '-' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesNegativeChargeCount) { return {}; }
DescriptorResult SmilesNegativeChargeCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, '-');
}

// Stereochemistry forward count
DECLARE_DESCRIPTOR(SmilesStereoForwardCount, RegexCounts, "Count of '/' stereochemistry in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesStereoForwardCount) { return {}; }
DescriptorResult SmilesStereoForwardCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, '/');
}

// Stereochemistry backward count
DECLARE_DESCRIPTOR(SmilesStereoBackwardCount, RegexCounts, "Count of '\\' stereochemistry in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesStereoBackwardCount) { return {}; }
DescriptorResult SmilesStereoBackwardCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, '\\');
}

// Dot disconnected count
DECLARE_DESCRIPTOR(SmilesDotCount, RegexCounts, "Count of '.' disconnected fragments in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesDotCount) { return {}; }
DescriptorResult SmilesDotCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    return countChar(smiles, '.');
}

// --- Start of complex patterns (using RE2) ---

// Methyl group: C, [CH3], C(H)(H)H, [C;H3], [C@@H3], [C@H3], [C;H3;!R], etc.
DECLARE_DESCRIPTOR(SmilesMethylCount, RegexCounts, "Count of methyl groups (C, [CH3], C(H)(H)H, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesMethylCount) { return {}; }
DescriptorResult SmilesMethylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C([^a-zA-Z0-9]|$)|"         // terminal C
        "\\[CH3\\]|"                 // explicit [CH3]
        "C\\(H\\)\\(H\\)H|"          // C(H)(H)H
        "\\[C;H3\\]|"                // [C;H3]
        "\\[C@@H3\\]|"               // [C@@H3]
        "\\[C@H3\\]|"                // [C@H3]
        "\\[C;H3;!R\\]"              // [C;H3;!R]
    );
    return countMatches(smiles, pattern);
}

// Ethyl group: CC, [CH2][CH3], C(C)C, [CH2CH3], [CH2][C;H3], C[CH2], [C;H2][C;H3], [CH2]C, C[CH3], [CH2][CH2][CH3]
DECLARE_DESCRIPTOR(SmilesEthylCount, RegexCounts, "Count of ethyl groups (CC, [CH2][CH3], C(C)C, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesEthylCount) { return {}; }
DescriptorResult SmilesEthylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "CC([^a-zA-Z0-9]|$)|"        // terminal CC
        "\\[CH2\\]\\[CH3\\]|"        // [CH2][CH3]
        "C\\(C\\)C|"                 // C(C)C
        "\\[CH2CH3\\]|"              // [CH2CH3]
        "\\[CH2\\]\\[C;H3\\]|"       // [CH2][C;H3]
        "C\\[CH2\\]|"                // C[CH2]
        "\\[C;H2\\]\\[C;H3\\]|"      // [C;H2][C;H3]
        "\\[CH2\\]C|"                // [CH2]C
        "C\\[CH3\\]|"                // C[CH3]
        "\\[CH2\\]\\[CH2\\]\\[CH3\\]"// [CH2][CH2][CH3]
    );
    return countMatches(smiles, pattern);
}

// Amino group: N, [NH2], N([H])[H], [N;H2], [NH2+], [NH2-], [N;H2;!R], N(H)H, [NH3], [NH3+], [NH], [N;H1], N(H), [N;H], [N;H3], [N;H0], [N;H1;R0], [N;H2;R1], N[H], N([H])
DECLARE_DESCRIPTOR(SmilesAminoCount, RegexCounts, "Count of amino groups (N, [NH2], N([H])[H], etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAminoCount) { return {}; }
DescriptorResult SmilesAminoCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "N([^a-zA-Z0-9]|$)|"         // terminal N
        "\\[NH2\\]|"                 // explicit [NH2]
        "N\\(\\[H\\]\\)\\[H\\]|"     // N([H])[H]
        "\\[N;H2\\]|"                // [N;H2]
        "\\[NH2\\+\\]|"              // [NH2+]
        "\\[NH2-\\]|"                // [NH2-]
        "\\[N;H2;!R\\]|"             // [N;H2;!R]
        "N\\(H\\)H|"                 // N(H)H
        "\\[NH3\\]|"                 // [NH3]
        "\\[NH3\\+\\]|"              // [NH3+]
        "\\[NH\\]|"                  // [NH]
        "\\[N;H1\\]|"                // [N;H1]
        "N\\(H\\)|"                  // N(H)
        "\\[N;H\\]|"                 // [N;H]
        "\\[N;H3\\]|"                // [N;H3]
        "\\[N;H0\\]|"                // [N;H0]
        "\\[N;H1;R0\\]|"             // [N;H1;R0]
        "\\[N;H2;R1\\]|"             // [N;H2;R1]
        "N\\[H\\]|"                  // N[H]
        "N\\(\\[H\\]\\)"             // N([H])
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesEsterCount, RegexCounts, "Count of ester groups (C(=O)OC) in SMILES")
DESCRIPTOR_DEPENDENCIES(SmilesEsterCount) { return {}; }
DescriptorResult SmilesEsterCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Simpler pattern to match common ester representations
    static const RE2 pattern("C=?[\\(]?=?O[\\)]?O[Cc]");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesEtherCount, RegexCounts, "Count of ether groups (COC, [CH2]O[CH3], etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesEtherCount) { return {}; }
DescriptorResult SmilesEtherCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "COC|"                        // COC
        "\\[CH2\\]O\\[CH3\\]|"        // [CH2]O[CH3]
        "\\[CH2\\]O\\[CH2\\]|"        // [CH2]O[CH2]
        "\\[COC\\]|"                  // [COC]
        "\\[C;H2\\]O\\[C;H3\\]|"      // [C;H2]O[C;H3]
        "CO\\[CH3\\]|"                // CO[CH3]
        "\\[CH2\\]OC|"                // [CH2]OC
        "\\[C;H2\\]O\\[C;H2\\]|"      // [C;H2]O[C;H2]
        "COC\\(C\\)|"                 // COC(C)
        "C\\(O\\)C|"                  // C(O)C
        "\\[CH2\\]O|"                 // [CH2]O
        "O\\[CH2\\]|"                 // O[CH2]
        "\\[C;H2\\]O|"                // [C;H2]O
        "O\\[C;H2\\]|"                // O[C;H2]
        "COEt|"                       // COEt
        "EtOC|"                       // EtOC
        "COC\\[CH3\\]|"               // COC[CH3]
        "\\[CH2\\]OEt|"               // [CH2]OEt
        "EtO\\[CH2\\]"                // EtO[CH2]
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesThiolCount, RegexCounts, "Count of thiol groups (S[H], [SH], etc) in SMILES")
DESCRIPTOR_DEPENDENCIES(SmilesThiolCount) { return {}; }
DescriptorResult SmilesThiolCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "S[H]|"                      // S followed by H
        "\\[SH\\]|"                  // explicit [SH]
        "\\[S;H1\\]|"                // [S;H1]
        "\\[SH-\\]|"                 // [SH-]
        "\\[S-\\]|"                  // [S-]
        "\\[S;H1;!R\\]|"             // [S;H1;!R]
        "S\\(H\\)|"                  // S(H)
        "\\[SH2\\]|"                 // [SH2]
        "\\[S;H2\\]|"                // [S;H2]
        "S[H2]|"                     // S[H2]
        "S\\(H\\)H"                  // S(H)H
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesThioetherCount, RegexCounts, "Count of thioether groups (CSC, [CH2]S[CH3], etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesThioetherCount) { return {}; }
DescriptorResult SmilesThioetherCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "CSC|"                       // CSC
        "\\[CH2\\]S\\[CH3\\]|"       // [CH2]S[CH3]
        "\\[CH2\\]S\\[CH2\\]|"       // [CH2]S[CH2]
        "\\[CSC\\]|"                 // [CSC]
        "\\[C;H2\\]S\\[C;H3\\]|"     // [C;H2]S[C;H3]
        "CS\\[CH3\\]|"               // CS[CH3]
        "\\[CH2\\]SC|"               // [CH2]SC
        "\\[C;H2\\]S\\[C;H2\\]|"     // [C;H2]S[C;H2]
        "S\\(C\\)C|"                 // S(C)C
        "C\\(S\\)C|"                 // C(S)C
        "\\[CH2\\]S|"                // [CH2]S
        "S\\[CH2\\]|"                // S[CH2]
        "\\[C;H2\\]S|"               // [C;H2]S
        "S\\[C;H2\\]|"               // S[C;H2]
        "CSS|"                       // CSS
        "SCS"                        // SCS
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAmideCount, RegexCounts, "Count of amide groups (C(=O)N, NC(=O), etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAmideCount) { return {}; }
DescriptorResult SmilesAmideCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C\\(=O\\)N|"                // C(=O)N
        "NC\\(=O\\)|"                // NC(=O)
        "\\[C\\(=O\\)N\\]|"          // [C(=O)N]
        "\\[NC\\(=O\\)\\]|"          // [NC(=O)]
        "\\[C;H0;X3\\]\\(=O\\)N|"    // [C;H0;X3](=O)N
        "\\[C\\]\\(=O\\)N|"          // [C](=O)N
        "C\\(=O\\)\\[NH2\\]|"        // C(=O)[NH2]
        "C\\(=O\\)\\[NH\\]|"         // C(=O)[NH]
        "C\\(=O\\)N\\(\\[H\\]\\)|"   // C(=O)N([H])
        "C\\(=O\\)N\\[H\\]"          // C(=O)N[H]
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesNitrileCount, RegexCounts, "Count of nitrile groups (C#N, [C]#[N], etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesNitrileCount) { return {}; }
DescriptorResult SmilesNitrileCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C#N|"                       // C#N
        "\\[C\\]#\\[N\\]|"           // [C]#[N]
        "\\[C;H0\\]#\\[N\\]|"        // [C;H0]#[N]
        "\\[C;H0\\]#\\[N;H0\\]|"     // [C;H0]#[N;H0]
        "CN#N|"                      // CN#N
        "C#N\\[H\\]"                 // C#N[H]
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesSulfonylCount, RegexCounts, "Count of sulfonyl groups (S(=O)(=O), [S](=O)(=O), etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSulfonylCount) { return {}; }
DescriptorResult SmilesSulfonylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "S\\(=O\\)\\(=O\\)|"         // S(=O)(=O)
        "\\[S\\]\\(=O\\)\\(=O\\)|"   // [S](=O)(=O)
        "\\[S;X4\\]\\(=O\\)\\(=O\\)|"// [S;X4](=O)(=O)
        "S\\(=O\\)2|"                // S(=O)2
        "S\\(=O\\)O|"                // S(=O)O
        "S\\(=O\\)\\(O\\)"           // S(=O)(O)
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesSulfoxideCount, RegexCounts, "Count of sulfoxide groups (S(=O), [S](=O), etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSulfoxideCount) { return {}; }
DescriptorResult SmilesSulfoxideCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "S\\(=O\\)[^\\(=]|S\\(=O\\)$|" // S(=O) not followed by (=O)
        "\\[S\\]\\(=O\\)[^\\(=]|\\[S\\]\\(=O\\)$|" // [S](=O)
        "\\[S;X3\\]\\(=O\\)[^\\(=]|\\[S;X3\\]\\(=O\\)$|" // [S;X3](=O)
        "\\[SO\\]|"                  // [SO]
        "S=O|"                       // S=O
        "S\\(=O\\)\\[H\\]|"          // S(=O)[H]
        "S\\(=O\\)C"                 // S(=O)C
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAlkyneCount, RegexCounts, "Count of alkynes (C#C, [C]#[C], etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAlkyneCount) { return {}; }
DescriptorResult SmilesAlkyneCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C#C|"                       // C#C
        "\\[C\\]#\\[C\\]|"           // [C]#[C]
        "\\[C;H\\]#\\[C;H\\]"        // [C;H]#[C;H]
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAlkeneCount, RegexCounts, "Count of alkenes (C=C, [C]=[C], etc) in SMILES")
DESCRIPTOR_DEPENDENCIES(SmilesAlkeneCount) { return {}; }
DescriptorResult SmilesAlkeneCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C=C|"                       // C=C
        "\\[C\\]=\\[C\\]|"           // [C]=[C]
        "\\[C;H\\]=\\[C;H\\]"        // [C;H]=[C;H]
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesFiveRingCount, RegexCounts, "Count of 5-membered rings in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesFiveRingCount) { return {}; }
DescriptorResult SmilesFiveRingCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C1CCCC1|"                   // cyclopentane
        "c1cccc1|"                   // cyclopentadiene/aromatic
        "N1CCCC1|"                   // pyrrolidine
        "O1CCCC1|"                   // oxolane
        "S1CCCC1"                    // thiolane
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesSixRingCount, RegexCounts, "Count of 6-membered rings in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSixRingCount) { return {}; }
DescriptorResult SmilesSixRingCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C1CCCCC1|"                  // cyclohexane
        "c1ccccc1|"                  // benzene
        "N1CCCCC1|"                  // piperidine
        "O1CCCCC1|"                  // oxane
        "S1CCCCC1"                   // thiane
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesSevenRingCount, RegexCounts, "Count of 7-membered rings in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSevenRingCount) { return {}; }
DescriptorResult SmilesSevenRingCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Regex for ring closures is unreliable for *size* counting
    static const RE2 pattern("1.{5,6}1");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesEightRingCount, RegexCounts, "Count of 8-membered rings in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesEightRingCount) { return {}; }
DescriptorResult SmilesEightRingCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Regex for ring closures is unreliable for *size* counting
    static const RE2 pattern("1.{6,7}1");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIsotopeCount, RegexCounts, "Count of isotopic atoms ([13C], [2H], etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIsotopeCount) { return {}; }
DescriptorResult SmilesIsotopeCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[[0-9]+[A-Za-z]{1,2}"); // Match isotope label inside brackets
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesTertiaryAmineCount, RegexCounts, "Count of tertiary amines (N(C)(C)C) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesTertiaryAmineCount) { return {}; }
DescriptorResult SmilesTertiaryAmineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
     // This is an approximation, real check needs graph traversal
    static const RE2 pattern("N\\([^\\)]+\\)\\(.+\\)."); // N followed by branches
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesQuaternaryAmmoniumCount, RegexCounts, "Count of quaternary ammonium ions ([N+](C)(C)(C)C) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesQuaternaryAmmoniumCount) { return {}; }
DescriptorResult SmilesQuaternaryAmmoniumCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[N\\+\\]\\("); // [N+] followed by branches
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesPrimaryAmineCount, RegexCounts, "Count of primary amines (N) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPrimaryAmineCount) { return {}; }
DescriptorResult SmilesPrimaryAmineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Match N[H2] or N not followed by '(', lowercase, or digit
    static const RE2 pattern("N\\[H2\\]|N([^a-z0-9\\(]|$)");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesSecondaryAmineCount, RegexCounts, "Count of secondary amines (N(C)C) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSecondaryAmineCount) { return {}; }
DescriptorResult SmilesSecondaryAmineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Match N[H](...) or N(...) followed by non-branch
     static const RE2 pattern("N\\[H\\]\\(|N\\([^\\)]+\\)([^\\(]|$)");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesSulfonamideCount, RegexCounts, "Count of sulfonamide groups (S(=O)2N) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSulfonamideCount) { return {}; }
DescriptorResult SmilesSulfonamideCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("S\\(=O\\)\\(=O\\)N"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesPhosphonateCount, RegexCounts, "Count of phosphonate groups (P(=O)(O)C) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPhosphonateCount) { return {}; }
DescriptorResult SmilesPhosphonateCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("P\\(=O\\)\\(O\\)C"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAlkylHalideCount, RegexCounts, "Count of alkyl halides (C[ClBrIF]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAlkylHalideCount) { return {}; }
DescriptorResult SmilesAlkylHalideCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("C(\\[Cl\\]|\\[Br\\]|\\[I\\]|F)"); // Match C connected to bracketed Cl/Br/I or F
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAliphaticNitroCount, RegexCounts, "Count of aliphatic nitro groups (C[N+](=O)[O-]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAliphaticNitroCount) { return {}; }
DescriptorResult SmilesAliphaticNitroCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("[CX4]\\[N\\+\\]\\(=O\\)\\[O-\\]"); // Match C not 'c'
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAzideCount, RegexCounts, "Count of azide groups (N=[N+]=[N-]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAzideCount) { return {}; }
DescriptorResult SmilesAzideCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("N=\\[N\\+\\]=\\[N-\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesDiazoCount, RegexCounts, "Count of diazo groups (N=N) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesDiazoCount) { return {}; }
DescriptorResult SmilesDiazoCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Replace negative lookahead with proper RE2 pattern
    static const RE2 pattern("N=N([^=]|$)"); // N=N not followed by =
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesThiocyanateCount, RegexCounts, "Count of thiocyanate groups (SCN) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesThiocyanateCount) { return {}; }
DescriptorResult SmilesThiocyanateCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("SC#N"); // Explicit triple bond
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIsothiocyanateLinearCount, RegexCounts, "Count of isothiocyanate groups (NCS) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIsothiocyanateLinearCount) { return {}; }
DescriptorResult SmilesIsothiocyanateLinearCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("N=C=S"); // Explicit double bonds
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesGuanidineCount, RegexCounts, "Count of guanidine groups (N=C(N)N) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesGuanidineCount) { return {}; }
DescriptorResult SmilesGuanidineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("N=C\\(N\\)N"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesUreaCount, RegexCounts, "Count of urea patterns (NC(=O)N) in SMILES")
DESCRIPTOR_DEPENDENCIES(SmilesUreaCount) { return {}; }
DescriptorResult SmilesUreaCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("NC\\(=O\\)N|NC=ON");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAmidineCount, RegexCounts, "Count of amidine groups (C(=NH)NH2) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAmidineCount) { return {}; }
DescriptorResult SmilesAmidineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("C\\(=N\\)N"); // General amidine C(=N)N
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesImidazoleCount, RegexCounts, "Count of imidazole rings (c1nccn1) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesImidazoleCount) { return {}; }
DescriptorResult SmilesImidazoleCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("c1nccn1"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesPyridineCount, RegexCounts, "Count of pyridine rings (c1ccncc1) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPyridineCount) { return {}; }
DescriptorResult SmilesPyridineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("c1ccncc1"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesPyrroleCount, RegexCounts, "Count of pyrrole rings (c1cc[nH]c1) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPyrroleCount) { return {}; }
DescriptorResult SmilesPyrroleCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("c1cc\\[nH\\]c1"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesFuranCount, RegexCounts, "Count of furan rings (c1ccoc1) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesFuranCount) { return {}; }
DescriptorResult SmilesFuranCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("c1ccoc1"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesThiopheneCount, RegexCounts, "Count of thiophene rings (c1ccsc1) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesThiopheneCount) { return {}; }
DescriptorResult SmilesThiopheneCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("c1ccsc1"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAcidicSulfurCount, RegexCounts, "Count of acidic sulfurs (S in SO3H, SO2H) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAcidicSulfurCount) { return {}; }
DescriptorResult SmilesAcidicSulfurCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
     // Match S(=O)OH or S(=O)(=O)OH
    static const RE2 pattern("S\\(=O\\)(\\(=O\\))?O(\\[H\\])?");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAcidicPhosphorusCount, RegexCounts, "Count of acidic phosphorus (P in PO4H2, PO3H) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAcidicPhosphorusCount) { return {}; }
DescriptorResult SmilesAcidicPhosphorusCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Match P(=O)(O)O(H) or P(=O)(O)OH
    static const RE2 pattern("P\\(=O\\)\\(O\\)O(\\[H\\])?");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesPhenolCount, RegexCounts, "Count of phenol groups (cO) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPhenolCount) { return {}; }
DescriptorResult SmilesPhenolCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("cO(\\[H\\])?"); // Aromatic carbon attached to OH
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesCarboxylateCount, RegexCounts, "Count of carboxylate anions (C(=O)[O-]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesCarboxylateCount) { return {}; }
DescriptorResult SmilesCarboxylateCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("C\\(=O\\)\\[O-\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesSulfonateCount, RegexCounts, "Count of sulfonate anions (S(=O)(=O)[O-]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSulfonateCount) { return {}; }
DescriptorResult SmilesSulfonateCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("S\\(=O\\)\\(=O\\)\\[O-\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesPhosphateAnionCount, RegexCounts, "Count of phosphate anions (P(=O)(O)[O-]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPhosphateAnionCount) { return {}; }
DescriptorResult SmilesPhosphateAnionCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("P\\(=O\\)\\(O\\)\\[O-\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesImidazoleBasicCount, RegexCounts, "Count of basic imidazole rings (c1nccn1) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesImidazoleBasicCount) { return {}; }
DescriptorResult SmilesImidazoleBasicCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("c1nccn1"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesPyridineBasicCount, RegexCounts, "Count of basic pyridine rings (c1ccncc1) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPyridineBasicCount) { return {}; }
DescriptorResult SmilesPyridineBasicCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("c1ccncc1"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesPrimaryAmineBasicCount, RegexCounts, "Count of primary amines (N) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPrimaryAmineBasicCount) { return {}; }
DescriptorResult SmilesPrimaryAmineBasicCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
     // Same as SmilesPrimaryAmineCount
    static const RE2 pattern("N\\[H2\\]|N([^a-z0-9\\(]|$)");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesSecondaryAmineBasicCount, RegexCounts, "Count of secondary amines (N(C)C) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSecondaryAmineBasicCount) { return {}; }
DescriptorResult SmilesSecondaryAmineBasicCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
     // Same as SmilesSecondaryAmineCount
    static const RE2 pattern("N\\[H\\]\\(|N\\([^\\)]+\\)([^\\(]|$)");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesTertiaryAmineBasicCount, RegexCounts, "Count of tertiary amines (N(C)(C)C) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesTertiaryAmineBasicCount) { return {}; }
DescriptorResult SmilesTertiaryAmineBasicCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
     // Same as SmilesTertiaryAmineCount
    static const RE2 pattern("N\\([^\\)]+\\)\\(.+\\).");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesGuanidineBasicCount, RegexCounts, "Count of guanidine groups (N=C(N)N) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesGuanidineBasicCount) { return {}; }
DescriptorResult SmilesGuanidineBasicCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("N=C\\(N\\)N"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAmidineBasicCount, RegexCounts, "Count of amidine groups (C(=NH)NH2) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAmidineBasicCount) { return {}; }
DescriptorResult SmilesAmidineBasicCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("C\\(=N\\)N"); // General amidine C(=N)N
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesQuaternaryAmmoniumBasicCount, RegexCounts, "Count of quaternary ammonium ions ([N+](C)(C)(C)C) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesQuaternaryAmmoniumBasicCount) { return {}; }
DescriptorResult SmilesQuaternaryAmmoniumBasicCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[N\\+\\]\\("); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesImidazoliumCount, RegexCounts, "Count of imidazolium cations ([nH+] in aromatic rings) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesImidazoliumCount) { return {}; }
DescriptorResult SmilesImidazoliumCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[nH\\+\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesPyridiniumCount, RegexCounts, "Count of pyridinium cations ([n+]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPyridiniumCount) { return {}; }
DescriptorResult SmilesPyridiniumCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[n\\+\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesAmmoniumCount, RegexCounts, "Count of ammonium ions ([NH4+]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAmmoniumCount) { return {}; }
DescriptorResult SmilesAmmoniumCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[NH4\\+\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesSulfoniumCount, RegexCounts, "Count of sulfonium ions ([S+]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSulfoniumCount) { return {}; }
DescriptorResult SmilesSulfoniumCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[S\\+\\]"); // Match any [S+]
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesPhosphoniumCount, RegexCounts, "Count of phosphonium ions ([P+]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPhosphoniumCount) { return {}; }
DescriptorResult SmilesPhosphoniumCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[P\\+\\]"); // Match any [P+]
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesZwitterionCount, RegexCounts, "Count of SMILES with both positive and negative charges")
DESCRIPTOR_DEPENDENCIES(SmilesZwitterionCount) { return {}; }
DescriptorResult SmilesZwitterionCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Simpler check for any + and - in the same molecule
    static const RE2 posPattern("\\+");
    static const RE2 negPattern("-");
    bool hasPositive = RE2::PartialMatch(smiles, posPattern);
    bool hasNegative = RE2::PartialMatch(smiles, negPattern);
    return (hasPositive && hasNegative) ? 1 : 0;
}

DECLARE_DESCRIPTOR(SmilesIonizableOxygenCount, RegexCounts, "Count of ionizable oxygens ([O-], [OH]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizableOxygenCount) { return {}; }
DescriptorResult SmilesIonizableOxygenCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[O-\\]|\\[OH\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizableNitrogenCount, RegexCounts, "Count of ionizable nitrogens ([N+], [NH2], [NH3+], [NH], [n+], [nH+], [N-], [NH-], etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizableNitrogenCount) { return {}; }
DescriptorResult SmilesIonizableNitrogenCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "\\[[nN].*\\+\\]|"      // [N+], [n+], [NH3+], [nH+], etc
        "\\[NH2\\]|"            // [NH2]
        "\\[NH\\]|"             // [NH]
        "\\[N-\\]|"             // [N-]
        "\\[NH-\\]|"            // [NH-]
        "\\[N;H2\\]|"           // [N;H2]
        "\\[N;H1\\]|"           // [N;H1]
        "\\[N;H0\\]|"           // [N;H0]
        "\\[N;H3\\]|"           // [N;H3]
        "\\[N;H2;!R\\]|"        // [N;H2;!R]
        "\\[N;H1;!R\\]|"        // [N;H1;!R]
        "\\[N;H2;R1\\]|"        // [N;H2;R1]
        "\\[N;H1;R0\\]|"        // [N;H1;R0]
        "\\[N;H1;R1\\]|"        // [N;H1;R1]
        "\\[N;H1;R2\\]|"        // [N;H1;R2]
        "\\[NH2-\\]|"           // [NH2-]
        "\\[NH2\\+\\]|"         // [NH2+]
        "\\[NH-\\]|"            // [NH-]
        "\\[NH3\\]|"            // [NH3]
        "\\[NH3\\+\\]"          // [NH3+]
    );
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizableSulfurCount, RegexCounts, "Count of ionizable sulfurs ([S-], [SH]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizableSulfurCount) { return {}; }
DescriptorResult SmilesIonizableSulfurCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[S-\\]|\\[SH\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizablePhosphorusCount, RegexCounts, "Count of ionizable phosphorus ([P-], [PH]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizablePhosphorusCount) { return {}; }
DescriptorResult SmilesIonizablePhosphorusCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[P-\\]|\\[PH\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizableCarboxylCount, RegexCounts, "Count of ionizable carboxyls (C(=O)[O-], C(=O)OH) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizableCarboxylCount) { return {}; }
DescriptorResult SmilesIonizableCarboxylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("C\\(=O\\)(\\[O-\\]|O\\[H\\]|O)"); // Matches COOH and COO-
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizablePhenolCount, RegexCounts, "Count of ionizable phenols (cO, c[O-]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizablePhenolCount) { return {}; }
DescriptorResult SmilesIonizablePhenolCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("c(\\[O-\\]|O\\[H\\]|O)"); // Matches cOH and cO-
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizableThiolCount, RegexCounts, "Count of ionizable thiols ([SH], [S-]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizableThiolCount) { return {}; }
DescriptorResult SmilesIonizableThiolCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[SH\\]|\\[S-\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizableAmineCount, RegexCounts, "Count of ionizable amines ([NH2], [NH3+], [N+]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizableAmineCount) { return {}; }
DescriptorResult SmilesIonizableAmineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
     // Matches various protonated or potentially protonatable N
    static const RE2 pattern("\\[[nN].*\\+\\]|\\[NH2\\]|\\[NH\\]");
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizableImidazoleCount, RegexCounts, "Count of ionizable imidazoles ([nH], [nH+]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizableImidazoleCount) { return {}; }
DescriptorResult SmilesIonizableImidazoleCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[nH\\]|\\[nH\\+\\]"); // Already optimized
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizablePyridineCount, RegexCounts, "Count of ionizable pyridines ([n+], [nH+]) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizablePyridineCount) { return {}; }
DescriptorResult SmilesIonizablePyridineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("\\[n\\+\\]|\\[nH\\+\\]"); // Match protonated pyridine N
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizableSulfonamideCount, RegexCounts, "Count of ionizable sulfonamides (S(=O)(=O)[N-], S(=O)(=O)NH) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizableSulfonamideCount) { return {}; }
DescriptorResult SmilesIonizableSulfonamideCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("S\\(=O\\)\\(=O\\)(\\[N-\\]|N\\[H\\])"); // Match SO2NH and SO2N-
    return countMatches(smiles, pattern);
}

DECLARE_DESCRIPTOR(SmilesIonizablePhosphonateCount, RegexCounts, "Count of ionizable phosphonates (P(=O)(O)[O-], P(=O)(O)OH) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIonizablePhosphonateCount) { return {}; }
DescriptorResult SmilesIonizablePhosphonateCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("P\\(=O\\)\\(O\\)(\\[O-\\]|O\\[H\\])"); // Match POOH and POO-
    return countMatches(smiles, pattern);
}

// Long aliphatic chains (3+ carbons)
DECLARE_DESCRIPTOR(SmilesAliphaticChainCount, RegexCounts, "Count of CCC sequences in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAliphaticChainCount) { return {}; }
DescriptorResult SmilesAliphaticChainCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("CCC");
    return countMatches(smiles, pattern);
}

// Long aromatic sequences
DECLARE_DESCRIPTOR(SmilesAromaticSequenceCount, RegexCounts, "Count of ccc sequences in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticSequenceCount) { return {}; }
DescriptorResult SmilesAromaticSequenceCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("[cnos][cnos][cnos]");
    return countMatches(smiles, pattern);
}

// Complex ring closures
DECLARE_DESCRIPTOR(SmilesComplexRingClosureCount, RegexCounts, "Count of double-digit ring closures in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesComplexRingClosureCount) { return {}; }
DescriptorResult SmilesComplexRingClosureCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("[0-9][0-9]");
    return countMatches(smiles, pattern);
}

// Tertiary carbon count
DECLARE_DESCRIPTOR(SmilesTertiaryCarbonCount, RegexCounts, "Count of tertiary carbons C(C)(C)(C) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesTertiaryCarbonCount) { return {}; }
DescriptorResult SmilesTertiaryCarbonCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)|" // C(C)(C)(C)
        "\\[C\\]\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)|" // [C](C)(C)(C)
        "C\\(C\\)\\(C\\)C|" // C(C)(C)C
        "C\\(C\\)C\\(C\\)|" // C(C)C(C)
        "C\\(C\\)C\\(C\\)|" // C(C)C(C)
        "C\\(C\\)\\(C\\)\\(C\\)" // C(C)(C)(C)
    );
    return countMatches(smiles, pattern);
}

// CNO sequences 
DECLARE_DESCRIPTOR(SmilesHeteroatomSequenceCount, RegexCounts, "Count of CNO neighbor sequences in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesHeteroatomSequenceCount) { return {}; }
DescriptorResult SmilesHeteroatomSequenceCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "[CNO][CNO]|"         // CNO, NOC, OCN, etc.
        "[CNSO][CNSO]|"       // also include S for sulfur
        "[CNP][CNP]|"         // also include P for phosphorus
        "[CNOFPS][CNOFPS]"    // include F, P, S for broader heteroatom coverage
    );
    return countMatches(smiles, pattern);
}

// Conjugated systems
DECLARE_DESCRIPTOR(SmilesConjugatedSystemCount, RegexCounts, "Count of C=CC=C patterns in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesConjugatedSystemCount) { return {}; }
DescriptorResult SmilesConjugatedSystemCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C=CC=C|"             // C=CC=C
        "C=CC=N|"             // C=CC=N
        "C=CN=C|"             // C=CN=C
        "N=CC=C|"             // N=CC=C
        "C=CC=O|"             // C=CC=O
        "C=CC=S|"             // C=CC=S
        "C=CC#N|"             // C=CC#N
        "C=CC=P|"             // C=CC=P
        "C=CC#C|"             // C=CC#C
        "C=CC=C|"             // C=CC=C (repeat for completeness)
        "c=cc=c"              // aromatic conjugation (lowercase for aromatic)
    );
    return countMatches(smiles, pattern);
}

// Methoxy groups
DECLARE_DESCRIPTOR(SmilesMethoxyCount, RegexCounts, "Count of OCH3 groups in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesMethoxyCount) { return {}; }
DescriptorResult SmilesMethoxyCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "OC([^a-zA-Z0-9]|$)|"   // OC
        "OCH3|"                 // OCH3
        "O\\[CH3\\]|"           // O[CH3]
        "\\[O\\]C|"             // [O]C
        "\\[O\\]\\[CH3\\]|"     // [O][CH3]
        "COC([^a-zA-Z0-9]|$)|"  // COC (terminal)
        "CO[CH3]|"              // CO[CH3]
        "CO\\(C\\)|"            // CO(C)
        "COC\\(C\\)"            // COC(C)
    );
    return countMatches(smiles, pattern);
}

// Ethers with branching
DECLARE_DESCRIPTOR(SmilesBranchedEtherCount, RegexCounts, "Count of O(C)C ether patterns in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBranchedEtherCount) { return {}; }
DescriptorResult SmilesBranchedEtherCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "O\\(C\\)C|"            // O(C)C
        "O(C)C|"                // O(C)C
        "O(C)[CH3]|"            // O(C)[CH3]
        "O(C)C(C)|"             // O(C)C(C)
        "COC(C)|"               // COC(C)
        "CO(C)C|"               // CO(C)C
        "CO(C)[CH3]|"           // CO(C)[CH3]
        "CO(C)C(C)"             // CO(C)C(C)
    );
    return countMatches(smiles, pattern);
}

// Amines with branching
DECLARE_DESCRIPTOR(SmilesBranchedAmineCount, RegexCounts, "Count of N(C)C amine patterns in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBranchedAmineCount) { return {}; }
DescriptorResult SmilesBranchedAmineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "N\\(C\\)C|"            // N(C)C
        "N(C)C|"                // N(C)C
        "N(C)[CH3]|"            // N(C)[CH3]
        "N(C)C(C)|"             // N(C)C(C)
        "NC(C)C|"               // NC(C)C
        "NC(C)[CH3]|"           // NC(C)[CH3]
        "NC(C)C(C)"             // NC(C)C(C)
    );
    return countMatches(smiles, pattern);
}

// Long carbon chains
DECLARE_DESCRIPTOR(SmilesLongCarbonChainCount, RegexCounts, "Count of CCCC sequences in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesLongCarbonChainCount) { return {}; }
DescriptorResult SmilesLongCarbonChainCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "CCCC|"                 // CCCC
        "C\\(C\\)CC|"           // C(C)CC
        "CC\\(C\\)C|"           // CC(C)C
        "C\\(C\\)C\\(C\\)|"     // C(C)C(C)
        "C\\(C\\)CC(C)|"        // C(C)CC(C)
        "CCCCC"                 // CCCCC (longer chain)
    );
    return countMatches(smiles, pattern);
}

// --- Registration Functions ---
// (No changes needed in registration functions themselves)

void register_SmilesCarbonCountDescriptor() {
    auto descriptor = std::make_shared<SmilesCarbonCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNitrogenCountDescriptor() {
    auto descriptor = std::make_shared<SmilesNitrogenCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesOxygenCountDescriptor() {
    auto descriptor = std::make_shared<SmilesOxygenCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesSulfurCountDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfurCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesPhosphorusCountDescriptor() {
    auto descriptor = std::make_shared<SmilesPhosphorusCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHalogenCountDescriptor() {
    auto descriptor = std::make_shared<SmilesHalogenCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesBoronCountDescriptor() {
    auto descriptor = std::make_shared<SmilesBoronCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesSiliconCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesSiliconCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesFluorineCountDescriptor() {
    auto descriptor = std::make_shared<SmilesFluorineCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesChlorineCountDescriptor() {
    auto descriptor = std::make_shared<SmilesChlorineCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesBromineCountDescriptor() {
    auto descriptor = std::make_shared<SmilesBromineCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIodineCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIodineCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesMethylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesMethylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesEthylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesEthylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesAminoCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAminoCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesEsterCountDescriptor() {
    auto descriptor = std::make_shared<SmilesEsterCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesEtherCountDescriptor() {
    auto descriptor = std::make_shared<SmilesEtherCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesThiolCountDescriptor() {
    auto descriptor = std::make_shared<SmilesThiolCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesThioetherCountDescriptor() {
    auto descriptor = std::make_shared<SmilesThioetherCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAmideCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesAmideCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesNitrileCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesNitrileCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}


void register_SmilesSulfonylCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesSulfonylCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}


void register_SmilesAlkyneCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesAlkyneCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesAlkeneCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesAlkeneCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}


void register_SmilesFiveRingCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesFiveRingCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesSixRingCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesSixRingCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesSevenRingCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesSevenRingCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesEightRingCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesEightRingCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesIsotopeCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesIsotopeCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesTertiaryAmineCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesTertiaryAmineCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesQuaternaryAmmoniumCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesQuaternaryAmmoniumCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesPrimaryAmineCountDescriptor() {
    auto descriptor = std::make_shared<SmilesPrimaryAmineCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesSecondaryAmineCountDescriptor() {
    auto descriptor = std::make_shared<SmilesSecondaryAmineCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesSulfonamideCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesSulfonamideCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesPhosphonateCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesPhosphonateCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesAlkylHalideCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAlkylHalideCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAliphaticNitroCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesAliphaticNitroCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesAzideCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesAzideCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesDiazoCountDescriptor() {
    auto descriptor = std::make_shared<SmilesDiazoCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesThiocyanateCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesThiocyanateCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesIsothiocyanateLinearCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesIsothiocyanateLinearCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesGuanidineCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesGuanidineCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesUreaCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesUreaCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesAmidineCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesAmidineCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesImidazoleCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesImidazoleCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesPyridineCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesPyridineCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesPyrroleCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesPyrroleCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesFuranCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesFuranCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}

void register_SmilesThiopheneCountDescriptor() {
    // auto descriptor = std::make_shared<SmilesThiopheneCountDescriptor>(); // Zero Variance
    // auto& registry = DescriptorRegistry::getInstance();
    // registry.registerDescriptor(descriptor);
}


void register_SmilesAcidicPhosphorusCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAcidicPhosphorusCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesPhenolCountDescriptor() {
    auto descriptor = std::make_shared<SmilesPhenolCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesCarboxylateCountDescriptor() {
    auto descriptor = std::make_shared<SmilesCarboxylateCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesSulfonateCountDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfonateCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesPhosphateAnionCountDescriptor() {
    auto descriptor = std::make_shared<SmilesPhosphateAnionCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesImidazoleBasicCountDescriptor() {
    auto descriptor = std::make_shared<SmilesImidazoleBasicCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesPyridineBasicCountDescriptor() {
    auto descriptor = std::make_shared<SmilesPyridineBasicCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesPrimaryAmineBasicCountDescriptor() {
    auto descriptor = std::make_shared<SmilesPrimaryAmineBasicCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesSecondaryAmineBasicCountDescriptor() {
    auto descriptor = std::make_shared<SmilesSecondaryAmineBasicCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesTertiaryAmineBasicCountDescriptor() {
    auto descriptor = std::make_shared<SmilesTertiaryAmineBasicCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesGuanidineBasicCountDescriptor() {
    auto descriptor = std::make_shared<SmilesGuanidineBasicCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesAmidineBasicCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAmidineBasicCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesQuaternaryAmmoniumBasicCountDescriptor() {
    auto descriptor = std::make_shared<SmilesQuaternaryAmmoniumBasicCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesImidazoliumCountDescriptor() {
    auto descriptor = std::make_shared<SmilesImidazoliumCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesPyridiniumCountDescriptor() {
    auto descriptor = std::make_shared<SmilesPyridiniumCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesAmmoniumCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAmmoniumCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesSulfoniumCountDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfoniumCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesPhosphoniumCountDescriptor() {
    auto descriptor = std::make_shared<SmilesPhosphoniumCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesZwitterionCountDescriptor() {
    auto descriptor = std::make_shared<SmilesZwitterionCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesIonizableOxygenCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizableOxygenCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizableNitrogenCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizableNitrogenCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizableSulfurCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizableSulfurCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizablePhosphorusCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizablePhosphorusCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizableCarboxylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizableCarboxylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizablePhenolCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizablePhenolCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizableThiolCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizableThiolCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizableAmineCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizableAmineCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizableImidazoleCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizableImidazoleCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizablePyridineCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizablePyridineCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizableSulfonamideCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizableSulfonamideCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIonizablePhosphonateCountDescriptor() {
    auto descriptor = std::make_shared<SmilesIonizablePhosphonateCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticCarbonCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticCarbonCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticNitrogenCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticNitrogenCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticOxygenCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticOxygenCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticSulfurCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticSulfurCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticHalogenCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticHalogenCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesSingleBondCountDescriptor() {
    auto descriptor = std::make_shared<SmilesSingleBondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesDoubleBondCountDescriptor() {
    auto descriptor = std::make_shared<SmilesDoubleBondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesTripleBondCountDescriptor() {
    auto descriptor = std::make_shared<SmilesTripleBondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticBondCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticBondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesRingClosureCountDescriptor() {
    auto descriptor = std::make_shared<SmilesRingClosureCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesBranchOpenCountDescriptor() {
    auto descriptor = std::make_shared<SmilesBranchOpenCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesBranchCloseCountDescriptor() {
    auto descriptor = std::make_shared<SmilesBranchCloseCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesBracketAtomCountDescriptor() {
    auto descriptor = std::make_shared<SmilesBracketAtomCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesChiralCenterCountDescriptor() {
    auto descriptor = std::make_shared<SmilesChiralCenterCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesPositiveChargeCountDescriptor() {
    auto descriptor = std::make_shared<SmilesPositiveChargeCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNegativeChargeCountDescriptor() {
    auto descriptor = std::make_shared<SmilesNegativeChargeCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesStereoForwardCountDescriptor() {
    auto descriptor = std::make_shared<SmilesStereoForwardCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesStereoBackwardCountDescriptor() {
    auto descriptor = std::make_shared<SmilesStereoBackwardCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesDotCountDescriptor() {
    auto descriptor = std::make_shared<SmilesDotCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesAliphaticChainCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAliphaticChainCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesAromaticSequenceCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticSequenceCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesComplexRingClosureCountDescriptor() {
    auto descriptor = std::make_shared<SmilesComplexRingClosureCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesTertiaryCarbonCountDescriptor() {
    auto descriptor = std::make_shared<SmilesTertiaryCarbonCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesHeteroatomSequenceCountDescriptor() {
    auto descriptor = std::make_shared<SmilesHeteroatomSequenceCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesConjugatedSystemCountDescriptor() {
    auto descriptor = std::make_shared<SmilesConjugatedSystemCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesMethoxyCountDescriptor() {
    auto descriptor = std::make_shared<SmilesMethoxyCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesBranchedEtherCountDescriptor() {
    auto descriptor = std::make_shared<SmilesBranchedEtherCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesBranchedAmineCountDescriptor() {
    auto descriptor = std::make_shared<SmilesBranchedAmineCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesLongCarbonChainCountDescriptor() {
    auto descriptor = std::make_shared<SmilesLongCarbonChainCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

// Fluorinated fragments
DECLARE_DESCRIPTOR(SmilesFluorinatedFragmentCount, RegexCounts, "Count of CF patterns in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesFluorinatedFragmentCount) { return {}; }
DescriptorResult SmilesFluorinatedFragmentCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("CF");
    return countMatches(smiles, pattern);
}

// Mixed heteroatom rings
DECLARE_DESCRIPTOR(SmilesMixedHeteroatomRingCount, RegexCounts, "Count of rings with multiple heteroatoms in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesMixedHeteroatomRingCount) { return {}; }
DescriptorResult SmilesMixedHeteroatomRingCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("c1[nos]c[nos]c1");
    return countMatches(smiles, pattern);
}

// Carbonyl adjacent to heteroatom
DECLARE_DESCRIPTOR(SmilesCarbonylHeteroatomCount, RegexCounts, "Count of C(=O)[NOS] patterns in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesCarbonylHeteroatomCount) { return {}; }
DescriptorResult SmilesCarbonylHeteroatomCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("C\\(=O\\)[NOS]");
    return countMatches(smiles, pattern);
}

// Aromatic carbon adjacent to heteroatom
DECLARE_DESCRIPTOR(SmilesAromaticCarbonHeteroatomCount, RegexCounts, "Count of c[NOS] patterns in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticCarbonHeteroatomCount) { return {}; }
DescriptorResult SmilesAromaticCarbonHeteroatomCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("c[NOSPnospf]|[NOSPnospf]C");
    return countMatches(smiles, pattern);
}

// Trifluoromethyl groups
DECLARE_DESCRIPTOR(SmilesTrifluoromethylCount, RegexCounts, "Count of CF3 groups in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesTrifluoromethylCount) { return {}; }
DescriptorResult SmilesTrifluoromethylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("CF3|C\\(F\\)\\(F\\)F");
    return countMatches(smiles, pattern);
}

// Terminal Alkynes
DECLARE_DESCRIPTOR(SmilesTerminalAlkyneCount, RegexCounts, "Count of terminal C#C patterns in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesTerminalAlkyneCount) { return {}; }
DescriptorResult SmilesTerminalAlkyneCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("C#C([^a-zA-Z]|$)");
    return countMatches(smiles, pattern);
}

// Dimethyl patterns
DECLARE_DESCRIPTOR(SmilesDimethylCount, RegexCounts, "Count of CC(C)C patterns in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesDimethylCount) { return {}; }
DescriptorResult SmilesDimethylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("CC\\(C\\)C");
    return countMatches(smiles, pattern);
}

// Aromatic nitrogen adjacent to carbon
DECLARE_DESCRIPTOR(SmilesAromaticNitrogenCarbonCount, RegexCounts, "Count of nc patterns in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticNitrogenCarbonCount) { return {}; }
DescriptorResult SmilesAromaticNitrogenCarbonCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern("nc");
    return countMatches(smiles, pattern);
}

void register_SmilesFluorinatedFragmentCountDescriptor() {
    auto descriptor = std::make_shared<SmilesFluorinatedFragmentCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesMixedHeteroatomRingCountDescriptor() {
    auto descriptor = std::make_shared<SmilesMixedHeteroatomRingCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesCarbonylHeteroatomCountDescriptor() {
    auto descriptor = std::make_shared<SmilesCarbonylHeteroatomCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}



void register_SmilesTrifluoromethylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesTrifluoromethylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesTerminalAlkyneCountDescriptor() {
    auto descriptor = std::make_shared<SmilesTerminalAlkyneCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesDimethylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesDimethylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesAromaticNitrogenCarbonCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticNitrogenCarbonCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

// Cyclopropyl groups: C1CC1, [CH2]1[CH2][CH2]1, [C;H2]1[C;H2][C;H2]1, etc.
DECLARE_DESCRIPTOR(SmilesCyclopropylCount, RegexCounts, "Count of cyclopropyl groups (C1CC1, [CH2]1[CH2][CH2]1, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesCyclopropylCount) { return {}; }
DescriptorResult SmilesCyclopropylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C1CC1|"                       // C1CC1
        "\\[CH2\\]1\\[CH2\\]\\[CH2\\]1|" // [CH2]1[CH2][CH2]1
        "\\[C;H2\\]1\\[C;H2\\]\\[C;H2\\]1" // [C;H2]1[C;H2][C;H2]1
    );
    return countMatches(smiles, pattern);
}

// Cyclobutyl groups: C1CCC1, [CH2]1[CH2][CH2][CH2]1, [C;H2]1[C;H2][C;H2][C;H2]1, etc.
DECLARE_DESCRIPTOR(SmilesCyclobutylCount, RegexCounts, "Count of cyclobutyl groups (C1CCC1, [CH2]1[CH2][CH2][CH2]1, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesCyclobutylCount) { return {}; }
DescriptorResult SmilesCyclobutylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C1CCC1|"                       // C1CCC1
        "\\[CH2\\]1\\[CH2\\]\\[CH2\\]\\[CH2\\]1|" // [CH2]1[CH2][CH2][CH2]1
        "\\[C;H2\\]1\\[C;H2\\]\\[C;H2\\]\\[C;H2\\]1" // [C;H2]1[C;H2][C;H2][C;H2]1
    );
    return countMatches(smiles, pattern);
}

// Cyclopentyl groups: C1CCCC1, [CH2]1[CH2][CH2][CH2][CH2]1, etc.
DECLARE_DESCRIPTOR(SmilesCyclopentylCount, RegexCounts, "Count of cyclopentyl groups (C1CCCC1, [CH2]1[CH2][CH2][CH2][CH2]1, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesCyclopentylCount) { return {}; }
DescriptorResult SmilesCyclopentylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C1CCCC1|"                       // C1CCCC1
        "\\[CH2\\]1\\[CH2\\]\\[CH2\\]\\[CH2\\]\\[CH2\\]1" // [CH2]1[CH2][CH2][CH2][CH2]1
    );
    return countMatches(smiles, pattern);
}

// Cyclohexyl groups: C1CCCCC1, [CH2]1[CH2][CH2][CH2][CH2][CH2]1, etc.
DECLARE_DESCRIPTOR(SmilesCyclohexylCount, RegexCounts, "Count of cyclohexyl groups (C1CCCCC1, [CH2]1[CH2][CH2][CH2][CH2][CH2]1, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesCyclohexylCount) { return {}; }
DescriptorResult SmilesCyclohexylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C1CCCCC1|"                       // C1CCCCC1
        "\\[CH2\\]1\\[CH2\\]\\[CH2\\]\\[CH2\\]\\[CH2\\]\\[CH2\\]1" // [CH2]1[CH2][CH2][CH2][CH2][CH2]1
    );
    return countMatches(smiles, pattern);
}

// Carbamate groups: NC(=O)O, [NH2]C(=O)O, [NH]C(=O)O, etc.
DECLARE_DESCRIPTOR(SmilesCarbamatePatternsCount, RegexCounts, "Count of NC(=O)O carbamate patterns (NC(=O)O, [NH2]C(=O)O, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesCarbamatePatternsCount) { return {}; }
DescriptorResult SmilesCarbamatePatternsCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "NC\\(=O\\)O|"                  // NC(=O)O
        "\\[NH2\\]C\\(=O\\)O|"         // [NH2]C(=O)O
        "\\[NH\\]C\\(=O\\)O"           // [NH]C(=O)O
    );
    return countMatches(smiles, pattern);
}

// N-alkyl patterns: NC, [NH2]C, [NH]C, etc.
DECLARE_DESCRIPTOR(SmilesNAlkylCount, RegexCounts, "Count of NC patterns (NC, [NH2]C, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesNAlkylCount) { return {}; }
DescriptorResult SmilesNAlkylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "NC|"                           // NC
        "\\[NH2\\]C|"                   // [NH2]C
        "\\[NH\\]C"                     // [NH]C
    );
    return countMatches(smiles, pattern);
}

// Quaternary carbon atoms: C(C)(C)(C)C, [C](C)(C)(C)C, etc.
DECLARE_DESCRIPTOR(SmilesQuaternaryCarbonCount, RegexCounts, "Count of quaternary carbon patterns (C(C)(C)(C)C, [C](C)(C)(C)C, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesQuaternaryCarbonCount) { return {}; }
DescriptorResult SmilesQuaternaryCarbonCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)|" // C(C)(C)(C)C
        "\\[C\\]\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)" // [C](C)(C)(C)C
    );
    return countMatches(smiles, pattern);
}

// Carbon-phosphorus bonds: CP, PC, [C]P, C[P], etc.
DECLARE_DESCRIPTOR(SmilesCarbonPhosphorusBondCount, RegexCounts, "Count of CP patterns (CP, PC, [C]P, C[P], etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesCarbonPhosphorusBondCount) { return {}; }
DescriptorResult SmilesCarbonPhosphorusBondCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "CP|PC|"                       // CP, PC
        "\\[C\\]P|C\\[P\\]"           // [C]P, C[P]
    );
    return countMatches(smiles, pattern);
}

// Double bonds adjacent to rings: C=C1, [C]=C1, C=[C]1, etc.
DECLARE_DESCRIPTOR(SmilesDoubleBondRingJunctionCount, RegexCounts, "Count of C=C1 patterns (C=C1, [C]=C1, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesDoubleBondRingJunctionCount) { return {}; }
DescriptorResult SmilesDoubleBondRingJunctionCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C=C[0-9]|"                    // C=C1
        "\\[C\\]=C[0-9]|C=\\[C\\][0-9]"// [C]=C1, C=[C]1
    );
    return countMatches(smiles, pattern);
}

// Oxygen-containing aromatic rings: c1cccoc1, c1ccco1, c1cocc1, etc.
DECLARE_DESCRIPTOR(SmilesOxygenAromaticRingCount, RegexCounts, "Count of c1cccoc1 patterns (c1cccoc1, c1ccco1, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesOxygenAromaticRingCount) { return {}; }
DescriptorResult SmilesOxygenAromaticRingCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "c1cccoc1|"                    // c1cccoc1
        "c1ccco1|"                     // c1ccco1
        "c1cocc1"                      // c1cocc1
    );
    return countMatches(smiles, pattern);
}

// Branched amide groups: C(=O)N(C), [C](=O)N(C), etc.
DECLARE_DESCRIPTOR(SmilesBranchedAmideCount, RegexCounts, "Count of C(=O)N(C) branched amide patterns (C(=O)N(C), [C](=O)N(C), etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBranchedAmideCount) { return {}; }
DescriptorResult SmilesBranchedAmideCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "C\\(=O\\)N\\(C\\)|"           // C(=O)N(C)
        "\\[C\\]\\(=O\\)N\\(C\\)"      // [C](=O)N(C)
    );
    return countMatches(smiles, pattern);
}

// Fluoromethyl group: CF, [CH2]F, [CFH2], etc.
DECLARE_DESCRIPTOR(SmilesFluoromethylCount, RegexCounts, "Count of fluoromethyl (CF, [CH2]F, etc) groups in SMILES")
DESCRIPTOR_DEPENDENCIES(SmilesFluoromethylCount) { return {}; }
DescriptorResult SmilesFluoromethylCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "CF([^a-zA-Z0-9]|$)|"          // CF
        "\\[CH2\\]F|"                  // [CH2]F
        "\\[CFH2\\]"                   // [CFH2]
    );
    return countMatches(smiles, pattern);
}


void register_SmilesCyclopropylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesCyclopropylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesCyclobutylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesCyclobutylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesCyclopentylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesCyclopentylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesCyclohexylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesCyclohexylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesCarbamatePatternsCountDescriptor() {
    auto descriptor = std::make_shared<SmilesCarbamatePatternsCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNAlkylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesNAlkylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesQuaternaryCarbonCountDescriptor() {
    auto descriptor = std::make_shared<SmilesQuaternaryCarbonCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesCarbonPhosphorusBondCountDescriptor() {
    auto descriptor = std::make_shared<SmilesCarbonPhosphorusBondCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesDoubleBondRingJunctionCountDescriptor() {
    auto descriptor = std::make_shared<SmilesDoubleBondRingJunctionCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesOxygenAromaticRingCountDescriptor() {
    auto descriptor = std::make_shared<SmilesOxygenAromaticRingCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesBranchedAmideCountDescriptor() {
    auto descriptor = std::make_shared<SmilesBranchedAmideCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesFluoromethylCountDescriptor() {
    auto descriptor = std::make_shared<SmilesFluoromethylCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

// Nitro group (general, including rare notations)
DECLARE_DESCRIPTOR(SmilesGeneralNitroCount, RegexCounts, "Count of nitro groups (NO2, [N+](=O)[O-], N(=O)=O, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesGeneralNitroCount) { return {}; }
DescriptorResult SmilesGeneralNitroCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "NO2|"                        // NO2
        "N\\(=O\\)=O|"                // N(=O)=O
        "\\[N\\+\\]\\(=O\\)\\[O-\\]|"// [N+](=O)[O-]
        "N\\(O\\)O"                   // N(O)O
    );
    return countMatches(smiles, pattern);
}

// Ammonium (general, including [NH4+], [N+], [N+](C)(C)C, etc)
DECLARE_DESCRIPTOR(SmilesGeneralAmmoniumCount, RegexCounts, "Count of ammonium ions ([NH4+], [N+], [N+](C)(C)C, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesGeneralAmmoniumCount) { return {}; }
DescriptorResult SmilesGeneralAmmoniumCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "\\[NH4\\+\\]|"               // [NH4+]
        "\\[N\\+\\]|"                 // [N+]
        "\\[N\\+\\]\\([^\\)]*\\)"     // [N+](...)
    );
    return countMatches(smiles, pattern);
}

// Guanidinium (protonated guanidine)
DECLARE_DESCRIPTOR(SmilesGuanidiniumCount, RegexCounts, "Count of guanidinium ions ([NH2]C(=NH2+)NH2, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesGuanidiniumCount) { return {}; }
DescriptorResult SmilesGuanidiniumCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "\\[NH2\\]C\\(=N\\[H2\\+\\]\\)\\[NH2\\]" // [NH2]C(=NH2+)NH2
    );
    return countMatches(smiles, pattern);
}

// Sulfonic acid group (SO3H, [SO3H], S(=O)(=O)O, etc)
DECLARE_DESCRIPTOR(SmilesSulfonicAcidCount, RegexCounts, "Count of sulfonic acid groups (SO3H, [SO3H], S(=O)(=O)O, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSulfonicAcidCount) { return {}; }
DescriptorResult SmilesSulfonicAcidCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "SO3H|"                       // SO3H
        "\\[SO3H\\]|"                 // [SO3H]
        "S\\(=O\\)\\(=O\\)O"          // S(=O)(=O)O
    );
    return countMatches(smiles, pattern);
}

// Phosphonic acid group (PO3H2, [PO3H2], P(=O)(O)O, etc)
DECLARE_DESCRIPTOR(SmilesPhosphonicAcidCount, RegexCounts, "Count of phosphonic acid groups (PO3H2, [PO3H2], P(=O)(O)O, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPhosphonicAcidCount) { return {}; }
DescriptorResult SmilesPhosphonicAcidCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "PO3H2|"                      // PO3H2
        "\\[PO3H2\\]|"                // [PO3H2]
        "P\\(=O\\)\\(O\\)O"           // P(=O)(O)O
    );
    return countMatches(smiles, pattern);
}

// Isocyanide group (NC, [NC], [N+]#C-, etc)
DECLARE_DESCRIPTOR(SmilesIsocyanideCount, RegexCounts, "Count of isocyanide groups (NC, [NC], [N+]#C-, etc) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIsocyanideCount) { return {}; }
DescriptorResult SmilesIsocyanideCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 pattern(
        "N#C|"                        // N#C
        "\\[NC\\]|"                   // [NC]
        "\\[N\\+\\]#C-"               // [N+]#C-
    );
    return countMatches(smiles, pattern);
}

void register_SmilesSulfonicAcidCountDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfonicAcidCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesPhosphonicAcidCountDescriptor() {
    auto descriptor = std::make_shared<SmilesPhosphonicAcidCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesGeneralAmmoniumCountDescriptor() {}

void register_SmilesGeneralNitroCountDescriptor() {}

void register_SmilesGuanidiniumCountDescriptor() {}

void register_SmilesIsocyanideCountDescriptor() {}


void register_SmilesAcidicSulfurCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAcidicSulfurCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAromaticCarbonHeteroatomCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticCarbonHeteroatomCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesSulfoxideCountDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfoxideCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

} // namespace desfact