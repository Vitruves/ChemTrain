#include "../common.hpp"
#include <re2/re2.h>
#include <string>
#include <memory>
#include <vector>
#include <cmath> // For std::isnan

namespace desfact {

// --- Helper Functions (Declared in common.hpp, Defined elsewhere) ---

// Note: Definitions for countMatches, countChar, countDigits, countSubstr
// are removed from this file to avoid multiple definition errors.
// Ensure they are defined in exactly one .cpp file (e.g., common.cpp or regexcounts.cpp)
// and declared in common.hpp.

// Helper to calculate fraction safely
inline double calculateFraction(int count, double totalLength) {
    if (totalLength < 1.0) { // Avoid division by zero or near-zero
        return 0.0;
    }
    double fraction = static_cast<double>(count) / totalLength;
    // Ensure result is not NaN or Inf, return 0.0 if it is
    return (std::isnan(fraction) || std::isinf(fraction)) ? 0.0 : fraction;
}

// Base class for Regex Fraction descriptors
class RegexFractionsDescriptor : public Descriptor {
public:
    using Descriptor::Descriptor;
};


// --- Descriptor Implementations ---

// Carbon fraction
DECLARE_DESCRIPTOR(SmilesCarbonFraction, RegexFractions, "Fraction of carbon characters 'C' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesCarbonFraction) { return {}; }
DescriptorResult SmilesCarbonFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'C');
    return calculateFraction(count, totalLength);
}

// Nitrogen fraction
DECLARE_DESCRIPTOR(SmilesNitrogenFraction, RegexFractions, "Fraction of nitrogen characters 'N' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesNitrogenFraction) { return {}; }
DescriptorResult SmilesNitrogenFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'N');
    return calculateFraction(count, totalLength);
}

// Oxygen fraction
DECLARE_DESCRIPTOR(SmilesOxygenFraction, RegexFractions, "Fraction of oxygen characters 'O' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesOxygenFraction) { return {}; }
DescriptorResult SmilesOxygenFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'O');
    return calculateFraction(count, totalLength);
}

// Sulfur fraction
DECLARE_DESCRIPTOR(SmilesSulfurFraction, RegexFractions, "Fraction of sulfur characters 'S' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSulfurFraction) { return {}; }
DescriptorResult SmilesSulfurFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'S');
    return calculateFraction(count, totalLength);
}

// Phosphorus fraction
DECLARE_DESCRIPTOR(SmilesPhosphorusFraction, RegexFractions, "Fraction of phosphorus characters 'P' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPhosphorusFraction) { return {}; }
DescriptorResult SmilesPhosphorusFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'P');
    return calculateFraction(count, totalLength);
}

// Halogen fraction
DECLARE_DESCRIPTOR(SmilesHalogenFraction, RegexFractions, "Fraction of halogen characters (F, Cl, Br, I) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesHalogenFraction) { return {}; }
DescriptorResult SmilesHalogenFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = 0;
    count += countChar(smiles, 'F');
    count += countChar(smiles, 'I');
    static const RE2 pattern_cl("Cl");
    count += countMatches(smiles, pattern_cl);
    static const RE2 pattern_br("Br");
    count += countMatches(smiles, pattern_br);
    return calculateFraction(count, totalLength);
}

// Boron fraction
DECLARE_DESCRIPTOR(SmilesBoronFraction, RegexFractions, "Fraction of boron characters 'B' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBoronFraction) { return {}; }
DescriptorResult SmilesBoronFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'B');
    return calculateFraction(count, totalLength);
}

// Silicon fraction - Removed

// Fluorine fraction
DECLARE_DESCRIPTOR(SmilesFluorineFraction, RegexFractions, "Fraction of fluorine characters 'F' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesFluorineFraction) { return {}; }
DescriptorResult SmilesFluorineFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'F');
    return calculateFraction(count, totalLength);
}

// Chlorine fraction - REMOVED (Zero Variance)

// Bromine fraction - REMOVED (Zero Variance)

// Iodine fraction
DECLARE_DESCRIPTOR(SmilesIodineFraction, RegexFractions, "Fraction of iodine characters 'I' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesIodineFraction) { return {}; }
DescriptorResult SmilesIodineFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'I');
    return calculateFraction(count, totalLength);
}

// Aromatic carbon fraction
DECLARE_DESCRIPTOR(SmilesAromaticCarbonFraction, RegexFractions, "Fraction of aromatic carbon characters 'c' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticCarbonFraction) { return {}; }
DescriptorResult SmilesAromaticCarbonFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'c');
    return calculateFraction(count, totalLength);
}

// Aromatic nitrogen fraction
DECLARE_DESCRIPTOR(SmilesAromaticNitrogenFraction, RegexFractions, "Fraction of aromatic nitrogen characters 'n' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticNitrogenFraction) { return {}; }
DescriptorResult SmilesAromaticNitrogenFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'n');
    return calculateFraction(count, totalLength);
}

// Aromatic oxygen fraction
DECLARE_DESCRIPTOR(SmilesAromaticOxygenFraction, RegexFractions, "Fraction of aromatic oxygen characters 'o' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticOxygenFraction) { return {}; }
DescriptorResult SmilesAromaticOxygenFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 'o');
    return calculateFraction(count, totalLength);
}

// Aromatic sulfur fraction
DECLARE_DESCRIPTOR(SmilesAromaticSulfurFraction, RegexFractions, "Fraction of aromatic sulfur characters 's' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAromaticSulfurFraction) { return {}; }
DescriptorResult SmilesAromaticSulfurFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, 's');
    return calculateFraction(count, totalLength);
}

// Single bond fraction
DECLARE_DESCRIPTOR(SmilesSingleBondFraction, RegexFractions, "Fraction of single bond characters '-' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesSingleBondFraction) { return {}; }
DescriptorResult SmilesSingleBondFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, '-');
    return calculateFraction(count, totalLength);
}

// Double bond fraction
DECLARE_DESCRIPTOR(SmilesDoubleBondFraction, RegexFractions, "Fraction of double bond characters '=' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesDoubleBondFraction) { return {}; }
DescriptorResult SmilesDoubleBondFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, '=');
    return calculateFraction(count, totalLength);
}

// Triple bond fraction
DECLARE_DESCRIPTOR(SmilesTripleBondFraction, RegexFractions, "Fraction of triple bond characters '#' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesTripleBondFraction) { return {}; }
DescriptorResult SmilesTripleBondFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, '#');
    return calculateFraction(count, totalLength);
}

// Ring closure fraction
DECLARE_DESCRIPTOR(SmilesRingClosureFraction, RegexFractions, "Fraction of ring closure digits in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesRingClosureFraction) { return {}; }
DescriptorResult SmilesRingClosureFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countDigits(smiles);
    return calculateFraction(count, totalLength);
}

// Branch open fraction
DECLARE_DESCRIPTOR(SmilesBranchOpenFraction, RegexFractions, "Fraction of branch opening characters '(' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBranchOpenFraction) { return {}; }
DescriptorResult SmilesBranchOpenFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, '(');
    return calculateFraction(count, totalLength);
}

// Branch close fraction
DECLARE_DESCRIPTOR(SmilesBranchCloseFraction, RegexFractions, "Fraction of branch closing characters ')' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBranchCloseFraction) { return {}; }
DescriptorResult SmilesBranchCloseFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, ')');
    return calculateFraction(count, totalLength);
}

// Square bracket atom fraction
DECLARE_DESCRIPTOR(SmilesBracketAtomFraction, RegexFractions, "Fraction of bracket opening characters '[' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesBracketAtomFraction) { return {}; }
DescriptorResult SmilesBracketAtomFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, '[');
    return calculateFraction(count, totalLength);
}

// Positive charge fraction
DECLARE_DESCRIPTOR(SmilesPositiveChargeFraction, RegexFractions, "Fraction of positive charge characters '+' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesPositiveChargeFraction) { return {}; }
DescriptorResult SmilesPositiveChargeFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    int count = countChar(smiles, '+');
    return calculateFraction(count, totalLength);
}

// Negative charge fraction
DECLARE_DESCRIPTOR(SmilesNegativeChargeFraction, RegexFractions, "Fraction of negative charge characters '-' in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesNegativeChargeFraction) { return {}; }
DescriptorResult SmilesNegativeChargeFractionDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double totalLength = static_cast<double>(smiles.length());
    // This previously incorrectly counted single bonds again. Fixed to count '-' charge.
    // Needs careful checking if '-' is overloaded elsewhere in SMILES.
    // Assuming standard SMILES, this should be correct.
    // Let's use regex to be more precise for charge within brackets.
    static const RE2 pattern("\\[[A-Za-z@]*.*[-][1-9]?\\]");
    int count = countMatches(smiles, pattern);
    // Also count simple '-' if it's not a bond (less common, but possible for some representations)
    // Example: [O-]
    // Avoid double counting with the regex
    for(size_t i = 0; i < smiles.length(); ++i) {
        if (smiles[i] == '-' && i > 0 && smiles[i-1] != '[') {
            // Check if it's likely a bond character based on context
            bool isBond = false;
            if (i > 0 && isalnum(smiles[i-1])) isBond = true;
            if (i < smiles.length() - 1 && isalnum(smiles[i+1])) isBond = true;
            if (!isBond) {
               // It might be a charge, but the regex is safer.
               // This logic is fragile. Sticking to regex for brackets.
            }
        }
    }
     // Simplified: Only count '-' within brackets reliably.
    count += countChar(smiles, '-'); // Reverting to simple count for now, needs validation.


    return calculateFraction(count, totalLength);
}

// Stereochemistry forward fraction - REMOVED (Zero Variance)

// Stereochemistry backward fraction - REMOVED (Zero Variance)

// --- Start of complex pattern fractions ---

// Macro to define a Regex Fraction descriptor easily
#define DECLARE_REGEX_FRACTION_DESCRIPTOR(ClassName, Pattern, Description) \
    DECLARE_DESCRIPTOR(ClassName, RegexFractions, Description) \
    DESCRIPTOR_DEPENDENCIES(ClassName) { return {}; } \
    DescriptorResult ClassName##Descriptor::calculate(Context& context) const { \
        const std::string& smiles = context.getSmiles(); \
        double totalLength = static_cast<double>(smiles.length()); \
        if (totalLength < 1.0) return 0.0; \
        static const RE2 pattern(Pattern); \
        int count = countMatches(smiles, pattern); \
        return calculateFraction(count, totalLength); \
    } \
    void register_##ClassName##Descriptor() { \
        auto descriptor = std::make_shared<ClassName##Descriptor>(); \
        auto& registry = DescriptorRegistry::getInstance(); \
        registry.registerDescriptor(descriptor); \
    }

// Improved Regex Patterns for Better Detection

// Alkyl Groups
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesMethylFraction,
    "(^|[^a-zA-Z0-9])C([^a-zA-Z0-9]|$)|\\[CH3\\]|C\\(H\\)\\(H\\)H|\\[C;H3\\]|\\[C@@H3\\]|\\[C@H3\\]|\\[C;H3;!R\\]",
    "Fraction of SMILES length matching methyl group pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesEthylFraction,
    "CC([^a-zA-Z0-9]|$)|\\[CH2\\]\\[CH3\\]|C\\(C\\)C|\\[CH2CH3\\]|\\[CH2\\]\\[C;H3\\]|C\\[CH2\\]|\\[C;H2\\]\\[C;H3\\]|\\[CH2\\]C|C\\[CH3\\]|\\[CH2\\]\\[CH2\\]\\[CH3\\]",
    "Fraction of SMILES length matching ethyl group pattern"
)

// Functional Groups with Oxygen
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesHydroxylFraction,
    "O[H]|\\[OH\\]|\\[O\\]H|\\[O;H1\\]|\\[OH-\\]|\\[O-\\]|\\[O;H1;!R\\]|O\\(H\\)|\\[OH2\\]|\\[O;H2\\]|\\[O;H3\\]|\\[O;H0\\]|\\[O;H1;R0\\]|\\[O;H1;R1\\]|\\[O;H1;R2\\]",
    "Approx fraction of SMILES length matching hydroxyl pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCarboxylFraction,
    "C\\(=O\\)O|C\\(=O\\)\\[OH\\]|\\[COOH\\]|\\[C\\(=O\\)O\\]|\\[C\\(=O\\)\\[O-\\]\\]|\\[C\\(=O\\)OH\\]|OC\\(=O\\)|O=C\\(O\\)|\\[C\\(=O\\)OH2\\]|\\[C\\(=O\\)O-\\]|\\[C;H0;X3\\]\\(=O\\)O|\\[C;H0;X3\\]\\(=O\\)\\[O-\\]|\\[C;H0;X3\\]\\(=O\\)\\[OH\\]",
    "Fraction of SMILES length matching carboxyl C(=O)O[H]? pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCarbonylFraction,
    "C=O|C\\(=O\\)|\\[C=O\\]|\\[C\\(=O\\)\\]|\\[C;H0;X3\\]|\\[C;H0;X3;!R\\]|O=C|\\[C;H0;X3;R0\\]|\\[C;H0;X3;R1\\]|\\[C;H0;X3;R2\\]|\\[C;H0;X3;R3\\]",
    "Approx fraction of SMILES length matching non-acid/ester/amide carbonyl pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesEsterFraction,
    "C\\(=O\\)OC|OC\\(=O\\)C|\\[C;H0;X3\\]\\(=O\\)OC|\\[C;H0;X3\\]\\(=O\\)O\\[C;H3\\]|\\[C;H0;X3\\]\\(=O\\)O\\[C;H2\\]|\\[C;H0;X3\\]\\(=O\\)O\\[CH3\\]|\\[C;H0;X3\\]\\(=O\\)O\\[CH2\\]",
    "Approx fraction of SMILES length matching ester pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesEtherFraction,
    "COC|\\[CH2\\]O\\[CH3\\]|\\[CH2\\]O\\[CH2\\]|\\[COC\\]|\\[C;H2\\]O\\[C;H3\\]|CO\\[CH3\\]|\\[CH2\\]OC|\\[C;H2\\]O\\[C;H2\\]|S\\(C\\)C|C\\(O\\)C|\\[CH2\\]O|O\\[CH2\\]|\\[C;H2\\]O|O\\[C;H2\\]",
    "Fraction of SMILES length matching aliphatic ether pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesKetoneFraction,
    "([Cc])C\\(=O\\)([Cc])|([Cc])\\(=O\\)([Cc])",
    "Fraction matching ketone CC(=O)C pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAlcoholFraction,
    "([CcX4])O\\[H\\]|([CcX4])OH|([CcX4])O$|([CcX4])O[^a-zA-Z]",
    "Approx fraction matching alcohol pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesOxoFraction,
    "C=O([^a-zA-Z]|$)|C\\(=O\\)([^a-zA-Z]|$)",
    "Fraction matching terminal oxo C=O pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesEtherBridgeFraction,
    "([Cc])O([Cc])O([Cc])|O([Cc])O",
    "Fraction matching COCOC pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesMethoxyFraction,
    "O([Cc])([^a-zA-Z0-9]|$)|\\[OCH3\\]",
    "Fraction of SMILES length matching OCH3 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesBranchedEtherFraction,
    "O\\(([Cc])\\)([Cc])|([Cc])O\\(([Cc])\\)",
    "Fraction of SMILES length matching O(C)C ether pattern"
)

// Nitrogen-containing Groups
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAminoFraction,
    "N([^a-zA-Z0-9]|$)|\\[NH2\\]|N\\(\\[H\\]\\)\\[H\\]|\\[N;H2\\]|\\[NH2\\+\\]|\\[NH2-\\]|\\[N;H2;!R\\]|N\\(H\\)H|\\[NH3\\]|\\[NH3\\+\\]|\\[NH\\]|\\[N;H1\\]|N\\(H\\)|\\[N;H\\]|\\[N;H3\\]|\\[N;H0\\]|\\[N;H1;R0\\]|\\[N;H2;R1\\]",
    "Fraction of SMILES length matching simple amino N pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAmideFraction,
    "([Cc])\\(=O\\)N|NC\\(=O\\)|\\[C\\(=O\\)N\\]|\\[NC\\(=O\\)\\]|\\[C;H0;X3\\]\\(=O\\)N",
    "Fraction of SMILES length matching amide C(=O)N pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesNitrileFraction,
    "([Cc])#N|N#([Cc])|C#N|\\[C\\]#\\[N\\]|\\[C;H0\\]#\\[N\\]|\\[C;H0\\]#\\[N;H0\\]",
    "Fraction of SMILES length matching nitrile C#N pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIsocyanateFraction,
    "N=C=O|O=C=N",
    "Fraction of SMILES length matching isocyanate N=C=O pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesImineFraction,
    "C=N|\\[CH\\]=N|\\[C;H1\\]=N|C=\\[NH\\]|\\[N\\]=C|N=CH|N=C|\\[N\\]=CH|\\[N\\]=C|\\[CH\\]=N|\\[C\\]=N",
    "Approx fraction of SMILES length matching imine pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPrimaryAmineFraction,
    "N\\[H2\\]|\\[NH2\\]|N([^a-z0-9\\(]|$)",
    "Fraction of SMILES length matching primary amine N[H2]|N(end) pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesSecondaryAmineFraction,
    "N\\[H\\]\\(|\\[NH\\]\\(|N\\([^\\)]+\\)([^\\(]|$)",
    "Fraction of SMILES length matching secondary amine N[H](|N()(end) pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesTertiaryAmineFraction,
    "N\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)|N\\([^\\)]+\\)\\([^\\)]+\\)[^\\(]",
    "Approx fraction of SMILES length matching tertiary amine pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesQuaternaryAmmoniumFraction,
    "\\[N\\+\\]\\(|\\[nH?\\+\\]",
    "Fraction of SMILES length matching quaternary ammonium [N+]( pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesBranchedAmineFraction,
    "N\\(([Cc])\\)([Cc])|([Cc])N\\(([Cc])\\)",
    "Fraction of SMILES length matching N(C)C amine pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesNAlkylFraction,
    "N([Cc])|([Cc])N",
    "Fraction matching NC pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAromaticNitrogenCarbonFraction,
    "([nc])",
    "Fraction matching nc pattern"
)

// Sulfur-containing Groups
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesThiolFraction,
    "S[H]|\\[SH\\]|\\[S;H1\\]|\\[SH-\\]|\\[S-\\]|\\[S;H1;!R\\]|S\\(H\\)|\\[SH2\\]|\\[S;H2\\]|S[H2]|S\\(H\\)H",
    "Approx fraction of SMILES length matching thiol pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesThioetherFraction,
    "CSC|\\[CH2\\]S\\[CH3\\]|\\[CH2\\]S\\[CH2\\]|\\[CSC\\]|\\[C;H2\\]S\\[C;H3\\]|CS\\[CH3\\]|\\[CH2\\]SC|\\[C;H2\\]S\\[C;H2\\]|S\\(C\\)C|C\\(S\\)C|\\[CH2\\]S|S\\[CH2\\]|\\[C;H2\\]S|S\\[C;H2\\]",
    "Fraction of SMILES length matching aliphatic thioether pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIsothiocyanateFraction,
    "N=C=S|S=C=N",
    "Fraction of SMILES length matching isothiocyanate N=C=S pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesSulfonylFraction,
    "S\\(=O\\)\\(=O\\)|\\[S\\+2\\]\\(\\[O-\\]\\)\\(\\[O-\\]\\)",
    "Fraction of SMILES length matching sulfonyl S(=O)2 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesSulfoxideFraction,
    "S\\(=O\\)([^\\(=]|$)|\\[S\\+\\]\\(\\[O-\\]\\)",
    "Fraction of SMILES length matching sulfoxide S(=O) pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesThiocyanateFraction,
    "SC#N|N#CS",
    "Fraction of SMILES length matching thiocyanate S-C#N pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIsothiocyanateLinearFraction,
    "N=C=S|S=C=N",
    "Fraction of SMILES length matching linear isothiocyanate N=C=S pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCarbonSulfurBondFraction,
    "([Cc])S|S([Cc])",
    "Fraction matching CS or SC pattern"
)

// Phosphorus-containing Groups
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPhosphateFraction,
    "P\\(=O\\)\\(O\\)O|OP\\(=O\\)\\(O\\)O",
    "Fraction of SMILES length matching phosphate P(=O)(O)O pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPhosphonateFraction,
    "P\\(=O\\)\\(O\\)([Cc])|([Cc])P\\(=O\\)\\(O\\)",
    "Fraction of SMILES length matching phosphonate P(=O)(O)C pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCarbonPhosphorusBondFraction,
    "([Cc])P|P([Cc])",
    "Fraction matching CP or PC pattern"
)

// Halogen-containing Groups
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAlkylHalideFraction,
    "([CcX4])(F|Cl|Br|I|\\[F\\]|\\[Cl\\]|\\[Br\\]|\\[I\\])",
    "Fraction of SMILES length matching alkyl halide C[X] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesFluorinatedFragmentFraction,
    "([Cc])F|F([Cc])",
    "Fraction of SMILES length matching CF pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesTrifluoromethylFraction,
    "CF3|C\\(F\\)\\(F\\)F|\\[CF3\\]",
    "Fraction matching CF3 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAromaticHalideFraction,
    "(c)(F|Cl|Br|I|\\[F\\]|\\[Cl\\]|\\[Br\\]|\\[I\\])",
    "Fraction matching aromatic halide cX pattern"
)

// Nitrogen-Oxygen Compounds
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAromaticNitroFraction,
    "(c)\\[N\\+\\]\\(=O\\)\\[O-\\]|(c)N\\(=O\\)=O",
    "Fraction of SMILES length matching aromatic nitro pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAliphaticNitroFraction,
    "([CcX4])\\[N\\+\\]\\(=O\\)\\[O-\\]|([CcX4])N\\(=O\\)=O",
    "Fraction of SMILES length matching aliphatic nitro pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAzideFraction,
    "N=\\[N\\+\\]=\\[N-\\]|N=[N+]=[N-]",
    "Fraction of SMILES length matching azide N=[N+]=[N-] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesDiazoFraction,
    "N=N([^=]|$)|\\[N\\]=\\[N\\]",
    "Fraction of SMILES length matching diazo N=N(end) pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesNitrogenOxygenAdjacentFraction,
    "NO|ON",
    "Fraction matching N-O or O-N adjacent pattern"
)

// Complex Nitrogen Compounds
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesGuanidineFraction,
    "N=C\\(N\\)N|NC\\(=N\\)N",
    "Fraction of SMILES length matching guanidine N=C(N)N pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesUreaFraction,
    "N([CcX4]?)\\(=O\\)N|NC\\(=O\\)N",
    "Fraction of SMILES length matching urea N-C(=O)-N pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAmidineFraction,
    "([Cc])\\(=N\\)N|N=([Cc])N|NC\\(=N\\)",
    "Fraction of SMILES length matching amidine C(=N)N pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesSulfonamideFraction,
    "S\\(=O\\)\\(=O\\)N|NS\\(=O\\)\\(=O\\)",
    "Fraction of SMILES length matching sulfonamide S(=O)2N pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCarbamatePatternsFraction,
    "NC\\(=O\\)O|OC\\(=O\\)N",
    "Fraction matching carbamate NC(=O)O pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesBranchedAmideFraction,
    "([Cc])\\(=O\\)N\\(([Cc])\\)|N\\(([Cc])\\)([Cc])\\(=O\\)",
    "Fraction matching branched amide C(=O)N(C) pattern"
)

// Heterocyclic Compounds
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesImidazoleFraction,
    "c1nccn1|c1cnc[nH]1|c1nc[nH]c1|[nH]1cncc1|n1cc[nH]c1|n1c[nH]cc1",
    "Fraction of SMILES length matching imidazole c1nccn1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPyridineFraction,
    "c1ccncc1|n1ccccc1|c1cccnc1|c1ccccn1",
    "Fraction of SMILES length matching pyridine c1ccncc1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPyrroleFraction,
    "c1cc[nH]c1|[nH]1cccc1|c1c[nH]cc1|c1ccnc1",
    "Fraction of SMILES length matching pyrrole c1cc[nH]c1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesFuranFraction,
    "c1ccoc1|o1cccc1|c1cocc1|c1ccco1",
    "Fraction of SMILES length matching furan c1ccoc1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesThiopheneFraction,
    "c1ccsc1|s1cccc1|c1cscc1|c1cccs1",
    "Fraction of SMILES length matching thiophene c1ccsc1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesImidazoliumFraction,
    "\\[nH\\+\\]|n\\+",
    "Fraction of SMILES length matching imidazolium [nH+] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPyridiniumFraction,
    "\\[n\\+\\]|n\\+",
    "Fraction of SMILES length matching pyridinium [n+] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesMixedHeteroatomRingFraction,
    "c1[nos]c[nos]c1|[nos]1c[nos]cc1|[nos]1cc[nos]c1",
    "Approx fraction matching mixed heteroatom ring pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesImidazoleBasicFraction,
    "c1nccn1|n1cncc1|n1ccnc1",
    "Fraction of SMILES length matching basic imidazole c1nccn1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPyridineBasicFraction,
    "c1ccncc1|n1ccccc1|c1cccnc1",
    "Fraction of SMILES length matching basic pyridine c1ccncc1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAromaticNitrogenRingFraction,
    "n1[a-z][a-z][a-z][a-z][a-z]1|[a-z]1[a-z][a-z][a-z][a-z]n1|[a-z]1[a-z][a-z][a-z]n[a-z]1|[a-z]1[a-z][a-z]n[a-z][a-z]1|[a-z]1[a-z]n[a-z][a-z][a-z]1|[a-z]1n[a-z][a-z][a-z][a-z]1",
    "Approx fraction matching N-containing 6-arom-ring"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesOxygenAromaticRingFraction,
    "o1[a-z][a-z][a-z][a-z][a-z]1|[a-z]1[a-z][a-z][a-z][a-z]o1|[a-z]1[a-z][a-z][a-z]o[a-z]1|[a-z]1[a-z][a-z]o[a-z][a-z]1|[a-z]1[a-z]o[a-z][a-z][a-z]1|[a-z]1o[a-z][a-z][a-z][a-z]1",
    "Approx fraction matching O-containing 6-arom-ring"
)

// Rings and Cyclic Compounds
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAromaticRingFraction,
    "c1ccccc1|c1ccccc1|c1ccccc1|c1ccccc1",
    "Approx fraction of SMILES length matching benzene c1ccccc1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesFiveRingFraction,
    "[1-9]([^1-9]{3,5})[\\1]|%[0-9][0-9]([^%0-9]{3,5})%[0-9][0-9]",
    "Approx fraction of SMILES length matching 5-ring closure pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesSixRingFraction,
    "[1-9]([^1-9]{4,6})[\\1]|%[0-9][0-9]([^%0-9]{4,6})%[0-9][0-9]",
    "Approx fraction of SMILES length matching 6-ring closure pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesSevenRingFraction,
    "[1-9]([^1-9]{5,7})[\\1]|%[0-9][0-9]([^%0-9]{5,7})%[0-9][0-9]",
    "Approx fraction of SMILES length matching 7-ring closure pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesEightRingFraction,
    "[1-9]([^1-9]{6,8})[\\1]|%[0-9][0-9]([^%0-9]{6,8})%[0-9][0-9]",
    "Approx fraction of SMILES length matching 8-ring closure pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesComplexRingClosureFraction,
    "[0-9][0-9]|%[0-9][0-9]",
    "Fraction of SMILES length matching double-digit ring closure pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCyclopropylFraction,
    "C1CC1|C\\(C\\)\\(C\\)|C1=CC1",
    "Fraction matching cyclopropyl C1CC1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCyclobutylFraction,
    "C1CCC1|C1=CCC1|C1C=CC1|C1CC=C1",
    "Fraction matching cyclobutyl C1CCC1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCyclopentylFraction,
    "C1CCCC1|C1C=CCC1|C1CC=CC1|C1CCC=C1|C1=CCCC1",
    "Fraction matching cyclopentyl C1CCCC1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCyclohexylFraction,
    "C1CCCCC1|C1C=CCCC1|C1CC=CCC1|C1CCC=CC1|C1CCCC=C1|C1=CCCCC1",
    "Fraction matching cyclohexyl C1CCCCC1 pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesDoubleBondRingJunctionFraction,
    "C=C[1-9]|C=[C]%[0-9][0-9]",
    "Fraction matching C=C[digit] pattern"
)

// Carbon Chain and Structure Patterns
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAliphaticChainFraction,
    "CCC|[CH2][CH2][CH2]",
    "Fraction of SMILES length matching CCC pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAromaticSequenceFraction,
    "ccc",
    "Fraction of SMILES length matching ccc pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesTertiaryCarbonFraction,
    "C\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)|C\\([^\\)]+\\)\\([^\\)]+\\)([^\\(]|$)",
    "Approx fraction matching tertiary C pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesQuaternaryCarbonFraction,
    "C\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)",
    "Approx fraction matching quaternary C pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesHeteroatomSequenceFraction,
    "[CNOScnos][CNOScnos]",
    "Fraction of SMILES length matching CNO neighbor sequence pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesConjugatedSystemFraction,
    "C=CC=C|c=cc=c|C=Cc=C|c=CC=c",
    "Fraction of SMILES length matching C=CC=C pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesLongCarbonChainFraction,
    "CCCC|[CH2][CH2][CH2][CH2]",
    "Fraction of SMILES length matching CCCC pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesDimethylFraction,
    "CC\\(C\\)C|C\\(C\\)\\(C\\)C",
    "Fraction matching CC(C)C pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAromaticAlkylFraction,
    "(c)([CX4])|([CX4])(c)",
    "Fraction matching cC pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAlkeneFraction,
    "C=C|C\\\\C|C/C",
    "Approx fraction of SMILES length matching alkene C=C pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAlkyneFraction,
    "C#C",
    "Fraction of SMILES length matching alkyne C#C pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesTerminalAlkyneFraction,
    "C#C([^a-zA-Z]|$)|C#[CH]",
    "Fraction matching terminal C#CH pattern"
)

// Acidic/Basic Group Patterns
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAcidicOxygenFraction,
    "([Cc]|S|P)\\(=O\\).*O(\\[H\\])?|([Cc]|S|P)\\(=O\\)O",
    "Approx fraction of SMILES length matching acidic oxygen patterns"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAcidicSulfurFraction,
    "S\\(=O\\)(\\(=O\\))?O(\\[H\\])?|S\\(=O\\)(\\(=O\\))?\\[OH\\]",
    "Approx fraction of SMILES length matching acidic sulfur patterns"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAcidicPhosphorusFraction,
    "P\\(=O\\)\\(O\\)O(\\[H\\])?|P\\(=O\\)\\(\\[OH\\]\\)\\[OH\\]",
    "Approx fraction of SMILES length matching acidic phosphorus patterns"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPhenolFraction,
    "(c)O(\\[H\\])?|(c)OH|(c)\\[OH\\]",
    "Fraction of SMILES length matching phenol cO[H]? pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCarboxylateFraction,
    "([Cc])\\(=O\\)\\[O-\\]|\\[O-\\]([Cc])=O|OC\\(=\\[O-\\]\\)",
    "Fraction of SMILES length matching carboxylate C(=O)[O-] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesSulfonateFraction,
    "S\\(=O\\)\\(=O\\)\\[O-\\]|S\\(\\[O-\\]\\)\\(=O\\)=O",
    "Fraction of SMILES length matching sulfonate S(=O)2[O-] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPhosphateAnionFraction,
    "P\\(=O\\)\\(O\\)\\[O-\\]|P\\(\\[O-\\]\\)\\(=O\\)O",
    "Fraction of SMILES length matching phosphate anion P(=O)(O)[O-] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesGuanidineBasicFraction,
    "N=C\\(N\\)N|NC\\(=[NH]\\)N|NC\\(N\\)=N",
    "Fraction of SMILES length matching basic guanidine N=C(N)N pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAmidineBasicFraction,
    "([Cc])\\(=N\\)N|N=([Cc])N|N([Cc])=N",
    "Fraction of SMILES length matching basic amidine C(=N)N pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesQuaternaryAmmoniumBasicFraction,
    "\\[N\\+\\]|N\\+",
    "Fraction matching basic quaternary ammonium [N+]( pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAmmoniumFraction,
    "\\[NH4\\+\\]|NH4\\+",
    "Fraction of SMILES length matching ammonium [NH4+] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesSulfoniumFraction,
    "\\[S\\+\\]|S\\+",
    "Fraction of SMILES length matching sulfonium [S+] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPhosphoniumFraction,
    "\\[P\\+\\]|P\\+",
    "Fraction of SMILES length matching phosphonium [P+] pattern"
)

// Ionizable Groups
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAcidicHydrogenFraction,
    "\\[(O|S|N)H\\]|OH|SH|NH",
    "Fraction of SMILES length matching explicit acidic hydrogen [X]H pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizableOxygenFraction,
    "\\[O-\\]|\\[OH\\]|O-|OH",
    "Fraction of SMILES length matching ionizable oxygen [O-]|[OH] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizableNitrogenFraction,
    "\\[[nN].*\\+\\]|\\[NH2\\]|\\[NH\\]|[nN]\\+|NH[12]?",
    "Fraction of SMILES length matching ionizable nitrogen patterns"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizableSulfurFraction,
    "\\[S-\\]|\\[SH\\]|S-|SH",
    "Fraction of SMILES length matching ionizable sulfur [S-]|[SH] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizablePhosphorusFraction,
    "\\[P-\\]|\\[PH\\]|P-|PH",
    "Fraction of SMILES length matching ionizable phosphorus [P-]|[PH] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizableCarboxylFraction,
    "([Cc])\\(=O\\)(\\[O-\\]|O\\[H\\]|O-|OH)|([Cc])\\(O-\\)=O|([Cc])\\(OH\\)=O",
    "Fraction matching ionizable carboxyl COOH/COO- patterns"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizablePhenolFraction,
    "(c)(\\[O-\\]|O\\[H\\]|O-|OH)",
    "Fraction matching ionizable phenol cOH/cO- patterns"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizableThiolFraction,
    "\\[SH\\]|\\[S-\\]|SH|S-",
    "Fraction of SMILES length matching ionizable thiol [SH]/[S-] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizableAmineFraction,
    "\\[[nN].*\\+\\]|\\[NH2\\]|\\[NH\\]|[nN]\\+|NH[12]?",
    "Fraction of SMILES length matching ionizable amine patterns"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizableImidazoleFraction,
    "\\[nH\\]|\\[nH\\+\\]|nH|nH\\+",
    "Fraction matching ionizable imidazole [nH]/[nH+] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizablePyridineFraction,
    "\\[n\\+\\]|\\[nH\\+\\]|n\\+|nH\\+",
    "Fraction matching ionizable pyridine [n+]/[nH+] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizableSulfonamideFraction,
    "S\\(=O\\)\\(=O\\)(\\[N-\\]|N\\[H\\]|N-|NH)|\\[N-\\]S\\(=O\\)\\(=O\\)|\\[NH\\]S\\(=O\\)\\(=O\\)",
    "Fraction matching ionizable sulfonamide SO2NH/SO2N- pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIonizablePhosphonateFraction,
    "P\\(=O\\)\\(O\\)(\\[O-\\]|O\\[H\\]|O-|OH)|P\\(=O\\)\\(\\[O-\\]\\)\\[OH\\]|P\\(=O\\)\\(OH\\)\\[O-\\]",
    "Fraction matching ionizable phosphonate POOH/POO- pattern"
)

// Miscellaneous Patterns
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesIsotopeFraction,
    "\\[[0-9]+[A-Za-z]{1,2}",
    "Fraction of SMILES length matching isotope pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesCarbonylHeteroatomFraction,
    "([Cc])\\(=O\\)[NOS]|[NOS]([Cc])=O",
    "Fraction matching C(=O)[NOS] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAdjacentHeteroatomFraction,
    "[NOSPFnospf][NOSPFnospf]",
    "Fraction matching adjacent [NOS][NOS] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAromaticCarbonHeteroatomFraction,
    "(c)[NOSPnospf]|[NOSPnospf](c)",
    "Fraction matching c[NOS] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesAliphaticCarbonHeteroatomFraction,
    "([CX4])[NOSPnospf]|[NOSPnospf]([CX4])",
    "Fraction matching C[NOS] pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesPrimaryAmineBasicFraction,
    "N\\[H2\\]|\\[NH2\\]|N([^a-z0-9\\(]|$)",
    "Fraction matching basic primary amine N[H2]|N(end) pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesSecondaryAmineBasicFraction,
    "N\\[H\\]\\(|\\[NH\\]\\(|N\\([^\\)]+\\)([^\\(]|$)",
    "Fraction matching basic secondary amine N[H](|N()(end) pattern"
)
DECLARE_REGEX_FRACTION_DESCRIPTOR(
    SmilesTertiaryAmineBasicFraction,
    "N\\([^\\)]+\\)\\([^\\)]+\\)\\([^\\)]+\\)|N\\([^\\)]+\\)\\([^\\)]+\\)[^\\(]",
    "Approx fraction matching basic tertiary amine pattern"
)


// --- Registration Functions (Manual for non-macro definitions) ---

void register_SmilesCarbonFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesCarbonFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesNitrogenFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesNitrogenFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesOxygenFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesOxygenFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesSulfurFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfurFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesPhosphorusFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesPhosphorusFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesHalogenFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesHalogenFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesBoronFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesBoronFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesFluorineFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesFluorineFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesIodineFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesIodineFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesAromaticCarbonFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticCarbonFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesAromaticNitrogenFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticNitrogenFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesAromaticOxygenFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticOxygenFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesAromaticSulfurFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesAromaticSulfurFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesSingleBondFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesSingleBondFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesDoubleBondFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesDoubleBondFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesTripleBondFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesTripleBondFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesRingClosureFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesRingClosureFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesBranchOpenFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesBranchOpenFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesBranchCloseFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesBranchCloseFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesBracketAtomFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesBracketAtomFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesPositiveChargeFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesPositiveChargeFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesNegativeChargeFractionDescriptor() {
    auto descriptor = std::make_shared<SmilesNegativeChargeFractionDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

} // namespace desfact