#include "../common.hpp"
#include <re2/re2.h>
#include <string>
#include <vector>
#include <cmath>
#include <GraphMol/ROMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace desfact {

// ============== pKa-BASED DESCRIPTORS ===============

// 2. Estimated Acidity of Carboxylic Acids
DECLARE_DESCRIPTOR(SmilesCarboxylicAcidAcidity, PkaEstimation, "Estimated pKa of carboxylic acid groups in the SMILES string (~ 4.2 for benzoic acid)")
DESCRIPTOR_DEPENDENCIES(SmilesCarboxylicAcidAcidity) { return {}; }
DescriptorResult SmilesCarboxylicAcidAcidityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    
    // Expanded: C(=O)O, C(O)=O, COOH, [C](=O)O, C(=O)[O-], [C](O)=O, [C](=O)[O-]
    static const RE2 pattern(
        "C\\(=O\\)O"                  // C(=O)O
        "|C\\(O\\)=O"                 // C(O)=O
        "|COOH"                       // COOH
        "|\\[C\\]\\(=O\\)O"           // [C](=O)O
        "|C\\(=O\\)\\[O-\\]"          // C(=O)[O-]
        "|\\[C\\]\\(O\\)=O"           // [C](O)=O
        "|\\[C\\]\\(=O\\)\\[O-\\]"    // [C](=O)[O-]
        "|C\\(=O\\)\\[OH\\]"          // C(=O)[OH]
        "|C\\(=O\\)O\\[H\\]"          // C(=O)O[H]
        "|C\\(=O\\)O\\[H\\]"          // C(=O)O[H]
        "|C\\(=O\\)O\\[H\\]"          // C(=O)O[H]
        "|C\\(=O\\)O\\[H\\]"          // C(=O)O[H]
        "|C\\(=O\\)O\\[H\\]"          // C(=O)O[H]
    );
    if (!RE2::PartialMatch(smiles, pattern)) {
        return 0.0;
    }
    
    // Expanded aromatic: c.*C(=O)O, c.*COOH, c.*C(=O)[O-]
    static const RE2 aromatic_pattern(
        "c.*(C\\(=O\\)O|COOH|C\\(=O\\)\\[O-\\]|C\\(=O\\)O\\[H\\])"
    );
    if (RE2::PartialMatch(smiles, aromatic_pattern)) {
        return 4.2;
    }
    
    return 4.8;
}


// 4. Estimated Acidity of Phenols
DECLARE_DESCRIPTOR(SmilesPhenolAcidity, PkaEstimation, "Estimated pKa of phenol groups in the SMILES string (~ 10.0)")
DESCRIPTOR_DEPENDENCIES(SmilesPhenolAcidity) { return {}; }
DescriptorResult SmilesPhenolAcidityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    
    // Expanded: cO, cOH, c[OH], c(O)
    static const RE2 pattern(
        "cO"               // cO
        "|cOH"             // cOH
        "|c\\[OH\\]"       // c[OH]
        "|c\\(O\\)"        // c(O)
        "|c\\(OH\\)"       // c(OH)
        "|c\\[O\\]H"       // c[O]H
        "|cO\\[H\\]"       // cO[H]
    );
    if (!RE2::PartialMatch(smiles, pattern)) {
        return 0.0;
    }
    
    static const RE2 nitro_pattern(
        "\\[N\\+\\]\\(=O\\)\\[O-\\]" // [N+](=O)[O-]
        "|NO2"                       // NO2
    );
    if (RE2::PartialMatch(smiles, nitro_pattern)) {
        return 7.2;
    }
    
    return 10.0;
}

// 7. Amine Count
DECLARE_DESCRIPTOR(SmilesAmineCount, PkaEstimation, "Count of amine groups (N) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesAmineCount) { return {}; }
DescriptorResult SmilesAmineCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Match N or NH not followed by =, +, or #
    static const RE2 pattern(
        "N"                 // N
        "|N(H)"             // N(H)
        "|N\\[H\\]"         // N[H]
        "|N\\(H\\)"         // N(H)
        "|N\\(H2\\)"        // N(H2)
        "|N\\[H2\\]"        // N[H2]
        "|N\\(H3\\)"        // N(H3)
        "|N\\[H3\\]"        // N[H3]
        "|N\\(H\\)\\[H\\]"  // N(H)[H]
        "|N\\(H\\)\\(H\\)"  // N(H)(H)
        "|N\\(H\\)\\(H\\)\\(H\\)" // N(H)(H)(H)
    );
    return countMatches(smiles, pattern);
}

// 8. Estimated Basicity of Amines
DECLARE_DESCRIPTOR(SmilesAmineBasicity, PkaEstimation, "Estimated pKa of amine groups in the SMILES string (~ 9-10)")
DESCRIPTOR_DEPENDENCIES(SmilesAmineBasicity) { return {}; }
DescriptorResult SmilesAmineBasicityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();

    // Check if the molecule contains amine
    static const RE2 pattern(
        "N"                 // N
        "|N(H)"             // N(H)
        "|N\\[H\\]"         // N[H]
        "|N\\(H\\)"         // N(H)
        "|N\\(H2\\)"        // N(H2)
        "|N\\[H2\\]"        // N[H2]
        "|N\\(H3\\)"        // N(H3)
        "|N\\[H3\\]"        // N[H3]
        "|N\\(H\\)\\[H\\]"  // N(H)[H]
        "|N\\(H\\)\\(H\\)"  // N(H)(H)
        "|N\\(H\\)\\(H\\)\\(H\\)" // N(H)(H)(H)
    );
    if (!RE2::PartialMatch(smiles, pattern)) {
        return 0.0; // No amine groups
    }

    // Check if it's aniline or derivative
    static const RE2 aniline_pattern(
        "cN"                // cN
        "|cN(H)"            // cN(H)
        "|cN\\[H\\]"        // cN[H]
        "|cN\\(H\\)"        // cN(H)
    );
    if (RE2::PartialMatch(smiles, aniline_pattern)) {
        return 4.6; // Approximate pKa for aniline
    }

    // Check if it's primary, secondary, or tertiary amine
    static const RE2 tertiary_pattern(
        "N\\([^)]*\\)\\([^)]*\\)\\([^)]*\\)" // N(R)(R)(R)
        "|N\\(C\\)\\(C\\)\\(C\\)"            // N(C)(C)(C)
        "|N\\(CH3\\)\\(CH3\\)\\(CH3\\)"      // N(CH3)(CH3)(CH3)
    );
    if (RE2::PartialMatch(smiles, tertiary_pattern)) {
        return 9.8; // Tertiary amines (e.g., Trimethylamine ~9.8)
    }

    static const RE2 secondary_pattern(
        "N\\([^)]*\\)\\([^)]*\\)"            // N(R)(R)
        "|N\\(C\\)\\(C\\)"                   // N(C)(C)
        "|N\\(CH3\\)\\(CH3\\)"               // N(CH3)(CH3)
    );
    if (RE2::PartialMatch(smiles, secondary_pattern)) {
        return 10.7; // Secondary amines
    }

    return 10.6; // Primary amines (default)
}

// 10. Estimated Acidity of Thiols
DECLARE_DESCRIPTOR(SmilesThiolAcidity, PkaEstimation, "Estimated pKa of thiol groups in the SMILES string (~ 8-9)")
DESCRIPTOR_DEPENDENCIES(SmilesThiolAcidity) { return {}; }
DescriptorResult SmilesThiolAcidityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();

    // Expanded: S, SH, S[H], S(H)
    static const RE2 pattern(
        "S"                 // S
        "|SH"               // SH
        "|S\\[H\\]"         // S[H]
        "|S\\(H\\)"         // S(H)
    );
    if (!RE2::PartialMatch(smiles, pattern)) {
        return 0.0;
    }

    // Aromatic thiol: cS, cSH, c[SH], c(S)
    static const RE2 aromatic_pattern(
        "cS"                // cS
        "|cSH"              // cSH
        "|cS\\[H\\]"        // cS[H]
        "|cS\\(H\\)"        // cS(H)
    );
    if (RE2::PartialMatch(smiles, aromatic_pattern)) {
        return 6.5;
    }

    // Mercaptoethanol: OC.*S(H|\\[H\\])?([^=]|$)
    static const RE2 hych_pattern(
        "OC.*S"             // OC...S
        "|OC.*SH"           // OC...SH
        "|OC.*S\\[H\\]"     // OC...S[H]
        "|OC.*S\\(H\\)"     // OC...S(H)
    );
    if (RE2::PartialMatch(smiles, hych_pattern)) {
        return 9.5;
    }

    return 10.5;
}

// 11. Alpha-Hydrogen Acidity Index
DECLARE_DESCRIPTOR(SmilesAlphaHAcidityIndex, PkaEstimation, "Estimated alpha-hydrogen acidity based on neighboring groups")
DESCRIPTOR_DEPENDENCIES(SmilesAlphaHAcidityIndex) { return {}; }
DescriptorResult SmilesAlphaHAcidityIndexDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    
    // Check if molecule contains ketones or esters
    bool has_carbonyl = RE2::PartialMatch(smiles, RE2("C=O"));
    bool has_ester = RE2::PartialMatch(smiles, RE2("C\\(=O\\)O"));
    bool has_nitro = RE2::PartialMatch(smiles, RE2("\\[N\\+\\]\\(=O\\)\\[O-\\]"));
    bool has_cyano = RE2::PartialMatch(smiles, RE2("C#N"));
    
    // Estimate acidity index (higher means more acidic, lower pKa)
    double acidity_index = 0.0;
    
    if (has_carbonyl) acidity_index += 5.0;
    if (has_ester) acidity_index += 3.0;
    if (has_nitro) acidity_index += 8.0;
    if (has_cyano) acidity_index += 6.0;
    
    // Check for diketones (acetylacetone, etc.)
    if (RE2::PartialMatch(smiles, RE2("C\\(=O\\)CC\\(=O\\)"))) {
        acidity_index += 10.0; // Beta-diketones are more acidic
    }
    
    return acidity_index;
}

// 12. Imide Acidity
DECLARE_DESCRIPTOR(SmilesImideAcidity, PkaEstimation, "Estimated pKa of imide groups in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesImideAcidity) { return {}; }
DescriptorResult SmilesImideAcidityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    
    // Check if molecule contains imide pattern
    static const RE2 pattern("C\\(=O\\)NC\\(=O\\)");
    if (RE2::PartialMatch(smiles, pattern)) {
        // Check for succinimide
        if (RE2::PartialMatch(smiles, RE2("C1CC\\(=O\\)NC1=O"))) {
            return 9.5; // Succinimide pKa
        }
        
        // Check for phthalimide
        if (RE2::PartialMatch(smiles, RE2("c1ccccc1C\\(=O\\)NC\\(=O\\)"))) {
            return 8.3; // Phthalimide pKa
        }
        
        return 10.0; // Generic imide pKa
    }
    
    return 0.0; // No imide groups
}

// 13. Hydrogen Bond Donor Count
DECLARE_DESCRIPTOR(SmilesHBondDonorCount, PkaEstimation, "Count of hydrogen bond donors (OH, NH) in the SMILES string")
DESCRIPTOR_DEPENDENCIES(SmilesHBondDonorCount) { return {}; }
DescriptorResult SmilesHBondDonorCountDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    static const RE2 oh_pattern("O(H)?([^=)|$])");
    int oh_count = countMatches(smiles, oh_pattern);
    static const RE2 nh_pattern("N(H)?([^=+#]|$)");
    int nh_count = countMatches(smiles, nh_pattern);
    return oh_count + nh_count;
}

// 14. pKa-based Hydrogen Bond Strength Index
DECLARE_DESCRIPTOR(SmilesPkaHBondStrength, PkaEstimation, "Strength of hydrogen bonds based on pKa differences")
DESCRIPTOR_DEPENDENCIES(SmilesPkaHBondStrength) { return {}; }
DescriptorResult SmilesPkaHBondStrengthDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    double strength_index = 0.0;
    static const RE2 carboxylic_pattern("C\\(=O\\)O");
    if (RE2::PartialMatch(smiles, carboxylic_pattern)) {
        strength_index += 5.0;
    }
    static const RE2 phenol_pattern("cO(H)?([^=)|$])");
    if (RE2::PartialMatch(smiles, phenol_pattern)) {
        strength_index += 3.0;
    }
    static const RE2 alcohol_pattern("C.*O(H)?([^=)|$])");
    if (RE2::PartialMatch(smiles, alcohol_pattern)) {
        strength_index += 1.0;
    }
    static const RE2 amine_pattern("N(H)?([^=+#]|$)");
    if (RE2::PartialMatch(smiles, amine_pattern)) {
        strength_index += 2.0;
    }
    return strength_index;
}

// 15. Amide Resonance Stability Index
DECLARE_DESCRIPTOR(SmilesAmideResonanceIndex, PkaEstimation, "Resonance stability index for amides based on pKa")
DESCRIPTOR_DEPENDENCIES(SmilesAmideResonanceIndex) { return {}; }
DescriptorResult SmilesAmideResonanceIndexDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    
    // Check if molecule contains amide
    static const RE2 pattern("C\\(=O\\)N");
    if (!RE2::PartialMatch(smiles, pattern)) {
        return 0.0; // No amide groups
    }
    
    double resonance_index = 1.0; // Base value for amides
    
    // Electron-withdrawing groups on nitrogen reduce resonance
    static const RE2 n_ew_pattern("C\\(=O\\)N\\([^)]*\\[.+\\]\\)");
    if (RE2::PartialMatch(smiles, n_ew_pattern)) {
        resonance_index -= 0.3;
    }
    
    // Electron-withdrawing groups on carbonyl enhance resonance
    static const RE2 c_ew_pattern("\\[.+\\].*C\\(=O\\)N");
    if (RE2::PartialMatch(smiles, c_ew_pattern)) {
        resonance_index += 0.3;
    }
    
    return resonance_index;
}

// 17. Estimated Acidity of Sulfonic Acids
DECLARE_DESCRIPTOR(SmilesSulfonicAcidAcidity, PkaEstimation, "Estimated pKa of sulfonic acid groups (usually around -1 to 2)")
DESCRIPTOR_DEPENDENCIES(SmilesSulfonicAcidAcidity) { return {}; }
DescriptorResult SmilesSulfonicAcidAcidityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    
    // Check if the molecule contains sulfonic acid
    static const RE2 pattern("S\\(=O\\)\\(=O\\)O[H]?");
    if (!RE2::PartialMatch(smiles, pattern)) {
        return 0.0; // No sulfonic acid groups
    }
    
    // Approximate pKa for sulfonic acids
    return 0.5; // Most sulfonic acids have pKa around -1 to 2
}


// 19. Estimated Acidity of Phosphonic Acids
DECLARE_DESCRIPTOR(SmilesPhosphonicAcidAcidity, PkaEstimation, "Estimated pKa of phosphonic acid groups (first ~2.0, second ~7.0)")
DESCRIPTOR_DEPENDENCIES(SmilesPhosphonicAcidAcidity) { return {}; }
DescriptorResult SmilesPhosphonicAcidAcidityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    
    // Check if the molecule contains phosphonic acid
    static const RE2 pattern("P\\(=O\\)\\(O[H]?\\)O[H]?");
    if (!RE2::PartialMatch(smiles, pattern)) {
        return 0.0; // No phosphonic acid groups
    }
    
    // Return first pKa
    return 2.0; // First ionization of phosphonic acids
}

// 21. Estimated Basicity of Nitrogen Heterocycles
DECLARE_DESCRIPTOR(SmilesNHeterocycleBasicity, PkaEstimation, "Estimated pKa of nitrogen-containing heterocycles")
DESCRIPTOR_DEPENDENCIES(SmilesNHeterocycleBasicity) { return {}; }
DescriptorResult SmilesNHeterocycleBasicityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    
    // Check for pyridine
    static const RE2 pyridine_pattern("n1ccccc1");
    if (RE2::PartialMatch(smiles, pyridine_pattern)) {
        return 5.2; // pKa of pyridine
    }
    
    // Check for imidazole
    static const RE2 imidazole_pattern("n1cncc1");
    if (RE2::PartialMatch(smiles, imidazole_pattern)) {
        return 7.0; // pKa of imidazole
    }
    
    // Check for piperidine
    static const RE2 piperidine_pattern("N1CCCCC1");
    if (RE2::PartialMatch(smiles, piperidine_pattern)) {
        return 11.2; // pKa of piperidine
    }
    
    // Check for pyrrolidine
    static const RE2 pyrrolidine_pattern("N1CCCC1");
    if (RE2::PartialMatch(smiles, pyrrolidine_pattern)) {
        return 11.3; // pKa of pyrrolidine
    }
    
    return 0.0; // No basic heterocycles found
}

// 30. Estimated Acidity of Sulfonamides
DECLARE_DESCRIPTOR(SmilesSulfonamideAcidity, PkaEstimation, "Estimated pKa of sulfonamide groups in the SMILES string (~ 10)")
DESCRIPTOR_DEPENDENCIES(SmilesSulfonamideAcidity) { return {}; }
DescriptorResult SmilesSulfonamideAcidityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    
    // Check if the molecule contains sulfonamide
    static const RE2 pattern(
        "S\\(=O\\)\\(=O\\)N"         // S(=O)(=O)N
        "|S\\(=O\\)\\(=O\\)N(H)"     // S(=O)(=O)N(H)
        "|S\\(=O\\)\\(=O\\)N\\[H\\]" // S(=O)(=O)N[H]
        "|S\\(=O\\)\\(=O\\)N\\(H\\)" // S(=O)(=O)N(H)
    );
    if (!RE2::PartialMatch(smiles, pattern)) {
        return 0.0; // No sulfonamide groups
    }
    
    // From Bordwell table: Sulfonamides generally have pKa values ~10
    // Specific examples from table: Sulphadiazine 6.48, Sulphapyridine 8.43, Sulphathiazole 7.12
    
    // Check if it contains heterocyclic component
    static const RE2 heterocycle_pattern("S\\(=O\\)\\(=O\\)N.*n");
    if (RE2::PartialMatch(smiles, heterocycle_pattern)) {
        return 7.0; // Average for heterocyclic sulfonamides
    }
    
    return 10.0; // Default pKa for sulfonamides
}

// === NEW DESCRIPTORS BASED ON BORDWELL TABLE ===

// Hydroxamic Acid Acidity (from Bordwell p10, p11)
DECLARE_DESCRIPTOR(SmilesHydroxamicAcidAcidity, PkaEstimation, "Estimated pKa of hydroxamic acid groups (R-C(=O)NHOH) (~ 9.0)")
DESCRIPTOR_DEPENDENCIES(SmilesHydroxamicAcidAcidity) { return {}; }
DescriptorResult SmilesHydroxamicAcidAcidityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Pattern looks for C(=O)N followed by OH or O connected to something else
    static const RE2 pattern("C\\(=O\\)N.*O[H]?"); 
    if (RE2::PartialMatch(smiles, pattern)) {
        // From Bordwell table: Benzo- 8.88, Aceto- 9.40
        return 9.0; // Approximate pKa for hydroxamic acids
    }
    return 0.0; // No hydroxamic acid groups
}

// Sulfinic Acid Acidity (from Bordwell p11)
DECLARE_DESCRIPTOR(SmilesSulfinicAcidAcidity, PkaEstimation, "Estimated pKa of sulfinic acid groups (R-S(=O)OH) (~ 2.0)")
DESCRIPTOR_DEPENDENCIES(SmilesSulfinicAcidAcidity) { return {}; }
DescriptorResult SmilesSulfinicAcidAcidityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Pattern looks for S(=O)OH or S(=O)O
    static const RE2 pattern("S\\(=O\\)O[H]?"); 
    if (RE2::PartialMatch(smiles, pattern)) {
        // From Bordwell table: Benzene-/Toluene-sulfinic acids ~1.8-2.2
        return 2.0; // Approximate pKa for sulfinic acids
    }
    return 0.0; // No sulfinic acid groups
}

// Hydrazine Basicity (from Bordwell p28)
DECLARE_DESCRIPTOR(SmilesHydrazineBasicity, PkaEstimation, "Estimated pKaH of hydrazine groups (R-NHNH-R') (~ 8.0)")
DESCRIPTOR_DEPENDENCIES(SmilesHydrazineBasicity) { return {}; }
DescriptorResult SmilesHydrazineBasicityDescriptor::calculate(Context& context) const {
    const std::string& smiles = context.getSmiles();
    // Pattern looks for N, optional bracketed content [H], [+], etc., then N.
    // Avoids unsupported lookahead/lookbehind. Matches NN, N[H]N, N[H2]N etc.
    static const RE2 pattern("N(\\[[^\\]]+])?N"); 
    if (RE2::PartialMatch(smiles, pattern)) {
        // From Bordwell table: Hydrazine 8.07, Methylhydrazine 7.87
        // This is a heuristic and might misidentify some complex cases.
        return 8.0; // Approximate pKaH for hydrazines
    }
    return 0.0; // No simple hydrazine groups detected
}

// === REGISTRATION FUNCTIONS ===

void register_SmilesCarboxylicAcidCountDescriptor() {
    auto descriptor = std::make_shared<SmilesCarboxylicAcidAcidityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesCarboxylicAcidAcidityDescriptor() {
    auto descriptor = std::make_shared<SmilesCarboxylicAcidAcidityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesPhenolAcidityDescriptor() {
    auto descriptor = std::make_shared<SmilesPhenolAcidityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}



void register_SmilesAmineCountDescriptor() {
    auto descriptor = std::make_shared<SmilesAmineCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAmineBasicityDescriptor() {
    auto descriptor = std::make_shared<SmilesAmineBasicityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesThiolAcidityDescriptor() {
    auto descriptor = std::make_shared<SmilesThiolAcidityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAlphaHAcidityIndexDescriptor() {
    auto descriptor = std::make_shared<SmilesAlphaHAcidityIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesImideAcidityDescriptor() {
    auto descriptor = std::make_shared<SmilesImideAcidityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesHBondDonorCountDescriptor() {
    auto descriptor = std::make_shared<SmilesHBondDonorCountDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesPkaHBondStrengthDescriptor() {
    auto descriptor = std::make_shared<SmilesPkaHBondStrengthDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesAmideResonanceIndexDescriptor() {
    auto descriptor = std::make_shared<SmilesAmideResonanceIndexDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesSulfonicAcidAcidityDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfonicAcidAcidityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesPhosphonicAcidAcidityDescriptor() {
    auto descriptor = std::make_shared<SmilesPhosphonicAcidAcidityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesNHeterocycleBasicityDescriptor() {
    auto descriptor = std::make_shared<SmilesNHeterocycleBasicityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}


void register_SmilesSulfonamideAcidityDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfonamideAcidityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

// New registration functions
void register_SmilesHydroxamicAcidAcidityDescriptor() {
    auto descriptor = std::make_shared<SmilesHydroxamicAcidAcidityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesSulfinicAcidAcidityDescriptor() {
    auto descriptor = std::make_shared<SmilesSulfinicAcidAcidityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_SmilesHydrazineBasicityDescriptor() {
    auto descriptor = std::make_shared<SmilesHydrazineBasicityDescriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

} // namespace desfact