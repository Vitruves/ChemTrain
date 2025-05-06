#include "../common.hpp"
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/MolOps.h> // For distance matrix
#include <GraphMol/AtomIterators.h>
#include <vector>
#include <cmath>

namespace desfact {

// --- Broto-Moreau Autocorrelation Calculation ---

// Structure to cache precomputed data
struct BrotoMoreauCache {
    std::vector<std::vector<int>> distMatrix;
    // Original properties
    std::vector<double> masses;
    std::vector<double> vdwVolumes;
    std::vector<double> electronegativities;
    std::vector<double> polarizabilities;
    // New properties
    std::vector<double> atomicNumbers;
    std::vector<double> covalentRadii;
    std::vector<double> ionizationEnergies;
    std::vector<double> degrees;
    std::vector<double> valenceElectrons;

    // Original autocorrelation vectors
    std::vector<double> atsMass;
    std::vector<double> atsVolume;
    std::vector<double> atsEN;
    std::vector<double> atsPolarizability;
    // New autocorrelation vectors
    std::vector<double> atsZ;
    std::vector<double> atsR;
    std::vector<double> atsI;
    std::vector<double> atsD;
    std::vector<double> atsS;

    const RDKit::ROMol* populatedForMol = nullptr;
    unsigned int maxLag = 8; // Default max lag

    void reset() {
        distMatrix.clear();
        // Clear original property vectors
        masses.clear();
        vdwVolumes.clear();
        electronegativities.clear();
        polarizabilities.clear();
        // Clear new property vectors
        atomicNumbers.clear();
        covalentRadii.clear();
        ionizationEnergies.clear();
        degrees.clear();
        valenceElectrons.clear();
        
        // Reset original autocorrelation vectors
        atsMass.assign(maxLag + 1, 0.0);
        atsVolume.assign(maxLag + 1, 0.0);
        atsEN.assign(maxLag + 1, 0.0);
        atsPolarizability.assign(maxLag + 1, 0.0);
        // Reset new autocorrelation vectors
        atsZ.assign(maxLag + 1, 0.0);
        atsR.assign(maxLag + 1, 0.0);
        atsI.assign(maxLag + 1, 0.0);
        atsD.assign(maxLag + 1, 0.0);
        atsS.assign(maxLag + 1, 0.0);
        
        populatedForMol = nullptr;
    }

    // Calculate atomic properties and distance matrix
    void computePropertiesAndDistances(const RDKit::ROMol* mol) {
        if (!mol) return;
        unsigned int nAtoms = mol->getNumAtoms();
        
        // Resize original property vectors
        masses.resize(nAtoms);
        vdwVolumes.resize(nAtoms);
        electronegativities.resize(nAtoms);
        polarizabilities.resize(nAtoms);
        // Resize new property vectors
        atomicNumbers.resize(nAtoms);
        covalentRadii.resize(nAtoms);
        ionizationEnergies.resize(nAtoms);
        degrees.resize(nAtoms);
        valenceElectrons.resize(nAtoms);

        for (unsigned int i = 0; i < nAtoms; ++i) {
            const RDKit::Atom* atom = mol->getAtomWithIdx(i);
            int atomicNum = atom->getAtomicNum();
            
            // Original properties
            masses[i] = element::getAtomicMass(atomicNum);
            double vdwRadius = element::getVanDerWaalsRadius(atomicNum);
            vdwVolumes[i] = (4.0/3.0) * M_PI * std::pow(vdwRadius * 1e10, 3); // Convert m to Å
            electronegativities[i] = element::getElectronegativity(atomicNum);
            polarizabilities[i] = element::getPolarizability(atomicNum);
            
            // New properties
            atomicNumbers[i] = static_cast<double>(atomicNum);
            covalentRadii[i] = element::getCovalentRadius(atomicNum) * 1e10; // Convert m to Å
            ionizationEnergies[i] = element::getIonizationEnergy(atomicNum);
            degrees[i] = static_cast<double>(atom->getDegree());
            valenceElectrons[i] = static_cast<double>(atom->getTotalValence());
        }

        distMatrix = mol::calculateDistanceMatrix(mol);
    }

    // Calculate autocorrelation for a specific property
    void computeAutocorrelation(const std::vector<double>& prop, std::vector<double>& atsVector) {
        unsigned int nAtoms = prop.size();
        atsVector.assign(maxLag + 1, 0.0); // Reset specific vector
        std::vector<unsigned int> pairCounts(maxLag + 1, 0); // Count of pairs at each distance

        // Calculate mean of property values
        double mean = 0.0;
        for (double val : prop) {
            mean += val;
        }
        mean /= nAtoms;

        // Center the property values
        std::vector<double> centeredProp(nAtoms);
        for (unsigned int i = 0; i < nAtoms; ++i) {
            centeredProp[i] = prop[i] - mean;
        }

        // Include self-correlations at distance 0
        for (unsigned int i = 0; i < nAtoms; ++i) {
            atsVector[0] += centeredProp[i] * centeredProp[i];
            pairCounts[0]++;
        }

        // Calculate correlations for all other distances
        for (unsigned int i = 0; i < nAtoms; ++i) {
            for (unsigned int j = i + 1; j < nAtoms; ++j) {
                int dist = distMatrix[i][j];
                if (dist > 0 && static_cast<unsigned int>(dist) <= maxLag) {
                    atsVector[dist] += centeredProp[i] * centeredProp[j];
                    pairCounts[dist]++;
                }
            }
        }

        // Normalize by number of pairs at each distance
        for (unsigned int d = 0; d <= maxLag; ++d) {
            if (pairCounts[d] > 0) {
                atsVector[d] /= pairCounts[d];
            }
        }
    }

    void computeAll(const RDKit::ROMol* mol) {
        if (!mol) return;
        if (populatedForMol == mol) return; // Already computed

        reset();
        populatedForMol = mol;

        computePropertiesAndDistances(mol);

        // Compute original autocorrelations
        computeAutocorrelation(masses, atsMass);
        computeAutocorrelation(vdwVolumes, atsVolume);
        computeAutocorrelation(electronegativities, atsEN);
        computeAutocorrelation(polarizabilities, atsPolarizability);
        
        // Compute new autocorrelations
        computeAutocorrelation(atomicNumbers, atsZ);
        computeAutocorrelation(covalentRadii, atsR);
        computeAutocorrelation(ionizationEnergies, atsI);
        computeAutocorrelation(degrees, atsD);
        computeAutocorrelation(valenceElectrons, atsS);
    }
};

static thread_local BrotoMoreauCache tlsBrotoMoreauCache;

// Helper to ensure cache is populated
inline void ensureCachePopulated(Context& context) {
    const RDKit::ROMol* mol = context.getMolecule();
    if (tlsBrotoMoreauCache.populatedForMol != mol || !mol) {
        tlsBrotoMoreauCache.computeAll(mol);
    }
}

// --- Base Class ---
class BrotoMoreauDescriptor : public Descriptor {
public:
    using Descriptor::Descriptor; // Inherit constructor
    std::string getCategory() const override { return "BrotoMoreauAutocorrelation"; }
    std::vector<std::string> getDependencies() const override { return {}; } // Self-contained calculation

protected:
    unsigned int lag;
    std::string propertyName;

    BrotoMoreauDescriptor(unsigned int d, const std::string& prop, const std::string& desc)
        : lag(d), propertyName(prop), baseName(prop + std::to_string(d)), description(desc) {}

    std::string getName() const override { return baseName + "Descriptor"; }
    std::string getDescription() const override { return description; }

    virtual const std::vector<double>& getATSVector(const BrotoMoreauCache& cache) const = 0;

    DescriptorResult calculate(Context& context) const override {
        ensureCachePopulated(context);
        const auto& atsVector = getATSVector(tlsBrotoMoreauCache);
        if (lag < atsVector.size()) {
            return atsVector[lag];
        }
        return 0.0; // Return 0 if lag is out of bounds
    }

private:
    std::string baseName;
    std::string description;
};

// --- Concrete Descriptor Classes for each property and lag ---

// Macro to define descriptors for a specific property
#define DEFINE_BROTO_MOREAU_DESCRIPTORS(PROP_NAME, PROP_IDENTIFIER, ATS_VECTOR_MEMBER) \
    class PROP_NAME##Descriptor : public BrotoMoreauDescriptor { \
    public: \
        PROP_NAME##Descriptor(unsigned int d) \
            : BrotoMoreauDescriptor(d, #PROP_IDENTIFIER, "Broto-Moreau Autocorrelation - " #PROP_NAME " - Lag " + std::to_string(d)) {} \
    protected: \
        const std::vector<double>& getATSVector(const BrotoMoreauCache& cache) const override { \
            return cache.ATS_VECTOR_MEMBER; \
        } \
    }; \
    void register_##PROP_NAME##Descriptors() { \
        auto& registry = DescriptorRegistry::getInstance(); \
        for (unsigned int d = 1; d <= 8; ++d) { \
            registry.registerDescriptor(std::make_shared<PROP_NAME##Descriptor>(d)); \
        } \
    }

// Define original descriptors for Mass (m), Volume (v), Electronegativity (e), Polarizability (p)
DEFINE_BROTO_MOREAU_DESCRIPTORS(ATSdMass, m, atsMass);
DEFINE_BROTO_MOREAU_DESCRIPTORS(ATSdVolume, v, atsVolume);
DEFINE_BROTO_MOREAU_DESCRIPTORS(ATSdEN, e, atsEN);
DEFINE_BROTO_MOREAU_DESCRIPTORS(ATSdPolarizability, p, atsPolarizability);

// Define new descriptors for Atomic Number (z), Covalent Radius (r), 
// Ionization Energy (i), Atom Degree (d), Valence Electrons (s)
DEFINE_BROTO_MOREAU_DESCRIPTORS(ATSdZ, z, atsZ);
DEFINE_BROTO_MOREAU_DESCRIPTORS(ATSdRadius, r, atsR);
DEFINE_BROTO_MOREAU_DESCRIPTORS(ATSdIonization, i, atsI);
DEFINE_BROTO_MOREAU_DESCRIPTORS(ATSdDegree, d, atsD);
DEFINE_BROTO_MOREAU_DESCRIPTORS(ATSdValence, s, atsS);

// --- Individual Registration Functions (Called by update_registry.py) ---
#define REGISTER_SINGLE_BROTO_MOREAU(PROP_NAME, PROP_IDENTIFIER, LAG) \
    void register_ATSd##PROP_IDENTIFIER##LAG##Descriptor() { \
        auto& registry = DescriptorRegistry::getInstance(); \
        registry.registerDescriptor(std::make_shared<PROP_NAME##Descriptor>(LAG)); \
    }

// Original properties - Mass
// REGISTER_SINGLE_BROTO_MOREAU(ATSdMass, m, 1) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdMass, m, 2) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdMass, m, 3) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdMass, m, 4) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdMass, m, 5) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdMass, m, 6) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdMass, m, 7) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdMass, m, 8) // Zero Variance

// Original properties - Volume
REGISTER_SINGLE_BROTO_MOREAU(ATSdVolume, v, 1)
REGISTER_SINGLE_BROTO_MOREAU(ATSdVolume, v, 2)
REGISTER_SINGLE_BROTO_MOREAU(ATSdVolume, v, 3)
REGISTER_SINGLE_BROTO_MOREAU(ATSdVolume, v, 4)
REGISTER_SINGLE_BROTO_MOREAU(ATSdVolume, v, 5)
REGISTER_SINGLE_BROTO_MOREAU(ATSdVolume, v, 6)
REGISTER_SINGLE_BROTO_MOREAU(ATSdVolume, v, 7)
REGISTER_SINGLE_BROTO_MOREAU(ATSdVolume, v, 8)

// Original properties - Electronegativity
REGISTER_SINGLE_BROTO_MOREAU(ATSdEN, e, 1)
REGISTER_SINGLE_BROTO_MOREAU(ATSdEN, e, 2)
REGISTER_SINGLE_BROTO_MOREAU(ATSdEN, e, 3)
REGISTER_SINGLE_BROTO_MOREAU(ATSdEN, e, 4)
REGISTER_SINGLE_BROTO_MOREAU(ATSdEN, e, 5)
REGISTER_SINGLE_BROTO_MOREAU(ATSdEN, e, 6)
REGISTER_SINGLE_BROTO_MOREAU(ATSdEN, e, 7)
REGISTER_SINGLE_BROTO_MOREAU(ATSdEN, e, 8)

// Original properties - Polarizability
REGISTER_SINGLE_BROTO_MOREAU(ATSdPolarizability, p, 1)
REGISTER_SINGLE_BROTO_MOREAU(ATSdPolarizability, p, 2)
REGISTER_SINGLE_BROTO_MOREAU(ATSdPolarizability, p, 3)
REGISTER_SINGLE_BROTO_MOREAU(ATSdPolarizability, p, 4)
REGISTER_SINGLE_BROTO_MOREAU(ATSdPolarizability, p, 5)
REGISTER_SINGLE_BROTO_MOREAU(ATSdPolarizability, p, 6)
REGISTER_SINGLE_BROTO_MOREAU(ATSdPolarizability, p, 7)
REGISTER_SINGLE_BROTO_MOREAU(ATSdPolarizability, p, 8)

// New properties - Atomic Number
REGISTER_SINGLE_BROTO_MOREAU(ATSdZ, z, 1)
REGISTER_SINGLE_BROTO_MOREAU(ATSdZ, z, 2)
REGISTER_SINGLE_BROTO_MOREAU(ATSdZ, z, 3)
REGISTER_SINGLE_BROTO_MOREAU(ATSdZ, z, 4)
REGISTER_SINGLE_BROTO_MOREAU(ATSdZ, z, 5)
REGISTER_SINGLE_BROTO_MOREAU(ATSdZ, z, 6)
REGISTER_SINGLE_BROTO_MOREAU(ATSdZ, z, 7)
REGISTER_SINGLE_BROTO_MOREAU(ATSdZ, z, 8)

// New properties - Covalent Radius
REGISTER_SINGLE_BROTO_MOREAU(ATSdRadius, r, 1)
REGISTER_SINGLE_BROTO_MOREAU(ATSdRadius, r, 2)
REGISTER_SINGLE_BROTO_MOREAU(ATSdRadius, r, 3)
REGISTER_SINGLE_BROTO_MOREAU(ATSdRadius, r, 4)
REGISTER_SINGLE_BROTO_MOREAU(ATSdRadius, r, 5)
REGISTER_SINGLE_BROTO_MOREAU(ATSdRadius, r, 6)
REGISTER_SINGLE_BROTO_MOREAU(ATSdRadius, r, 7)
REGISTER_SINGLE_BROTO_MOREAU(ATSdRadius, r, 8)

// New properties - Ionization Energy
// REGISTER_SINGLE_BROTO_MOREAU(ATSdIonization, i, 1) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdIonization, i, 2) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdIonization, i, 3) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdIonization, i, 4) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdIonization, i, 5) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdIonization, i, 6) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdIonization, i, 7) // Zero Variance
// REGISTER_SINGLE_BROTO_MOREAU(ATSdIonization, i, 8) // Zero Variance

// New properties - Atom Degree
REGISTER_SINGLE_BROTO_MOREAU(ATSdDegree, d, 1)
REGISTER_SINGLE_BROTO_MOREAU(ATSdDegree, d, 2)
REGISTER_SINGLE_BROTO_MOREAU(ATSdDegree, d, 3)
REGISTER_SINGLE_BROTO_MOREAU(ATSdDegree, d, 4)
REGISTER_SINGLE_BROTO_MOREAU(ATSdDegree, d, 5)
REGISTER_SINGLE_BROTO_MOREAU(ATSdDegree, d, 6)
REGISTER_SINGLE_BROTO_MOREAU(ATSdDegree, d, 7)
REGISTER_SINGLE_BROTO_MOREAU(ATSdDegree, d, 8)

// New properties - Valence Electrons
REGISTER_SINGLE_BROTO_MOREAU(ATSdValence, s, 1)
REGISTER_SINGLE_BROTO_MOREAU(ATSdValence, s, 2)
REGISTER_SINGLE_BROTO_MOREAU(ATSdValence, s, 3)
REGISTER_SINGLE_BROTO_MOREAU(ATSdValence, s, 4)
REGISTER_SINGLE_BROTO_MOREAU(ATSdValence, s, 5)
REGISTER_SINGLE_BROTO_MOREAU(ATSdValence, s, 6)
REGISTER_SINGLE_BROTO_MOREAU(ATSdValence, s, 7)
REGISTER_SINGLE_BROTO_MOREAU(ATSdValence, s, 8)

// --- Registration Functions ---
void register_BrotoMoreauDescriptors() {
    // register_ATSdMassDescriptors(); // Zero Variance
    register_ATSdVolumeDescriptors();
    register_ATSdENDescriptors();
    register_ATSdPolarizabilityDescriptors();
    
    register_ATSdZDescriptors();
    register_ATSdRadiusDescriptors();
    // register_ATSdIonizationDescriptors(); // Zero Variance
    register_ATSdDegreeDescriptors();
    register_ATSdValenceDescriptors();
}

// ATS (Autocorrelation of Topological Structure) Descriptor base class
class ATSDescriptor : public Descriptor {
public:
    using Descriptor::Descriptor;
    std::string getCategory() const override { return "AutocorrelationDescriptors"; }
};

// Implementations for each ATS identity descriptor
DECLARE_DESCRIPTOR(ATSdi1, ATSDescriptor, "ATS identity autocorrelation descriptor of lag 1")
DESCRIPTOR_DEPENDENCIES(ATSdi1) { return {}; }
DescriptorResult ATSdi1Descriptor::calculate(Context& context) const {
    const RDKit::ROMol* mol = context.getMolecule();
    if (!mol) return 0.0;
    
    unsigned int nAtoms = mol->getNumAtoms();
    auto distMatrix = mol::calculateDistanceMatrix(mol);
    
    double sum = 0.0;
    unsigned int pairCount = 0;
    
    for (unsigned int i = 0; i < nAtoms; ++i) {
        for (unsigned int j = i; j < nAtoms; ++j) {
            if (distMatrix[i][j] == 1) { // Lag 1
                sum += 1.0; // Identity function (always 1)
                pairCount++;
            }
        }
    }
    
    return pairCount > 0 ? sum / pairCount : 0.0;
}

DECLARE_DESCRIPTOR(ATSdi2, ATSDescriptor, "ATS identity autocorrelation descriptor of lag 2")
DESCRIPTOR_DEPENDENCIES(ATSdi2) { return {}; }
DescriptorResult ATSdi2Descriptor::calculate(Context& context) const {
    const RDKit::ROMol* mol = context.getMolecule();
    if (!mol) return 0.0;
    
    unsigned int nAtoms = mol->getNumAtoms();
    auto distMatrix = mol::calculateDistanceMatrix(mol);
    
    double sum = 0.0;
    unsigned int pairCount = 0;
    
    for (unsigned int i = 0; i < nAtoms; ++i) {
        for (unsigned int j = i; j < nAtoms; ++j) {
            if (distMatrix[i][j] == 2) { // Lag 2
                sum += 1.0;
                pairCount++;
            }
        }
    }
    
    return pairCount > 0 ? sum / pairCount : 0.0;
}

DECLARE_DESCRIPTOR(ATSdi3, ATSDescriptor, "ATS identity autocorrelation descriptor of lag 3")
DESCRIPTOR_DEPENDENCIES(ATSdi3) { return {}; }
DescriptorResult ATSdi3Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdi4, ATSDescriptor, "ATS identity autocorrelation descriptor of lag 4")
DESCRIPTOR_DEPENDENCIES(ATSdi4) { return {}; }
DescriptorResult ATSdi4Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdi5, ATSDescriptor, "ATS identity autocorrelation descriptor of lag 5")
DESCRIPTOR_DEPENDENCIES(ATSdi5) { return {}; }
DescriptorResult ATSdi5Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdi6, ATSDescriptor, "ATS identity autocorrelation descriptor of lag 6")
DESCRIPTOR_DEPENDENCIES(ATSdi6) { return {}; }
DescriptorResult ATSdi6Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdi7, ATSDescriptor, "ATS identity autocorrelation descriptor of lag 7")
DESCRIPTOR_DEPENDENCIES(ATSdi7) { return {}; }
DescriptorResult ATSdi7Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdi8, ATSDescriptor, "ATS identity autocorrelation descriptor of lag 8")
DESCRIPTOR_DEPENDENCIES(ATSdi8) { return {}; }
DescriptorResult ATSdi8Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

// Implementations for each ATS mass descriptor
DECLARE_DESCRIPTOR(ATSdm1, ATSDescriptor, "ATS mass autocorrelation descriptor of lag 1")
DESCRIPTOR_DEPENDENCIES(ATSdm1) { return {}; }
DescriptorResult ATSdm1Descriptor::calculate(Context& context) const {
    const RDKit::ROMol* mol = context.getMolecule();
    if (!mol) return 0.0;
    
    unsigned int nAtoms = mol->getNumAtoms();
    auto distMatrix = mol::calculateDistanceMatrix(mol);
    
    // Calculate atomic masses
    std::vector<double> masses(nAtoms);
    for (unsigned int i = 0; i < nAtoms; ++i) {
        const RDKit::Atom* atom = mol->getAtomWithIdx(i);
        masses[i] = element::getAtomicMass(atom->getAtomicNum());
    }
    
    // Calculate mean of masses
    double mean = 0.0;
    for (double mass : masses) {
        mean += mass;
    }
    mean /= nAtoms;
    
    // Center the masses
    std::vector<double> centeredMasses(nAtoms);
    for (unsigned int i = 0; i < nAtoms; ++i) {
        centeredMasses[i] = masses[i] - mean;
    }
    
    double sum = 0.0;
    unsigned int pairCount = 0;
    
    for (unsigned int i = 0; i < nAtoms; ++i) {
        for (unsigned int j = i; j < nAtoms; ++j) {
            if (distMatrix[i][j] == 1) { // Lag 1
                sum += centeredMasses[i] * centeredMasses[j];
                pairCount++;
            }
        }
    }
    
    return pairCount > 0 ? sum / pairCount : 0.0;
}

DECLARE_DESCRIPTOR(ATSdm2, ATSDescriptor, "ATS mass autocorrelation descriptor of lag 2")
DESCRIPTOR_DEPENDENCIES(ATSdm2) { return {}; }
DescriptorResult ATSdm2Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdm3, ATSDescriptor, "ATS mass autocorrelation descriptor of lag 3")
DESCRIPTOR_DEPENDENCIES(ATSdm3) { return {}; }
DescriptorResult ATSdm3Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdm4, ATSDescriptor, "ATS mass autocorrelation descriptor of lag 4")
DESCRIPTOR_DEPENDENCIES(ATSdm4) { return {}; }
DescriptorResult ATSdm4Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdm5, ATSDescriptor, "ATS mass autocorrelation descriptor of lag 5")
DESCRIPTOR_DEPENDENCIES(ATSdm5) { return {}; }
DescriptorResult ATSdm5Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdm6, ATSDescriptor, "ATS mass autocorrelation descriptor of lag 6")
DESCRIPTOR_DEPENDENCIES(ATSdm6) { return {}; }
DescriptorResult ATSdm6Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdm7, ATSDescriptor, "ATS mass autocorrelation descriptor of lag 7")
DESCRIPTOR_DEPENDENCIES(ATSdm7) { return {}; }
DescriptorResult ATSdm7Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

DECLARE_DESCRIPTOR(ATSdm8, ATSDescriptor, "ATS mass autocorrelation descriptor of lag 8")
DESCRIPTOR_DEPENDENCIES(ATSdm8) { return {}; }
DescriptorResult ATSdm8Descriptor::calculate(Context& context) const {
    return 0.0; // Placeholder
}

// Registration functions
void register_ATSdi1Descriptor() {
    auto descriptor = std::make_shared<ATSdi1Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdi2Descriptor() {
    auto descriptor = std::make_shared<ATSdi2Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdi3Descriptor() {
    auto descriptor = std::make_shared<ATSdi3Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdi4Descriptor() {
    auto descriptor = std::make_shared<ATSdi4Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdi5Descriptor() {
    auto descriptor = std::make_shared<ATSdi5Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdi6Descriptor() {
    auto descriptor = std::make_shared<ATSdi6Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdi7Descriptor() {
    auto descriptor = std::make_shared<ATSdi7Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdi8Descriptor() {
    auto descriptor = std::make_shared<ATSdi8Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdm1Descriptor() {
    auto descriptor = std::make_shared<ATSdm1Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdm2Descriptor() {
    auto descriptor = std::make_shared<ATSdm2Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdm3Descriptor() {
    auto descriptor = std::make_shared<ATSdm3Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdm4Descriptor() {
    auto descriptor = std::make_shared<ATSdm4Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdm5Descriptor() {
    auto descriptor = std::make_shared<ATSdm5Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdm6Descriptor() {
    auto descriptor = std::make_shared<ATSdm6Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdm7Descriptor() {
    auto descriptor = std::make_shared<ATSdm7Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

void register_ATSdm8Descriptor() {
    auto descriptor = std::make_shared<ATSdm8Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}

} // namespace desfact
