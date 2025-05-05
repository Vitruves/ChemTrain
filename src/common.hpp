// src/common.hpp
#pragma once
#include <GraphMol/ROMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/SmilesParse/SmilesWrite.h> // Needed for MolToSmiles
#include <variant>
#include <map>
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <any>
#include <unordered_map>
#include <re2/re2.h>
#include <optional>
#include <stdexcept>
#include <cmath>    // Needed for std::sqrt, std::log, std::pow, std::floor, std::log2, std::abs

namespace re2 { class RE2; } // Forward declaration for RE2

namespace desfact {

//==============================================================================
// CORE TYPES AND FORWARD DECLARATIONS
//==============================================================================
// Descriptor result (can hold any supported return type)
using DescriptorResult = std::variant<double, int, std::string, std::vector<double>>;

// Forward declarations
class Context;
class Observer;
class Descriptor;
class MoleculeTraverser;
class DescriptorPipeline;

//==============================================================================
// CONTEXT - SHARED CALCULATION ENVIRONMENT
//==============================================================================
class Context {
public:
    // Store data for reuse across calculations
    template<typename T>
    void setProperty(const std::string& key, const T& value) {
        properties_[key] = value;
    }

    // Retrieve stored data
    template<typename T>
    const T& getProperty(const std::string& key) const {
        auto it = properties_.find(key);
        if (it == properties_.end()) {
            throw std::runtime_error("Property not found: " + key);
        }
        try {
            return std::any_cast<const T&>(it->second);
        } catch (const std::bad_any_cast& e) {
             throw std::runtime_error("Property type mismatch for key: " + key + ". Error: " + e.what());
        }
    }

    // Check if property exists
    bool hasProperty(const std::string& key) const {
        return properties_.find(key) != properties_.end();
    }

    // Get or compute a property
    template<typename T>
    const T& getOrCompute(const std::string& key, std::function<T()> computeFunc) {
        if (!hasProperty(key)) {
            setProperty(key, computeFunc());
        }
        return getProperty<T>(key);
    }

    // Access current molecule
    const RDKit::ROMol* getMolecule() const { return molecule_; }

    // Reset context for next molecule
    void reset() {
        properties_.clear();
        molecule_ = nullptr;
        cachedSmiles_.reset();
    }

    const std::string& getSmiles() const {
        if (!molecule_) {
            throw std::runtime_error("Molecule not set in Context before calling getSmiles()");
        }
        if (!cachedSmiles_) {
            cachedSmiles_ = RDKit::MolToSmiles(*molecule_);
        }
        return *cachedSmiles_;
    }

private:
    friend class MoleculeTraverser;

    const RDKit::ROMol* molecule_ = nullptr;
    std::map<std::string, std::any> properties_;

    mutable std::optional<std::string> cachedSmiles_;

    // Set current molecule (called by traverser)
    void setMolecule(const RDKit::ROMol* mol) {
        molecule_ = mol;
        properties_.clear();
        cachedSmiles_.reset();
    }
};

//==============================================================================
// OBSERVER INTERFACE - COLLECTS PROPERTIES DURING TRAVERSAL
//==============================================================================
class Observer {
public:
    virtual ~Observer() = default;

    // Initialize before traversal
    virtual void init(Context& context) = 0;

    // Process an atom
    virtual void observeAtom(const RDKit::Atom* atom, Context& context) = 0;

    // Process a bond
    virtual void observeBond(const RDKit::Bond* bond, Context& context) = 0;

    // Finalize after traversal
    virtual void finalize(Context& context) = 0;

    // Get observer name
    virtual std::string getName() const = 0;

    // Get observer dependencies
    virtual std::vector<std::string> getDependencies() const = 0;
};

//==============================================================================
// DESCRIPTOR INTERFACE - CALCULATES MOLECULAR DESCRIPTORS
//==============================================================================
class Descriptor {
public:
    virtual ~Descriptor() = default;

    // Calculate the descriptor
    virtual DescriptorResult calculate(Context& context) const = 0;

    // Get descriptor name
    virtual std::string getName() const = 0;

    // Get descriptor category
    virtual std::string getCategory() const = 0;

    // Get descriptor dependencies (observers or other descriptors)
    virtual std::vector<std::string> getDependencies() const = 0;

    // Get description
    virtual std::string getDescription() const = 0;
};

//==============================================================================
// MOLECULE TRAVERSER - SINGLE-PASS MOLECULE ANALYSIS
//==============================================================================
class MoleculeTraverser {
public:
    // Add an observer
    void addObserver(std::shared_ptr<Observer> observer);

    // Traverse a molecule with all registered observers
    void traverse(const RDKit::ROMol* mol, Context& context);

    // Reset traverser
    void reset();

    const std::vector<std::shared_ptr<Observer>>& getObservers() const { return observers_; }

private:
    std::vector<std::shared_ptr<Observer>> observers_;

    // Resolve dependencies between observers
    void resolveObserverDependencies();
};

//==============================================================================
// DESCRIPTOR PIPELINE - COMPLETE CALCULATION WORKFLOW
//==============================================================================
class DescriptorPipeline {
public:
    // Add a descriptor to calculate
    void addDescriptor(std::shared_ptr<Descriptor> descriptor);

    // Add a descriptor by name
    void addDescriptor(const std::string& name);

    // Process a single molecule
    std::map<std::string, DescriptorResult> process(const RDKit::ROMol* mol);

    // Process a stream of molecules
    void processStream(std::istream& input, std::ostream& output);

    // Get the names of the descriptors in the pipeline
    std::vector<std::string> getDescriptorNames() const;

private:
    MoleculeTraverser traverser_;
    std::vector<std::shared_ptr<Descriptor>> descriptors_;

    // Resolve all dependencies
    void resolveDependencies();
};

//==============================================================================
// ELEMENT PROPERTIES - ATOMIC PROPERTY DATABASE
//==============================================================================
namespace element {
    // Get Pauling electronegativity for element
    double getElectronegativity(int atomicNum);

    // Get Allred-Rochow electronegativity for element
    double getAllredRochowEN(int atomicNum);

    // Get atomic mass in kg
    double getAtomicMass(int atomicNum);

    // Get ionization energy in J
    double getIonizationEnergy(int atomicNum);

    // Get electron affinity in J
    double getElectronAffinity(int atomicNum);

    // Get covalent radius in m
    double getCovalentRadius(int atomicNum);

    // Get van der Waals radius in m
    double getVanDerWaalsRadius(int atomicNum);

    // Get polarizability in Å³
    double getPolarizability(int atomicNum);

    // Get density in kg/m³
    double getDensity(int atomicNum);

    // Get thermal conductivity in W/m·K
    double getThermalConductivity(int atomicNum);

    // Get specific heat in J/g·K
    double getSpecificHeat(int atomicNum);

    // Get molar volume in cm³/mol
    double getMolarVolume(int atomicNum);

    // Get atomic number from element symbol
    int getAtomicNumber(const std::string& symbol);

    // Get element symbol from atomic number
    std::string getElementSymbol(int atomicNum);

    // Get periodic table group for element
    int getPeriodicGroup(int atomicNum);

    // Get periodic table period for element
    int getPeriodicPeriod(int atomicNum);

    // Check element classification
    bool isMainGroupElement(int atomicNum);
    bool isTransitionMetal(int atomicNum);
    bool isAlkaliMetal(int atomicNum);
    bool isAlkalineEarthMetal(int atomicNum);
    bool isNobleGas(int atomicNum);
    bool isHalogen(int atomicNum);
    bool isMetal(int atomicNum);

    // Check if element is a heteroatom (not C or H) - Defined based on atomic number
    inline bool isHeteroatom(int atomicNum) {
        return atomicNum != 6 && atomicNum != 1;
    }
}

//==============================================================================
// MOLECULAR UTILITIES - COMMON MOLECULE ANALYSIS FUNCTIONS
//==============================================================================
namespace mol {
    // Initialize ring info if needed
    void ensureRingInfo(const RDKit::ROMol* mol);

    // Count atoms/bonds matching a predicate
    int countAtoms(const RDKit::ROMol* mol, std::function<bool(const RDKit::Atom*)> predicate);
    int countBonds(const RDKit::ROMol* mol, std::function<bool(const RDKit::Bond*)> predicate);

    // Get atom/bond fraction matching a predicate
    double atomFraction(const RDKit::ROMol* mol, std::function<bool(const RDKit::Atom*)> predicate);
    double bondFraction(const RDKit::ROMol* mol, std::function<bool(const RDKit::Bond*)> predicate);

    // Bond type checks
    bool isSingleBond(const RDKit::Bond* bond);
    bool isDoubleBond(const RDKit::Bond* bond);
    bool isTripleBond(const RDKit::Bond* bond);
    bool isAromaticBond(const RDKit::Bond* bond);
    bool isRotatableBond(const RDKit::Bond* bond);

    // Atom type checks
    bool isHeteroatom(const RDKit::Atom* atom); // Overload for Atom*
    bool isHalogen(const RDKit::Atom* atom);
    bool isMetal(const RDKit::Atom* atom);
    bool isInRing(const RDKit::Atom* atom);
    bool isAromatic(const RDKit::Atom* atom);
    bool isBasicNitrogen(const RDKit::Atom* atom);
    bool isAcidicOxygen(const RDKit::Atom* atom);

    // Distance and connectivity
    int getAtomDistance(const RDKit::ROMol* mol, unsigned int idx1, unsigned int idx2);
    int getBondCount(const RDKit::ROMol* mol);
    int getHeavyAtomCount(const RDKit::ROMol* mol);
    int getRotatableBondCount(const RDKit::ROMol* mol);

    // Ring properties
    int getRingCount(const RDKit::ROMol* mol);
    int getAromaticRingCount(const RDKit::ROMol* mol);

    // Stereo chemistry
    bool hasStereoChemistry(const RDKit::ROMol* mol);

    // Statistical properties
    double getAverageAtomicMass(const RDKit::ROMol* mol);
    double getAverageElectronegativity(const RDKit::ROMol* mol);
    double getMaxElectronegativityDifference(const RDKit::ROMol* mol);
    double sumAtomProperty(const RDKit::ROMol* mol, std::function<double(const RDKit::Atom*)> propertyFn);
    double sumElectronegativity(const RDKit::ROMol* mol);
    double sumPolarizability(const RDKit::ROMol* mol);
    double sumCovalentRadii(const RDKit::ROMol* mol);

    // Property vector generation and statistics
    std::vector<double> getAtomPropertyVector(const RDKit::ROMol* mol, std::function<double(const RDKit::Atom*)> propertyFn);
    std::vector<double> getElectronegativityVector(const RDKit::ROMol* mol);
    std::vector<double> getPolarizabilityVector(const RDKit::ROMol* mol);
    double vectorMean(const std::vector<double>& vec);
    double vectorVariance(const std::vector<double>& vec);
    double vectorStdDev(const std::vector<double>& vec);
    double vectorMin(const std::vector<double>& vec);
    double vectorMax(const std::vector<double>& vec);
    double vectorRange(const std::vector<double>& vec);

    // Autocorrelation
    std::vector<double> calculateAutocorrelation(
        const RDKit::ROMol* mol,
        std::function<double(const RDKit::Atom*)> propertyFn,
        unsigned int maxDistance = 5
    );
    std::vector<double> getElectronegativityAutocorrelation(const RDKit::ROMol* mol, unsigned int maxDistance = 5);
    std::vector<double> getPolarizabilityAutocorrelation(const RDKit::ROMol* mol, unsigned int maxDistance = 5);

    // Distance matrix and related functions
    std::vector<std::vector<int>> calculateDistanceMatrix(const RDKit::ROMol* mol);
    std::vector<std::pair<const RDKit::Atom*, const RDKit::Atom*>>
    getAtomPairsAtDistance(const RDKit::ROMol* mol, unsigned int distance);
    int getAtomsAtDistanceCount(const RDKit::ROMol* mol, unsigned int atomIdx, unsigned int distance);

    // Graph indices
    int getMolecularEccentricity(const RDKit::ROMol* mol);
    int getMolecularDiameter(const RDKit::ROMol* mol);
    int getMolecularRadius(const RDKit::ROMol* mol);
    int getWienerIndex(const RDKit::ROMol* mol);
    double getBalabanJIndex(const RDKit::ROMol* mol);
    double getRandicConnectivityIndex(const RDKit::ROMol* mol);

    // Approximated properties
    double getTopologicalPolarSurfaceArea(const RDKit::ROMol* mol);
    double getFsp3(const RDKit::ROMol* mol);
    double getBertzComplexityIndex(const RDKit::ROMol* mol);
    double calculateLogP(const RDKit::ROMol* mol);
    double calculateMolarRefractivity(const RDKit::ROMol* mol);
    double calculateVanDerWaalsVolume(const RDKit::ROMol* mol);
    double calculateApproxSurfaceArea(const RDKit::ROMol* mol);

    // Complexity and structure
    int getFragmentCount(const RDKit::ROMol* mol);
    int getRingRotatableBondCount(const RDKit::ROMol* mol);

    // Vector/Histogram utilities
    std::vector<double> combineAtomPropertyVectors(
        const RDKit::ROMol* mol,
        std::function<double(const RDKit::Atom*)> propertyFn1,
        std::function<double(const RDKit::Atom*)> propertyFn2,
        std::function<double(double, double)> combiner
    );
    double harmonicMean(const std::vector<double>& vec);
    double geometricMean(const std::vector<double>& vec);
    std::map<int, int> generateAtomPropertyHistogram(
        const RDKit::ROMol* mol,
        std::function<double(const RDKit::Atom*)> propertyFn,
        double binWidth = 1.0
    );
    double calculatePropertyEntropy(
        const RDKit::ROMol* mol,
        std::function<double(const RDKit::Atom*)> propertyFn,
        double binWidth = 1.0
    );

}

//==============================================================================
// STANDARD OBSERVERS - PRE-BUILT PROPERTY COLLECTORS
//==============================================================================

// Element count observer
class ElementCountObserver : public Observer {
public:
    std::string getName() const override { return "ElementCountObserver"; }
    std::vector<std::string> getDependencies() const override { return {}; }

    void init(Context& context) override;
    void observeAtom(const RDKit::Atom* atom, Context& context) override;
    void observeBond(const RDKit::Bond* bond, Context& context) override;
    void finalize(Context& context) override;
};

// Ring information observer
class RingInfoObserver : public Observer {
public:
    std::string getName() const override { return "RingInfoObserver"; }
    std::vector<std::string> getDependencies() const override { return {}; }

    void init(Context& context) override;
    void observeAtom(const RDKit::Atom* atom, Context& context) override;
    void observeBond(const RDKit::Bond* bond, Context& context) override;
    void finalize(Context& context) override;
};

// Electronegativity observer
class ElectronegativityObserver : public Observer {
public:
    std::string getName() const override { return "ElectronegativityObserver"; }
    std::vector<std::string> getDependencies() const override { return {}; }

    void init(Context& context) override;
    void observeAtom(const RDKit::Atom* atom, Context& context) override;
    void observeBond(const RDKit::Bond* bond, Context& context) override;
    void finalize(Context& context) override;
};

//==============================================================================
// DESCRIPTOR MACROS - FOR EASY DESCRIPTOR CREATION
//==============================================================================

// For descriptor definitions
#ifndef REGISTRY_IMPLEMENTATION_ONLY
#define DECLARE_DESCRIPTOR(NAME, CATEGORY, DESCRIPTION) \
    class NAME##Descriptor : public Descriptor { \
    public: \
        std::string getName() const override { return #NAME "Descriptor"; } \
        std::string getCategory() const override { return #CATEGORY; } \
        std::string getDescription() const override { return DESCRIPTION; } \
        std::vector<std::string> getDependencies() const override; \
        DescriptorResult calculate(Context& context) const override; \
    };
#else
#define DECLARE_DESCRIPTOR(NAME, CATEGORY, DESCRIPTION)
#endif

// Implement descriptor dependencies
#define DESCRIPTOR_DEPENDENCIES(NAME) \
    std::vector<std::string> NAME##Descriptor::getDependencies() const

// Implement descriptor calculation
#define DESCRIPTOR_CALCULATE(NAME) \
    DescriptorResult NAME##Descriptor::calculate(Context& context) const

//==============================================================================
// DESCRIPTOR REGISTRY - GLOBAL REGISTRY (NOW FOR OBSERVERS TOO)
//==============================================================================

class DescriptorRegistry {
private:
    std::unordered_map<std::string, std::shared_ptr<Descriptor>> descriptors_;
    std::unordered_map<std::string, std::shared_ptr<Observer>> observers_;

    DescriptorRegistry() = default;

public:
    DescriptorRegistry(const DescriptorRegistry&) = delete;
    DescriptorRegistry& operator=(const DescriptorRegistry&) = delete;

    static DescriptorRegistry& getInstance() {
        static DescriptorRegistry instance;
        return instance;
    }

    // --- Descriptor Methods ---
    void registerDescriptor(std::shared_ptr<Descriptor> descriptor);

    std::shared_ptr<Descriptor> getDescriptor(const std::string& name) const;

    std::vector<std::string> getDescriptors() const;

    // --- Observer Methods ---
    void registerObserver(std::shared_ptr<Observer> observer);

    std::shared_ptr<Observer> getObserver(const std::string& name) const;

    const std::unordered_map<std::string, std::shared_ptr<Observer>>& getObservers() const;
};

// --- Helper Function Declarations ---
int countMatches(const std::string& smiles, const re2::RE2& pattern);
int countChar(const std::string& s, char c);
int countDigits(const std::string& s);
int countSubstr(const std::string& text, const std::string& pattern);

// Observer Registration Function Declarations (Optional, depends on usage)
// void registerElementCountObserver();
// void registerRingInfoObserver();
// void registerElectronegativityObserver();

} // namespace desfact