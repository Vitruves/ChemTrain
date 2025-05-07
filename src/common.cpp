#include "common.hpp"
#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <queue>
#include <cmath>
#include <string>
#include <stdexcept>
#include <vector>
#include <map>
#include <memory> // For std::shared_ptr
#include <GraphMol/MolOps.h> // Ensure MolOps is included
#include <RDGeneral/types.h> // For RDKit::MAX_DIST_VAL

namespace desfact {

//==============================================================================
// MOLECULE TRAVERSER IMPLEMENTATION
//==============================================================================
void MoleculeTraverser::addObserver(std::shared_ptr<Observer> observer) {
    bool found = false;
    for(const auto& obs : observers_) {
        if (obs->getName() == observer->getName()) {
            found = true;
            break;
        }
    }
    if (!found) {
        observers_.push_back(observer);
    }
}

void MoleculeTraverser::traverse(const RDKit::ROMol* mol, Context& context) {
    if (!mol) return;
    
    context.setMolecule(mol);
    resolveObserverDependencies();
    
    for (auto& observer : observers_) {
        observer->init(context);
    }
    
    for (const RDKit::Atom* atom : mol->atoms()) {
        for (auto& observer : observers_) {
            observer->observeAtom(atom, context);
        }
    }
    
    for (const RDKit::Bond* bond : mol->bonds()) {
        for (auto& observer : observers_) {
            observer->observeBond(bond, context);
        }
    }
    
    for (auto& observer : observers_) {
        observer->finalize(context);
    }
}

void MoleculeTraverser::reset() {
    observers_.clear();
}

void MoleculeTraverser::resolveObserverDependencies() {
    if (observers_.size() < 2) return;

    std::vector<std::shared_ptr<Observer>> sortedObservers;
    std::unordered_map<std::string, std::shared_ptr<Observer>> observerMap;
    std::unordered_map<std::string, std::unordered_set<std::string>> adj;
    std::unordered_map<std::string, int> inDegree;

    for (const auto& obs : observers_) {
        observerMap[obs->getName()] = obs;
        inDegree[obs->getName()] = 0;
    }

    for (const auto& obs : observers_) {
        std::string dependentName = obs->getName();
        for (const auto& depName : obs->getDependencies()) {
            if (observerMap.count(depName)) {
                adj[depName].insert(dependentName);
                inDegree[dependentName]++;
            }
        }
    }

    std::queue<std::string> q;
    for (const auto& pair : inDegree) {
        if (pair.second == 0) {
            q.push(pair.first);
        }
    }

    while (!q.empty()) {
        std::string u = q.front();
        q.pop();
        sortedObservers.push_back(observerMap[u]);

        if (adj.count(u)) {
            for (const std::string& v : adj[u]) {
                if (--inDegree[v] == 0) {
                    q.push(v);
                }
            }
        }
    }

    if (sortedObservers.size() != observers_.size()) {
        return;
    }

    observers_ = std::move(sortedObservers);
}

//==============================================================================
// DESCRIPTOR PIPELINE IMPLEMENTATION
//==============================================================================
void DescriptorPipeline::addDescriptor(const std::string& name) {
    auto& registry = DescriptorRegistry::getInstance();
    auto descriptor = registry.getDescriptor(name);
    if (!descriptor) {
        throw std::runtime_error("Descriptor not found: " + name);
    }
    bool found = false;
    for(const auto& desc : descriptors_) {
        if (desc->getName() == descriptor->getName()) {
            found = true;
            break;
        }
    }
    if (!found) {
        descriptors_.push_back(descriptor);
    }
}

std::map<std::string, DescriptorResult> DescriptorPipeline::process(const RDKit::ROMol* mol) {
    std::map<std::string, DescriptorResult> results;
    if (!mol) return results;
    
    Context context;
    resolveDependencies();
    traverser_.traverse(mol, context);
    
    for (auto& descriptor : descriptors_) {
        try {
            results[descriptor->getName()] = descriptor->calculate(context);
        } catch (const std::exception& e) {
            results[descriptor->getName()] = "Error: " + std::string(e.what());
        }
    }
    
    return results;
}

std::vector<std::string> DescriptorPipeline::getDescriptorNames() const {
    std::vector<std::string> names;
    names.reserve(descriptors_.size());
    for (const auto& desc : descriptors_) {
        names.push_back(desc->getName());
    }
    return names;
}

void DescriptorPipeline::resolveDependencies() {
    traverser_.reset();
    auto& registry = DescriptorRegistry::getInstance();

    std::unordered_set<std::string> requiredObserverNames;
    std::vector<std::string> q = {};

    for(const auto& desc : descriptors_) {
        for(const auto& depName : desc->getDependencies()) {
            if(registry.getObserver(depName)) {
                 if (requiredObserverNames.find(depName) == requiredObserverNames.end()) {
                    requiredObserverNames.insert(depName);
                    q.push_back(depName);
                 }
            }
        }
    }

    size_t head = 0;
    while(head < q.size()) {
        std::string currentObserverName = q[head++];
        auto observer = registry.getObserver(currentObserverName);
        if(observer) {
            for(const auto& subDepName : observer->getDependencies()) {
                if(registry.getObserver(subDepName)) {
                    if (requiredObserverNames.find(subDepName) == requiredObserverNames.end()) {
                       requiredObserverNames.insert(subDepName);
                       q.push_back(subDepName);
                    }
                }
            }
        }
    }

    for (const auto& name : requiredObserverNames) {
        auto observer = registry.getObserver(name);
        if (observer) {
            traverser_.addObserver(observer);
        } else {
             std::cerr << "Warning: Could not find required observer dependency: " << name << std::endl;
        }
    }
}

//==============================================================================
// ELEMENT PROPERTIES IMPLEMENTATION
//==============================================================================
namespace element {

// Constants for all elements (H-Sn, atomic numbers 1-50)
// We use unordered_map for faster lookups
static const std::unordered_map<int, double> atomicMass = {
    {1, 1.6735e-27}, {2, 6.6465e-27}, {3, 1.1526e-26}, {4, 1.4966e-26}, {5, 1.796e-26},
    {6, 1.9944e-26}, {7, 2.3259e-26}, {8, 2.656e-26}, {9, 3.1548e-26}, {10, 3.3509e-26},
    {11, 3.8175e-26}, {12, 4.0438e-26}, {13, 4.4804e-26}, {14, 4.6637e-26}, {15, 5.1441e-26},
    {16, 5.3126e-26}, {17, 5.887e-26}, {18, 6.6335e-26}, {19, 6.4924e-26}, {20, 6.656e-26},
    {21, 7.946e-26}, {22, 7.946e-26}, {23, 8.493e-26}, {24, 8.634e-26}, {25, 9.127e-26},
    {26, 9.288e-26}, {27, 9.787e-26}, {28, 9.744e-26}, {29, 1.055e-25}, {30, 1.088e-25},
    {31, 1.151e-25}, {32, 1.206e-25}, {33, 1.244e-25}, {34, 1.31e-25}, {35, 1.328e-25},
    {36, 1.392e-25}, {37, 1.419e-25}, {38, 1.487e-25}, {39, 1.51e-25}, {40, 1.644e-25},
    {41, 1.745e-25}, {42, 1.593e-25}, {43, 1.811e-25}, {44, 1.774e-25}, {45, 1.867e-25},
    {46, 1.77e-25}, {47, 1.791e-25}, {48, 1.865e-25}, {49, 2.1e-25}, {50, 1.97e-25}
};

// Pauling electronegativity values
static const std::unordered_map<int, double> paulingEN = {
    {1, 2.20}, {3, 0.98}, {4, 1.57}, {5, 2.04}, {6, 2.55}, {7, 3.04}, {8, 3.44}, {9, 3.98},
    {11, 0.93}, {12, 1.31}, {13, 1.61}, {14, 1.90}, {15, 2.19}, {16, 2.58}, {17, 3.16},
    {19, 0.82}, {20, 1.00}, {21, 1.36}, {22, 1.54}, {23, 1.63}, {24, 1.66}, {25, 1.55},
    {26, 1.83}, {27, 1.88}, {28, 1.91}, {29, 1.90}, {30, 1.65}, {31, 1.81}, {32, 2.01},
    {33, 2.18}, {34, 2.55}, {35, 2.96}, {37, 0.82}, {38, 0.95}, {39, 1.22}, {40, 1.33},
    {41, 1.60}, {42, 2.16}, {43, 1.90}, {44, 2.20}, {45, 2.28}, {46, 2.20}, {47, 1.93},
    {48, 1.69}, {49, 1.78}, {50, 1.96}
};

// Allred-Rochow electronegativity values
static const std::unordered_map<int, double> allredRochowEN = {
    {1, 2.20}, {3, 0.97}, {4, 1.47}, {5, 2.01}, {6, 2.50}, {7, 3.07}, {8, 3.50}, {9, 4.10},
    {11, 1.01}, {12, 1.23}, {13, 1.47}, {14, 1.74}, {15, 2.06}, {16, 2.44}, {17, 2.83}, 
    {19, 0.91}, {20, 1.00}, {21, 1.36}, {22, 1.32}, {23, 1.45}, {24, 1.56}, {25, 1.60},
    {26, 1.64}, {27, 1.70}, {28, 1.75}, {29, 1.75}, {30, 1.65}, {31, 1.81}, {32, 2.02},
    {33, 2.20}, {34, 2.48}, {35, 2.74}, {37, 0.89}, {38, 0.96}, {39, 1.22}, {40, 1.32},
    {41, 1.60}, {42, 2.16}, {43, 1.90}, {44, 2.20}, {45, 2.28}, {46, 2.20}, {47, 1.93},
    {48, 1.69}, {49, 1.78}, {50, 1.96}
};

// Ionization energy (J)
static const std::unordered_map<int, double> ionizationEnergy = {
    {1, 2.179e-18}, {2, 3.94e-18}, {3, 5.39e-19}, {4, 9.32e-19}, {5, 8.3e-19},
    {6, 1.805e-18}, {7, 2.327e-18}, {8, 2.183e-18}, {9, 2.791e-18}, {10, 3.09e-18},
    {11, 8.236e-19}, {12, 1.226e-18}, {13, 9.8e-19}, {14, 1.079e-18}, {15, 1.681e-18},
    {16, 1.661e-18}, {17, 2.078e-18}, {18, 2.273e-18}, {19, 6.98e-19}, {20, 9.796e-19},
    {21, 6.56e-19}, {22, 6.82e-19}, {23, 6.74e-19}, {24, 6.77e-19}, {25, 7.43e-19},
    {26, 7.9e-19}, {27, 7.86e-19}, {28, 7.64e-19}, {29, 7.73e-19}, {30, 9.39e-19},
    {31, 6e-19}, {32, 7.9e-19}, {33, 9.81e-19}, {34, 9.75e-19}, {35, 1.18e-18},
    {36, 1.35e-18}, {37, 4.18e-19}, {38, 5.7e-19}, {39, 6.17e-19}, {40, 6.33e-19},
    {41, 6.76e-19}, {42, 6.88e-19}, {43, 7.28e-19}, {44, 7.36e-19}, {45, 7.46e-19},
    {46, 8.34e-19}, {47, 7.58e-19}, {48, 8.99e-19}, {49, 5.79e-19}, {50, 7.34e-19}
};

// Electron affinity (J)
static const std::unordered_map<int, double> electronAffinity = {
    {1, 1.209e-19}, {3, 6.47e-20}, {5, 5.84e-20}, {6, 2.024e-19}, {7, -1.12e-20},
    {8, 2.341e-19}, {9, 5.451e-19}, {11, 8.78e-20}, {12, 0.0}, {13, 7.08e-20},
    {14, 2.23e-19}, {15, 1.197e-19}, {16, 2.077e-19}, {17, 3.612e-19}, {19, 8.03e-20},
    {20, 2.39e-20}, {21, 1.92e-20}, {22, 3.69e-20}, {23, 5.15e-20}, {24, 6.46e-20},
    {25, 1.24e-19}, {26, 2.42e-20}, {27, 6.3e-20}, {28, 1.15e-19}, {29, 1.24e-19},
    {30, 0.0}, {31, 3.02e-20}, {32, 1.2e-19}, {33, 7.77e-20}, {34, 2.02e-19},
    {35, 3.365e-19}, {37, 4e-20}, {38, 4.4e-20}, {39, 3.5e-20}, {40, 4.15e-20},
    {41, 8.91e-20}, {42, 7.38e-20}, {43, 5.49e-20}, {44, 1.05e-19}, {45, 1.1e-19},
    {46, 5.6e-20}, {47, 1.3e-19}, {48, 0.0}, {49, 3e-20}, {50, 1.17e-19}
};

// Covalent radius (m)
static const std::unordered_map<int, double> covalentRadius = {
    {1, 3.1e-11}, {3, 1.28e-10}, {4, 9.6e-11}, {5, 8.4e-11}, {6, 7.6e-11},
    {7, 7.1e-11}, {8, 6.6e-11}, {9, 5.7e-11}, {11, 1.66e-10}, {12, 1.41e-10},
    {13, 1.21e-10}, {14, 1.11e-10}, {15, 1.07e-10}, {16, 1.05e-10}, {17, 1.02e-10},
    {19, 2.03e-10}, {20, 1.76e-10}, {21, 1.7e-10}, {22, 1.6e-10}, {23, 1.53e-10},
    {24, 1.39e-10}, {25, 1.39e-10}, {26, 1.32e-10}, {27, 1.26e-10}, {28, 1.24e-10},
    {29, 1.32e-10}, {30, 1.22e-10}, {31, 1.22e-10}, {32, 1.2e-10}, {33, 1.19e-10},
    {34, 1.2e-10}, {35, 1.2e-10}, {37, 2.2e-10}, {38, 1.95e-10}, {39, 1.8e-10},
    {40, 1.75e-10}, {41, 1.64e-10}, {42, 1.54e-10}, {43, 1.47e-10}, {44, 1.46e-10},
    {45, 1.42e-10}, {46, 1.39e-10}, {47, 1.45e-10}, {48, 1.44e-10}, {49, 1.42e-10},
    {50, 1.39e-10}
};

// Van der Waals radius (m)
static const std::unordered_map<int, double> vanDerWaalsRadius = {
    {1, 1.2e-10}, {2, 1.4e-10}, {3, 1.82e-10}, {6, 1.7e-10}, {7, 1.55e-10},
    {8, 1.52e-10}, {9, 1.47e-10}, {10, 1.54e-10}, {11, 2.27e-10}, {12, 1.73e-10},
    {13, 1.84e-10}, {14, 2.1e-10}, {15, 1.8e-10}, {16, 1.8e-10}, {17, 1.75e-10},
    {18, 1.88e-10}, {19, 2.75e-10}, {20, 2.31e-10}, {21, 2.11e-10}, {22, 2e-10},
    {23, 2e-10}, {24, 2e-10}, {25, 2e-10}, {26, 1.94e-10}, {27, 2e-10},
    {28, 1.63e-10}, {29, 1.4e-10}, {30, 1.39e-10}, {31, 1.87e-10}, {32, 2.11e-10},
    {33, 1.85e-10}, {34, 1.9e-10}, {35, 1.85e-10}, {36, 2.02e-10}, {37, 2.98e-10},
    {38, 2.49e-10}, {39, 2.27e-10}, {40, 2.16e-10}, {41, 2.08e-10}, {42, 2.01e-10},
    {43, 2e-10}, {44, 2e-10}, {45, 2e-10}, {46, 1.63e-10}, {47, 1.72e-10},
    {48, 1.58e-10}, {49, 1.93e-10}, {50, 2.17e-10}
};

// Polarizability (Å³)
static const std::unordered_map<int, double> polarizability = {
    {1, 0.666}, {2, 0.204}, {3, 24.3}, {4, 5.6}, {5, 3.03}, {6, 1.76}, {7, 1.1},
    {8, 0.802}, {9, 0.557}, {10, 0.395}, {11, 24.11}, {12, 10.6}, {13, 6.8},
    {14, 5.38}, {15, 7.43}, {16, 2.9}, {17, 2.18}, {18, 1.64}, {19, 43.4},
    {20, 22.8}, {21, 18.1}, {22, 10.6}, {23, 8.4}, {24, 6.2}, {25, 7.4},
    {26, 7.0}, {27, 6.4}, {28, 6.2}, {29, 6.4}, {30, 6.2}, {31, 7.7},
    {32, 6.1}, {33, 4.4}, {34, 3.8}, {35, 3.05}, {36, 2.48}, {37, 47.3},
    {38, 28.9}, {39, 11.6}, {40, 9.5}, {41, 7.3}, {42, 5.7}, {43, 5.3},
    {44, 5.3}, {45, 6.1}, {46, 6.1}, {47, 7.2}, {48, 7.1}, {49, 8.4}, {50, 7.4}
};

// Density (kg/m³)
static const std::unordered_map<int, double> density = {
    {1, 0.08988}, {2, 0.1786}, {3, 534}, {4, 1848}, {5, 2460}, {6, 2267},
    {7, 1.2506}, {8, 1.429}, {9, 1.696}, {10, 0.9002}, {11, 968}, {12, 1738},
    {13, 2700}, {14, 2330}, {15, 1823}, {16, 2070}, {17, 3.2}, {18, 1.784},
    {19, 862}, {20, 1550}, {21, 2989}, {22, 4507}, {23, 6110}, {24, 7190},
    {25, 7470}, {26, 7874}, {27, 8900}, {28, 8908}, {29, 8960}, {30, 7140},
    {31, 5907}, {32, 5323}, {33, 5727}, {34, 4819}, {35, 3120}, {36, 3.75},
    {37, 1532}, {38, 2630}, {39, 4472}, {40, 6511}, {41, 8570}, {42, 10280},
    {43, 11500}, {44, 12370}, {45, 12410}, {46, 12020}, {47, 10490}, {48, 8650},
    {49, 7310}, {50, 7287}
};

// Thermal conductivity (W/m·K)
static const std::unordered_map<int, double> thermalConductivity = {
    {1, 0.1815}, {2, 0.1513}, {3, 84.7}, {4, 200}, {5, 27}, {6, 140},
    {7, 0.026}, {8, 0.0266}, {9, 0.0277}, {10, 0.0491}, {11, 142}, {12, 156},
    {13, 235}, {14, 148}, {15, 0.236}, {16, 0.269}, {17, 0.0098}, {18, 0.0177},
    {19, 102}, {20, 201}, {21, 15.8}, {22, 21.9}, {23, 30.7}, {24, 93.7},
    {25, 7.8}, {26, 80.4}, {27, 100}, {28, 90.9}, {29, 401}, {30, 116},
    {31, 29}, {32, 60.2}, {33, 50}, {34, 0.52}, {35, 0.12}, {36, 0.0095},
    {37, 58}, {38, 35}, {39, 17}, {40, 22.6}, {41, 53.7}, {42, 138},
    {43, 50.6}, {44, 117}, {45, 150}, {46, 71.8}, {47, 429}, {48, 96.6},
    {49, 81.8}, {50, 66.8}
};

// Specific heat (J/g·K)
static const std::unordered_map<int, double> specificHeat = {
    {1, 14.304}, {2, 5.193}, {3, 3.582}, {4, 1.825}, {5, 1.026}, {6, 0.709},
    {7, 1.04}, {8, 0.918}, {9, 0.824}, {10, 1.03}, {11, 1.228}, {12, 1.023},
    {13, 0.897}, {14, 0.705}, {15, 0.769}, {16, 0.71}, {17, 0.479}, {18, 0.52},
    {19, 0.753}, {20, 0.647}, {21, 0.568}, {22, 0.523}, {23, 0.489}, {24, 0.449},
    {25, 0.479}, {26, 0.449}, {27, 0.421}, {28, 0.444}, {29, 0.385}, {30, 0.388},
    {31, 0.371}, {32, 0.32}, {33, 0.329}, {34, 0.321}, {35, 0.474}, {36, 0.248},
    {37, 0.363}, {38, 0.301}, {39, 0.298}, {40, 0.278}, {41, 0.265}, {42, 0.251},
    {43, 0.21}, {44, 0.238}, {45, 0.243}, {46, 0.244}, {47, 0.235}, {48, 0.232},
    {49, 0.233}, {50, 0.228}
};

// Molar volume (cm³/mol)
static const std::unordered_map<int, double> molarVolume = {
    {1, 14.1}, {2, 27.2}, {3, 13.1}, {4, 4.9}, {5, 4.6}, {6, 5.3},
    {7, 17.3}, {8, 14}, {9, 17.1}, {10, 16.5}, {11, 23.7}, {12, 13.97},
    {13, 10}, {14, 12.1}, {15, 17}, {16, 15.5}, {17, 17.9}, {18, 22.4},
    {19, 45.3}, {20, 29.9}, {21, 15}, {22, 10.6}, {23, 8.78}, {24, 7.23},
    {25, 7.39}, {26, 7.09}, {27, 6.7}, {28, 6.59}, {29, 7.11}, {30, 9.16},
    {31, 11.8}, {32, 13.6}, {33, 13.1}, {34, 16.5}, {35, 23.5}, {36, 38.9},
    {37, 55.9}, {38, 33.4}, {39, 19.9}, {40, 14.1}, {41, 10.9}, {42, 9.4},
    {43, 8.5}, {44, 8.3}, {45, 8.2}, {46, 8.9}, {47, 10.3}, {48, 13.1},
    {49, 15.7}, {50, 16.3}
};

// Replace existing getElectronegativity function with a more complete implementation
double getElectronegativity(int atomicNum) {
    auto it = paulingEN.find(atomicNum);
    return (it != paulingEN.end()) ? it->second : 0.0;
}

// Get Allred-Rochow electronegativity
double getAllredRochowEN(int atomicNum) {
    auto it = allredRochowEN.find(atomicNum);
    return (it != allredRochowEN.end()) ? it->second : 0.0;
}

// Get atomic mass in kg
double getAtomicMass(int atomicNum) {
    auto it = atomicMass.find(atomicNum);
    return (it != atomicMass.end()) ? it->second : 0.0;
}

// Get ionization energy in J
double getIonizationEnergy(int atomicNum) {
    auto it = ionizationEnergy.find(atomicNum);
    return (it != ionizationEnergy.end()) ? it->second : 0.0;
}

// Get electron affinity in J
double getElectronAffinity(int atomicNum) {
    auto it = electronAffinity.find(atomicNum);
    return (it != electronAffinity.end()) ? it->second : 0.0;
}

// Get covalent radius in m
double getCovalentRadius(int atomicNum) {
    auto it = covalentRadius.find(atomicNum);
    return (it != covalentRadius.end()) ? it->second : 0.0;
}

// Get van der Waals radius in m
double getVanDerWaalsRadius(int atomicNum) {
    auto it = vanDerWaalsRadius.find(atomicNum);
    return (it != vanDerWaalsRadius.end()) ? it->second : 0.0;
}

// Get polarizability in Å³
double getPolarizability(int atomicNum) {
    auto it = polarizability.find(atomicNum);
    return (it != polarizability.end()) ? it->second : 0.0;
}

// Get density in kg/m³
double getDensity(int atomicNum) {
    auto it = density.find(atomicNum);
    return (it != density.end()) ? it->second : 0.0;
}

// Get thermal conductivity in W/m·K
double getThermalConductivity(int atomicNum) {
    auto it = thermalConductivity.find(atomicNum);
    return (it != thermalConductivity.end()) ? it->second : 0.0;
}

// Get specific heat in J/g·K
double getSpecificHeat(int atomicNum) {
    auto it = specificHeat.find(atomicNum);
    return (it != specificHeat.end()) ? it->second : 0.0;
}

// Get molar volume in cm³/mol
double getMolarVolume(int atomicNum) {
    auto it = molarVolume.find(atomicNum);
    return (it != molarVolume.end()) ? it->second : 0.0;
}

// Get atomic number from element symbol
int getAtomicNumber(const std::string& symbol) {
    static const std::unordered_map<std::string, int> symbolToZ = {
        {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
        {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18},
        {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26},
        {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30}, {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34},
        {"Br", 35}, {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40}, {"Nb", 41}, {"Mo", 42},
        {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50}
    };
    
    auto it = symbolToZ.find(symbol);
    return (it != symbolToZ.end()) ? it->second : 0;
}

// Get element symbol from atomic number
std::string getElementSymbol(int atomicNum) {
    static const std::unordered_map<int, std::string> zToSymbol = {
        {1, "H"}, {2, "He"}, {3, "Li"}, {4, "Be"}, {5, "B"}, {6, "C"}, {7, "N"}, {8, "O"}, {9, "F"}, {10, "Ne"},
        {11, "Na"}, {12, "Mg"}, {13, "Al"}, {14, "Si"}, {15, "P"}, {16, "S"}, {17, "Cl"}, {18, "Ar"},
        {19, "K"}, {20, "Ca"}, {21, "Sc"}, {22, "Ti"}, {23, "V"}, {24, "Cr"}, {25, "Mn"}, {26, "Fe"},
        {27, "Co"}, {28, "Ni"}, {29, "Cu"}, {30, "Zn"}, {31, "Ga"}, {32, "Ge"}, {33, "As"}, {34, "Se"},
        {35, "Br"}, {36, "Kr"}, {37, "Rb"}, {38, "Sr"}, {39, "Y"}, {40, "Zr"}, {41, "Nb"}, {42, "Mo"},
        {43, "Tc"}, {44, "Ru"}, {45, "Rh"}, {46, "Pd"}, {47, "Ag"}, {48, "Cd"}, {49, "In"}, {50, "Sn"}
    };
    
    auto it = zToSymbol.find(atomicNum);
    return (it != zToSymbol.end()) ? it->second : "";
}

// Calculate periodic table group for an element
int getPeriodicGroup(int atomicNum) {
    static const std::unordered_map<int, int> elementGroups = {
        {1, 1}, {2, 18}, {3, 1}, {4, 2}, {5, 13}, {6, 14}, {7, 15}, {8, 16}, {9, 17}, {10, 18},
        {11, 1}, {12, 2}, {13, 13}, {14, 14}, {15, 15}, {16, 16}, {17, 17}, {18, 18},
        {19, 1}, {20, 2}, {21, 3}, {22, 4}, {23, 5}, {24, 6}, {25, 7}, {26, 8},
        {27, 9}, {28, 10}, {29, 11}, {30, 12}, {31, 13}, {32, 14}, {33, 15}, {34, 16},
        {35, 17}, {36, 18}, {37, 1}, {38, 2}, {39, 3}, {40, 4}, {41, 5}, {42, 6},
        {43, 7}, {44, 8}, {45, 9}, {46, 10}, {47, 11}, {48, 12}, {49, 13}, {50, 14}
    };
    
    auto it = elementGroups.find(atomicNum);
    return (it != elementGroups.end()) ? it->second : 0;
}

// Calculate periodic table period for an element
int getPeriodicPeriod(int atomicNum) {
    static const std::unordered_map<int, int> elementPeriods = {
        {1, 1}, {2, 1}, 
        {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}, {9, 2}, {10, 2},
        {11, 3}, {12, 3}, {13, 3}, {14, 3}, {15, 3}, {16, 3}, {17, 3}, {18, 3},
        {19, 4}, {20, 4}, {21, 4}, {22, 4}, {23, 4}, {24, 4}, {25, 4}, {26, 4},
        {27, 4}, {28, 4}, {29, 4}, {30, 4}, {31, 4}, {32, 4}, {33, 4}, {34, 4},
        {35, 4}, {36, 4}, {37, 5}, {38, 5}, {39, 5}, {40, 5}, {41, 5}, {42, 5},
        {43, 5}, {44, 5}, {45, 5}, {46, 5}, {47, 5}, {48, 5}, {49, 5}, {50, 5}
    };
    
    auto it = elementPeriods.find(atomicNum);
    return (it != elementPeriods.end()) ? it->second : 0;
}

// Functions to determine element classification
bool isMainGroupElement(int atomicNum) {
    int group = getPeriodicGroup(atomicNum);
    return (group >= 1 && group <= 2) || (group >= 13 && group <= 18);
}

bool isTransitionMetal(int atomicNum) {
    int group = getPeriodicGroup(atomicNum);
    int period = getPeriodicPeriod(atomicNum);
    return (period >= 4) && (group >= 3 && group <= 12);
}

bool isAlkaliMetal(int atomicNum) {
    return getPeriodicGroup(atomicNum) == 1 && atomicNum > 2;
}

bool isAlkalineEarthMetal(int atomicNum) {
    return getPeriodicGroup(atomicNum) == 2 && atomicNum > 2;
}

bool isNobleGas(int atomicNum) {
    return getPeriodicGroup(atomicNum) == 18 && atomicNum > 1;
}

// Check if an element is a halogen (group 17)
bool isHalogen(int atomicNum) {
    // Halogens: F(9), Cl(17), Br(35), I(53), At(85)
    return atomicNum == 9 || atomicNum == 17 || atomicNum == 35 || 
           atomicNum == 53 || atomicNum == 85;
}

// Check if an element is a metal
bool isMetal(int atomicNum) {
    // Everything except for these elements is considered a metal
    // Non-metals: H(1), He(2), B(5), C(6), N(7), O(8), F(9), Ne(10),
    //             Si(14), P(15), S(16), Cl(17), Ar(18), 
    //             Ge(32) [metalloid], As(33) [metalloid], Se(34), Br(35), Kr(36),
    //             Sb(51) [metalloid], Te(52) [metalloid], I(53), Xe(54)
    static const std::unordered_set<int> nonMetals = {
        1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 
        32, 33, 34, 35, 36, 51, 52, 53, 54
    };
    
    return nonMetals.find(atomicNum) == nonMetals.end() && atomicNum > 0;
}

} // namespace element

//==============================================================================
// MOLECULAR UTILITIES IMPLEMENTATION
//==============================================================================
namespace mol {

void ensureRingInfo(const RDKit::ROMol* mol) {
    if (!mol) return;
    
    if (!mol->getRingInfo()->isInitialized()) {
        RDKit::MolOps::findSSSR(*const_cast<RDKit::ROMol*>(mol));
    }
}

int countAtoms(const RDKit::ROMol* mol, std::function<bool(const RDKit::Atom*)> predicate) {
    if (!mol) return 0;
    
    int count = 0;
    for (const RDKit::Atom* atom : mol->atoms()) {
        if (predicate(atom)) {
            count++;
        }
    }
    return count;
}

int countBonds(const RDKit::ROMol* mol, std::function<bool(const RDKit::Bond*)> predicate) {
    if (!mol) return 0;
    
    int count = 0;
    for (const RDKit::Bond* bond : mol->bonds()) {
        if (predicate(bond)) {
            count++;
        }
    }
    return count;
}

double atomFraction(const RDKit::ROMol* mol, std::function<bool(const RDKit::Atom*)> predicate) {
    if (!mol) return 0.0;
    int total = mol->getNumAtoms();
    if (total == 0) return 0.0;
    
    int matching = countAtoms(mol, predicate);
    return static_cast<double>(matching) / total;
}

double bondFraction(const RDKit::ROMol* mol, std::function<bool(const RDKit::Bond*)> predicate) {
    if (!mol) return 0.0;
    int total = mol->getNumBonds();
    if (total == 0) return 0.0;
    
    int matching = countBonds(mol, predicate);
    return static_cast<double>(matching) / total;
}

// Bond type utilities
bool isSingleBond(const RDKit::Bond* bond) {
    return bond && bond->getBondType() == RDKit::Bond::SINGLE;
}

bool isDoubleBond(const RDKit::Bond* bond) {
    return bond && bond->getBondType() == RDKit::Bond::DOUBLE;
}

bool isTripleBond(const RDKit::Bond* bond) {
    return bond && bond->getBondType() == RDKit::Bond::TRIPLE;
}

bool isAromaticBond(const RDKit::Bond* bond) {
    return bond && bond->getBondType() == RDKit::Bond::AROMATIC;
}

bool isRotatableBond(const RDKit::Bond* bond) {
    if (!bond || bond->getBondType() != RDKit::Bond::SINGLE) return false;
    if (bond->getIsAromatic() || bond->getIsConjugated()) return false;
    
    const RDKit::ROMol* mol = &bond->getOwningMol();
    if (mol->getRingInfo()->numBondRings(bond->getIdx()) > 0) return false;
    
    const RDKit::Atom* a1 = bond->getBeginAtom();
    const RDKit::Atom* a2 = bond->getEndAtom();
    
    // Skip terminal bonds
    if (a1->getDegree() == 1 || a2->getDegree() == 1) return false;
    
    // Skip C-H bonds
    if ((a1->getAtomicNum() == 1 || a2->getAtomicNum() == 1) &&
        (a1->getAtomicNum() == 6 || a2->getAtomicNum() == 6)) return false;
    
    return true;
}

// Atom type utilities
bool isHeteroatom(const RDKit::Atom* atom) {
    return atom && atom->getAtomicNum() != 6 && atom->getAtomicNum() != 1;
}

bool isHalogen(const RDKit::Atom* atom) {
    return atom && element::isHalogen(atom->getAtomicNum());
}

bool isMetal(const RDKit::Atom* atom) {
    return atom && element::isMetal(atom->getAtomicNum());
}

bool isInRing(const RDKit::Atom* atom) {
    return atom && atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx()) > 0;
}

bool isAromatic(const RDKit::Atom* atom) {
    return atom && atom->getIsAromatic();
}

bool isBasicNitrogen(const RDKit::Atom* atom) {
    if (!atom || atom->getAtomicNum() != 7) return false;
    
    // Check for aromatic nitrogen (usually not basic)
    if (atom->getIsAromatic()) return false;
    
    // Check for basic nitrogen configurations
    // This is a simplified model, real pKa prediction is complex
    int valence = atom->getTotalValence();
    int formalCharge = atom->getFormalCharge();
    
    // Typical basic nitrogen has 3 valence (like in amines)
    if (valence == 3 && formalCharge == 0) return true;
    
    // Positively charged nitrogens are not basic
    if (formalCharge > 0) return false;
    
    return false; // Default case
}

bool isAcidicOxygen(const RDKit::Atom* atom) {
    if (!atom || atom->getAtomicNum() != 8) return false;
    
    // Check for carboxylic acid oxygens
    if (atom->getDegree() == 1) {
        // Get the atom's neighbors using the correct RDKit API
        for (const auto& bondIdx : boost::make_iterator_range(atom->getOwningMol().getAtomBonds(atom))) {
            const RDKit::Bond* bond = atom->getOwningMol()[bondIdx];
            if (bond->getBondType() == RDKit::Bond::DOUBLE) {
                // Check if it's connected to a carbon with another oxygen
                const RDKit::Atom* neighbor = bond->getOtherAtom(atom);
                if (neighbor->getAtomicNum() == 6) {
                    for (const auto& nbri : boost::make_iterator_range(neighbor->getOwningMol().getAtomNeighbors(neighbor))) {
                        const RDKit::Atom* nb = neighbor->getOwningMol()[nbri];
                        if (nb->getAtomicNum() == 8 && nb != atom) {
                            return true;
                        }
                    }
                }
            }
        }
    }
    
    return false;
}

// Distance and connectivity utilities
int getAtomDistance(const RDKit::ROMol* mol, unsigned int idx1, unsigned int idx2) {
    if (!mol) return -1;
    
    if (idx1 >= mol->getNumAtoms() || idx2 >= mol->getNumAtoms()) return -1;
    
    // Use the function that returns the path length directly
    auto path = RDKit::MolOps::getShortestPath(*mol, idx1, idx2);
    // Return the size of the path minus 1 (which is the number of bonds traversed)
    return path.empty() ? -1 : static_cast<int>(path.size() - 1);
}

int getBondCount(const RDKit::ROMol* mol) {
    return mol ? mol->getNumBonds() : 0;
}

int getHeavyAtomCount(const RDKit::ROMol* mol) {
    return mol ? mol->getNumHeavyAtoms() : 0;
}

int getRotatableBondCount(const RDKit::ROMol* mol) {
    if (!mol) return 0;
    
    int count = 0;
    for (const RDKit::Bond* bond : mol->bonds()) {
        if (isRotatableBond(bond)) {
            count++;
        }
    }
    return count;
}

// Ring utilities
int getRingCount(const RDKit::ROMol* mol) {
    if (!mol) return 0;
    ensureRingInfo(mol);
    return mol->getRingInfo()->numRings();
}

int getAromaticRingCount(const RDKit::ROMol* mol) {
    if (!mol) return 0;
    
    ensureRingInfo(mol);
    int count = 0;
    
    for (const auto& ring : mol->getRingInfo()->atomRings()) {
        bool isAromatic = true;
        for (auto atomIdx : ring) {
            if (!mol->getAtomWithIdx(atomIdx)->getIsAromatic()) {
                isAromatic = false;
                break;
            }
        }
        if (isAromatic) count++;
    }
    
    return count;
}

// Descriptive state utilities
bool hasStereoChemistry(const RDKit::ROMol* mol) {
    if (!mol) return false;
    
    for (const RDKit::Atom* atom : mol->atoms()) {
        if (atom->getChiralTag() != RDKit::Atom::CHI_UNSPECIFIED) {
            return true;
        }
    }
    
    for (const RDKit::Bond* bond : mol->bonds()) {
        if (bond->getStereo() != RDKit::Bond::STEREONONE) {
            return true;
        }
    }
    
    return false;
}

// Statistical utilities for atom properties
double getAverageAtomicMass(const RDKit::ROMol* mol) {
    if (!mol || mol->getNumAtoms() == 0) return 0.0;
    
    double sum = 0.0;
    for (const RDKit::Atom* atom : mol->atoms()) {
        sum += element::getAtomicMass(atom->getAtomicNum()) * 6.022e23; // Convert to g/mol
    }
    return sum / mol->getNumAtoms();
}

double getAverageElectronegativity(const RDKit::ROMol* mol) {
    if (!mol || mol->getNumAtoms() == 0) return 0.0;
    
    double sum = 0.0;
    int count = 0;
    for (const RDKit::Atom* atom : mol->atoms()) {
        double en = element::getElectronegativity(atom->getAtomicNum());
        if (en > 0) { // Skip atoms with unknown electronegativity
            sum += en;
            count++;
        }
    }
    return count > 0 ? sum / count : 0.0;
}

// Electronegativity difference between connected atoms
double getMaxElectronegativityDifference(const RDKit::ROMol* mol) {
    if (!mol || mol->getNumBonds() == 0) return 0.0;
    
    double maxDiff = 0.0;
    for (const RDKit::Bond* bond : mol->bonds()) {
        const RDKit::Atom* a1 = bond->getBeginAtom();
        const RDKit::Atom* a2 = bond->getEndAtom();
        double en1 = element::getElectronegativity(a1->getAtomicNum());
        double en2 = element::getElectronegativity(a2->getAtomicNum());
        double diff = std::abs(en1 - en2);
        maxDiff = std::max(maxDiff, diff);
    }
    return maxDiff;
}

// Calculate sum of a specific atom property
double sumAtomProperty(const RDKit::ROMol* mol, std::function<double(const RDKit::Atom*)> propertyFn) {
    if (!mol) return 0.0;
    
    double sum = 0.0;
    for (const RDKit::Atom* atom : mol->atoms()) {
        sum += propertyFn(atom);
    }
    return sum;
}

// Sum of atomic electronegativities
double sumElectronegativity(const RDKit::ROMol* mol) {
    return sumAtomProperty(mol, [](const RDKit::Atom* a) { 
        return element::getElectronegativity(a->getAtomicNum()); 
    });
}

// Sum of atomic polarizabilities
double sumPolarizability(const RDKit::ROMol* mol) {
    return sumAtomProperty(mol, [](const RDKit::Atom* a) { 
        return element::getPolarizability(a->getAtomicNum()); 
    });
}

// Sum of covalent radii
double sumCovalentRadii(const RDKit::ROMol* mol) {
    return sumAtomProperty(mol, [](const RDKit::Atom* a) { 
        return element::getCovalentRadius(a->getAtomicNum()); 
    });
}

// Get property vectors for all atoms
std::vector<double> getAtomPropertyVector(const RDKit::ROMol* mol, std::function<double(const RDKit::Atom*)> propertyFn) {
    std::vector<double> result;
    if (!mol) return result;
    
    result.reserve(mol->getNumAtoms());
    for (const RDKit::Atom* atom : mol->atoms()) {
        result.push_back(propertyFn(atom));
    }
    return result;
}

// Vector of atom electronegativities
std::vector<double> getElectronegativityVector(const RDKit::ROMol* mol) {
    return getAtomPropertyVector(mol, [](const RDKit::Atom* a) { 
        return element::getElectronegativity(a->getAtomicNum()); 
    });
}

// Vector of atom polarizabilities
std::vector<double> getPolarizabilityVector(const RDKit::ROMol* mol) {
    return getAtomPropertyVector(mol, [](const RDKit::Atom* a) { 
        return element::getPolarizability(a->getAtomicNum()); 
    });
}

// Simple vector statistics helpers
double vectorMean(const std::vector<double>& vec) {
    if (vec.empty()) return 0.0;
    return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

double vectorVariance(const std::vector<double>& vec) {
    if (vec.size() <= 1) return 0.0;
    double mean = vectorMean(vec);
    double sumSq = 0.0;
    for (double val : vec) {
        double diff = val - mean;
        sumSq += diff * diff;
    }
    return sumSq / vec.size();
}

double vectorStdDev(const std::vector<double>& vec) {
    return std::sqrt(vectorVariance(vec));
}

double vectorMin(const std::vector<double>& vec) {
    if (vec.empty()) return 0.0;
    return *std::min_element(vec.begin(), vec.end());
}

double vectorMax(const std::vector<double>& vec) {
    if (vec.empty()) return 0.0;
    return *std::max_element(vec.begin(), vec.end());
}

double vectorRange(const std::vector<double>& vec) {
    if (vec.empty()) return 0.0;
    return vectorMax(vec) - vectorMin(vec);
}

// Autocorrelation calculations
std::vector<double> calculateAutocorrelation(
    const RDKit::ROMol* mol, 
    std::function<double(const RDKit::Atom*)> propertyFn,
    unsigned int maxDistance
) {
    std::vector<double> result(maxDistance + 1, 0.0);
    if (!mol) return result;
    
    // Get atom properties
    std::vector<double> properties;
    properties.reserve(mol->getNumAtoms());
    for (const RDKit::Atom* atom : mol->atoms()) {
        properties.push_back(propertyFn(atom));
    }
    
    // Calculate paths between all atoms up to maxDistance
    unsigned int nAtoms = mol->getNumAtoms();
    
    for (unsigned int i = 0; i < nAtoms; ++i) {
        for (unsigned int j = i; j < nAtoms; ++j) {
            auto path = RDKit::MolOps::getShortestPath(*mol, i, j);
            int dist = path.empty() ? -1 : static_cast<int>(path.size() - 1);
            if (dist >= 0 && static_cast<unsigned int>(dist) <= maxDistance) {
                result[dist] += properties[i] * properties[j];
            }
        }
    }
    
    return result;
}

// Specific autocorrelation functions
std::vector<double> getElectronegativityAutocorrelation(const RDKit::ROMol* mol, unsigned int maxDistance) {
    return calculateAutocorrelation(
        mol,
        [](const RDKit::Atom* a) { return element::getElectronegativity(a->getAtomicNum()); },
        maxDistance
    );
}

std::vector<double> getPolarizabilityAutocorrelation(const RDKit::ROMol* mol, unsigned int maxDistance) {
    return calculateAutocorrelation(
        mol,
        [](const RDKit::Atom* a) { return element::getPolarizability(a->getAtomicNum()); },
        maxDistance
    );
}

// Additional general utility functions

// Calculate path distance matrix using RDKit's shortest path algorithm
std::vector<std::vector<int>> calculateDistanceMatrix(const RDKit::ROMol* mol) {
    if (!mol) return {};
    unsigned int nAtoms = mol->getNumAtoms();
    std::vector<std::vector<int>> distMatrix(nAtoms, std::vector<int>(nAtoms, 0));

    if (nAtoms == 0) return distMatrix;

    // Use RDKit's optimized getDistanceMat for topological distances.
    // useBO=false, useAtomWts=false results in topological distances.
    // The returned pointer is to a 1D array (row-major order, nAtoms x nAtoms).
    // RDKit manages the memory for dmat; do not delete it.
    double *dmat = RDKit::MolOps::getDistanceMat(*mol, false, false, true);

    // Define MAX_DIST_VAL locally as a workaround if RDKit::MAX_DIST_VAL is not found.
    // The typical RDKit value is 1.0e8.
    const double LOCAL_MAX_DIST_VAL = 1.0e8;
    
    for (unsigned int i = 0; i < nAtoms; ++i) {
        for (unsigned int j = 0; j < nAtoms; ++j) {
            double val = dmat[i * nAtoms + j];
            // RDKit::MAX_DIST_VAL (approx 1.0e8) is used for disconnected pairs.
            // For consistency with previous -1 for no path:
            if (val >= LOCAL_MAX_DIST_VAL) {
                distMatrix[i][j] = -1;
            } else {
                // Round to nearest int for topological distances, handles potential float inaccuracies.
                distMatrix[i][j] = static_cast<int>(std::round(val));
            }
        }
    }
    
    return distMatrix;
}

// Get all atom pairs at a specific topological distance
std::vector<std::pair<const RDKit::Atom*, const RDKit::Atom*>> 
getAtomPairsAtDistance(const RDKit::ROMol* mol, unsigned int distance) {
    std::vector<std::pair<const RDKit::Atom*, const RDKit::Atom*>> pairs;
    if (!mol) return pairs;
    
    unsigned int nAtoms = mol->getNumAtoms();
    for (unsigned int i = 0; i < nAtoms; ++i) {
        for (unsigned int j = i+1; j < nAtoms; ++j) {
            auto path = RDKit::MolOps::getShortestPath(*mol, i, j);
            int dist = path.empty() ? -1 : static_cast<int>(path.size() - 1);
            if (dist == static_cast<int>(distance)) {
                pairs.emplace_back(mol->getAtomWithIdx(i), mol->getAtomWithIdx(j));
            }
        }
    }
    
    return pairs;
}

// Calculate number of atoms at specific topological distance from an atom
int getAtomsAtDistanceCount(const RDKit::ROMol* mol, unsigned int atomIdx, unsigned int distance) {
    if (!mol || atomIdx >= mol->getNumAtoms()) return 0;
    
    int count = 0;
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
        if (i == atomIdx) continue;
        
        auto path = RDKit::MolOps::getShortestPath(*mol, atomIdx, i);
        int dist = path.empty() ? -1 : static_cast<int>(path.size() - 1);
        if (dist == static_cast<int>(distance)) {
            count++;
        }
    }
    
    return count;
}

// Molecular graph eccentricity (maximum topological distance from any atom)
int getMolecularEccentricity(const RDKit::ROMol* mol) {
    if (!mol || mol->getNumAtoms() == 0) return 0;
    
    int maxDist = 0;
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
        int atomMaxDist = 0;
        for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
            if (i == j) continue;
            
            auto path = RDKit::MolOps::getShortestPath(*mol, i, j);
            int dist = path.empty() ? -1 : static_cast<int>(path.size() - 1);
            atomMaxDist = std::max(atomMaxDist, dist);
        }
        maxDist = std::max(maxDist, atomMaxDist);
    }
    
    return maxDist;
}

// Calculate molecular diameter (maximum shortest path between any two atoms)
int getMolecularDiameter(const RDKit::ROMol* mol) {
    return getMolecularEccentricity(mol);
}

// Calculate molecular radius (minimum eccentricity of any atom)
int getMolecularRadius(const RDKit::ROMol* mol) {
    if (!mol || mol->getNumAtoms() == 0) return 0;
    
    int minEccentricity = std::numeric_limits<int>::max();
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
        int atomMaxDist = 0;
        for (unsigned int j = 0; j < mol->getNumAtoms(); ++j) {
            if (i == j) continue;
            
            auto path = RDKit::MolOps::getShortestPath(*mol, i, j);
            int dist = path.empty() ? -1 : static_cast<int>(path.size() - 1);
            atomMaxDist = std::max(atomMaxDist, dist);
        }
        minEccentricity = std::min(minEccentricity, atomMaxDist);
    }
    
    return minEccentricity == std::numeric_limits<int>::max() ? 0 : minEccentricity;
}

// Calculate Wiener index (sum of all shortest paths between atoms)
int getWienerIndex(const RDKit::ROMol* mol) {
    if (!mol || mol->getNumAtoms() <= 1) return 0;
    
    int sum = 0;
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
        for (unsigned int j = i+1; j < mol->getNumAtoms(); ++j) {
            auto path = RDKit::MolOps::getShortestPath(*mol, i, j);
            int dist = path.empty() ? -1 : static_cast<int>(path.size() - 1);
            if (dist > 0) {
                sum += dist;
            }
        }
    }
    
    return sum;
}

// Get Balaban J index - a topological index related to molecular connectivity
double getBalabanJIndex(const RDKit::ROMol* mol) {
    if (!mol || mol->getNumAtoms() <= 2) return 0.0;
    
    // Calculate distance sums for each atom (si)
    std::vector<int> distanceSums(mol->getNumAtoms(), 0);
    int pathSum = 0;
    
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
        for (unsigned int j = i+1; j < mol->getNumAtoms(); ++j) {
            auto path = RDKit::MolOps::getShortestPath(*mol, i, j);
            int dist = path.empty() ? -1 : static_cast<int>(path.size() - 1);
            if (dist > 0) {
                distanceSums[i] += dist;
                distanceSums[j] += dist;
                pathSum += dist;
            }
        }
    }
    
    // Calculate edge contribution
    double sum = 0.0;
    for (const RDKit::Bond* bond : mol->bonds()) {
        unsigned int i = bond->getBeginAtomIdx();
        unsigned int j = bond->getEndAtomIdx();
        
        if (distanceSums[i] > 0 && distanceSums[j] > 0) {
            sum += 1.0 / std::sqrt(distanceSums[i] * distanceSums[j]);
        }
    }
    
    int m = mol->getNumBonds();
    int n = mol->getNumAtoms();
    int cycles = m - n + 1;  // Number of cycles in the graph
    
    if (cycles < 1) return 0.0;  // Acyclic molecules have J=0 by definition
    
    return (double)(m) / (double)(cycles) * sum;
}

// Calculates Randic connectivity index (simple version)
double getRandicConnectivityIndex(const RDKit::ROMol* mol) {
    if (!mol || mol->getNumBonds() == 0) return 0.0;
    
    double sum = 0.0;
    for (const RDKit::Bond* bond : mol->bonds()) {
        const RDKit::Atom* a1 = bond->getBeginAtom();
        const RDKit::Atom* a2 = bond->getEndAtom();
        
        int deg1 = a1->getDegree();
        int deg2 = a2->getDegree();
        
        if (deg1 > 0 && deg2 > 0) {
            sum += 1.0 / std::sqrt(deg1 * deg2);
        }
    }
    
    return sum;
}

// Topological polar surface area approximation (fast version)
double getTopologicalPolarSurfaceArea(const RDKit::ROMol* mol) {
    if (!mol) return 0.0;
    
    // Simple and fast approximation based on atom contributions
    // This is much simpler than RDKit's full implementation
    double tpsa = 0.0;
    
    for (const RDKit::Atom* atom : mol->atoms()) {
        int atomNum = atom->getAtomicNum();
        
        if (atomNum == 7) {  // Nitrogen
            if (atom->getIsAromatic()) {
                tpsa += 15.6;  // Aromatic N
            } else {
                tpsa += 17.5;  // Aliphatic N
            }
        } else if (atomNum == 8) {  // Oxygen
            if (atom->getIsAromatic()) {
                tpsa += 13.1;  // Aromatic O
            } else {
                tpsa += 17.0;  // Aliphatic O
            }
        } else if (atomNum == 15) {  // Phosphorus
            tpsa += 13.6;
        } else if (atomNum == 16) {  // Sulfur
            tpsa += 25.3;
        }
    }
    
    return tpsa;
}

// Complexity metrics

// Calculate fragment complexity (number of unique fragments)
int getFragmentCount(const RDKit::ROMol* mol) {
    if (!mol) return 0;
    
    // Use proper RDKit type for fragments
    std::vector<std::vector<int>> fragments;
    RDKit::MolOps::getMolFrags(*mol, fragments);
    
    return fragments.size();
}

// Calculate fsp3 character (fraction of sp3 hybridized carbons)
double getFsp3(const RDKit::ROMol* mol) {
    if (!mol) return 0.0;
    
    int carbonCount = 0;
    int sp3Count = 0;
    
    for (const RDKit::Atom* atom : mol->atoms()) {
        if (atom->getAtomicNum() == 6) {  // Carbon
            carbonCount++;
            
            // Count neighbors and check for planarity
            if (atom->getDegree() <= 2 || // Linear sp or sp carbons
                atom->getIsAromatic() ||   // Aromatic carbons
                atom->getHybridization() == RDKit::Atom::SP2 ||
                atom->getHybridization() == RDKit::Atom::SP) {
                // Not sp3
            } else {
                sp3Count++;
            }
        }
    }
    
    if (carbonCount == 0) return 0.0;
    return static_cast<double>(sp3Count) / carbonCount;
}

// Calculate Bertz complexity index (simplified version)
double getBertzComplexityIndex(const RDKit::ROMol* mol) {
    if (!mol) return 0.0;
    
    // This is a simplified version - the full algorithm is more complex
    double n = mol->getNumAtoms();
    double e = mol->getNumBonds();
    double cyclomatic = e - n + 1;  // For connected graphs
    
    // Get atom diversity factor (unique element count)
    std::set<int> uniqueAtoms;
    for (const RDKit::Atom* atom : mol->atoms()) {
        uniqueAtoms.insert(atom->getAtomicNum());
    }
    double uniqueAtomCount = uniqueAtoms.size();
    
    // Get bond diversity factor
    std::set<RDKit::Bond::BondType> uniqueBonds;
    for (const RDKit::Bond* bond : mol->bonds()) {
        uniqueBonds.insert(bond->getBondType());
    }
    double uniqueBondCount = uniqueBonds.size();
    
    // Combine factors
    return (n + e) * std::log(n + e) * (1 + cyclomatic) * uniqueAtomCount * uniqueBondCount;
}

// Get number of rotatable bonds in rings
int getRingRotatableBondCount(const RDKit::ROMol* mol) {
    if (!mol) return 0;
    
    int count = 0;
    
    mol::ensureRingInfo(mol);
    for (const RDKit::Bond* bond : mol->bonds()) {
        if (bond->getBondType() == RDKit::Bond::SINGLE && 
            mol->getRingInfo()->numBondRings(bond->getIdx()) > 0 &&
            !bond->getIsAromatic()) {
            
            count++;
        }
    }
    
    return count;
}

// Molecular property calculations

// Calculate logP using atomic contributions (simplified Crippen method)
double calculateLogP(const RDKit::ROMol* mol) {
    if (!mol) return 0.0;
    
    // Simple contribution model (very simplified version of Crippen method)
    // Real implementations would use a more sophisticated approach
    double logP = 0.0;
    
    for (const RDKit::Atom* atom : mol->atoms()) {
        int atomNum = atom->getAtomicNum();
        
        if (atomNum == 6) {  // Carbon
            if (atom->getIsAromatic()) {
                logP += 0.13;  // Aromatic carbon
            } else {
                logP += 0.08;  // Aliphatic carbon
            }
        } else if (atomNum == 7) {  // Nitrogen
            logP -= 0.2;
        } else if (atomNum == 8) {  // Oxygen
            logP -= 0.5;
        } else if (atomNum == 9) {  // Fluorine
            logP += 0.14;
        } else if (atomNum == 17) {  // Chlorine
            logP += 0.5;
        } else if (atomNum == 35) {  // Bromine
            logP += 0.6;
        } else if (atomNum == 53) {  // Iodine
            logP += 0.73;
        } else if (atomNum == 16) {  // Sulfur
            logP += 0.25;
        } else if (atomNum == 15) {  // Phosphorus
            logP += 0.2;
        }
    }
    
    return logP;
}

// Calculate molecular refractivity (simplified approach)
double calculateMolarRefractivity(const RDKit::ROMol* mol) {
    if (!mol) return 0.0;
    
    // Simplified atomic contribution model
    double mr = 0.0;
    
    for (const RDKit::Atom* atom : mol->atoms()) {
        int atomNum = atom->getAtomicNum();
        
        if (atomNum == 6) {  // Carbon
            if (atom->getIsAromatic()) {
                mr += 3.1;  // Aromatic carbon
            } else {
                mr += 2.4;  // Aliphatic carbon
            }
        } else if (atomNum == 7) {  // Nitrogen
            mr += 1.7;
        } else if (atomNum == 8) {  // Oxygen
            mr += 1.3;
        } else if (atomNum == 9) {  // Fluorine
            mr += 0.5;
        } else if (atomNum == 17) {  // Chlorine
            mr += 5.5;
        } else if (atomNum == 35) {  // Bromine
            mr += 8.7;
        } else if (atomNum == 53) {  // Iodine
            mr += 14.0;
        } else if (atomNum == 16) {  // Sulfur
            mr += 7.3;
        }
    }
    
    // Add bond contributions
    for (const RDKit::Bond* bond : mol->bonds()) {
        if (bond->getBondType() == RDKit::Bond::DOUBLE) {
            mr += 1.7;
        } else if (bond->getBondType() == RDKit::Bond::TRIPLE) {
            mr += 2.2;
        }
    }
    
    return mr;
}

// Calculate Van der Waals volume of molecule in cubic angstroms
double calculateVanDerWaalsVolume(const RDKit::ROMol* mol) {
    if (!mol) return 0.0;
    
    double volume = 0.0;
    
    // Simple atom-based calculation (ignores overlaps)
    for (const RDKit::Atom* atom : mol->atoms()) {
        int atomNum = atom->getAtomicNum();
        double radius = element::getVanDerWaalsRadius(atomNum);
        
        // Convert from meters to angstroms and calculate volume
        double radiusAng = radius * 1e10; 
        volume += (4.0/3.0) * M_PI * std::pow(radiusAng, 3);
    }
    
    return volume;
}

// Calculate surface area (very simple approximation)
double calculateApproxSurfaceArea(const RDKit::ROMol* mol) {
    if (!mol) return 0.0;
    
    double area = 0.0;
    
    // Simple atom-based calculation (ignores overlaps)
    for (const RDKit::Atom* atom : mol->atoms()) {
        int atomNum = atom->getAtomicNum();
        double radius = element::getVanDerWaalsRadius(atomNum);
        
        // Convert from meters to angstroms
        double radiusAng = radius * 1e10;
        area += 4.0 * M_PI * std::pow(radiusAng, 2);
    }
    
    // Apply empirical correction factor to account for overlaps
    // This is a very rough approximation
    return area * 0.5;
}

// Utility for vector operations with atom properties

// Apply a mathematical operation between two atom property vectors
std::vector<double> combineAtomPropertyVectors(
    const RDKit::ROMol* mol,
    std::function<double(const RDKit::Atom*)> propertyFn1,
    std::function<double(const RDKit::Atom*)> propertyFn2,
    std::function<double(double, double)> combiner
) {
    std::vector<double> result;
    if (!mol) return result;
    
    result.reserve(mol->getNumAtoms());
    for (const RDKit::Atom* atom : mol->atoms()) {
        double prop1 = propertyFn1(atom);
        double prop2 = propertyFn2(atom);
        result.push_back(combiner(prop1, prop2));
    }
    
    return result;
}

// Calculate harmonic mean for a property vector
double harmonicMean(const std::vector<double>& vec) {
    if (vec.empty()) return 0.0;
    
    double sum = 0.0;
    int count = 0;
    
    for (double val : vec) {
        if (val != 0.0) {
            sum += 1.0 / val;
            count++;
        }
    }
    
    return (count > 0) ? count / sum : 0.0;
}

// Calculate geometric mean for a property vector
double geometricMean(const std::vector<double>& vec) {
    if (vec.empty()) return 0.0;
    
    double product = 1.0;
    int count = 0;
    
    for (double val : vec) {
        if (val > 0.0) {
            product *= val;
            count++;
        }
    }
    
    return (count > 0) ? std::pow(product, 1.0 / count) : 0.0;
}

// Histogram utilities

// Generate a histogram for an atom property
std::map<int, int> generateAtomPropertyHistogram(
    const RDKit::ROMol* mol,
    std::function<double(const RDKit::Atom*)> propertyFn,
    double binWidth
) {
    std::map<int, int> histogram;
    if (!mol) return histogram;
    
    for (const RDKit::Atom* atom : mol->atoms()) {
        double value = propertyFn(atom);
        int bin = static_cast<int>(std::floor(value / binWidth));
        histogram[bin]++;
    }
    
    return histogram;
}

// Entropy of atom property distribution
double calculatePropertyEntropy(
    const RDKit::ROMol* mol,
    std::function<double(const RDKit::Atom*)> propertyFn,
    double binWidth
) {
    if (!mol || mol->getNumAtoms() == 0) return 0.0;
    
    // Generate histogram
    auto histogram = generateAtomPropertyHistogram(mol, propertyFn, binWidth);
    
    // Calculate Shannon entropy
    double entropy = 0.0;
    double totalCount = mol->getNumAtoms();
    
    for (const auto& bin : histogram) {
        double p = bin.second / totalCount;
        entropy -= p * std::log2(p);
    }
    
    return entropy;
}

} // namespace mol

//==============================================================================
// STANDARD OBSERVERS IMPLEMENTATION
//==============================================================================

void ElementCountObserver::init(Context& context) {
    // Clear previous data
}

void ElementCountObserver::observeAtom(const RDKit::Atom* atom, Context& context) {
    // Count elements by atomic number
    std::map<int, int> counts;
    if (context.hasProperty("element_counts")) {
        counts = context.getProperty<std::map<int, int>>("element_counts");
    }
    
    counts[atom->getAtomicNum()]++;
    context.setProperty("element_counts", counts);
}

void ElementCountObserver::observeBond(const RDKit::Bond* bond, Context& context) {
    // Nothing to do
}

void ElementCountObserver::finalize(Context& context) {
    // Calculate total atoms
    const auto& counts = context.getProperty<std::map<int, int>>("element_counts");
    int total = 0;
    for (const auto& pair : counts) {
        total += pair.second;
    }
    context.setProperty("total_atoms", total);
}

void RingInfoObserver::init(Context& context) {
    mol::ensureRingInfo(context.getMolecule());
}

void RingInfoObserver::observeAtom(const RDKit::Atom* atom, Context& context) {
    // Track which atoms are in rings
    std::map<unsigned int, bool> atomInRing;
    if (context.hasProperty("atom_in_ring")) {
        atomInRing = context.getProperty<std::map<unsigned int, bool>>("atom_in_ring");
    }
    
    const RDKit::ROMol* mol = context.getMolecule();
    atomInRing[atom->getIdx()] = mol->getRingInfo()->numAtomRings(atom->getIdx()) > 0;
    context.setProperty("atom_in_ring", atomInRing);
}

void RingInfoObserver::observeBond(const RDKit::Bond* bond, Context& context) {
    // Track which bonds are in rings
    std::map<unsigned int, bool> bondInRing;
    if (context.hasProperty("bond_in_ring")) {
        bondInRing = context.getProperty<std::map<unsigned int, bool>>("bond_in_ring");
    }
    
    const RDKit::ROMol* mol = context.getMolecule();
    bondInRing[bond->getIdx()] = mol->getRingInfo()->numBondRings(bond->getIdx()) > 0;
    context.setProperty("bond_in_ring", bondInRing);
}

void RingInfoObserver::finalize(Context& context) {
    // Calculate ring statistics
    const RDKit::ROMol* mol = context.getMolecule();
    
    // Store ring counts by size
    std::map<int, int> ringSizeCounts;
    for (const auto& ring : mol->getRingInfo()->atomRings()) {
        ringSizeCounts[ring.size()]++;
    }
    context.setProperty("ring_size_counts", ringSizeCounts);
    
    // Store total ring count
    context.setProperty("total_rings", mol->getRingInfo()->numRings());
}

void ElectronegativityObserver::init(Context& context) {
    // Clear data
}

void ElectronegativityObserver::observeAtom(const RDKit::Atom* atom, Context& context) {
    // Store electronegativity for each atom
    std::map<unsigned int, double> enValues;
    if (context.hasProperty("atom_en_values")) {
        enValues = context.getProperty<std::map<unsigned int, double>>("atom_en_values");
    }
    
    double en = element::getElectronegativity(atom->getAtomicNum());
    enValues[atom->getIdx()] = en;
    context.setProperty("atom_en_values", enValues);
}

void ElectronegativityObserver::observeBond(const RDKit::Bond* bond, Context& context) {
    // Nothing to do
}

void ElectronegativityObserver::finalize(Context& context) {
    // Calculate electronegativity statistics
    const auto& enValues = context.getProperty<std::map<unsigned int, double>>("atom_en_values");
    
    double sum = 0.0;
    double maxVal = -1.0;
    double minVal = 999.0;
    
    for (const auto& pair : enValues) {
        sum += pair.second;
        maxVal = std::max(maxVal, pair.second);
        minVal = std::min(minVal, pair.second);
    }
    
    // Store statistics
    int count = enValues.size();
    context.setProperty("en_sum", sum);
    context.setProperty("en_max", maxVal);
    context.setProperty("en_min", minVal);
    context.setProperty("en_mean", count > 0 ? sum / count : 0.0);
}

//==============================================================================
// OBSERVER REGISTRATION FUNCTIONS
//==============================================================================
// These functions create and register the standard observers
void registerElementCountObserver() {
    auto observer = std::make_shared<ElementCountObserver>();
    DescriptorRegistry::getInstance().registerObserver(observer);
}

void registerRingInfoObserver() {
    auto observer = std::make_shared<RingInfoObserver>();
    DescriptorRegistry::getInstance().registerObserver(observer);
}

void registerElectronegativityObserver() {
    auto observer = std::make_shared<ElectronegativityObserver>();
    DescriptorRegistry::getInstance().registerObserver(observer);
}

//==============================================================================
// DESCRIPTOR REGISTRY IMPLEMENTATION
//==============================================================================

// Instead of trying to implement all methods directly in common.cpp, create a custom implementation 
// class that delegates to the real registry. This avoids all the redefinition issues.

// Define a registry implementation that won't clash with the auto-generated one
class RegistryImplementation {
private:
    RegistryImplementation() = default;
    std::mutex mutex_;
    std::map<std::string, std::shared_ptr<Descriptor>> descriptors_;
    std::map<std::string, std::shared_ptr<Observer>> observers_;

public:
    static RegistryImplementation& getInstance() {
        static RegistryImplementation instance;
        return instance;
    }

    void registerDescriptor(std::shared_ptr<Descriptor> descriptor) {
        std::lock_guard<std::mutex> lock(mutex_);
        descriptors_[descriptor->getName()] = descriptor;
    }

    std::shared_ptr<Descriptor> getDescriptor(const std::string& name) const {
        auto it = descriptors_.find(name);
        if (it != descriptors_.end()) return it->second;
        return nullptr;
    }

    std::vector<std::shared_ptr<Descriptor>> getDescriptors() const {
        std::vector<std::shared_ptr<Descriptor>> result;
        for (const auto& kv : descriptors_) {
            result.push_back(kv.second);
        }
        return result;
    }

    void registerObserver(std::shared_ptr<Observer> observer) {
        std::lock_guard<std::mutex> lock(mutex_);
        observers_[observer->getName()] = observer;
    }

    std::shared_ptr<Observer> getObserver(const std::string& name) const {
        auto it = observers_.find(name);
        if (it != observers_.end()) return it->second;
        return nullptr;
    }
};

void DescriptorRegistry::registerDescriptor(std::shared_ptr<Descriptor> descriptor) {
    RegistryImplementation::getInstance().registerDescriptor(descriptor);
}

std::shared_ptr<Descriptor> DescriptorRegistry::getDescriptor(const std::string& name) const {
    return RegistryImplementation::getInstance().getDescriptor(name);
}

std::vector<std::string> DescriptorRegistry::getDescriptors() const {
    std::vector<std::string> result;
    for (const auto& kv : RegistryImplementation::getInstance().getDescriptors()) {
        result.push_back(kv->getName());
    }
    return result;
}

void DescriptorRegistry::registerObserver(std::shared_ptr<Observer> observer) {
    RegistryImplementation::getInstance().registerObserver(observer);
}

std::shared_ptr<Observer> DescriptorRegistry::getObserver(const std::string& name) const {
    return RegistryImplementation::getInstance().getObserver(name);
}

} // namespace desfact