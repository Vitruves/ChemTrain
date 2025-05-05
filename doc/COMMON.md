# DescriptorFactory Utility Functions Reference

## Element Properties (namespace desfact::element)

### Basic Element Information

| Function | Description | Usage |
|----------|-------------|-------|
| `getAtomicNumber(const std::string& symbol)` | Get atomic number from element symbol | `int atomicNum = element::getAtomicNumber("C");` |
| `getElementSymbol(int atomicNum)` | Get element symbol from atomic number | `std::string symbol = element::getElementSymbol(6);` |
| `getPeriodicGroup(int atomicNum)` | Get periodic table group | `int group = element::getPeriodicGroup(atom->getAtomicNum());` |
| `getPeriodicPeriod(int atomicNum)` | Get periodic table period | `int period = element::getPeriodicPeriod(atom->getAtomicNum());` |

### Element Properties

| Function | Description | Usage |
|----------|-------------|-------|
| `getElectronegativity(int atomicNum)` | Get Pauling electronegativity | `double en = element::getElectronegativity(atom->getAtomicNum());` |
| `getAllredRochowEN(int atomicNum)` | Get Allred-Rochow electronegativity | `double en = element::getAllredRochowEN(atom->getAtomicNum());` |
| `getAtomicMass(int atomicNum)` | Get atomic mass in kg | `double mass = element::getAtomicMass(atom->getAtomicNum());` |
| `getIonizationEnergy(int atomicNum)` | Get ionization energy in J | `double ie = element::getIonizationEnergy(atom->getAtomicNum());` |
| `getElectronAffinity(int atomicNum)` | Get electron affinity in J | `double ea = element::getElectronAffinity(atom->getAtomicNum());` |
| `getCovalentRadius(int atomicNum)` | Get covalent radius in m | `double rad = element::getCovalentRadius(atom->getAtomicNum());` |
| `getVanDerWaalsRadius(int atomicNum)` | Get van der Waals radius in m | `double rad = element::getVanDerWaalsRadius(atom->getAtomicNum());` |
| `getPolarizability(int atomicNum)` | Get polarizability in Å³ | `double pol = element::getPolarizability(atom->getAtomicNum());` |
| `getDensity(int atomicNum)` | Get density in kg/m³ | `double den = element::getDensity(atom->getAtomicNum());` |
| `getThermalConductivity(int atomicNum)` | Get thermal conductivity in W/m·K | `double tc = element::getThermalConductivity(atom->getAtomicNum());` |
| `getSpecificHeat(int atomicNum)` | Get specific heat in J/g·K | `double sh = element::getSpecificHeat(atom->getAtomicNum());` |
| `getMolarVolume(int atomicNum)` | Get molar volume in cm³/mol | `double mv = element::getMolarVolume(atom->getAtomicNum());` |

### Element Classification

| Function | Description | Usage |
|----------|-------------|-------|
| `isMainGroupElement(int atomicNum)` | Check if element is a main group element | `bool isMain = element::isMainGroupElement(atom->getAtomicNum());` |
| `isTransitionMetal(int atomicNum)` | Check if element is a transition metal | `bool isTM = element::isTransitionMetal(atom->getAtomicNum());` |
| `isAlkaliMetal(int atomicNum)` | Check if element is an alkali metal | `bool isAlkali = element::isAlkaliMetal(atom->getAtomicNum());` |
| `isAlkalineEarthMetal(int atomicNum)` | Check if element is an alkaline earth metal | `bool isAEM = element::isAlkalineEarthMetal(atom->getAtomicNum());` |
| `isNobleGas(int atomicNum)` | Check if element is a noble gas | `bool isNoble = element::isNobleGas(atom->getAtomicNum());` |
| `isHalogen(int atomicNum)` | Check if element is a halogen | `bool isHal = element::isHalogen(atom->getAtomicNum());` |
| `isMetal(int atomicNum)` | Check if element is a metal | `bool isMet = element::isMetal(atom->getAtomicNum());` |
| `isHeteroatom(int atomicNum)` | Check if element is not C or H | `bool isHet = element::isHeteroatom(atom->getAtomicNum());` |

## Molecule Utilities (namespace desfact::mol)

### Ring Utilities

| Function | Description | Usage |
|----------|-------------|-------|
| `ensureRingInfo(const RDKit::ROMol* mol)` | Initialize ring info if needed | `mol::ensureRingInfo(mol);` |
| `getRingCount(const RDKit::ROMol* mol)` | Get total number of rings | `int rings = mol::getRingCount(mol);` |
| `getAromaticRingCount(const RDKit::ROMol* mol)` | Get number of aromatic rings | `int arom = mol::getAromaticRingCount(mol);` |

### Atom Classification

| Function | Description | Usage |
|----------|-------------|-------|
| `isHeteroatom(const RDKit::Atom* atom)` | Check if atom is not C or H | `bool isHet = mol::isHeteroatom(atom);` |
| `isHalogen(const RDKit::Atom* atom)` | Check if atom is a halogen | `bool isHal = mol::isHalogen(atom);` |
| `isMetal(const RDKit::Atom* atom)` | Check if atom is a metal | `bool isMet = mol::isMetal(atom);` |
| `isInRing(const RDKit::Atom* atom)` | Check if atom is in a ring | `bool inRing = mol::isInRing(atom);` |
| `isAromatic(const RDKit::Atom* atom)` | Check if atom is aromatic | `bool isArom = mol::isAromatic(atom);` |
| `isBasicNitrogen(const RDKit::Atom* atom)` | Check if N atom is basic | `bool isBasic = mol::isBasicNitrogen(atom);` |
| `isAcidicOxygen(const RDKit::Atom* atom)` | Check if O atom is acidic | `bool isAcidic = mol::isAcidicOxygen(atom);` |

### Bond Classification

| Function | Description | Usage |
|----------|-------------|-------|
| `isSingleBond(const RDKit::Bond* bond)` | Check if bond is single | `bool isSingle = mol::isSingleBond(bond);` |
| `isDoubleBond(const RDKit::Bond* bond)` | Check if bond is double | `bool isDouble = mol::isDoubleBond(bond);` |
| `isTripleBond(const RDKit::Bond* bond)` | Check if bond is triple | `bool isTriple = mol::isTripleBond(bond);` |
| `isAromaticBond(const RDKit::Bond* bond)` | Check if bond is aromatic | `bool isArom = mol::isAromaticBond(bond);` |
| `isRotatableBond(const RDKit::Bond* bond)` | Check if bond is rotatable | `bool isRot = mol::isRotatableBond(bond);` |

### Counting Functions

| Function | Description | Usage |
|----------|-------------|-------|
| `countAtoms(const RDKit::ROMol* mol, predicate)` | Count atoms matching a predicate | `int n = mol::countAtoms(mol, [](const RDKit::Atom* a) { return a->getAtomicNum() == 7; });` |
| `countBonds(const RDKit::ROMol* mol, predicate)` | Count bonds matching a predicate | `int n = mol::countBonds(mol, [](const RDKit::Bond* b) { return mol::isDoubleBond(b); });` |
| `getBondCount(const RDKit::ROMol* mol)` | Get total number of bonds | `int bonds = mol::getBondCount(mol);` |
| `getHeavyAtomCount(const RDKit::ROMol* mol)` | Get number of non-hydrogen atoms | `int heavy = mol::getHeavyAtomCount(mol);` |
| `getRotatableBondCount(const RDKit::ROMol* mol)` | Get number of rotatable bonds | `int rot = mol::getRotatableBondCount(mol);` |

### Ratio Functions

| Function | Description | Usage |
|----------|-------------|-------|
| `atomFraction(const RDKit::ROMol* mol, predicate)` | Get fraction of atoms matching predicate | `double frac = mol::atomFraction(mol, [](const RDKit::Atom* a) { return mol::isHeteroatom(a); });` |
| `bondFraction(const RDKit::ROMol* mol, predicate)` | Get fraction of bonds matching predicate | `double frac = mol::bondFraction(mol, [](const RDKit::Bond* b) { return mol::isAromaticBond(b); });` |

### Property Calculation

| Function | Description | Usage |
|----------|-------------|-------|
| `getAverageAtomicMass(const RDKit::ROMol* mol)` | Calculate avg atomic mass | `double avgMass = mol::getAverageAtomicMass(mol);` |
| `getAverageElectronegativity(const RDKit::ROMol* mol)` | Calculate avg electronegativity | `double avgEN = mol::getAverageElectronegativity(mol);` |
| `getMaxElectronegativityDifference(const RDKit::ROMol* mol)` | Max electronegativity difference | `double maxDiff = mol::getMaxElectronegativityDifference(mol);` |
| `sumAtomProperty(const RDKit::ROMol* mol, propertyFn)` | Sum custom atom property | `double sum = mol::sumAtomProperty(mol, [](const RDKit::Atom* a) { return a->getTotalDegree(); });` |
| `sumElectronegativity(const RDKit::ROMol* mol)` | Sum all atom electronegativities | `double sumEN = mol::sumElectronegativity(mol);` |
| `sumPolarizability(const RDKit::ROMol* mol)` | Sum all atom polarizabilities | `double sumPol = mol::sumPolarizability(mol);` |
| `sumCovalentRadii(const RDKit::ROMol* mol)` | Sum all covalent radii | `double sumRad = mol::sumCovalentRadii(mol);` |

### Vector Operations

| Function | Description | Usage |
|----------|-------------|-------|
| `getAtomPropertyVector(const RDKit::ROMol* mol, propertyFn)` | Get vector of atom properties | `auto vec = mol::getAtomPropertyVector(mol, [](const RDKit::Atom* a) { return a->getTotalDegree(); });` |
| `getElectronegativityVector(const RDKit::ROMol* mol)` | Get vector of atom electronegativities | `auto enVec = mol::getElectronegativityVector(mol);` |
| `getPolarizabilityVector(const RDKit::ROMol* mol)` | Get vector of atom polarizabilities | `auto polVec = mol::getPolarizabilityVector(mol);` |
| `vectorMean(const std::vector<double>& vec)` | Calculate mean of vector | `double mean = mol::vectorMean(enVec);` |
| `vectorVariance(const std::vector<double>& vec)` | Calculate variance of vector | `double var = mol::vectorVariance(enVec);` |
| `vectorStdDev(const std::vector<double>& vec)` | Calculate std deviation of vector | `double sd = mol::vectorStdDev(enVec);` |
| `vectorMin(const std::vector<double>& vec)` | Find minimum in vector | `double min = mol::vectorMin(enVec);` |
| `vectorMax(const std::vector<double>& vec)` | Find maximum in vector | `double max = mol::vectorMax(enVec);` |
| `vectorRange(const std::vector<double>& vec)` | Calculate range of vector | `double range = mol::vectorRange(enVec);` |

### Autocorrelation Functions

| Function | Description | Usage |
|----------|-------------|-------|
| `calculateAutocorrelation(const RDKit::ROMol* mol, propertyFn, maxDistance)` | Calculate atom property autocorrelation | `auto ac = mol::calculateAutocorrelation(mol, [](const RDKit::Atom* a) { return a->getDegree(); }, 5);` |
| `getElectronegativityAutocorrelation(const RDKit::ROMol* mol, maxDistance)` | Get electronegativity autocorrelation | `auto enAC = mol::getElectronegativityAutocorrelation(mol, 5);` |
| `getPolarizabilityAutocorrelation(const RDKit::ROMol* mol, maxDistance)` | Get polarizability autocorrelation | `auto polAC = mol::getPolarizabilityAutocorrelation(mol, 5);` |

### Distance and Connectivity

| Function | Description | Usage |
|----------|-------------|-------|
| `getAtomDistance(const RDKit::ROMol* mol, idx1, idx2)` | Get distance between atoms | `int dist = mol::getAtomDistance(mol, atom1->getIdx(), atom2->getIdx());` |
| `hasStereoChemistry(const RDKit::ROMol* mol)` | Check if molecule has stereo | `bool hasStereo = mol::hasStereoChemistry(mol);` |


## Topological Descriptors

| Function | Description | Usage |
|----------|-------------|-------|
| `calculateDistanceMatrix(const RDKit::ROMol* mol)` | Get full topological distance matrix | `auto distMat = mol::calculateDistanceMatrix(mol);` |
| `getAtomPairsAtDistance(const RDKit::ROMol* mol, unsigned int distance)` | Get all atom pairs at specific distance | `auto pairs = mol::getAtomPairsAtDistance(mol, 3);` |
| `getAtomsAtDistanceCount(const RDKit::ROMol* mol, unsigned int atomIdx, unsigned int distance)` | Count atoms at given distance from an atom | `int count = mol::getAtomsAtDistanceCount(mol, 0, 2);` |
| `getMolecularEccentricity(const RDKit::ROMol* mol)` | Get maximum topological distance | `int ecc = mol::getMolecularEccentricity(mol);` |
| `getMolecularDiameter(const RDKit::ROMol* mol)` | Get molecular diameter | `int diam = mol::getMolecularDiameter(mol);` |
| `getMolecularRadius(const RDKit::ROMol* mol)` | Get molecular radius | `int rad = mol::getMolecularRadius(mol);` |
| `getWienerIndex(const RDKit::ROMol* mol)` | Calculate Wiener index | `int wiener = mol::getWienerIndex(mol);` |
| `getBalabanJIndex(const RDKit::ROMol* mol)` | Calculate Balaban J index | `double j = mol::getBalabanJIndex(mol);` |
| `getRandicConnectivityIndex(const RDKit::ROMol* mol)` | Calculate Randić connectivity index | `double chi = mol::getRandicConnectivityIndex(mol);` |

## Physicochemical Properties

| Function | Description | Usage |
|----------|-------------|-------|
| `getTopologicalPolarSurfaceArea(const RDKit::ROMol* mol)` | Calculate approximated TPSA | `double tpsa = mol::getTopologicalPolarSurfaceArea(mol);` |
| `calculateLogP(const RDKit::ROMol* mol)` | Calculate approximated logP | `double logp = mol::calculateLogP(mol);` |
| `calculateMolarRefractivity(const RDKit::ROMol* mol)` | Calculate molar refractivity | `double mr = mol::calculateMolarRefractivity(mol);` |
| `calculateVanDerWaalsVolume(const RDKit::ROMol* mol)` | Calculate VdW volume | `double vol = mol::calculateVanDerWaalsVolume(mol);` |
| `calculateApproxSurfaceArea(const RDKit::ROMol* mol)` | Calculate approximated surface area | `double area = mol::calculateApproxSurfaceArea(mol);` |

## Molecular Complexity Metrics

| Function | Description | Usage |
|----------|-------------|-------|
| `getFragmentCount(const RDKit::ROMol* mol)` | Count molecular fragments | `int frags = mol::getFragmentCount(mol);` |
| `getFsp3(const RDKit::ROMol* mol)` | Calculate fraction of sp3 carbons | `double fsp3 = mol::getFsp3(mol);` |
| `getBertzComplexityIndex(const RDKit::ROMol* mol)` | Calculate Bertz complexity index | `double bertz = mol::getBertzComplexityIndex(mol);` |
| `getRingRotatableBondCount(const RDKit::ROMol* mol)` | Count rotatable bonds in rings | `int rbc = mol::getRingRotatableBondCount(mol);` |

## Statistical Utilities

| Function | Description | Usage |
|----------|-------------|-------|
| `harmonicMean(const std::vector<double>& vec)` | Calculate harmonic mean | `double hm = mol::harmonicMean(values);` |
| `geometricMean(const std::vector<double>& vec)` | Calculate geometric mean | `double gm = mol::geometricMean(values);` |
| `generateAtomPropertyHistogram(const RDKit::ROMol* mol, propertyFn, binWidth)` | Create histogram of atom properties | `auto hist = mol::generateAtomPropertyHistogram(mol, [](const RDKit::Atom* a) { return a->getDegree(); });` |
| `calculatePropertyEntropy(const RDKit::ROMol* mol, propertyFn, binWidth)` | Calculate entropy of property distribution | `double entropy = mol::calculatePropertyEntropy(mol, [](const RDKit::Atom* a) { return a->getDegree(); });` |

## Vector Utilities

| Function | Description | Usage |
|----------|-------------|-------|
| `combineAtomPropertyVectors(const RDKit::ROMol* mol, propertyFn1, propertyFn2, combiner)` | Combine atom properties with custom operation | `auto combined = mol::combineAtomPropertyVectors(mol, fn1, fn2, [](double a, double b) { return a*b; });` |
