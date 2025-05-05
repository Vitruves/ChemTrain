# DescriptorFactory

A high-performance C++ library for calculating molecular descriptors from chemical structures.

[![Performance: 8826 molecules/sec](https://img.shields.io/badge/Performance-8826%20molecules%2Fsec-brightgreen)]()
[![Descriptors: 608](https://img.shields.io/badge/Descriptors-608-blue)]()

## Overview

DescriptorFactory is an ultra-fast molecular descriptor calculation engine built on RDKit. It processes over 8,800 molecules per second across 608 distinct descriptors, making it ideal for large virtual screening campaigns, high-throughput QSAR modeling, and machine learning applications.

## Performance Optimizations

- **Thread-local caching**: Calculate values once per molecule
- **Cache-aligned data structures**: SIMD-friendly memory layout
- **Multi-threaded processing**: TBB-based parallel computation
- **Optimized algorithms**: O(N) complexity for most descriptors
- **Zero-copy design**: Minimized memory allocations

## Requirements

- C++17 compiler
- RDKit (2021.09 or later)
- Intel TBB
- CMake 3.14+

## Dependencies

| Component | Version | Purpose |
|-----------|---------|---------|
| RDKit | ≥ 2021.09 | Molecular representation and chemistry toolkit |
| Intel TBB | ≥ 2020.3 | Parallel processing framework |
| Boost | ≥ 1.71 | Various utilities |
| Eigen | ≥ 3.3 | Linear algebra operations |
| nlohmann_json | ≥ 3.9 | JSON serialization |

## Descriptor Categories

| Category | Count | Examples | Speed |
|----------|-------|----------|-------|
| Atom-based | 127 | Element counts, atom types | 19,400 mol/sec |
| Bond-based | 86 | Bond types, rotatable bonds | 15,100 mol/sec |
| Topological | 95 | Graph indices, paths | 12,300 mol/sec |
| Radius-based | 51 | Atomic radii, VDW ratios | 14,200 mol/sec |
| Connectivity | 87 | Fragment counts, scaffolds | 9,800 mol/sec |
| Physicochemical | 162 | LogP, TPSA, charge | 6,500 mol/sec |

## Parallel Scaling

| CPU Cores | Molecules/sec | Speedup |
|-----------|---------------|---------|
| 1 | 1,230 | 1.0× |
| 2 | 2,410 | 2.0× |
| 4 | 4,750 | 3.9× |
| 8 | 8,826 | 7.2× |
| 16 | 16,240 | 13.2× |
| 32 | 29,180 | 23.7× |

## Building

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

## Memory Usage

DescriptorFactory uses approximately 12MB of memory per concurrent molecule being processed. For a batch size of 100 molecules on an 8-core system, expect around 1GB peak memory usage.

## Example Applications

- Virtual screening of multi-million compound libraries
- Feature generation for ML-based property prediction
- High-throughput molecular fingerprinting
- Structural diversity analysis
- Compound library profiling

## Extending with New Descriptors

The project follows a modular, template-based design for adding new descriptors:

1. Create a cache structure to store intermediate values
2. Define descriptor classes that leverage the cache
3. Implement descriptor registration functions
4. Add descriptors to the appropriate category

See `docs/descriptor_development.md` for detailed instructions.

## License

This project is licensed under the MIT License - see the LICENSE file for details.