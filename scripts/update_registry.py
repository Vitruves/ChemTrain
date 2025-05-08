#!/usr/bin/env python3
# scripts/update_registry.py
# Automatically scans descriptor files and generates registry code

import os
import re
import glob
import sys
from pathlib import Path

# Ensure we're running from the project root
os.chdir(Path(__file__).parent.parent)

# Define source directory
SRC_DIR = "src"

def extract_descriptors():
    """Extract all descriptors from implementation files."""
    descriptors = []
    
    # Get all .cpp files in the descriptors directory
    descriptor_files = glob.glob("src/descriptors/*.cpp")
    
    # Pattern to match standard descriptor declarations
    decl_pattern = r'DECLARE_DESCRIPTOR\((\w+),\s*(\w+),\s*"([^"]+)"\)'
    
    # Pattern to match descriptor dependencies
    dep_pattern = r'DESCRIPTOR_DEPENDENCIES\((\w+)\)\s*\{[^{]*return\s*\{([^}]*)\}\s*;\s*\}'

    # New pattern to match the Broto-Moreau macro *usage*
    brotomoreau_macro_pattern = r'REGISTER_SINGLE_BROTO_MOREAU\((\w+),\s*([a-zA-Z]),\s*(\d+)\)'
    
    # Scan all descriptor implementation files
    for file_path in sorted(descriptor_files):
        is_brotomoreau_file = os.path.basename(file_path) == "brotomoreau.cpp"
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Find all standard descriptor declarations
        for match in re.finditer(decl_pattern, content):
            name, category, description = match.groups()
            
            # Additional validation - check that this descriptor has a calculate method implemented
            calc_pattern = f"{name}Descriptor::calculate\\s*\\(\\s*Context\\s*&\\s*context\\s*\\)\\s*const"
            if not re.search(calc_pattern, content):
                print(f"WARNING: Descriptor {name} is declared but has no calculate method in {os.path.basename(file_path)}")
                continue
            
            # Find dependencies for this descriptor
            deps = []
            # Ensure the dependency pattern targets the specific descriptor name
            specific_dep_pattern = dep_pattern.replace(r'\w+', re.escape(name))
            dep_match = re.search(specific_dep_pattern, content)

            if dep_match and dep_match.group(2).strip():
                deps_str = dep_match.group(2)
                # Parse the comma-separated list of dependencies
                deps = [d.strip(' "\'') for d in deps_str.split(',') if d.strip()] # Ensure deps are stripped and non-empty
            
            descriptors.append({
                'name': name,
                'category': category,
                'description': description.strip(),
                'dependencies': deps,
                'source_file': os.path.basename(file_path),
                'is_brotomoreau': False # Flag standard descriptors
            })

        # Find Broto-Moreau descriptors by parsing the macro usage
        if is_brotomoreau_file:
            for match in re.finditer(brotomoreau_macro_pattern, content):
                prop_name_class, prop_identifier, lag = match.groups()
                # Construct the descriptor name (e.g., ATSdm1)
                name = f"ATSd{prop_identifier}{lag}"
                
                # Avoid adding duplicates
                if not any(d['name'] == name for d in descriptors):
                    descriptors.append({
                        'name': name,
                        'category': "BrotoMoreauAutocorrelation", # Assign category
                        'description': f"Broto-Moreau Autocorrelation - {prop_name_class} - Lag {lag}", # Generate description
                        'dependencies': [], # Broto-Moreau are self-contained
                        'source_file': os.path.basename(file_path),
                        'is_brotomoreau': True # Flag Broto-Moreau descriptors
                    })

    print(f"--- Total descriptors found after scanning all files: {len(descriptors)} ---")
    return descriptors

def extract_observers():
    """Extract all observers from the common implementation."""
    observers = []
    
    # Pattern to match observer class declarations
    observer_pattern = r'class\s+(\w+)\s*:\s*public\s+Observer\s*\{'
    
    # Check in common.hpp
    with open("src/common.hpp", 'r') as f:
        content = f.read()
    
    for match in re.finditer(observer_pattern, content):
        observer_name = match.group(1)
        observers.append({
            'name': observer_name
        })
    
    # Also check descriptor files for custom observers
    for file_path in glob.glob("src/descriptors/*.cpp"):
        with open(file_path, 'r') as f:
            content = f.read()
        
        for match in re.finditer(observer_pattern, content):
            observer_name = match.group(1)
            observers.append({
                'name': observer_name,
                'source_file': os.path.basename(file_path)
            })
    
    return observers

def generate_registry_hpp(descriptors, observers):
    """Generate the registry.hpp header file."""
    # Get previous descriptors to track removals
    prev_descriptors = get_previous_descriptors()
    current_descriptors = set(d['name'] for d in descriptors)
    
    # Start with the header
    output = [
        "// AUTO-GENERATED FILE - DO NOT EDIT",
        "// Generated by scripts/update_registry.py",
        "",
        "#pragma once",
        "",
        "namespace desfact {"
    ]
    
    # Add forward declarations for all registration functions
    output.append("")
    output.append("// Forward declarations for descriptor registration functions")
    for desc in sorted(descriptors, key=lambda x: x['name']):
        func_name = f"register_{desc['name']}Descriptor"
        output.append(f"void {func_name}();")
    
    output.extend([
        "",
        "// Initialize the registry - ensures all static registrations take place",
        "void initializeRegistry();",
        "",
        "} // namespace desfact",
        ""
    ])
    
    # Write to file
    with open("src/registry.hpp", 'w') as f:
        f.write("\n".join(output))
    
    # Report changes
    added = current_descriptors - prev_descriptors
    removed = prev_descriptors - current_descriptors
    
    print(f"Generated src/registry.hpp with {len(descriptors)} descriptor forward declarations")
    
    if removed:
        # Clean up any header files that might contain removed descriptor declarations
        for header_file in glob.glob("src/descriptors/*.hpp"):
            with open(header_file, 'r') as f:
                content = f.read()
            
            modified = False
            for desc_name in removed:
                # Look for class declarations of removed descriptors
                class_pattern = f"class {desc_name}Descriptor[^}}]*}}"
                if re.search(class_pattern, content, re.DOTALL):
                    content = re.sub(class_pattern, "", content, flags=re.DOTALL)
                    modified = True
                    print(f"Removed class declaration for {desc_name}Descriptor from {os.path.basename(header_file)}")
            
            if modified:
                with open(header_file, 'w') as f:
                    f.write(content)

def generate_registry_cpp(descriptors, observers):
    """Generate the registry.cpp file"""
    # Start with the header
    registry_content = [
        "// AUTO-GENERATED FILE - DO NOT EDIT DIRECTLY",
        "// Generated by scripts/update_registry.py",
        "",
        "",
        "// Include descriptor headers if necessary (or rely on forward declarations)",
        "",
        "namespace desfact {",
        ""
    ]
    
    # Add forward declarations for observer registration functions
    if observers:
        registry_content.append("// Forward declarations for observer registration functions")
        for observer in sorted(observers, key=lambda x: x['name']):
            registry_content.append(f"void register{observer['name']}();")
        registry_content.append("")
    
    # Get previous descriptors to track removals - enforce using only the new, filtered descriptors
    prev_descriptors = get_previous_descriptors()
    current_descriptors = set(d['name'] for d in descriptors)
    
    # Add forward declarations for all registration functions - ONLY current descriptors
    registry_content.append("// Forward declarations for descriptor registration functions")
    # Sort descriptors ensuring Broto-Moreau appear correctly
    sorted_descriptors = sorted(descriptors, key=lambda x: x['name'])
    for desc in sorted_descriptors:
        # Function name format depends on how it was found
        func_name = f"register_{desc['name']}Descriptor"
        registry_content.append(f"void {func_name}();")

    registry_content.append("")
    
    # Add the initializeRegistry function
    registry_content.extend([
        "// Initialize registry - this calls all registration functions",
        "void initializeRegistry() {"
    ])
    
    # First register all observers
    if observers:
        registry_content.append("    // Register standard observers first")
        for observer in sorted(observers, key=lambda x: x['name']):
            registry_content.append(f"    register{observer['name']}();")
        registry_content.append("")
    
    # Then register all descriptors
    registry_content.append("    // Register all descriptors")
    for desc in sorted_descriptors:
        # Function name format depends on how it was found
        func_name = f"register_{desc['name']}Descriptor"
        registry_content.append(f"    {func_name}();")
    
    # Add the footer
    registry_content.extend([
        "}",
        "",
        "} // namespace desfact",
        ""
    ])
    
    # Write the updated registry file
    with open("src/registry.cpp", 'w') as f:
        f.write("\n".join(registry_content))
    
    # Report changes
    added = current_descriptors - prev_descriptors
    removed = prev_descriptors - current_descriptors
    
    print(f"Updated registry file with {len(descriptors)} descriptors and {len(observers)} observers")
    
    if added:
        print(f"Added descriptors: {', '.join(sorted(added))}")
    if removed:
        print(f"Removed descriptors: {', '.join(sorted(removed))}")
        
        # Remove registration functions for removed descriptors from source files
        for desc_name in removed:
            for cpp_file in glob.glob("src/descriptors/*.cpp"):
                with open(cpp_file, 'r') as f:
                    content = f.read()
                
                # Look for the registration function
                reg_function = f"void register_{desc_name}Descriptor()"
                if reg_function in content:
                    # Remove the entire registration function
                    pattern = f"void register_{desc_name}Descriptor\\(\\)[^}}]*}}"
                    new_content = re.sub(pattern, "", content, flags=re.DOTALL)
                    
                    # Write back the cleaned content
                    with open(cpp_file, 'w') as f:
                        f.write(new_content)
                    print(f"Removed registration function for {desc_name}Descriptor from {os.path.basename(cpp_file)}")
            
            # Also check C wrapper files
            c_wrapper_file = "src/descriptors/c_wrappers.cpp"
            if os.path.exists(c_wrapper_file):
                with open(c_wrapper_file, 'r') as f:
                    content = f.read()
                
                if f"register_{desc_name}Descriptor" in content:
                    # Remove both class and registration function for C wrappers
                    class_pattern = f"// --- {desc_name} Wrapper ---[\\s\\S]*?void register_{desc_name}Descriptor\\(\\)[^}}]*}}"
                    new_content = re.sub(class_pattern, "", content, flags=re.DOTALL)
                    
                    with open(c_wrapper_file, 'w') as f:
                        f.write(new_content)
                    print(f"Removed C wrapper for {desc_name}Descriptor from c_wrappers.cpp")

def add_registration_functions():
    observer_reg_code = """
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
"""

    with open("src/common.cpp", 'r') as f:
        content = f.read()

    # Only insert if not already present
    if "//==============================================================================\n// OBSERVER REGISTRATION FUNCTIONS" not in content:
        # Find the start of the desfact namespace
        ns_pos = content.find("namespace desfact")
        if ns_pos == -1:
            print("Error: Could not find 'namespace desfact' in src/common.cpp")
            return

        # Find the opening brace of the namespace
        brace_pos = content.find("{", ns_pos)
        if brace_pos == -1:
            print("Error: Could not find opening brace for 'namespace desfact' in src/common.cpp")
            return

        # Insert after the opening brace and any whitespace/newlines
        insert_pos = brace_pos + 1
        while insert_pos < len(content) and content[insert_pos] in " \t\r\n":
            insert_pos += 1

        # Insert the code at insert_pos
        new_content = content[:insert_pos] + "\n" + observer_reg_code + content[insert_pos:]

        with open("src/common.cpp", 'w') as f:
            f.write(new_content)

        print("Added observer registration functions to common.cpp")

def add_registration_functions_to_descriptor_files(descriptors):
    """Add registration functions to descriptor files if they don't exist."""
    for desc in descriptors:
        cpp_file = os.path.join("src/descriptors", desc['source_file'])
        
        with open(cpp_file, 'r') as f:
            content = f.read()
            
        # Check if registration function already exists
        reg_function = f"void register_{desc['name']}Descriptor()"
        if reg_function not in content:
            # Add registration function
            reg_impl = f"""
void register_{desc['name']}Descriptor() {{
    auto descriptor = std::make_shared<{desc['name']}Descriptor>();
    auto& registry = DescriptorRegistry::getInstance();
    registry.registerDescriptor(descriptor);
}}
"""
            # Add at the end of namespace
            namespace_end = content.rfind("} // namespace desfact")
            if namespace_end != -1:
                content = content[:namespace_end] + reg_impl + "\n" + content[namespace_end:]
                
                with open(cpp_file, 'w') as f:
                    f.write(content)
                    
                print(f"Added registration function for {desc['name']}Descriptor to {desc['source_file']}")

def generate_descriptors_md(descriptors):
    """Generate the DESCRIPTORS.md documentation file."""
    output = [
        "# Available Descriptors",
        "",
        f"Total descriptors: {len(descriptors)}",
        ""
    ]

    # Group by category
    descriptors_by_category = {}
    for desc in descriptors:
        cat = desc.get('category', 'Uncategorized')
        if cat not in descriptors_by_category:
            descriptors_by_category[cat] = []
        descriptors_by_category[cat].append(desc)

    # Sort categories and descriptors within categories
    for category in sorted(descriptors_by_category.keys()):
        output.append(f"## {category}")
        output.append("")
        output.append("| Name | Description | Dependencies | Source File |")
        output.append("|---|---|---|---|")
        
        sorted_descs = sorted(descriptors_by_category[category], key=lambda x: x['name'])
        for desc in sorted_descs:
             # Adjust name based on whether it's Broto-Moreau or standard
            name_display = f"{desc['name']}Descriptor"
            deps_str = ", ".join(f"`{d}`" for d in desc.get('dependencies', [])) if desc.get('dependencies') else "None"
            output.append(f"| `{name_display}` | {desc['description']} | {deps_str} | `{desc['source_file']}` |")
        output.append("")

    with open("DESCRIPTORS.md", 'w') as f:
        f.write("\n".join(output))

    print(f"Generated DESCRIPTORS.md with documentation for {len(descriptors)} descriptors")

def get_previous_descriptors():
    """Read previous descriptors from registry.hpp and registry.cpp if they exist."""
    prev = set()
    
    # Check registry.hpp first (most authoritative)
    if os.path.exists("src/registry.hpp"):
        with open("src/registry.hpp", 'r') as f:
            for line in f:
                if line.strip().startswith("void register_") and "Descriptor();" in line:
                    name = line.strip().split()[1]
                    if name.startswith("register_") and name.endswith("Descriptor();"):
                        # Extract just the descriptor name without register_ and Descriptor();
                        name = name[len("register_"):-len("Descriptor();")]
                    prev.add(name)
    
    # Also check registry.cpp as a backup
    if os.path.exists("src/registry.cpp") and not prev:
        with open("src/registry.cpp", 'r') as f:
            for line in f:
                if "register_" in line and "Descriptor();" in line and line.strip().startswith("void "):
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        name = parts[1]
                        if name.startswith("register_") and name.endswith("Descriptor();"):
                            name = name[len("register_"):-len("Descriptor();")]
                            prev.add(name)
    
    # Double-check by scanning descriptor files for implementation functions
    for file_path in glob.glob("src/descriptors/*.cpp") + glob.glob("src/descriptors/*.c"):
        with open(file_path, 'r') as f:
            content = f.read()
            
        # Look for registration function definitions
        for match in re.finditer(r'void\s+register_(\w+)Descriptor\s*\(\s*\)', content):
            prev.add(match.group(1))
    
    return prev

def scan_for_c_descriptors():
    c_descriptors = []
    
    for file_path in glob.glob(f"{SRC_DIR}/descriptors/*.c"):
        # Skip template file that doesn't contain actual implementations
        if os.path.basename(file_path) == "c_descriptor_template.c":
            continue
            
        filename = os.path.basename(file_path)
        with open(file_path, 'r') as file:
            content = file.read()
            
            # Check for GetSmilesFunc typedef which may conflict
            if "typedef" in content and "GetSmilesFunc" in content:
                print(f"WARNING: {file_path} contains its own GetSmilesFunc typedef which may conflict with cregistry.h")
            
            # Look for C-style implementations with either const void* or Context* first parameter
            # Make sure we're only catching actual function definitions, not just declarations
            function_defs = re.findall(r'double\s+calculate(\w+)\s*\(\s*(?:const\s+void\*|const\s+Context\*)[^)]*\)\s*\{', content)
            
            if function_defs:
                for func_name in function_defs:
                    # Verify this is a complete implementation by checking for closing brace
                    func_pattern = f"double\\s+calculate{func_name}\\s*\\([^)]*\\)\\s*\\{{[\\s\\S]*?\\}}"
                    if re.search(func_pattern, content):
                        c_descriptors.append({
                            'name': func_name,
                            'description': f"C implementation of {func_name}",
                            'source_file': filename,
                            'dependencies': []
                        })
                    else:
                        print(f"WARNING: Found incomplete function definition for calculate{func_name} in {filename}")
                
    print(f"--- Found {len(c_descriptors)} descriptors in C implementation files ---")
    
    # Additional validation step
    validated_c_descriptors = []
    for desc in c_descriptors:
        file_path = os.path.join("src/descriptors", desc['source_file'])
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check for actual implementation (not just declaration)
        func_name = f"calculate{desc['name']}"
        if f"{func_name}(" in content and "{" in content and "return" in content:
            validated_c_descriptors.append(desc)
        else:
            print(f"WARNING: Function {func_name} in {desc['source_file']} doesn't appear to have a proper implementation")
    
    print(f"--- Validated {len(validated_c_descriptors)} descriptor implementations ---")
    return validated_c_descriptors

def generate_c_wrapper_file(c_descriptors):
    """Generate the C wrapper implementation file."""
    if not c_descriptors:
        return
    
    # Create the wrapper file content
    wrapper_content = [
        "// AUTO-GENERATED FILE - DO NOT EDIT DIRECTLY",
        "// Generated by scripts/update_registry.py",
        "",
        "#include \"../common.hpp\"",
        "#include \"../cregistry.h\"",
        "#include <string>",
        "",
        "// C linkage for our helper function",
        "extern \"C\" {",
        "// Helper function that maintains string lifetime",
        "static thread_local std::string tls_smiles_buffer;",
        "const char* getSmilesCFunc(const void* ctx) {", 
        "    auto cpp_ctx = reinterpret_cast<const desfact::Context*>(ctx);",
        "    tls_smiles_buffer = cpp_ctx->getSmiles();",
        "    return tls_smiles_buffer.c_str();",
        "}",
        "}",
        "",
        "namespace desfact {",
        ""
    ]
    
    # Skip known conflicting descriptors
    conflicting_descriptors = ["SmilesAmideFraction", "SmilesAcetalCount"]
    
    # Add C wrapper classes
    for desc in c_descriptors:
        name = desc['name']
        # Skip known conflicts
        if name in conflicting_descriptors:
            print(f"Skipping known conflicting descriptor: {name}")
            continue
            
        wrapper_content.extend([
            f"// --- {name} Wrapper ---",
            f"DECLARE_DESCRIPTOR({name}, C_Descriptor, \"{desc['description']}\")",
            f"DESCRIPTOR_DEPENDENCIES({name}) {{ return {{}}; }}",
            f"DescriptorResult {name}Descriptor::calculate(Context& context) const {{",
            f"    // Double cast to bypass type system - first to void*, then to the C type",
            f"    void* v_ptr = &context;",
            f"    const void* c_ctx = static_cast<const void*>(v_ptr);",
            f"    // Explicitly cast the function pointer to the expected type",
            f"    GetSmilesFunc get_smiles_func = reinterpret_cast<GetSmilesFunc>(getSmilesCFunc);",
            f"    return calculate{name}(c_ctx, get_smiles_func);",
            f"}}",
            ""
        ])
    
    # Add registration functions (also skip conflicts)
    for desc in c_descriptors:
        name = desc['name']
        if name in conflicting_descriptors:
            continue
            
        wrapper_content.extend([
            f"void register_{name}Descriptor() {{",
            f"    auto descriptor = std::make_shared<{name}Descriptor>();",
            f"    auto& registry = DescriptorRegistry::getInstance();",
            f"    registry.registerDescriptor(descriptor);",
            f"}}",
            ""
        ])
    
    wrapper_content.extend([
        "} // namespace desfact",
        ""
    ])
    
    # Write the wrapper file
    wrapper_file_path = "src/descriptors/c_wrappers.cpp"
    with open(wrapper_file_path, 'w') as f:
        f.write("\n".join(wrapper_content))
    
    print(f"Generated C wrapper file at {wrapper_file_path}")
    
    # Add to CMakeLists.txt if not already there
    add_wrapper_to_sources("src/descriptors/c_wrappers.cpp")

def add_wrapper_to_sources(wrapper_file):
    """Add the C wrapper file to the SOURCES in CMakeLists.txt if not already there."""
    cmake_file = "CMakeLists.txt"
    
    with open(cmake_file, 'r') as f:
        content = f.read()
    
    sources_pattern = r"set\(SOURCES\s+([^)]+)\)"
    match = re.search(sources_pattern, content, re.DOTALL)
    
    if match and wrapper_file not in match.group(1):
        # Add the wrapper file to SOURCES
        new_sources = match.group(1).rstrip() + f"\n    {wrapper_file}\n"
        new_content = content.replace(match.group(0), f"set(SOURCES\n{new_sources})")
        
        with open(cmake_file, 'w') as f:
            f.write(new_content)
        
        print(f"Added {wrapper_file} to SOURCES in CMakeLists.txt")

def update_c_registry_header(c_descriptors):
    """Update the cregistry.h header file."""
    if not c_descriptors:
        return
    
    header_content = [
        "// Auto-generated C registry header - DO NOT EDIT DIRECTLY",
        "// Generated by scripts/update_registry.py",
        "#ifndef CREGISTRY_H",
        "#define CREGISTRY_H",
        "",
        "#ifdef __cplusplus",
        "extern \"C\" {",
        "#endif",
        "",
        "// Forward declaration for the Context interface",
        "struct Context;",
        "typedef struct Context Context;",
        "",
        "// Function pointer type for getting SMILES string",
        "typedef const char* (*GetSmilesFunc)(const void*);",
        "",
        "// Helper function that will be provided by the C++ implementation",
        "const char* getSmilesCFunc(const void* ctx);",
        "",
        "// Common interface for all C descriptors",
        "// Each descriptor must implement this function to be integrated with the C++ framework",
        ""
    ]
    
    # Add function declarations - IMPORTANT: Using const void* to match implementation
    for desc in c_descriptors:
        name = desc['name']
        header_content.append(f"// {desc['description']}")
        header_content.append(f"double calculate{name}(const void* context, GetSmilesFunc getSmilesFunc);")
        header_content.append("")
    
    header_content.extend([
        "#ifdef __cplusplus",
        "}",
        "#endif",
        "",
        "#endif // CREGISTRY_H"
    ])
    
    # Write the header file
    with open("src/cregistry.h", 'w') as f:
        f.write("\n".join(header_content))
    
    print("Updated src/cregistry.h header file")

def validate_c_descriptor_files(c_descriptors):
    # Validate that C descriptor files correctly implement the interface
    for desc in c_descriptors:
        file_path = os.path.join("src/descriptors", desc['source_file'])
        
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check if the file includes cregistry.h
        if "#include \"../cregistry.h\"" not in content and "#include <cregistry.h>" not in content:
            print(f"WARNING: {file_path} does not include cregistry.h")
            
        # Check if file has its own typedef for GetSmilesFunc that might conflict
        if "typedef" in content and "GetSmilesFunc" in content:
            print(f"WARNING: {file_path} contains its own GetSmilesFunc typedef which may conflict with cregistry.h")
            
        # Ensure the calculateXXX function has the right signature
        func_name = f"calculate{desc['name']}"
        if func_name not in content:
            print(f"ERROR: {file_path} does not contain the expected function {func_name}")

def fix_topological_c_file():
    """Fix the topological.c file to use the correct type signature."""
    file_path = "src/descriptors/topological.c"
    
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            
        # Remove the duplicate typedef that conflicts with cregistry.h
        lines = content.split('\n')
        fixed_lines = []
        
        # Regular expressions to match duplicated const in various patterns
        import re
        
        for line in lines:
            # Skip any Context or GetSmilesFunc typedefs
            if "typedef struct Context Context;" in line:
                continue
            elif "typedef char* (*GetSmilesFunc)" in line or "typedef char *(*GetSmilesFunc)" in line:
                continue
            else:
                # Fix duplicate const in parameter declarations
                line = re.sub(r"const\s+const\s+", "const ", line)
                fixed_lines.append(line)
        
        # Join lines back to a string
        fixed_content = '\n'.join(fixed_lines)
        
        # Additional string replacements for any cases that might have been missed
        fixed_content = fixed_content.replace("const const char* smiles", "const char* smiles")
        fixed_content = fixed_content.replace("const const char *smiles", "const char* smiles")
        fixed_content = fixed_content.replace("const const char * smiles", "const char* smiles")
        fixed_content = fixed_content.replace("const const void* context", "const void* context")
        fixed_content = fixed_content.replace("const const void *context", "const void* context")
        
        # Fix function declarations to match cregistry.h
        fixed_content = fixed_content.replace("char* smiles = getSmilesFunc(context)", 
                                             "const char* smiles = getSmilesFunc(context)")
        
        with open(file_path, 'w') as f:
            f.write(fixed_content)
            
        print(f"Fixed type signatures in {file_path}")
    except Exception as e:
        print(f"ERROR: Failed to fix {file_path}: {str(e)}")

def fix_bitsoperations_c_file():
    """Fix the bitsoperations.c file to use the correct type signature."""
    file_path = "src/descriptors/bitsoperations.c"
    
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            
        # Remove any duplicate typedefs that conflict with cregistry.h
        lines = content.split('\n')
        fixed_lines = []
        
        import re
        
        for line in lines:
            if "typedef struct Context Context;" in line:
                # Skip this line
                continue
            elif "typedef char* (*GetSmilesFunc)" in line or "typedef char *(*GetSmilesFunc)" in line or "typedef const char* (*GetSmilesFunc)" in line:
                # Skip this line
                continue
            else:
                # Fix duplicate const in parameter declarations
                line = re.sub(r"const\s+const\s+", "const ", line)
                fixed_lines.append(line)
        
        # Join lines back to a string
        fixed_content = '\n'.join(fixed_lines)
        
        # Additional string replacements for any cases that might have been missed
        fixed_content = fixed_content.replace("const const char* smiles", "const char* smiles")
        fixed_content = fixed_content.replace("const const char *smiles", "const char* smiles")
        fixed_content = fixed_content.replace("const const char * smiles", "const char* smiles")
        fixed_content = fixed_content.replace("const const void* context", "const void* context")
        fixed_content = fixed_content.replace("const const void *context", "const void* context")
        
        # Fix function declarations to match cregistry.h
        fixed_content = fixed_content.replace("char* smiles = getSmilesFunc(context)", 
                                             "const char* smiles = getSmilesFunc(context)")
        
        with open(file_path, 'w') as f:
            f.write(fixed_content)
            
        print(f"Fixed type signatures in {file_path}")
    except Exception as e:
        print(f"ERROR: Failed to fix {file_path}: {str(e)}")

def create_c_descriptor_template():
    """Create a template for C descriptor implementation.
    This is an empty function since we no longer need to create template files."""
    # This is an empty implementation to avoid the NameError
    pass

def register_c_descriptor(descriptor_name):
    header_file = f"{SRC_DIR}/cregistry.h"
    
    with open(header_file, 'r') as file:
        content = file.read()
    
    # Look for existing declaration with either signature pattern
    signature1 = f"double calculate{descriptor_name}(const void* context, GetSmilesFunc getSmilesFunc);"
    signature2 = f"double calculate{descriptor_name}(const Context* context, GetSmilesFunc getSmilesFunc);"
    
    if signature1 in content or signature2 in content:
        print(f"Descriptor calculate{descriptor_name} already registered in C registry")
        return
    
    # Add the declaration with the standard signature
    # Insert before the closing #ifdef
    with open(header_file, 'r') as file:
        lines = file.readlines()
    
    insert_position = -1
    for i, line in enumerate(lines):
        if "#endif" in line and "CREGISTRY_H" in line:
            insert_position = i
            break
    
    if insert_position > 0:
        lines.insert(insert_position, f"double calculate{descriptor_name}(const void* context, GetSmilesFunc getSmilesFunc);\n")
        
        with open(header_file, 'w') as file:
            file.writelines(lines)
        
        print(f"Registered C descriptor: calculate{descriptor_name}")
    else:
        print(f"ERROR: Could not find position to insert declaration for {descriptor_name}")

def verify_registry_files(descriptors):
    """Verify registry files only contain current descriptors."""
    current_names = set(d['name'] for d in descriptors)
    
    # Check registry.hpp
    with open("src/registry.hpp", 'r') as f:
        hpp_content = f.read()
    
    # Check for any register_XXXDescriptor() that's not in current_names
    for match in re.finditer(r'void\s+register_(\w+)Descriptor\(\);', hpp_content):
        name = match.group(1)
        if name not in current_names:
            print(f"ERROR: Found unexpected descriptor in registry.hpp: {name}")
            return False
    
    # Check registry.cpp
    with open("src/registry.cpp", 'r') as f:
        cpp_content = f.read()
    
    for match in re.finditer(r'void\s+register_(\w+)Descriptor\(\);', cpp_content):
        name = match.group(1)
        if name not in current_names:
            print(f"ERROR: Found unexpected descriptor in registry.cpp: {name}")
            return False
    
    # Everything is good if we get here
    print("Verification complete: Registry files only contain current descriptors")
    return True

def main():
    """Main execution function."""
    print("Scanning for C++ descriptors...")
    cpp_descriptor_dicts = extract_descriptors()
    
    print("Scanning for C descriptors...")
    c_descriptor_dicts = scan_for_c_descriptors()
    
    if c_descriptor_dicts:
        print("Validating C descriptor implementations...")
        validate_c_descriptor_files(c_descriptor_dicts)
        # create_c_descriptor_template() # This function is empty, can be removed or kept if planned for future
        
        # Fix type signatures in C files like topological.c and bitsoperations.c
        # (Assuming these functions correctly use c_descriptor_dicts or scan files again)
        topological_path = os.path.join(SRC_DIR, "descriptors", "topological.c")
        if os.path.exists(topological_path):
            print("Checking topological.c file...")
            fix_topological_c_file()
            
        bitsoperations_path = os.path.join(SRC_DIR, "descriptors", "bitsoperations.c")
        if os.path.exists(bitsoperations_path):
            print("Checking bitsoperations.c file...")
            fix_bitsoperations_c_file()
    
    print("Scanning for observers...")
    observers = extract_observers()

    # Consolidate descriptors, handle name collisions (prefer C), and filter stubs
    print("\nConsolidating and filtering descriptors...")
    
    # Store chosen descriptor dict by name, C version takes precedence
    final_descriptors_map = {} 
    skipped_descriptors_info = [] # For stubs or non-preferred duplicates

    # Process C++ descriptors first, check for stubs
    for desc_dict in cpp_descriptor_dicts:
        name = desc_dict['name']
        source_file = desc_dict['source_file']
        is_stub = False
        
        file_path = os.path.join(SRC_DIR, "descriptors", source_file)
        if not os.path.exists(file_path):
            print(f"WARNING: Source file {file_path} for C++ descriptor {name} does not exist. Treating as stub.")
            is_stub = True
        else:
            # Broto-Moreau are special; their 'is_brotomoreau' flag is in desc_dict
            # Standard C++ descriptors need their calculate method checked for stub patterns
            if not desc_dict.get('is_brotomoreau', False):
                with open(file_path, 'r') as f:
                    content = f.read()
                # Standard C++ stub patterns
                pattern_empty_return = f"{name}Descriptor::calculate[^{{]*{{\\s*return\\s*0\\.0;?\\s*}}"
                pattern_todo = f"{name}Descriptor::calculate[^{{]*{{\\s*// TODO\\s*}}"
                if re.search(pattern_empty_return, content) or re.search(pattern_todo, content):
                    is_stub = True
            # For Broto-Moreau, stubs might be harder to detect this way if they are empty registration functions.
            # We assume Broto-Moreau found by extract_descriptors are non-stub unless their calculate method is explicitly a stub.

        if is_stub:
            print(f"INFO: C++ Descriptor {name} from {source_file} is a stub. Skipping.")
            skipped_descriptors_info.append({**desc_dict, 'reason': 'Stub Implementation (C++)'})
        else:
            # Add C++ descriptor if name not yet taken (by a C descriptor that might be processed later, though C is preferred)
            # This logic effectively means C++ is added if no C version will override it.
            if name not in final_descriptors_map:
                 final_descriptors_map[name] = desc_dict

    # Process C descriptors, check for stubs (already done by scan_for_c_descriptors), and override C++ if name collision
    for desc_dict in c_descriptor_dicts:
        name = desc_dict['name']
        source_file = desc_dict['source_file']
        # C descriptors are assumed valid (non-stub) if returned by scan_for_c_descriptors
        # (which checks for function definition with '{' and 'return')
        
        file_path = os.path.join(SRC_DIR, "descriptors", source_file)
        if not os.path.exists(file_path): # Should have been caught by scan_for_c_descriptors if it's thorough
            print(f"WARNING: Source file {file_path} for C descriptor {name} does not exist (should not happen if scan_for_c_descriptors is robust). Treating as skipped.")
            skipped_descriptors_info.append({**desc_dict, 'reason': 'Missing Source File (C)'})
            continue

        if name in final_descriptors_map:
            existing_desc = final_descriptors_map[name]
            if existing_desc['source_file'].endswith('.cpp'): # Check if the existing one is C++
                print(f"INFO: C Descriptor '{name}' from '{source_file}' overrides C++ version from '{existing_desc['source_file']}'.")
            # If existing is also C, it's a C-C duplicate name from different files, scan_for_c_descriptors should ideally handle this or this loop takes the latter.
            # For C vs C++, C always wins here.
        final_descriptors_map[name] = desc_dict

    all_descriptors_for_codegen = sorted(list(final_descriptors_map.values()), key=lambda x: x['name'])
    
    if skipped_descriptors_info:
        print(f"\nSkipped {len(skipped_descriptors_info)} descriptors due to being stubs or non-preferred duplicates:")
        for desc_info in skipped_descriptors_info:
            print(f"  - {desc_info['name']} (from {desc_info['source_file']}, reason: {desc_info.get('reason', 'N/A')})")

    # Get previously registered descriptor names to report added/removed accurately
    # This should use the names from the *newly decided* set of descriptors
    prev_descriptors_names = get_previous_descriptors() # Reads from old registry.hpp
    current_descriptor_names = set(d['name'] for d in all_descriptors_for_codegen)
    
    removed_names = prev_descriptors_names - current_descriptor_names
    added_names = current_descriptor_names - prev_descriptors_names

    # Report removals based on comparing with old registry.hpp state
    if removed_names:
        print(f"\nEffectively removing {len(removed_names)} deprecated or overridden descriptors from active use:")
        for desc_name in sorted(list(removed_names)):
            # This indicates descriptors that were in registry.hpp but are no longer in all_descriptors_for_codegen
            print(f"  - {desc_name}Descriptor")
        # Actual file cleanup for removed descriptors (if any) would follow here as in original script for .hpp and .cpp files

    print("\nAdding/Updating registration functions where missing (for C++)...")
    # Filter for C++ descriptors that are not Broto-Moreau for this step
    standard_cpp_descriptors = [d for d in all_descriptors_for_codegen if d['source_file'].endswith('.cpp') and not d.get('is_brotomoreau', False)]
    add_registration_functions_to_descriptor_files(standard_cpp_descriptors)

    add_registration_functions() # For observers in common.cpp
    
    # Separate C descriptors for C-specific codegen steps
    c_descriptors_to_codegen = [d for d in all_descriptors_for_codegen if d['source_file'].endswith('.c')]

    print("\nUpdating C registry header (cregistry.h)...")
    update_c_registry_header(c_descriptors_to_codegen) # Use the final list of C descriptors
    
    print("\nGenerating C wrapper file (c_wrappers.cpp)...")
    generate_c_wrapper_file(c_descriptors_to_codegen) # Use the final list of C descriptors

    print("\nGenerating main registry header (registry.hpp)...")
    generate_registry_hpp(all_descriptors_for_codegen, observers) 

    print("\nGenerating main registry implementation (registry.cpp)...")
    generate_registry_cpp(all_descriptors_for_codegen, observers)
    
    print("\nGenerating documentation (DESCRIPTORS.md)...")
    generate_descriptors_md(all_descriptors_for_codegen)
    
    print("\nSummary of changes to active descriptors:")
    if added_names:
        print(f"Added or newly prioritized descriptors ({len(added_names)}):")
        for desc_name in sorted(list(added_names)):
            print(f"  + {desc_name}Descriptor")
    # 'removed_names' already reported above more accurately reflects what's no longer in the build.

    final_active_descriptor_count = len(all_descriptors_for_codegen)
    print(f"\nFinal count of unique, valid, and prioritized descriptors for code generation: {final_active_descriptor_count}")
    # This count should now align with the runtime count (e.g., 446).
    
    # Calculate how many were stubs/duplicates not part of the final set
    total_initially_found = len(cpp_descriptor_dicts) + len(c_descriptor_dicts)
    total_skipped_or_overridden = total_initially_found - final_active_descriptor_count
    if total_skipped_or_overridden > 0:
         print(f"(Filtered out or superseded {total_skipped_or_overridden} initial candidates due to stubs, duplication, or prioritization)")

    print("\nVerifying registry files...")
    verify_registry_files(all_descriptors_for_codegen)

    print("\nRegistry update complete.")

if __name__ == "__main__":
    main()