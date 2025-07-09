#!/usr/bin/env python3
"""Comprehensive script to fix all remaining flake8 issues in parallel."""

import re
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


def fix_unused_imports_comprehensive(file_path):
    """Remove all unused imports."""
    print(f"Fixing unused imports in {file_path}")
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    for line in lines:
        # Remove specific unused imports
        if (
            line.strip() == "import gzip" or
            line.strip() == "from collections import namedtuple" or
            line.strip() == "from typing import Set, Union" or
            line.strip() == "from typing import Set" or
            line.strip() == "from typing import Union" or
            line.strip() == "from concurrent.futures import as_completed"
        ):
            continue
        new_lines.append(line)
    
    with open(file_path, 'w') as f:
        f.writelines(new_lines)


def fix_fstring_comprehensive(file_path):
    """Fix all f-strings without placeholders."""
    print(f"Fixing f-strings in {file_path}")
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Replace f-strings without placeholders
    content = re.sub(r'f"([^"]*?)"', lambda m: f'"{m.group(1)}"' if '{' not in m.group(1) else m.group(0), content)
    content = re.sub(r"f'([^']*?)'", lambda m: f"'{m.group(1)}'" if '{' not in m.group(1) else m.group(0), content)
    
    with open(file_path, 'w') as f:
        f.write(content)


def fix_line_length_manual(file_path):
    """Fix some line length issues manually."""
    print(f"Fixing line length in {file_path}")
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    for line in lines:
        # Remove trailing whitespace
        line = line.rstrip() + '\n' if line.endswith('\n') else line.rstrip()
        
        # Convert blank lines with whitespace to empty lines
        if line.strip() == '':
            line = '\n' if line.endswith('\n') else ''
        
        new_lines.append(line)
    
    with open(file_path, 'w') as f:
        f.writelines(new_lines)


def process_file(file_path):
    """Process a single file to fix issues."""
    if not Path(file_path).exists():
        return f"File not found: {file_path}"
    
    try:
        fix_unused_imports_comprehensive(file_path)
        fix_fstring_comprehensive(file_path)
        fix_line_length_manual(file_path)
        return f"Fixed: {file_path}"
    except Exception as e:
        return f"Error processing {file_path}: {e}"


def main():
    """Fix all files in parallel."""
    files_to_fix = [
        "./tests/fixtures/geneburden/generate_enhanced_test_data.py",
        "./tests/fixtures/giab/call_variants_gatk.py", 
        "./tests/fixtures/giab/download_giab_bams.py",
        "./tests/integration/test_parallel_processing_integration.py",
        "./tests/test_gene_burden_comprehensive.py",
        "./tests/test_gene_burden_integration.py",
        "./variantcentrifuge/gene_burden.py",
    ]
    
    with ThreadPoolExecutor(max_workers=len(files_to_fix)) as executor:
        futures = {executor.submit(process_file, file_path): file_path for file_path in files_to_fix}
        
        for future in as_completed(futures):
            result = future.result()
            print(result)


if __name__ == "__main__":
    main()