#!/usr/bin/env python3

import sys
import os
import re
from pathlib import Path

def load_mapping(mapping_file):
“”“Load the mapping from long IDs to short names”””
mapping = {}
with open(mapping_file, ‘r’) as f:
for line in f:
line = line.strip()
if line:
parts = line.split(’\t’)
if len(parts) >= 2:
long_id = parts[0]
short_name = parts[1]
mapping[long_id] = short_name
return mapping

def extract_sample_id(file_path):
“”“Extract sample ID from various file path formats”””
# Get the filename without path
filename = os.path.basename(file_path)

```
# Remove common suffixes like _sorted_dedup.bam, .bam, .vcf, etc.
filename = re.sub(r'(_sorted)?(_dedup)?(_final)?\.(bam|vcf|gz)$', '', filename)

# Look for patterns like SAMEA followed by numbers
match = re.search(r'(SAMEA\d+)', filename)
if match:
    return match.group(1)

# If no SAMEA pattern found, return the cleaned filename
return filename
```

def main():
if len(sys.argv) != 3:
print(“Usage: python script.py <mapping_file> <sample_list_file>”)
print(“mapping_file: tab-separated file with long_id\tshort_name”)
print(“sample_list_file: file with current sample names/paths”)
sys.exit(1)

```
mapping_file = sys.argv[1]
sample_list_file = sys.argv[2]

# Load the ID to short name mapping
id_mapping = load_mapping(mapping_file)

# Process the sample list
with open(sample_list_file, 'r') as f:
    for line in f:
        current_name = line.strip()
        if not current_name:
            continue
            
        # Extract the sample ID from the current name
        sample_id = extract_sample_id(current_name)
        
        # Find the corresponding short name
        if sample_id in id_mapping:
            short_name = id_mapping[sample_id]
            print(f"{current_name}\t{short_name}")
        else:
            print(f"Warning: No mapping found for {sample_id} (from {current_name})", file=sys.stderr)
```

if **name** == “**main**”:
main()s