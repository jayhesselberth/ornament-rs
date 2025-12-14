#!/bin/bash
# Download MODOMICS modification data
# Run from the data/ directory or provide output path

set -e

OUTPUT_DIR="${1:-.}"

echo "Downloading MODOMICS modification database..."

# Download modifications list (JSON format)
curl -s "https://genesilico.pl/modomics/api/modifications" \
    -o "${OUTPUT_DIR}/modomics_modifications.json"

echo "Downloaded modifications to ${OUTPUT_DIR}/modomics_modifications.json"

# Show summary
if command -v python3 &> /dev/null; then
    python3 -c "
import json
with open('${OUTPUT_DIR}/modomics_modifications.json') as f:
    data = json.load(f)
print(f'Total modifications: {len(data)}')

# Group by reference base
by_base = {}
for k, v in data.items():
    bases = tuple(v.get('reference_moiety', ['unknown']))
    by_base[bases] = by_base.get(bases, 0) + 1
print('By parent base:')
for base, count in sorted(by_base.items()):
    print(f'  {base}: {count}')
"
fi

echo ""
echo "MODOMICS data downloaded successfully!"
echo ""
echo "The JSON file contains modification definitions with:"
echo "  - name, short_name: Modification names"
echo "  - reference_moiety: Parent base (A, C, G, U)"
echo "  - formula, mass_*: Chemical properties"
echo "  - smile: Chemical structure"
echo "  - new_abbrev: MODOMICS unicode character"
echo ""
echo "Note: Position-specific expectations (which modifications at which"
echo "tRNA positions) are currently hardcoded in the Rust database based"
echo "on literature. To update those, edit:"
echo "  crates/ornament-core/src/modification/database.rs"
