#!/usr/bin/env python3
import sys

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <input_file>", file=sys.stderr)
    sys.exit(1)

input_file = sys.argv[1]
seen = {}

with open(input_file, 'r') as f:
    for line in f:
        line = line.rstrip('\n')
        if not line:
            continue
        fields = line.split('\t')
        original_id = fields[0]

        if original_id not in seen:
            seen[original_id] = 1
            print(line)
        else:
            new_id = f"{original_id}_{seen[original_id]}"
            seen[original_id] += 1
            fields[0] = new_id
            print('\t'.join(fields))