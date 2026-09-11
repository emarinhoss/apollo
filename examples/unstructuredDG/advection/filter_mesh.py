#!/usr/bin/env python3
"""
Filter Gmsh file to remove boundary line elements, keeping only triangular cells.
This prevents DMPlex from importing line elements as cells.
"""

import sys

def filter_gmsh_mesh(input_file, output_file):
    """Remove type-1 (line) elements from Gmsh file, keep only type-2 (triangles)."""

    with open(input_file, 'r') as f:
        lines = f.readlines()

    output_lines = []
    in_elements = False
    triangle_count = 0
    line_count = 0

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line == '$Elements':
            in_elements = True
            output_lines.append(lines[i])
            i += 1
            # Read element count
            element_count = int(lines[i].strip())
            # We'll update this later
            count_line_index = len(output_lines)
            output_lines.append('')  # Placeholder for count
            i += 1
            continue

        if line == '$EndElements':
            in_elements = False
            # Update element count
            output_lines[count_line_index] = f"{triangle_count}\n"
            output_lines.append(lines[i])
            i += 1
            continue

        if in_elements:
            parts = line.split()
            if len(parts) >= 2:
                element_type = int(parts[1])
                if element_type == 1:  # Line element
                    line_count += 1
                    i += 1
                    continue
                elif element_type == 2:  # Triangle element
                    triangle_count += 1
                    output_lines.append(lines[i])
                    i += 1
                    continue
            # Other element types - keep them
            output_lines.append(lines[i])
            i += 1
            continue

        # Not in elements section - keep as is
        output_lines.append(lines[i])
        i += 1

    with open(output_file, 'w') as f:
        f.writelines(output_lines)

    print(f"Filtered mesh: removed {line_count} line elements, kept {triangle_count} triangular elements")
    return triangle_count

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python filter_mesh.py <input.msh> [output.msh]")
        print("If output is not specified, creates <input>_filtered.msh")
        sys.exit(1)

    input_mesh = sys.argv[1]
    if len(sys.argv) > 2:
        output_mesh = sys.argv[2]
    else:
        output_mesh = input_mesh.replace('.msh', '_filtered.msh')

    filter_gmsh_mesh(input_mesh, output_mesh)
    print(f"Filtered mesh written to: {output_mesh}")
