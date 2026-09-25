#!/usr/bin/env python3
"""
Batch convert Python 2 scripts to Python 3 compatibility.
Handles common conversion patterns needed for Apollo scripts.
"""

import re
import sys
import os
from pathlib import Path

def fix_print_statements(content):
    """Convert print statements to print functions."""
    lines = content.split('\n')
    fixed_lines = []

    for line in lines:
        # Skip if already a print function or commented
        if re.match(r'^\s*#', line) or re.search(r'\bprint\s*\(', line):
            fixed_lines.append(line)
            continue

        # Match print statement: print <anything until end of line or comment>
        match = re.match(r'^(\s*)print\s+(.+?)(\s*#.*)?$', line)
        if match:
            indent = match.group(1)
            expr = match.group(2).rstrip()
            comment = match.group(3) or ''
            # Convert to print function
            fixed_line = f'{indent}print({expr}){comment}'
            fixed_lines.append(fixed_line)
        else:
            fixed_lines.append(line)

    return '\n'.join(fixed_lines)

def fix_exceptions(content):
    """Fix exception syntax."""
    # except Type, var: → except Type as var:
    content = re.sub(r'except\s+(\w+(?:\.\w+)*)\s*,\s*(\w+)\s*:', r'except \1 as \2:', content)

    # raise Exception, message → raise Exception(message)
    content = re.sub(r'raise\s+(\w+)\s*,\s*(.+)', r'raise \1(\2)', content)

    # Fix the specific error in visualize_data.py
    content = re.sub(
        r'raise\s+Exception\("(.+?)"\)\s*%\s*\(([^)]+)\)',
        r'raise Exception("\1" % (\2))',
        content
    )

    return content

def fix_string_module(content):
    """Replace string module functions with string methods."""
    replacements = {
        r'string\.split\(([^,]+)\)': r'\1.split()',
        r'string\.split\(([^,]+),\s*([^)]+)\)': r'\1.split(\2)',
        r'string\.strip\(([^)]+)\)': r'\1.strip()',
        r'string\.find\(([^,]+),\s*([^)]+)\)': r'\1.find(\2)',
        r'string\.replace\(([^,]+),\s*([^,]+),\s*([^)]+)\)': r'\1.replace(\2, \3)',
        r'string\.join\(([^,]+),\s*([^)]+)\)': r'\2.join(\1)',
        r'string\.lower\(([^)]+)\)': r'\1.lower()',
        r'string\.upper\(([^)]+)\)': r'\1.upper()',
    }

    for pattern, replacement in replacements.items():
        content = re.sub(pattern, replacement, content)

    return content

def fix_imports(content):
    """Fix deprecated imports."""
    # Remove compiler module imports (removed in Python 3.9)
    content = re.sub(r'import\s+compiler\s*\n', '', content)
    content = re.sub(r'from\s+compiler\s+import\s+.*\n', '', content)

    # Replace compiler.compile with built-in compile
    content = re.sub(r'compiler\.compile\(', 'compile(', content)

    return content

def fix_exec_syntax(content):
    """Fix exec statement syntax."""
    # exec code in env, env → exec(code, env, env)
    content = re.sub(r'exec\s+(.+?)\s+in\s+(.+?),\s*(.+?)(?:\n|$)',
                     r'exec(\1, \2, \3)\n', content)

    return content

def fix_time_module(content):
    """Fix deprecated time.clock()."""
    content = re.sub(r'time\.clock\(\)', 'time.perf_counter()', content)

    return content

def fix_escape_sequences(content):
    """Fix invalid escape sequences by using raw strings."""
    # Find regex patterns that should be raw strings
    lines = content.split('\n')
    fixed_lines = []

    for line in lines:
        # If line contains a string with \s, \d, etc. and looks like a regex
        if re.search(r'''["'].*?\\[sdwSDW]''', line) and not line.strip().startswith('#'):
            # Make it a raw string if not already
            line = re.sub(r'''(["'])(.*?\\[sdwSDW].*?)\1''', r'r\1\2\1', line)
        fixed_lines.append(line)

    return '\n'.join(fixed_lines)

def convert_file(filepath):
    """Convert a single Python file to Python 3."""
    print(f"Converting {filepath}...")

    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except Exception as e:
        print(f"  Error reading {filepath}: {e}")
        return False

    original_content = content

    # Apply all fixes
    content = fix_print_statements(content)
    content = fix_exceptions(content)
    content = fix_string_module(content)
    content = fix_imports(content)
    content = fix_exec_syntax(content)
    content = fix_time_module(content)
    content = fix_escape_sequences(content)

    # Only write if changed
    if content != original_content:
        try:
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"  ✓ Updated {filepath}")
            return True
        except Exception as e:
            print(f"  Error writing {filepath}: {e}")
            return False
    else:
        print(f"  - No changes needed for {filepath}")
        return None

def main():
    scripts_dir = Path(__file__).parent
    py_files = sorted(scripts_dir.glob('*.py'))

    # Exclude this script itself
    py_files = [f for f in py_files if f.name != 'convert_to_python3.py']

    print(f"Found {len(py_files)} Python files to convert\n")

    updated = 0
    unchanged = 0
    errors = 0

    for py_file in py_files:
        result = convert_file(py_file)
        if result is True:
            updated += 1
        elif result is False:
            errors += 1
        else:
            unchanged += 1

    print(f"\n{'='*60}")
    print(f"Conversion complete!")
    print(f"  Updated:   {updated} files")
    print(f"  Unchanged: {unchanged} files")
    print(f"  Errors:    {errors} files")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
