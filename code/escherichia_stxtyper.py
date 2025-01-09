"""
Copyright 2025 Mary Maranga (gathonimaranga@gmailcom), Kat Holt
https://github.com/klebgenomics/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""





import os
import sys
import shutil
from pathlib import Path
import subprocess


def description():
    return 'Shiga toxin (stx) gene typing using STXTyper'


def prerequisite_modules():
    return []


def get_headers():
    """
    Define the headers for STXTyper results.
    """
    full_headers = [
        'Assembly', 'target_contig', 'stx_type', 'operon', 'identity', 
        'target_start', 'target_stop', 'target_strand', 
        'A_reference', 'A_identity', 'A_reference_subtype', 'A_coverage',
        'B_reference', 'B_reference_subtype', 'B_identity', 'B_coverage'
    ]
    return full_headers


def add_cli_options(parser):
    """
    Add command-line options for this module.
    """
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('-t', '--threads', type=int, default=8, metavar='',
                       help="Number of threads for STXTyper (default: %(default)s).")
    group.add_argument('-q', '--quiet', action='store_true', default=False,
                       help="Suppress additional STXTyper output (default: %(default)s).")
    return group


def check_cli_options(args):
    """
    Validate the command-line arguments.
    """
    if args.threads < 1:
        raise ValueError("The number of threads must be at least 1.")
    if not shutil.which('stxtyper'):
        sys.exit('Error: STXTyper is not installed or not in PATH.')


def check_external_programs():
    """
    Ensure the required external programs are available.
    """
    if not shutil.which('stxtyper'):
        sys.exit('Error: could not find STXTyper executable.')
    return ['stxtyper']


def run_stxtyper(nucleotide_file, threads, quiet=False):
    """
    Run STXTyper and return the output directly as a string.
    """
    cmd = [
        'stxtyper',
        '-n', nucleotide_file,
        '-o', '/dev/stdout',  # Redirect output to stdout for capture
        '--threads', str(threads)
    ]
    if quiet:
        cmd.append('-q')

    try:
        # Capture the output of STXTyper directly
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error running STXTyper: {e}")


def get_results(assemblies, args, previous_results):
    """
    Execute STXTyper for multiple assemblies and combine results into a single dictionary.
    Returns a dictionary of concatenated results.
    """
    full_headers = get_headers()
    combined_results = {}

    for assembly_path in assemblies:
        assembly_path = Path(assembly_path)

        # Run STXTyper and capture the output
        output = run_stxtyper(
            nucleotide_file=assembly_path,
            threads=args.threads,
            quiet=args.quiet
        )

        # Parse output into rows and add the Assembly identifier
        for line in output.splitlines():
            if line.startswith('#') or not line.strip():  # Skip header or empty lines
                continue
            parts = line.strip().split('\t')
            row_dict = {key: value for key, value in zip(full_headers[1:], parts)}  # Exclude 'Assembly'
            row_dict['Assembly'] = assembly_path.stem  # Add assembly identifier
            
            # Use a unique key for each row (e.g., Assembly + Contig + StxType)
            key = f"{assembly_path.stem}_{row_dict['target_contig']}_{row_dict['stx_type']}"
            combined_results[key] = row_dict

    return combined_results
