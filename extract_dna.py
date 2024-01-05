import re
import sys

def extract_dna_sequence(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Use regular expression to extract the DNA sequence
            match = re.search(r'([ACGT]+),', line)
            if match:
                dna_sequence = match.group(1)
                outfile.write(f'{dna_sequence}\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_dna.py input_file.txt output_file.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    extract_dna_sequence(input_file, output_file)
    print(f"DNA sequences extracted from {input_file} and saved to {output_file}.")
