import argparse
import sys
import textwrap

CODON_TABLE = {
    "UUU": "Phe", "UCU": "Ser", "UAU": "Tyr",  "UGU": "Cys",
    "UUC": "Phe", "UCC": "Ser", "UAC": "Tyr",  "UGC": "Cys",
    "UUA": "Leu", "UCA": "Ser", "UAA": "Stop", "UGA": "Stop",
    "UUG": "Leu", "UCG": "Ser", "UAG": "Stop", "UGG": "Trp",

    "CUU": "Leu", "CCU": "Pro", "CAU": "His",  "CGU": "Arg",
    "CUC": "Leu", "CCC": "Pro", "CAC": "His",  "CGC": "Arg",
    "CUA": "Leu", "CCA": "Pro", "CAA": "Gln",  "CGA": "Arg",
    "CUG": "Leu", "CCG": "Pro", "CAG": "Gln",  "CGG": "Arg",

    "AUU": "Ile", "ACU": "Thr", "AAU": "Asn",  "AGU": "Ser",
    "AUC": "Ile", "ACC": "Thr", "AAC": "Asn",  "AGC": "Ser",
    "AUA": "Ile", "ACA": "Thr", "AAA": "Lys",  "AGA": "Arg",
    "AUG": "Met", "ACG": "Thr", "AAG": "Lys",  "AGG": "Arg",

    "GUU": "Val", "GCU": "Ala", "GAU": "Asp",  "GGU": "Gly",
    "GUC": "Val", "GCC": "Ala", "GAC": "Asp",  "GGC": "Gly",
    "GUA": "Val", "GCA": "Ala", "GAA": "Glu",  "GGA": "Gly",
    "GUG": "Val", "GCG": "Ala", "GAG": "Glu",  "GGG": "Gly",
}

VALID_BASES = set("ACGU")
LABEL_WIDTH = 12
INDENT = " " * (2 + LABEL_WIDTH + 3)


def read_input(args) -> str:
    if args.file:
        with open(args.file, encoding="utf-8") as f:
            return f.read()
    if args.sequence:
        return " ".join(args.sequence)
    if not sys.stdin.isatty():
        return sys.stdin.read()
    print("Paste the mRNA sequence, then press Enter on an empty line:")
    lines = []
    while True:
        try:
            line = input()
        except EOFError:
            break
        if not line.strip():
            break
        lines.append(line)
    return "\n".join(lines)


def clean_sequence(raw: str) -> str:
    lines = [line for line in raw.splitlines() if not line.startswith(">")]
    text = "".join(lines).upper()
    return "".join(ch for ch in text if not ch.isspace() and not ch.isdigit())


def validate(seq: str) -> None:
    if not seq:
        raise ValueError("The sequence is empty.")
    if "T" in seq:
        raise ValueError(
            "Found 'T' in the sequence. This app expects mRNA, which uses U instead of T. "
            "If this is DNA from the coding strand, replace every T with U first."
        )
    bad = sorted(set(seq) - VALID_BASES)
    if bad:
        raise ValueError(
            f"Invalid character(s): {', '.join(bad)}. mRNA may only contain A, C, G and U."
        )
    if len(seq) < 3:
        raise ValueError("The sequence is shorter than one codon (3 bases).")


def read_from_start(seq: str, start: int) -> tuple[list[str], bool]:
    codons = []
    for i in range(start, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        codons.append(codon)
        if CODON_TABLE[codon] == "Stop":
            return codons, True
    return codons, False


def field(label: str, value: str) -> None:
    print(f"  {label:<{LABEL_WIDTH}} : {value}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find the closest AUG in an mRNA sequence and translate it "
                    "up to the stop codon."
    )
    parser.add_argument("sequence", nargs="*",
                        help="mRNA sequence, e.g. GGAUGUUUGGCUAA (spaces allowed)")
    parser.add_argument("-f", "--file", help="read the sequence from a text or FASTA file")
    parser.add_argument("-c", "--codons", action="store_true",
                        help="also print the codon-by-codon translation")
    args = parser.parse_args()

    if args.file and args.sequence:
        parser.error("give either a sequence or --file, not both")

    try:
        seq = clean_sequence(read_input(args))
        validate(seq)
    except (OSError, ValueError) as err:
        print(f"Error: {err}", file=sys.stderr)
        sys.exit(1)

    print()
    field("mRNA length", f"{len(seq)} bases")

    start = seq.find("AUG")
    if start == -1:
        field("AUG", "none - there is no start codon in the sequence")
        return

    codons, has_stop = read_from_start(seq, start)
    protein = [CODON_TABLE[codon] for codon in codons]
    amino_acids = len(protein) - 1 if has_stop else len(protein)

    field("AUG at base", f"{start + 1} (frame {start % 3 + 1})")
    if has_stop:
        stop_base = start + 3 * (len(codons) - 1) + 1
        field("Stop", f"{codons[-1]} at base {stop_base}")
    else:
        field("Stop", "none - the sequence ends before a stop codon")
    field("Amino acids", str(amino_acids))

    print(textwrap.fill("-".join(protein), width=80,
                        initial_indent=f"  {'Protein':<{LABEL_WIDTH}} : ",
                        subsequent_indent=INDENT))

    if args.codons:
        print("\n     #   Base  Codon  Amino acid")
        for i, (codon, aa) in enumerate(zip(codons, protein)):
            print(f"  {i + 1:4}  {start + 3 * i + 1:5}  {codon:5}  {aa}")


if __name__ == "__main__":
    main()
