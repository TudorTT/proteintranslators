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


def split_codons(seq: str) -> tuple[list[str], str]:
    usable = len(seq) - len(seq) % 3
    codons = [seq[i:i + 3] for i in range(0, usable, 3)]
    return codons, seq[usable:]


def translate(codons: list[str]) -> list[str]:
    return [CODON_TABLE[codon] for codon in codons]


def print_result(seq: str, codons: list[str], leftover: str,
                 protein: list[str], show_codons: bool) -> None:
    print()
    print(f"mRNA length : {len(seq)} bases")
    print(f"Codons read : {len(codons)}")
    if leftover:
        print(f"Warning     : length is not a multiple of 3 - the last "
              f"{len(leftover)} base(s) '{leftover}' were ignored.")

    stop_positions = [i + 1 for i, aa in enumerate(protein) if aa == "Stop"]
    if stop_positions:
        print(f"Stop codons : at codon # {', '.join(map(str, stop_positions))}")
    else:
        print("Stop codons : none found")

    print("\nProtein :")
    print(textwrap.fill("-".join(protein), width=80))



    if show_codons:
        print("\n   #   Base  Codon  Amino acid")
        for i, (codon, aa) in enumerate(zip(codons, protein)):
            print(f"{i + 1:4}  {3 * i + 1:5}  {codon:5}  {aa}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Translate an mRNA coding sequence into an amino acid sequence."
    )
    parser.add_argument("sequence", nargs="*",
                        help="mRNA sequence, e.g. AUGUUUGGCUAA (spaces allowed)")
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

    codons, leftover = split_codons(seq)
    protein = translate(codons)
    print_result(seq, codons, leftover, protein, args.codons)


if __name__ == "__main__":
    main()