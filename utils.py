"""A few miscellaneous scripts needed for primer map to function"""

# TODO this is only roughly translated, need to validate the conversion
import re


def check_genetic_code(array_of_patterns):
    z = 0
    codon = ""
    one_match = False
    test_sequence = "gggggaggtggcgaggaagatgacgtggtagttgtcgcggcagctgccaggagaagtagcaagaaaaataacatgataattatcacgacaactacctggtgatgttgctagtaatattacttgttatttttctcgtcatcttcccggcgacgtcgccagcaacatcacctgctacttctcccgccacctccc"
    while z < len(array_of_patterns):
        if not re.search(r"^\s*/[a-zA-Z\|\[\]]+/=[a-zA-Z\*]", array_of_patterns[z]):
            print("Genetic code error: one or more patterns have been entered incorrectly.")
            return False
        if not more_expression_check(array_of_patterns[z]):
            print("Genetic code error: one or more patterns have been entered incorrectly.")
            return False
        z += 1

    genetic_code_match_result = [None] * len(array_of_patterns)
    genetic_code_match_exp = [None] * len(array_of_patterns)
    for j in range(len(array_of_patterns)):
        genetic_code_match_exp[j] = re.compile(array_of_patterns[j].split("=")[0][1:-1], re.IGNORECASE)
        genetic_code_match_result[j] = array_of_patterns[j].split("=")[1]

    for i in range(0, len(test_sequence) - 2, 3):
        codon = test_sequence[i : i + 3]
        for j in range(len(genetic_code_match_exp)):
            if genetic_code_match_exp[j].search(codon):
                if one_match:
                    print(f"Genetic code error: more than one amino acid is coded by the codon: {codon}.")
                    return False
                one_match = True
        if not one_match:
            print("The genetic code expressions are missing a codon.")
            return False
        one_match = False

    return True


def check_rest_patterns(array_of_patterns):
    z = 0
    while z < len(array_of_patterns):
        if re.search(r"^\s*/[acgturyswkmbdhvn\[\]]+/\s+\([^/]+\)\d+", array_of_patterns[z], re.IGNORECASE) is None:
            print("One or more patterns have been entered incorrectly.")
            return False
        if not more_expression_check(array_of_patterns[z]):
            print("One or more patterns have been entered incorrectly.")
            return False
        z += 1
    return True


def convert_degenerates(sequence: str) -> str:
    sequence = sequence.lower()
    sequence = sequence.replace("t", "[TU]")
    sequence = sequence.replace("r", "[AGR]")
    sequence = sequence.replace("y", "[CTUY]")
    sequence = sequence.replace("s", "[GCS]")
    sequence = sequence.replace("w", "[ATUW]")
    sequence = sequence.replace("k", "[GTUK]")
    sequence = sequence.replace("m", "[ACM]")
    sequence = sequence.replace("b", "[CGTUBSKY]")
    sequence = sequence.replace("d", "[AGTUDRKW]")
    sequence = sequence.replace("h", "[ACTUHMYW]")
    sequence = sequence.replace("v", "[ACGVSMR]")
    sequence = sequence.replace("n", "[ACGTURYSWKMBDHVN]")
    return sequence


def more_expression_check(pattern):
    if (
        re.search(r"\[[A-Za-z\|]*\[", pattern)
        or re.search(r"\][A-Za-z\|]*\]", pattern)
        or re.search(r"\[\]", pattern)
        or re.search(r"/[A-Za-z\|]*\]", pattern)
        or re.search(r"\[[A-Za-z\|]*/", pattern)
        or re.search(r"\|\|", pattern)
        or re.search(r"/\|", pattern)
        or re.search(r"\|/", pattern)
        or re.search(r"\[.\]", pattern)
        or re.search(r"<", pattern)
        or re.search(r">", pattern)
    ):
        return False
    return True
