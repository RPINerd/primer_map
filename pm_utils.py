"""A few miscellaneous scripts needed for primer map to function"""

import re

from pm_classes import RestrictionSite, RestrictionSiteCollection


def check_genetic_code(patterns: list[str]) -> bool:

    codon = ""
    one_match = False
    test_sequence = "gggggaggtggcgaggaagatgacgtggtagttgtcgcggcagctgccaggagaagtagcaagaaaaataacatgataattatcacgacaactacctggtgatgttgctagtaatattacttgttatttttctcgtcatcttcccggcgacgtcgccagcaacatcacctgctacttctcccgccacctccc"
    for pattern in patterns:
        if not re.search(r"^\s*/[a-zA-Z\|\[\]]+/=[a-zA-Z\*]", pattern):
            print("Genetic code error: one or more patterns have been entered incorrectly.")
            return False

        if not more_expression_check(pattern):
            print("Genetic code error: one or more patterns have been entered incorrectly.")
            return False

    genetic_code_match_result = [None] * len(patterns)
    genetic_code_match_exp = [None] * len(patterns)
    for j in range(len(patterns)):
        genetic_code_match_exp[j] = re.compile(patterns[j].split("=")[0][1:-1], re.IGNORECASE)
        genetic_code_match_result[j] = patterns[j].split("=")[1]

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


def check_rest_patterns(patterns: list[str]) -> bool:

    for pattern in patterns:

        if re.search(r"^\s*/[acgturyswkmbdhvn\[\]]+/\s+\([^/]+\)\d+", pattern, re.IGNORECASE) is None:
            print("One or more patterns have been entered incorrectly.")
            return False

        if not more_expression_check(pattern):
            print("One or more patterns have been entered incorrectly.")
            return False

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


def find_restriction_sites(
    sequence: str, array_of_items: list[str], dna_conformation: str
) -> RestrictionSiteCollection:
    look_ahead: int = 50
    lower_limit: int = 0
    upper_limit: int = len(sequence)
    shift_value: int = 0
    cut_distance: int = 0
    match_exp: re.Pattern = None
    match_position: int = 0
    label: str = ""
    times_found: int = 0
    temp_array: list[RestrictionSite] = []

    rs_collection = RestrictionSiteCollection()

    if dna_conformation == "circular":
        shift_value = len(sequence[:look_ahead])
        sequence = sequence[-look_ahead:] + sequence + sequence[:look_ahead]
        lower_limit = 0 + shift_value
        upper_limit += shift_value

    for item in array_of_items:
        match_exp = re.compile(item.split("/")[1], re.IGNORECASE)
        cut_distance = int(re.search(r"\)\D*(\d+)", item).group(1))
        label = re.search(r"\(([^\(]+)\)", item).group(1)

        for match in match_exp.finditer(sequence):
            match_position = match.start() - cut_distance
            if lower_limit <= match_position < upper_limit:
                times_found += 1
                temp_array.append(
                    RestrictionSite(f"{label} {match_position - shift_value + 1}", match_position - shift_value)
                )

        for site in temp_array:
            site.set_cut_count(times_found)
            rs_collection.add_site(site)
        times_found = 0
        temp_array = []

    return rs_collection


def more_expression_check(pattern: str) -> bool:
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
