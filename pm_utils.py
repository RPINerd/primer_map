"""A few miscellaneous scripts needed for primer map to function"""

import re

from pm_classes import Match, MatchCollection, Primer, RestrictionSite, RestrictionSiteCollection

TRANSTABLE = str.maketrans("ATGCatcguUryRYkmKMbvBVdhDH", "TACGtagcaAyrYRmkMKvbVBhdHD")


def check_genetic_code(patterns: list[str]) -> bool:

    codon = ""
    one_match = False
    test_sequence = """
        gggggaggtggcgaggaagatgacgtggtagttgtcgcggcagctgccaggagaag\
        tagcaagaaaaataacatgataattatcacgacaactacctggtgatgttgctagt\
        aatattacttgttatttttctcgtcatcttcccggcgacgtcgccagcaacatcac\
        ctgctacttctcccgccacctccc\
        """
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


def complement(sequence: str) -> str:
    return sequence.translate(TRANSTABLE)


# TODO switch this to use a maketrans dict
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


def find_matches(primers: list[Primer], sequence: str, topology: str, is_reverse: bool) -> MatchCollection:
    match_collection = MatchCollection()
    original_length = len(sequence)

    if topology == "circular":
        look_ahead = 50
        shift_value = len(sequence[:look_ahead])
        upper_limit = len(sequence) + shift_value
        sequence = sequence[-look_ahead:] + sequence + sequence[:look_ahead]
        lower_limit = shift_value

        for primer in primers:
            re_pattern = primer.regex
            for match in re.finditer(re_pattern, sequence):
                match_position = match.end()
                if lower_limit <= match_position < upper_limit:
                    match_position -= shift_value
                    if match_position == 0:
                        match_position = original_length
                    match_collection.add_match(Match(primer["name"], match.group(0), match_position))
                    if is_reverse:
                        primer["has_reverse_match"] = True
                    else:
                        primer["has_forward_match"] = True
    else:
        for primer in primers:
            re_pattern = primer.regex
            for match in re.finditer(re_pattern, sequence):
                match_position = match.end()
                match_collection.add_match(Match(primer["name"], match.group(0), match_position))
                if is_reverse:
                    primer["has_reverse_match"] = True
                else:
                    primer["has_forward_match"] = True

    return match_collection


def find_restriction_sites(sequence: str, array_of_items: list[str], topology: str) -> RestrictionSiteCollection:
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

    if topology == "circular":
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


def reverse(sequence: str) -> str:
    return sequence[::-1]


def right_num(num: int, sequence: str, col_length: int, tab: str) -> str:
    temp_string = ""
    num_str = str(num)
    for _ in range(len(num_str), col_length):
        temp_string += " "
    num_str = temp_string + num_str + " "
    sequence += num_str + tab
    return sequence


def validate_sequence(sequence: str) -> str:
    """
    Validates that the provided target sequence is actually DNA, and less than 200,000,000 characters long.

    :param sequence: The target sequence to validate.
    :type sequence: str
    :return sequence: The validated sequence.
    :rtype: str
    :raises ValueError: If the sequence is empty, too long, or contains invalid characters.
    """

    sequence = sequence.replace("\n", "").replace("\r", "").replace(" ", "")

    if len(sequence) == 0:
        raise ValueError("Sequence is empty!")
    if len(sequence) > 200000000:
        raise ValueError(f"Sequence is too long! Limit is 200000000 characters (recieved {len(sequence)}).")
    if re.search(r"[^gatucryswkmbdhvnxGATUCRYSWKMBDHVNX]", sequence):
        raise ValueError("Sequence contains invalid characters!")

    return sequence


def write_restriction_sites(sequence: str, array_of_items: list, topology: str) -> None:
    # result_array = []
    look_ahead = 50
    lower_limit = 0
    upper_limit = len(sequence)
    shift_value = 0
    cut_distance = 0
    match_exp = None
    match_position = 0
    temp_string = ""
    background_class = ""
    # match_array = None
    times_found = 0

    if topology == "circular":
        shift_value = len(sequence[:look_ahead])
        sequence = sequence[-look_ahead:] + sequence + sequence[:look_ahead]
        lower_limit = 0 + shift_value
        upper_limit = upper_limit + shift_value

    print('<table border="1" width="100%" cellspacing="0" cellpadding="2"><tbody>')
    print('<tr><td class="title" width="200px">Site:</td><td class="title">Positions:</td></tr>')

    for item in array_of_items:
        temp_string = "none"
        background_class = "many"
        match_exp = re.compile(item.split("/")[1], re.IGNORECASE)
        cut_distance = int(re.search(r"\)\D*(\d+)", item).group(1))

        for match in match_exp.finditer(sequence):
            match_position = match.end() - cut_distance
            if lower_limit <= match_position < upper_limit:
                times_found += 1
                temp_string += ", " + str(match_position - shift_value + 1)

            if re.search(r"\d", temp_string):
                temp_string = temp_string.replace("none, ", "")

            if times_found == 0:
                background_class = "none"
            elif times_found == 1:
                background_class = "one"
            elif times_found == 2:
                background_class = "two"
            elif times_found == 3:
                background_class = "three"
            else:
                background_class = "many"

        site_name = re.search(r"\([^\(]+\)", item).group(0).replace("(", "").replace(")", "")
        print(
            f'<tr><td class="{background_class}">{site_name}</td><td class="{background_class}">{temp_string}</td></tr>'
        )

        times_found = 0

    print("</tbody></table>")
