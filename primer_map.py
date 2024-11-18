"""
    Primer Map | RPINerd, 11/18/24

    Primer Map is a sub-tool of the Sequence Manipulation Suite: A collection of simple JavaScript programs
    for generating, formatting, and analyzing short DNA and protein sequences.
    Original Source Author: Paul Stothard stothard@ualberta.ca

    Translated here is just the primer map function from the original source code to be able to run
    standalone or as a module.
"""

import argparse
import os
import re
import sys

from sms_codes import get_genetic_code_string
from sms_restriction import get_restriction_sites
from utils import check_genetic_code, check_rest_patterns, convert_degenerates


class Primer:
    def __init__(self, sequence, regex, name):
        self.sequence = sequence
        self.regex = regex
        self.name = name
        self.hasForwardMatch = False
        self.hasReverseMatch = False


def parse_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--target", type=str, required=True, help="Target sequence")
    parser.add_argument("-g", "--genetic_code", type=str, default="standard", required=False, help="Genetic code")
    parser.add_argument("-c", "--color", action="store_true", help="Color output")
    args = parser.parse_args()
    return args


def primer_map(args: argparse.Namespace) -> None:
    fasta_seq = ""
    title = ""
    restriction_site_collection = None
    forward_matches = None
    reverse_matches = None
    is_color = args.color

    genetic_code = get_genetic_code_string(args.genetic_code)
    restriction_sites = get_restriction_sites()

    # TODO raise as error
    if not check_genetic_code(genetic_code):
        return
    if not check_rest_patterns(restriction_sites):
        return

    primers = args.primers.split(",")
    new_primers = []
    re_pattern = re.compile(r"\(([^\(]+)\)\s*([A-Za-z]+)")
    for primer in primers:
        match_array = re_pattern.search(primer)
        if match_array:
            primer_name = match_array.group(1)
            if len(match_array.group(2)) < 10:
                print("Please enter primer sequences that are at least 10 bases long.")
                return
            primer_seq = match_array.group(2)
            primer_re = re.compile(convert_degenerates(primer_seq), re.IGNORECASE)
            new_primers.append(Primer(primer_seq, primer_re, primer_name))

    # Launch an output window to display mappings on
    """
    open_window("Primer Map", is_color)
    if args.forms[0].elements[8].options[args.forms[0].elements[8].selectedIndex].value == "shown":
        output_window.args.write('<span class="one">cuts once</span><br />\n')
        output_window.args.write('<span class="two">cuts twice</span><br />\n')
        output_window.args.write("\n")
        output_window.args.write('<span class="forward_primer">&gt;&gt;&gt;forward primer</span><br />\n')
        output_window.args.write('<span class="reverse_primer">&lt;&lt;&lt;reverse primer</span><br />\n')
        output_window.args.write("\n")
    """

    fasta = args.target
    # fasta_seq = remove_non_dna(fasta)
    open_pre()
    output_window.args.write(
        get_info_from_title_and_sequence_and_topology(
            title, fasta_seq, args.forms[0].elements[9].options[args.forms[0].elements[9].selectedIndex].value
        )
    )

    if args.forms[0].elements[8].options[args.forms[0].elements[8].selectedIndex].value == "shown":
        restriction_site_collection = find_restriction_sites(
            fasta_seq,
            restriction_sites,
            args.forms[0].elements[9].options[args.forms[0].elements[9].selectedIndex].value,
        )
        restriction_site_collection.sort_restriction_sites()

    forward_matches = find_matches(
        new_primers,
        fasta_seq,
        args.forms[0].elements[9].options[args.forms[0].elements[9].selectedIndex].value,
        False,
    )
    reverse_matches = find_matches(
        new_primers,
        reverse(complement(fasta_seq)),
        args.forms[0].elements[9].options[args.forms[0].elements[9].selectedIndex].value,
        True,
    )

    for match in forward_matches.matches:
        match.position -= len(match.matching_text)
        match.end = match.position + len(match.matching_text)
        if match.position < 0:
            match.position += len(fasta_seq)
        if match.end > len(fasta_seq):
            match.end -= len(fasta_seq)

    for match in reverse_matches.matches:
        match.position = len(fasta_seq) - match.position
        match.end = match.position + len(match.matching_text)
        if match.position < 0:
            match.position += len(fasta_seq)
        if match.end > len(fasta_seq):
            match.end -= len(fasta_seq)

    forward_matches.sort_matches()
    reverse_matches.sort_matches()

    layout_primer_map(
        fasta_seq,
        genetic_code,
        restriction_site_collection,
        forward_matches,
        reverse_matches,
        args.forms[0].elements[5].options[args.forms[0].elements[5].selectedIndex].value,
        args.forms[0].elements[6].options[args.forms[0].elements[6].selectedIndex].value,
    )

    output_window.args.write("\n")
    close_pre()

    if args.forms[0].elements[8].options[args.forms[0].elements[8].selectedIndex].value == "shown":
        write_restriction_sites(
            fasta_seq,
            restriction_sites,
            args.forms[0].elements[9].options[args.forms[0].elements[9].selectedIndex].value,
        )
        output_window.args.write("<br />\n")

    write_primer_sites(new_primers)

    output_window.args.write("<br />\n<br />\n")

    for primer in new_primers:
        primer.has_forward_match = False
        primer.has_reverse_match = False

    close_window()
    return


if __name__ == "__main__":
    args = parse_args()

    assert len(args.target) <= 200000000, f"Input sequence is too long, max length is 200000000, got {len(args.target)}"

    assert [x for x in args.target if x not in "ACGTatcg"], "Input sequence contains invalid characters"

    primer_map(args)
