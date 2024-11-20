"""
    Primer Map | RPINerd, 11/18/24

    Primer Map is a sub-tool of the Sequence Manipulation Suite: A collection of simple JavaScript programs
    for generating, formatting, and analyzing short DNA and protein sequences.
    Original Source Author: Paul Stothard stothard@ualberta.ca

    Translated here is just the primer map function from the original source code to be able to run
    standalone or as a module.
"""

import argparse
import re

from pm_classes import Primer
from pm_utils import (
    check_genetic_code,
    check_rest_patterns,
    complement,
    convert_degenerates,
    find_matches,
    find_restriction_sites,
    reverse,
    validate_sequence,
    write_restriction_sites,
)
from sms_codes import get_genetic_code_string
from sms_restriction import get_restriction_sites


def parse_args() -> argparse.ArgumentParser:
    """
    Simple argument parser for the primer map function.

    :param None:
    :return: args
    :rtype: argparse.Namespace
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--primers", type=str, required=True, help="Comma separated list of primers to map.")
    parser.add_argument(
        "-t", "--target", type=str, required=True, help="Reference sequence to map the primers against."
    )
    parser.add_argument(
        "-g",
        "--genetic-code",
        type=str,
        default="standard",
        required=False,
        help="!NOT IMPLEMENTED! Organism to use as the reference for codon translation.",
    )
    parser.add_argument(
        "-r",
        "--rs-sites",
        action="store_true",
        required=False,
        help="!NOT IMPLEMENTED! Whether to show restriction sites or not.",
    )
    parser.add_argument(
        "--topology", type=str, default="linear", required=False, help="Whether the sequence is linear or circular"
    )
    parser.add_argument(
        "--bp-per-line", type=int, default=100, required=False, help="Number of base pairs to display per line."
    )
    parser.add_argument(
        "--reading-frame", type=int, default=0, required=False, help="!NOT IMPLEMENTED! Reading frame(s) to display."
    )
    parser.add_argument("-c", "--color", action="store_true", help="!NOT IMPLEMENTED! Color output")

    return parser.parse_args()


def primer_map(args: argparse.Namespace) -> None:
    """
    Primer Map driver function

    :param args: Arguments provided by command line
    :type args: argparse.Namespace
    :return: None
    """
    try:
        fasta_seq = validate_sequence(args.target)
    except ValueError as e:
        raise ValueError(f"Invalid sequence: {e}")
    fasta_title = ""
    rs_collection = None
    forward_matches = None
    reverse_matches = None
    # is_color = args.color

    genetic_code = get_genetic_code_string(args.genetic_code)
    restriction_sites = get_restriction_sites()

    # TODO raise as error
    if not check_genetic_code(genetic_code):
        return
    if not check_rest_patterns(restriction_sites):
        return

    new_primers = []
    re_pattern = re.compile(r"\(([^\(]+)\)\s*([A-Za-z]+)")
    for primer in args.primers.split(","):
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

    print(
        f'Results for {args.topology} {len(fasta_seq)} residue sequence "{fasta_title}" starting "{fasta_seq[:10]}"\n'
    )

    if args.rs_sites:
        rs_collection = find_restriction_sites(
            fasta_seq,
            restriction_sites,
            args.topology,
        )
        rs_collection.sort_sites()

    forward_matches = find_matches(
        new_primers,
        fasta_seq,
        args.topology,
        False,
    )
    reverse_matches = find_matches(
        new_primers,
        reverse(complement(fasta_seq)),
        args.topology,
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
        rs_collection,
        forward_matches,
        reverse_matches,
        args.bp_per_line,
        args.reading_frame,
    )

    output_window.args.write("\n")

    if args.rs_sites:
        write_restriction_sites(
            fasta_seq,
            restriction_sites,
            args.topology,
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
    primer_map(args)
