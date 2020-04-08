#! /usr/bin/env python3

import sys
import re

# Crate main function
# import argparse, find_orf and translate
def main():
    import argparse
    import find_orf
    import translate

    # Create a command-line parser object
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    # Tell the parser what command-line arguments this script can receive
    parser.add_argument('sequence',
        metavar='SEQUENCE',
        type=str,
        help= ('The sequence to search for an open-reading frame. '
               'If the path flag (\'-p\'/\'--path\') is specified, '
               'then this should be a path to a file containing the '
               'sequence to be searched.'))

    parser.add_argument('-p', '--path',
        action='store_true',
        help= ('The sequence argument should be treated as a  path to a'
                'containing the sequence to be searched.'))


    parser.add_argument('-s', '--start-codons',
        type=str,
        nargs= '+', # one or more arguments
        default= ['AUG'],
        help=('One or more possible start codons.'))

    parser.add_argument('-x', '--stop-codons',
        type=str,
        nargs='+',
        default= ['UAA', 'UAG', 'UGA'],
        help= ('one or more possible stop codons.'))

     # Parse the command-line arguments into a 'dict'-like container
    args = parser.parse_args()

    # Check to see if the path option was set to True by the caller. If so, parse
    # the sequence from the path
    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    orf = find_orf.find_first_orf(sequence = sequence,
        start_codons = args.start_codons,
        stop_codons = args.stop_codons)
