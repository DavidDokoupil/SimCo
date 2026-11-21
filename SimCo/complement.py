from argparse import ArgumentParser

import spot

spot.setup()

from SimCo import SimCo

if __name__ == "__main__":
    arg_parser = ArgumentParser(description="Complements transition-based generalized Büchi automata (HOA format) "
                                            "with the construction of the LPAR 2026 paper; the complements are "
                                            "printed in HOA format.")
    arg_parser.add_argument('file',
                            type=str,
                            nargs='*',
                            help='files with automata to complement (default: standard input)',
                            default=['-'])
    arg_parser.add_argument('--jump',
                            action='store_true',
                            help='level jumping: jump to the highest possible level at once')
    arg_parser.add_argument('--merge',
                            action='store_true',
                            help='component merging: merge tainted components of the same level')
    arg_parser.add_argument('--postprocess', '--post-process', action='store_true',
                            help="simplify the complement with Spot's postprocessor (preference Small; level Medium "
                                 "for complements with fewer than 48 states, Low otherwise)")

    arguments = arg_parser.parse_args()

    flags = {
        'jumps': arguments.jump,
        'merge': arguments.merge,
        'postprocess': arguments.postprocess,
    }

    for f in arguments.file:
        for a in spot.automata(f):
            print(SimCo(a, flags).complement().to_str())
