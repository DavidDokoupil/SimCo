from argparse import ArgumentParser

import spot

spot.setup()

from SimCo import SimCo

if __name__ == "__main__":
    arg_parser = ArgumentParser(description="Complements an omega-automaton using Simple, General and Efficient construction.")
    arg_parser.add_argument('file', 
                            type=str, 
                            nargs='*', 
                            help='automata to process', 
                            default='')
    arg_parser.add_argument('--jump', 
                            action='store_true',
                            help='enable jumping to the highest possible level at once')
    arg_parser.add_argument('--merge',
                            action='store_true',
                            help='use the tainted component-merging optimization')
    arg_parser.add_argument('--scc-based', 
                            action='store_true',
                            help='use the SCC decomposition optimization')
    arg_parser.add_argument("--skip", 
                            action='store_true', 
                            help='skip calculation of successors for states in the lower automaton with tainted rightmost component')
    arg_parser.add_argument("--early-taint", 
                            action='store_true', 
                            help='forbid any level-up in the lower automaton')
    arg_parser.add_argument("--round-robin",
                            action='store_true',
                            help='enable the round-robin version of the lower automata for which the optimality was proven, also forces \'early-taint\'')

    arg_parser.add_argument('--postprocess', action='store_true', default='',
                            help='use spot postprocessing with custom configuration')
    
    arguments = arg_parser.parse_args()

    flags = {
        'jumps': arguments.jump,
        'merge': arguments.merge,
        'scc-based': arguments.scc_based,
        'skip': arguments.skip,
        'early-taint': arguments.early_taint or arguments.round_robin,
        'round-robin': arguments.round_robin,
        'postprocess': arguments.postprocess
    }
    
    if not str(*arguments.file):
        arg_parser.print_help()
    
    for a in spot.automata(*arguments.file):
        a = spot.complete(a)
        algorithm = SimCo(a, flags)

        complement = algorithm.complement()
        
        print(complement.to_str())