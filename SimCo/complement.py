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
                            default=['-'])
    arg_parser.add_argument('--jump', 
                            action='store_true',
                            help='enable jumping to the highest possible level at once')
    arg_parser.add_argument('--couple-jump', 
                            action='store_true',
                            help='enable jumping to the highest possible level at once such that all successors can have at least that level')
    arg_parser.add_argument('--merge',
                            action='store_true',
                            help='use the tainted component-merging optimization')
    arg_parser.add_argument('--scc-based', 
                            action='store_true',
                            help='use the SCC decomposition optimization')
    arg_parser.add_argument("--sim", 
                            action='store_true', 
                            help='use simulation-based greediness [Experimental]',
                            )
    arg_parser.add_argument("--skip", 
                            action='store_true', 
                            help='skip calculation of successors for states in the lower automaton with tainted rightmost component')
    arg_parser.add_argument("--early-taint", 
                            action='store_true', 
                            help='forbid any level-up in the lower automaton')
    arg_parser.add_argument("--test", 
                            action='store_true', 
                            )
    arg_parser.add_argument("--round-robin",
                            action='store_true',
                            help='enable the round-robin version of the lower automata for which the optimality was proven, also forces \'early-taint\'')
    arg_parser.add_argument('--postprocess', action='store_true', default='',
                            help='use spot postprocessing with custom configuration')
    arg_parser.add_argument('--scc_filt', action='store_true', default='',
                            help='use spot postprocessing with custom configuration')
    arg_parser.add_argument('--sim_useful', action='store_true', default='',
                            help='use spot postprocessing with custom configuration')
    arg_parser.add_argument('--sim_input', action='store_true', default='',
                            help='use spot postprocessing with custom configuration')
    
    arguments = arg_parser.parse_args()

    flags = {
        'couple-jumps': arguments.couple_jump,
        'jumps': arguments.jump,
        'merge': arguments.merge,
        'scc-based': arguments.scc_based,
        'sim': arguments.sim,
        'skip': arguments.skip,
        'early-taint': arguments.early_taint or arguments.round_robin,
        'test': arguments.test,
        'round-robin': arguments.round_robin,
        'postprocess': arguments.postprocess
    }
    
    if not str(*arguments.file):
        arg_parser.print_help()
        
    if flags['couple-jumps'] and flags['jumps']:
        print("The optimizations --couple-jump and --jump are exclusive. Only one can be enabled!")
        exit()
        
    
    for a in spot.automata(*arguments.file):
        if arguments.scc_filt:
            print(spot.scc_filter(a).to_str())
            continue
        if arguments.sim_useful:
            if spot.get_sim_matrix(a)[0].num_states() != a.num_states():
                print(a.to_str())
                # print(spot.get_sim_matrix(a)[0].to_str())
            continue 
        if arguments.sim_input:
            a = spot.get_sim_matrix(a)[0]
        # a = spot.complete(a)
        algorithm = SimCo(a, flags)

        complement = algorithm.complement()
            
        print(complement.to_str())