from abc import ABC, abstractmethod

from utils import unflatten, getStateSimulators

import spot
from buddy import bddtrue, bddfalse, bdd_support, bdd_satoneset


class Base(ABC):
    @abstractmethod
    def __init__(self, input, flags):
        self.input = spot.scc_filter(input)

        self.flags = flags
        
        if self.flags['skip']:
            self.input = spot.complete(self.input)
        
        if self.flags['sim']: # TODO: After it works, make compatible with other optimizations
            # self.input = spot.scc_filter(input)
            n = self.input.num_states()
            # sim_rel[i][j] => state "i" simulates "j" i>j
            self.sim_rel = unflatten(spot.get_sim_matrix(input)[1], n)
            # print(self.sim_rel)
            # print(self.input.to_str())
            # print(getStateSimulators(2,self.sim_rel))
        
        #TODO: Explore permutating the acc marks
        
        # Input Acceptance Processing
        self.maxLevels = []
        self.levelMaps = []
        self.parseAcceptance()
        
        self.start = self.input.get_init_state_number()
        
        self.aSCC = [] # Accepting SCCs (different for each tracker) ((Does not need to be saved))
        self.baSCC = [] # Bottom Accepting SCCs (different for each tracker)
        self.pSCC2SCC = [] # partition indexing -> Spot indexing (different for each tracker)
        self.pSCC2States = [] # aSCC index -> States of SCC (different for each tracker)
                
        self.sccInit() 
        
        self.APInit()
        
        self.complementInit()
        
        if not self.flags['scc-based']:
            self.pSCC2SCC = [{0:0} for _ in range(self.tCount)]
            self.pCount = 1
        
        self.cache = {}
        
        self.postprocessor = None if not self.flags['postprocess'] else spot.postprocessor()
        
        
        
    def minterms(self, label):
        cached = self.cache.setdefault('minterms', dict())
        
        if label in cached:
            return cached[label]
        
        minterms = []

        all_ = label

        while all_ != bddfalse:
            one = bdd_satoneset(all_, self.AP, bddfalse)
            all_ -= one
            minterms.append(one)
            
        cached[label] = minterms
        return minterms
        
    
    def parseAcceptance(self):
        self.acc = self.input.acc()

        if self.acc.uses_fin_acceptance():
            raise ValueError("Acceptance with \'Fin\' not allowed.")
        
        self.isGenBuchi = self.acc.is_generalized_buchi()
                                                           
        dnf_acc = self.acc.get_acceptance().to_dnf()
        top_disjuncts = dnf_acc.top_disjuncts()
        for conjunct in top_disjuncts:
            marks = set(conjunct.used_sets().sets())
            self.levelMaps.append({v:k for k,v in dict(enumerate(marks)).items()})
            self.maxLevels.append(len(self.levelMaps[-1].keys()))
        
        self.tCount = len(self.levelMaps)
        
        return
    

    def sccInit(self):
        subs = [] # Sub-Acceptance-Formulas
        for mapping in self.levelMaps:
            marks = list(mapping.keys())
            sub = ""
            for mark in marks:
                sub += f"Inf({mark})&"
            sub= sub[:-1]
            subs.append((len(marks), sub))

        for sub in subs:
            self.subAcceptanceSccInit(sub[0], sub[1]) # Init aSCC, bSCC, dicts for subformulas
        
        self.pCount = max([len(partitions) for partitions in self.pSCC2SCC])
                 
        return
    
    
    def subAcceptanceSccInit(self, count, formula):
        if not formula:
            formula = 'all'
        self.input.set_acceptance(count, formula)
        
        sccInfo = spot.scc_info(self.input)
        
        self.aSCC.append([i for i in range(sccInfo.scc_count()) if sccInfo.is_maximally_accepting_scc(i)])
        self.baSCC.append([])
        
        self.pSCC2SCC.append(dict(enumerate(self.aSCC[-1])))
        self.pSCC2SCC[-1][-1] = -1 # Denotes partition of states in non-accepting SCCs
        
        for pIdx, sIdx in self.pSCC2SCC[-1].items(): # aIdx -> aSCC index; sIdx -> spot index
            if pIdx == -1:
                continue

            inEdges = [e for e in sccInfo.inner_edges_of(sIdx)]
            allEdges = [e for e in sccInfo.edges_of(sIdx)]
            
            if len(inEdges) == len(allEdges):
                self.baSCC[-1].append(pIdx)
            
        SCC2States = dict((aIdx, sccInfo.states_of(aIdx)) for aIdx in self.aSCC[-1])
            
        allStates = set([e.src for e in self.input.edges()])
        assignedStates = set([state for SCC in SCC2States.values() for state in SCC])
            
        SCC2States[-1] = tuple(allStates - assignedStates)
        pSCC2States = {pIdx: SCC2States[sIdx] for pIdx, sIdx in self.pSCC2SCC[-1].items()}
        

        self.pSCC2States.append(pSCC2States)
        
        return
                
            
    def APInit(self):
        self.AP = bddtrue
        
        conds = {edge.cond for edge in self.input.edges()}
        
        for cond in conds:
            self.AP &= bdd_support(cond)
    
    
    def state2aSCC(self, tIdx, state):
        if not self.flags['scc-based']:
            return 0
        pSCC2States = self.pSCC2States[tIdx]
        
        for pIdx, states in pSCC2States.items():
            if state in states:
                return pIdx

        raise Exception("State not part of the input automaton")


    def complementInit(self): # Initialize variables tied to the output complement
        bdd_dict = spot.make_bdd_dict()
        self.out = spot.make_twa_graph(bdd_dict)
        self.out.copy_ap_of(self.input)

        if self.flags['scc-based']:
            self.out.set_generalized_buchi(sum([len(x) for x in self.aSCC]))
            return
        self.out.set_generalized_buchi((self.tCount))
        
        
    def pTemplate(self, tIdx, useSets=False):
        if useSets:
            return [set() for _ in range(len(self.pSCC2SCC[tIdx]))]
        return [[] for _ in range(len(self.pSCC2SCC[tIdx]))]
    
    
    def pMarksTemplate(self): # Template for marks of a single partition
        if self.flags['scc-based']:
            return [[False for _ in range(len(self.pSCC2SCC[tIdx]) - 1)] for tIdx in range(self.tCount)]
        return [[False for _ in range(len(self.pSCC2SCC[tIdx]))] for tIdx in range(self.tCount)]
  
    
    def makeGreedy(self, L, maxL, g, l, t, s):
        
        if self.flags['couple-jumps']:
            most:set = l[(L+1) % maxL][0]
            for i in range(len(l)):
                if not most.issubset(l[(L+ 1 + i) % maxL][0]):
                    l[(L+ 1 + i) % maxL][0] = set() # TODO: Verify correctness [Seems like correct]
         
        t -= g
        
        if self.flags['sim']:
            t = self.simClean(t, t)
            t = self.simClean(t, g)
        
        g |= t
        
            
        for i in range(maxL):
            # print("I:", i) # TODO: Check functionality
            l[(L - i) % maxL][0] -= g

            if self.flags['sim']:
                l[(L - i) % maxL][0] = self.simClean(l[(L - i) % maxL][0], l[(L - i) % maxL][0])
                l[(L - i) % maxL][0] = self.simClean(l[(L - i) % maxL][0], g)
            
            g |= l[(L - i) % maxL][0]
            
        s -= g
        
        if self.flags['sim']:
            # s = self.simClean(s, s)
            s = self.simClean(s, g)
        
        g |= s
        
        return
    
    
    def lClean(self, L, l):
        lShift = l[L + 1:] + l[:L + 1]
        lClean = []
        
        for i in range(len(lShift)):
            if lShift[i][0]:
                lClean.append(lShift[i])
        return lClean
    
    
    def pClean(self, pSuccessor):
        pClean = []
        for i in range(len(pSuccessor)):
            if pSuccessor[i][0]:
                pClean.append(pSuccessor[i])
        return pClean
    
    
    def lTransform(self, L, maxL, l, terminal): # Transforms leveled up components based on optimization
        if self.flags['early-taint']:
            al = set()
            for c in l:
                al |= set(c[0])
            newL = - 2 if terminal else - 3
            if self.flags['round-robin']:
                newL = -1
            return [[al, newL]] if len(al) else []
        
        lTransformed = [] 
        
        for c in l[::-1]:
            reset = not (maxL - 1 >= c[1] > L)
            newL = c[1] if not reset else - 2 if terminal else - 3
            if self.flags['round-robin']:
                newL = -1
            lTransformed.append([c[0], newL])

        return lTransformed
    
    def simClean(self, component, g):
        toRemove = set()
        for state in list(component):
            if not getStateSimulators(state, self.sim_rel).isdisjoint(g):
                toRemove.add(state)
            component -= toRemove
        return component
    
    
    def resolveLevels(self, dSuccessor, terminalCount):
        if self.flags['round-robin']:
            terminalIdx = self.terminalIdx(dSuccessor)
            if terminalIdx >= 0:
                if len(dSuccessor[terminalIdx][0]):
                    return
                nextNonterminalIdx = terminalIdx + 1

                for offset in range(len(dSuccessor)):
                    if dSuccessor[(nextNonterminalIdx + offset) % len(dSuccessor)][1] < 0:
                        dSuccessor[(nextNonterminalIdx + offset) % len(dSuccessor)][1] = - 2
                        return
                    
            for c in dSuccessor[::-1]:
                if c[1] == - 1:
                    c[1] = - 2
                    return

        for c in dSuccessor:
            if c[1] == - 3:
                c[1] = - 2 if terminalCount == 0 else - 1
                
                
    def stateAlwaysTerminal(self, state): # Rightmost component of any tracker and any partition is terminal
        if not self.flags['scc-based']:
            if not state.trackers[0][0]:
                return False
            return state.trackers[0][0][-1][1] == -2

        for tIdx in range(self.tCount):
            for pIdx in range(len(state.trackers[tIdx]) - 1):
                if (not state.trackers[tIdx][pIdx] or
                    pIdx not in self.baSCC[tIdx]):
                    continue

                if state.trackers[tIdx][pIdx][-1][1] == - 2:
                    return True
        return False
    
    
    def mergeTainted(self, state):
        mergedTrackers = []
        
        pEndOffset = 1 if self.flags['scc-based'] else 0 
        
        for tIdx in range(self.tCount):
            mergedTrackers.append([])
            for pIdx in range(len(state.trackers[tIdx]) - pEndOffset):
                mergedTrackers[-1].append([])
                if len(state.trackers[tIdx][pIdx]):
                    mergedTrackers[-1][-1].append(state.trackers[tIdx][pIdx][0])
                for cIdx in range(1, len(state.trackers[tIdx][pIdx])):
                    if (mergedTrackers[-1][-1][-1][1] in [-1,-2] and 
                        mergedTrackers[-1][-1][-1][1] == state.trackers[tIdx][pIdx][cIdx][1]): # Component has the same level  
                        mergedTrackers[-1][-1][-1][0] |= state.trackers[tIdx][pIdx][cIdx][0]
                        continue
                    mergedTrackers[-1][-1].append(state.trackers[tIdx][pIdx][cIdx])
            if self.flags['scc-based']:
                mergedTrackers[-1].append(state.trackers[tIdx][-1])
        
        state.trackers = mergedTrackers
    
    
    def propagateTerminal(self, state):
        propagatedTrackers = []
        
        pEndOffset = 1 if self.flags['scc-based'] else 0 
        
        for tIdx in range(self.tCount):
            propagatedTrackers.append([])
            for pIdx in range(len(state.trackers[tIdx]) - pEndOffset):
                propagatedTrackers[-1].append([])
                if len(state.trackers[tIdx][pIdx]):
                    propagatedTrackers[-1][-1].append(state.trackers[tIdx][pIdx][0])
                for cIdx in range(1, len(state.trackers[tIdx][pIdx])):
                    if (propagatedTrackers[-1][-1][-1][1] == -2 and 
                        state.trackers[tIdx][pIdx][cIdx][1] == -1): 
                        propagatedTrackers[-1][-1][-1][0] |= state.trackers[tIdx][pIdx][cIdx][0]
                        continue
                    propagatedTrackers[-1][-1].append(state.trackers[tIdx][pIdx][cIdx])
            if self.flags['scc-based']:
                propagatedTrackers[-1][-1].append(state.trackers[tIdx][-1])
        
        state.trackers = propagatedTrackers
        
    
    def terminalIdx(self, partition):
        for cIdx in range(len(partition)):
            if partition[cIdx][1] == -2:
                return cIdx
            
        return -1
    
    def simulatorContained(self, state, partition):
        simulators: set = getStateSimulators(state, self.sim_rel)
        for cl in partition:
            if not simulators.isdisjoint(cl[0]):
                return True
        return False
    
    def stateNames(self, stateMap):
        return [str(state) for state in stateMap]

    def postprocess(self):
        self.out = spot.scc_filter(self.out)
        
        self.postprocessor.set_pref(self.postprocessor.Small)

        # if self.out.num_states() < 128:
        #     self.postprocessor.set_level(self.postprocessor.High)
        if self.out.num_states() < 48:
            self.postprocessor.set_level(self.postprocessor.Medium)
        else:
            self.postprocessor.set_level(self.postprocessor.Low)
        
        self.out = self.postprocessor.run(self.out)

