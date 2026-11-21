from abc import ABC, abstractmethod

import spot
from buddy import bddtrue, bddfalse, bdd_support, bdd_satoneset


class Base(ABC):
    @abstractmethod
    def __init__(self, input, flags):
        self.input = spot.scc_filter(input)

        self.flags = flags

        # Input Acceptance Processing
        self.maxLevels = []
        self.levelMaps = []
        self.parseAcceptance()

        self.start = self.input.get_init_state_number()

        self.pSCC2SCC = [{0: 0} for _ in range(self.tCount)] # a single partition per tracker

        self.APInit()

        self.complementInit()

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


    def APInit(self):
        self.AP = bddtrue

        conds = {edge.cond for edge in self.input.edges()}

        for cond in conds:
            self.AP &= bdd_support(cond)


    def complementInit(self): # Initialize variables tied to the output complement
        bdd_dict = spot.make_bdd_dict()
        self.out = spot.make_twa_graph(bdd_dict)
        self.out.copy_ap_of(self.input)
        self.out.set_generalized_buchi((self.tCount))


    def pTemplate(self, tIdx, useSets=False):
        if useSets:
            return [set() for _ in range(len(self.pSCC2SCC[tIdx]))]
        return [[] for _ in range(len(self.pSCC2SCC[tIdx]))]


    def pMarksTemplate(self): # Template for marks of a single partition
        return [[False for _ in range(len(self.pSCC2SCC[tIdx]))] for tIdx in range(self.tCount)]


    def makeGreedy(self, L, maxL, g, l, t, s):
        t -= g
        g |= t

        for i in range(maxL):
            l[(L - i) % maxL][0] -= g
            g |= l[(L - i) % maxL][0]

        s -= g
        g |= s

        return


    def pClean(self, pSuccessor):
        pClean = []
        for i in range(len(pSuccessor)):
            if pSuccessor[i][0]:
                pClean.append(pSuccessor[i])
        return pClean


    def lTransform(self, L, maxL, l, terminal): # Transforms leveled up components based on optimization
        lTransformed = []

        for c in l[::-1]:
            reset = not (maxL - 1 >= c[1] > L)
            newL = c[1] if not reset else - 2 if terminal else - 3
            lTransformed.append([c[0], newL])

        return lTransformed


    def resolveLevels(self, dSuccessor, terminalCount):
        for c in dSuccessor:
            if c[1] == - 3:
                c[1] = - 2 if terminalCount == 0 else - 1


    def mergeTainted(self, state):
        mergedTrackers = []

        for tIdx in range(self.tCount):
            mergedTrackers.append([])
            for pIdx in range(len(state.trackers[tIdx])):
                mergedTrackers[-1].append([])
                if len(state.trackers[tIdx][pIdx]):
                    mergedTrackers[-1][-1].append(state.trackers[tIdx][pIdx][0])
                for cIdx in range(1, len(state.trackers[tIdx][pIdx])):
                    if (mergedTrackers[-1][-1][-1][1] in [-1,-2] and
                        mergedTrackers[-1][-1][-1][1] == state.trackers[tIdx][pIdx][cIdx][1]): # Component has the same level
                        mergedTrackers[-1][-1][-1][0] |= state.trackers[tIdx][pIdx][cIdx][0]
                        continue
                    mergedTrackers[-1][-1].append(state.trackers[tIdx][pIdx][cIdx])

        state.trackers = mergedTrackers


    def stateNames(self, stateMap):
        return [str(state) for state in stateMap]


    def postprocess(self): # Small automaton; level Medium below 48 states, Low otherwise
        self.out = spot.scc_filter(self.out)

        self.postprocessor.set_pref(self.postprocessor.Small)

        if self.out.num_states() < 48:
            self.postprocessor.set_level(self.postprocessor.Medium)
        else:
            self.postprocessor.set_level(self.postprocessor.Low)

        self.out = self.postprocessor.run(self.out)
