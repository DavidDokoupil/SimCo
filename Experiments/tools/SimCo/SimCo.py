from buddy import bddtrue

from base import Base
from state import State
 
class SimCo(Base):
    def __init__(self, input, flags):
        super().__init__(input, flags)
        
    def complement(self):
        start = self.initStart()
        
        stateMap = {start : self.out.new_state()}
        
        self.out.set_init_state(stateMap[start])
        
        # Construct in 2 parts to achieve better readability (Easier to go over the output automaton)
                
        # Upper states
        upQueue= [start]
        
        while upQueue:
            current = upQueue.pop(0)
            
            for minterm in self.minterms(bddtrue):
                successor, _ = self.successors(current, minterm, True)

                if successor not in stateMap:
                    stateMap[successor] = self.out.new_state()
                    upQueue.append(successor)
                    
                self.out.new_edge(stateMap[current], stateMap[successor], minterm, [])
        
        # Lower states 
        downQueue = [macroState for macroState in stateMap.keys()]

        while downQueue:
            current = downQueue.pop(0)
            
            for minterm in self.minterms(bddtrue):
                successor, marks = self.successors(current, minterm, False)
                marks = [mark for pMarks in marks for mark in pMarks]
                marks = [mIdx for mIdx in range(len(marks)) if marks[mIdx]] if not current.up else []
                
                if self.flags['skip'] and self.stateAlwaysTerminal(successor):
                    continue
                
                if successor not in stateMap:
                    stateMap[successor] = self.out.new_state()
                    downQueue.append(successor)

                self.out.new_edge(stateMap[current], stateMap[successor], minterm, marks)
        
        self.out.set_state_names(self.stateNames(stateMap))
        
        if self.postprocessor:
            self.postprocess()
        
        return self.out 
    
    
    def initStart(self): # Initializes the initial state 
        start = State(True, self.tCount, [len(self.pSCC2SCC[tIdx]) - 1 for tIdx in range(self.tCount)], self.isGenBuchi, self.flags['scc-based']) # TODO: Change after flags implemented
        for tIdx in range(self.tCount):
            if not self.flags['scc-based']:
                start.trackers[tIdx][0].append((set([self.start]), 0))
                continue
            pIdx = self.state2aSCC(tIdx, self.start)
            
            if pIdx == -1:
                start.trackers[tIdx][pIdx].add(self.start)
                continue
            
            start.trackers[tIdx][pIdx].append((set([self.start]), 0))
        return start
    
    
    def successors(self, state, minterm, up):
        successor = State(up, self.tCount, [len(self.pSCC2SCC[tIdx]) - 1 for tIdx in range(self.tCount)], self.isGenBuchi, True) # TODO: Change after flags implemented
        marks = self.pMarksTemplate()
        for tIdx in range(state.tCount):
            allExits = self.pTemplate(tIdx, True)
            contained = self.pTemplate(tIdx)

            for pIdx in self.pSCC2SCC[tIdx]:
                if self.flags['scc-based']:
                    exits = self.exitingSuccessors(state, minterm, tIdx, pIdx)
                    for exitTarget, exitPartition in zip(self.pSCC2SCC[tIdx].keys(), exits):
                        allExits[exitTarget] |= exitPartition
                     
                if pIdx == -1 or pIdx == len(state.trackers[tIdx]):
                    continue

                contained[pIdx], mark = self.containedSuccessors(state, minterm, tIdx, pIdx, up)
                marks[tIdx][pIdx] = mark
                

            successor.trackers[tIdx] = self.combinePartitions(allExits, contained, tIdx)
        
        if self.flags['merge']:
            self.mergeTainted(successor)
        
        return successor, marks
    
    def exitingSuccessors(self, state, minterm, tIdx, pIdx):
        partition = state.trackers[tIdx][pIdx]
        exited = self.pTemplate(tIdx, True)
        targetSCCs = [k for k in self.pSCC2SCC[tIdx].keys()]
        
        states = partition if pIdx == -1 else [state for component in partition for state in component[0]]
        for s in states:
            for e in self.input.out(s):
                if minterm not in self.minterms(e.cond):
                    continue
                
                for targetSCC in targetSCCs:
                    if e.dst in self.pSCC2States[tIdx][targetSCC]:
                        exited[targetSCC] |= {e.dst}
        
        return exited
    
    
    def containedSuccessors(self, state, minterm, tIdx, pIdx, up):
        uSuccessorPartition = []
        dSuccessorPartition = []
        
        partition = state.trackers[tIdx][pIdx]
        
        currentLevels = state.pLevels(tIdx, pIdx)
        terminalCount = state.pTerminalCount(tIdx, pIdx)
        
        g = set()
        
        nextNonterminalIdx = 0
        
        for idx in range(len(partition) - 1, -1 , -1):
            l, t, s = self.componentSuccessors(partition[idx], minterm, tIdx, pIdx)
            
            L = partition[idx][1]
            
            self.makeGreedy(L, self.maxLevels[tIdx], g, l, t, s)
            self.lClean(L, l)
            
            terminal = (L == -2 or -2 not in currentLevels or terminalCount == 0)
             
            if L == -2 and len(t) == 0: # Terminal died
                terminalCount -= 1
                if self.flags['round-robin']:
                    dSuccessorPartition.append([set(), - 2])
            
            if L < 0 and len(t) > 0:
                newL = -2 if terminal else - 3
                dSuccessorPartition.append([t, newL])

            
            for component in l[::-1]:
                uSuccessorPartition.append(component)

            dSuccessorPartition += self.lTransform(L, self.maxLevels[tIdx], l, terminal)                     

            
            uSuccessorPartition.append([s, L])
            dSuccessorPartition.append([s, L])
            
                
        self.resolveLevels(dSuccessorPartition, terminalCount)
        uSuccessorPartition = self.pClean(uSuccessorPartition)
        dSuccessorPartition = self.pClean(dSuccessorPartition)
        return (uSuccessorPartition[::-1], False) if up else (dSuccessorPartition[::-1], terminalCount == 0)

                        
    def componentSuccessors(self, component, minterm, tIdx, pIdx):
        l = [[set(), 0] for _ in range(self.maxLevels[tIdx])] # Leveled-up
        t = set() # Tainted
        s = set() # Stagnant
        
        S, L = component # States, Level
        isTainted = L < 0
        
        
        for state in S:
            for e in self.input.out(state):                
                if (minterm not in self.minterms(e.cond) or 
                    (self.flags['scc-based'] and e.dst not in self.pSCC2States[tIdx][pIdx])):
                    continue
                    
                if isTainted:
                    t |= {e.dst}
                
                levelMappedMarks = [self.levelMaps[tIdx][m] for m in e.acc.sets() if m in self.levelMaps[tIdx].keys()]
                
                if L in levelMappedMarks:
                    successorL = (L + 1) % self.maxLevels[tIdx]
                    if self.flags['jumps']: 
                        while successorL != L and successorL in levelMappedMarks:
                            successorL += 1
                            successorL %= self.maxLevels[tIdx]
                    l[successorL][0] |= {e.dst}
                    l[successorL][1] = successorL
                        
                else:
                    s |= {e.dst}
        return l,t,s
    
    
    def combinePartitions(self, exited, contained, tIdx):
        if not self.flags['scc-based']:
            return contained
        partitionStates = self.pTemplate(tIdx, True)
        tracker = self.pTemplate(tIdx, False)
        for partition, pIdx in zip(contained, self.pSCC2SCC[tIdx].keys()):
            if pIdx == -1:
                continue
            for S, _ in partition:
                partitionStates[pIdx] |= set(S)
                
        for pIdx in range(len(exited) - 1):
            exited[pIdx] -= partitionStates[pIdx]

        for pIdx, (incoming, partition) in enumerate(zip(exited, contained)):
            if len(incoming):
                tracker[pIdx].append([incoming, 0])
            if len(partition):
                tracker[pIdx] += list(partition)
        
        tracker[-1] = set(exited[-1]) if self.flags['scc-based'] else set()

        return tracker
    
    