class State():
    def __init__(self, up, tCount = 1, pCounts = [1], genBuchi = False, sccBased = False):
        self.up = up
        self.fTela = genBuchi
        self.sccBased = sccBased
        
        # Structure:
        #   Trackers (for conjuncts)
        #       Segments (for SCCs)
        #           Components (with levels)
        
        self.trackers = [[] for _ in range(tCount)]
        
        for tIdx, tracker in enumerate(self.trackers):
            tracker[:] = [[] for _ in range(pCounts[tIdx])]
            if sccBased:
                tracker.append(set())
                continue
            tracker.append([])

        self.tCount = len(self.trackers)

        return
    
    
    def __hash__(self):
        return hash(str(self))
    
    
    def __eq__(self, value):
        return hash(self) == hash(value)
    
    
    def __str__(self):
        flag = "▴" if self.up else "▾"
        
        result = "("
        for tracker in self.trackers:
            result += "⟨"
            
            for partition in tracker:
                if type(partition) == type(set()):
                    result += str(partition) if len(partition) else "{}"
                    continue
                
                result += "["
                for component in partition:
                    result += ("{" + ",".join([str(S) for S in component[0]]) + "}"
                               if len(component[0]) > 1
                               else str(component[0])) 
                    result += ":" + str(component[1])
                    result += ","
                    
                result = result[:-1] if result[-1] == "," else result
                result += "]|"
                
            result += "⟩"
            
        result += flag + ")"
        
        return result
    
    def pLevels(self, tIdx, pIdx):
        return set([level for _, level in self.trackers[tIdx][pIdx]])
    
    def pTerminalCount(self, tIdx, pIdx):
        return len([_ for _, level in self.trackers[tIdx][pIdx] if level == -2])
        