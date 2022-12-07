# python3

import sys
from IPython import embed
'''
Implement MaximalNonBranchingPaths.
     Input: The adjacency list of a graph whose nodes are integers.
     Output: The collection of all maximal nonbranching paths in this graph.
Sample Input:
1 -> 2
2 -> 3
3 -> 4,5
6 -> 7
7 -> 6
Sample Output:
1 -> 2 -> 3
3 -> 4
3 -> 5
7 -> 6 -> 7
A node v in a directed graph Graph is called a 1-in-1-out node if its indegree and outdegree are both equal to 1, 
i.e., in(v) = out(v) = 1.  We can rephrase the definition of a "maximal non-branching path" from the main text as
a path whose internal nodes are 1-in-1-out nodes and whose initial and final nodes are not 1-in-1-out nodes. Also, 
note that the definition from the main text does not handle the special case when Graph has a connected component 
that is an isolated cycle, in which all nodes are 1-in-1-out nodes.
The MaximalNonBranchingPaths pseudocode below generates all non-branching paths in a graph. It iterates through all
nodes of the graph that are not 1-in-1-out nodes and generates all non-branching paths starting at each such node. 
In a final step, MaximalNonBranchingPaths finds all isolated cycles in the graph.
    MaximalNonBranchingPaths(Graph)
        Paths ← empty list
        for each node v in Graph
            if v is not a 1-in-1-out node
                if out(v) > 0
                    for each outgoing edge (v, w) from v
                        NonBranchingPath ← the path consisting of the single edge (v, w)
                        while w is a 1-in-1-out node
                            extend NonBranchingPath by the outgoing edge (w, u) from w 
                            w ← u
                        add NonBranchingPath to the set Paths
        for each isolated cycle Cycle in Graph
            add Cycle to Paths
        return Paths
'''


class MaximalNonBrachingPath:
    def __init__(self):
        self.adj = self._input()
        self.updateAdj(self.adj)
        self.paths = self.findMaximalNonBranchingPaths()
        self.printPaths(self.paths)

    def _input(self):
        data = list(sys.stdin.read().strip().split())
        adj = dict()
        for i in range(len(data) // 3):
            nl = data[i * 3]
            nrList = data[i * 3 + 2].split(',')
            adj[nl] = adj.get(nl, []) + nrList
            for nr in nrList:
                if nr not in adj:
                    adj[nr] = []
        return adj

    def updateAdj(self, adj):
        self.n = len(adj)
        self.inDeg = dict()
        self.outDeg = dict()
        for w, vList in adj.items():
            self.inDeg[w] = self.inDeg.get(w, 0)
            for v in vList:
                self.inDeg[v] = self.inDeg.get(v, 0) + 1
            self.outDeg[w] = len(vList)
        return

    def findMaximalNonBranchingPaths(self):
        paths = []
        nodes1in1out = set()  # 1-in-1-out nodes
        nExplored = set()  # 1-in-1-out nodes which were explored
        for v in self.adj.keys():
            if not (1 == self.in_degree[v] and 1 == self.out_degree[v]):
                if self.out_degree[v] > 0:
                    for w in self.adj[v]:
                        nbPath = [v, w]  # NonBrachingPath
                        while 1 == self.in_degree[w] and 1 == self.out_degree[w]:
                            nExplored.add(w)
                            u = self.adj[w][0]
                            nbPath.append(u)
                            w = u
                        paths.append(nbPath)
            else:
                nodes1in1out.add(v)
        for v in nodes1in1out:
            if v not in nExplored:
                w = v
                nbPath = []
                while w in nodes1in1out:
                    nbPath.append(w)
                    if w == v and len(nbPath) > 1:
                        paths.append(nbPath)
                        for node in nbPath:
                            nExplored.add(node)
                        break
                    w = self.adj[w][0]
        return paths

    def printPaths(self, paths):
        for p in paths:
            print(' -> '.join(p))


if __name__ == "__main__":
    MaximalNonBrachingPath()