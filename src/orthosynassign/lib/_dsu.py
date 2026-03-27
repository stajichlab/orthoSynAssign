import numpy as np


class DSU:
    def __init__(self, n: int):
        # Each node is its own parent initially
        self.parent = np.arange(n, dtype=np.int32)
        # Rank is used to keep the tree flat (Union by Rank)
        self.rank = np.zeros(n, dtype=np.int32)

    def find(self, i: int) -> int:
        """Finds the root of the set containing i, with path compression."""
        if self.parent[i] == i:
            return i
        # Path compression: point node directly to the root
        self.parent[i] = self.find(self.parent[i])
        return self.parent[i]

    def union(self, i: int, j: int) -> bool:
        """Merges the sets containing i and j."""
        root_i = self.find(i)
        root_j = self.find(j)

        if root_i != root_j:
            # Attach the smaller tree under the larger tree
            if self.rank[root_i] < self.rank[root_j]:
                self.parent[root_i] = root_j
            elif self.rank[root_i] > self.rank[root_j]:
                self.parent[root_j] = root_i
            else:
                self.parent[root_i] = root_j
                self.rank[root_j] += 1
            return True
        return False
