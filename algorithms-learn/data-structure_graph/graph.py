from typing import Any, ClassVar, override
import numpy as np

class EmptyGraphError(Exception):
    def __init__(self, msg:str="Empty graph", *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.msg = msg

    def __str__(self) -> str:
        return self.msg


class Vertex:
    def __init__(self, idx:int, val:int|float|None=None):
        self.idx = idx # self.elem:str = elem
        self.val = val # self.mass:float = mass
                       # self.idx:int|str # 'C1', 'H2', 'N3' or 0,1,2... 


class Graph:
    '''基于邻接矩阵的实现'''
    def __init__(self, nvertices:int):
        self.nvertices = nvertices
        self.adj_mat = np.zeros((self.nvertices, self.nvertices), dtype=int)

    @override
    def __len__(self) -> int:
        return self.nvertices

    def add_edges(self, v_row_idx:int, v_col_idx:int, weight:int|float=1):
        self.adj_mat[v_row_idx, v_col_idx] = weight
        self.adj_mat[v_col_idx, v_row_idx] = weight

    def remove_edges(self, v_row_idx:int, v_col_idx:int):
        self.adj_mat[v_row_idx, v_col_idx] = 0
        self.adj_mat[v_col_idx, v_row_idx] = 0

    def add_vertex(self):
        self.nvertices += 1
        old = self.adj_mat.copy()
        self.adj_mat = np.zeros((self.nvertices, self.nvertices), dtype=int)
        self.adj_mat[:-1, :-1] = old

    def remove_vertex(self, v_idx):
        if v_idx == - 1:
            self.adj_mat = self.adj_mat[:-1, :-1]
            self.nvertices -= 1
            return
        if v_idx >= self.nvertices or v_idx < -1 or not isinstance(v_idx, int):
            raise IndexError("index out of bound")

        new_mat = np.zeros((self.nvertices-1, self.nvertices-1))
        new_mat[:v_idx, :v_idx] = self.adj_mat[:v_idx, :v_idx]
        new_mat[:v_idx, v_idx:] = self.adj_mat[:v_idx, v_idx+1:]
        new_mat[v_idx:, :v_idx] = self.adj_mat[v_idx+1:, :v_idx]
        new_mat[v_idx:, v_idx:] = self.adj_mat[v_idx+1:, v_idx+1:]
        self.adj_mat = new_mat
        self.nvertices -= 1
        

class Graph2:
    '''基于邻接表的实现
        邻接表似乎更普遍
        Todo: 或许用Vertex对象作键或/和值更好; 处理非连通图如:
            A --B   D--E
             \ /
              C
        以及:
            A --B   D
             \ /
              C
        vertex:{(neighbour1, weight1),(neighbour2, weight2)}
        除非需要明确矩阵操作，否则基于邻接表的图结构更普遍
    '''
    def __init__(self, nvertices:int):
        self.nvertices = nvertices
        self.adj_list = {i:set() for i in range(self.nvertices)}

    @override
    def __len__(self) -> int:
        return self.nvertices

    def add_edges(self, head:int, target:int, weight:int|float=1):
        '''adj_list[head].append(target)'''
        if (head in self.adj_list.keys()) and (target in self.adj_list.keys()):
            self.adj_list[head].add(target)
            self.adj_list[target].add(head)
        else:
            raise IndexError("index out of bound")

    def remove_edges(self, head:int, target:int):
        if (head in self.adj_list.keys()) and (target in self.adj_list.keys()):
            self.adj_list[head].remove(target)
            self.adj_list[target].remove(head)
        else:
            raise IndexError("index out of bound")

    def add_vertex(self):
        self.adj_list[self.nvertices] = set()
        self.nvertices += 1

    def remove_vertex(self, head:int):
        if head not in self.adj_list.keys():
            raise KeyError("target vertex not in adj_list")
        del self.adj_list[head]
        for key, val in self.adj_list.items():
            if head in val:
                self.adj_list[key].remove(head)
        self.nvertices -= 1

    def __BFS(self, v_start:int) -> set:
        '''set for some reason unexpectedly order the result'''
        if v_start not in self.adj_list.keys():
            raise KeyError("no such vertex")
        queue = [v_start]
        visited = set()
        visited.add(v_start)
        while queue != []:
            key = queue.pop(0)
            visited.update(self.adj_list[key])
            queue = [ii for ii in self.adj_list[key] if ii not in visited]
        return visited

    def BFS(self, v_start:int) -> list:
        if v_start not in self.adj_list.keys():
            raise KeyError("no such vertex")
        queue = [v_start]
        visited_ = []
        while queue != []:
            key = queue.pop(0)
            print(key, queue, end=' ')
            if key not in visited_:
                # remove duplicate, set() do this automatically
                visited_.append(key)
                print(visited_, end='')
                queue.extend(self.adj_list[key])
            print()
        return visited_

    def DFS(self, v_start:int, visited:list|None=None) -> list|None:
        if v_start not in self.adj_list.keys():
            raise KeyError("no such vertex")
        if visited is None:
            visited = []
        if v_start not in visited:
            visited.append(v_start)
            for ii in self.adj_list[v_start]:
                if ii not in visited:
                    self.DFS(ii, visited)
        
        return visited


    def display(self):
        for head, target in self.adj_list.items():
            print(f"Vertex{head} -> ", end='')
            for tt in target:
                print(f"{tt} -> ", end='')
            print()




if __name__ == '__main__':
    pass



