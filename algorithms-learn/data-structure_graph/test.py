from graph import EmptyGraphError, Vertex, Graph, Graph2

def main():
    nvertices = 5
    vertices = [Vertex(ii) for ii in range(nvertices)]
    graph = Graph(nvertices=nvertices)
    graph.add_edges(0,1)
    graph.add_edges(0,2)
    graph.add_edges(0,4)
    graph.add_edges(1,2)
    graph.add_edges(1,3)
    graph.add_edges(1,4)
    graph.add_edges(3,4)

    print(graph.adj_mat); print()
    graph.remove_vertex(2)

    print(graph.adj_mat)
    print(graph.adj_mat.shape)
    print(len(graph))

def main2():
    nvertices = 5
    graph = Graph2(nvertices)
    graph.add_edges(0,1)
    graph.add_edges(0,2)
    graph.add_edges(0,4)
    graph.add_edges(1,2)
    graph.add_edges(1,3)
    graph.add_edges(1,4)
    graph.add_edges(3,4)

    # graph.add_vertex()
    # graph.remove_vertex(1)

    graph.display()
    print(len(graph))

    # bfs_ = graph.BFS(0)
    # print(bfs_)

    dfs = graph.DFS(0)
    print(dfs)

if __name__ == '__main__':
    main2()









