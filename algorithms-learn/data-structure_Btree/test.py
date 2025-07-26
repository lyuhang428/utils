from btree import TreeNode, BTree

def main():
    btree = BTree()
    btree.add(TreeNode(42))

    for i in range(37,47):
        btree._add(TreeNode(i))

   

    btree.DFS_preorder(btree.root); print()
    btree.DFS_inorder(btree.root); print()
    btree.DFS_postorder(btree.root); print()

    print(len(btree))
    btree.BFS(btree.root)


if __name__ == '__main__':
    main()




