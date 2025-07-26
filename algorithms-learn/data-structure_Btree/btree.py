# -*- encoding: utf-8 -*-
# ref: https://menghaoguo.com/Data_Structure_with_Python_book/chapter7/section2.html
# https://menghaoguo.com/Data_Structure_with_Python_book/chapter7/section1.html
# https://python-data-structures-and-algorithms.readthedocs.io/zh/latest/14_%E6%A0%91%E4%B8%8E%E4%BA%8C%E5%8F%89%E6%A0%91/tree/
# https://cloud.tencent.com/developer/article/1678925


class EmptyTreeError(Exception):
    def __init__(self, msg="this is empty tree", *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.msg = msg

    def __str__(self) -> str:
        return self.msg


class TreeNode:
    def __init__(self, val:int|float):
        self.val = val
        self.left:TreeNode|None = None
        self.right:TreeNode|None = None

    def __repr__(self) -> str:
        return str(self.__dict__)




class BTree:
    def __init__(self, root:TreeNode|None=None):
        self.root = root
        self.__nodes_counter:int = 0

    def __len__(self) -> int:
        return self.__nodes_counter

    def add(self, tnode:TreeNode):
        '''bst style, left node is always smaller than right node'''
        if self.root is None:
            self.root = tnode
            self.__nodes_counter += 1
            return
        
        cur_tnode = self.root
        while True: #cur_tnode is not None:
            if tnode.val <= cur_tnode.val:
                if cur_tnode.left is None:
                    cur_tnode.left = tnode
                    self.__nodes_counter += 1
                    return 
                cur_tnode = cur_tnode.left
            else:
                if cur_tnode.right is None:
                    cur_tnode.right = tnode
                    self.__nodes_counter += 1
                    return 
                cur_tnode = cur_tnode.right

    def _add(self, tnode:TreeNode):
        '''full tree style, left to right, layer by layer'''
        if self.root is None:
            self.root = tnode
            self.__nodes_counter += 1
            return
        
        queue = [self.root]
        while queue != []:
            cur_tnode = queue.pop(0)
            if cur_tnode.left is None:
                cur_tnode.left = tnode
                self.__nodes_counter += 1
                return
            elif cur_tnode.right is None:
                cur_tnode.right = tnode
                self.__nodes_counter += 1
                return
            else:
                queue.append(cur_tnode.left)
                queue.append(cur_tnode.right)
 
    def DFS_preorder(self, root:TreeNode|None):
        '''depth first traverse, preorder'''
        if self.root is None:
            print("this is empty tree")
            return
        if root is None:
            return
        print(root.val, end=' ')
        self.DFS_preorder(root.left)
        self.DFS_preorder(root.right)

    def DFS_inorder(self, root:TreeNode|None):
        '''depth first traverse, inorder'''
        if self.root is None:
            print("this is empty tree")
            return
        if root is None:
            return 
        self.DFS_inorder(root.left)
        print(root.val, end=' ')
        self.DFS_inorder(root.right)

    def DFS_postorder(self, root:TreeNode|None):
        '''depth first traverse, postorder'''
        if self.root is None:
            print("this is empty tree")
            return 
        if root is None:
            return
        self.DFS_postorder(root.left)
        self.DFS_postorder(root.right)
        print(root.val, end=' ')

    def BFS(self, root:TreeNode|None):
        '''bredth first traverse '''
        if self.root is None:
            print("this is empty tree")
            return
        if root is None:
            return
        
        queue = [self.root]
        while queue != []:
            cur_tnode = queue.pop(0)
            print(cur_tnode.val, end=' ')
            if cur_tnode.left is not None:
                queue.append(cur_tnode.left)
            if cur_tnode.right is not None:
                queue.append(cur_tnode.right)


if __name__ == '__main__':
    pass






