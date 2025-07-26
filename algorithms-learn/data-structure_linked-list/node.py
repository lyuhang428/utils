# -*- encoding: utf-8 -*-
# None type in Python is similar to NULL ptr in C



class EmptyLinkedListError(Exception):
    '''used in LinkedList class, check if insert or delete in an empty linked list'''
    def __init__(self, msg="Empty linked list cannot be inserted or deleted"):
        super().__init__()
        self.msg = msg

    def __str__(self):
        return self.msg


class Node:
    def __init__(self, data:int|float):
        self.data = data
        self.next:Node|None = None
        self.prev:Node|None = None

    def __repr__(self) -> str:
        return str(self.__dict__)

    def __str__(self) -> str:
        return str(self.__dict__)


class LinkedList:
    '''doubly linked list'''

    def __init__(self, node:Node|None=None):
        self.head: Node | None = node
        self.tail: Node | None = node
        self.__idx_counter = 0

    def __len__(self) -> int:
        return self.__idx_counter
    
    def __iter__(self):
        '''generator, so ll can be converted to list(ll), and return iterator 
            for ii in ll:
                ...
            An iterable must have __iter__ but not necessary __next__
            An iterator must have both __iter__ and __next__
            A generator must have yield, and __iter__ and __next__ are generated automatically
            An iterable is an object that can return an iterator
            An iterator that is its own iterator and knows how to get the next value
            A generator is a function that automatically creates an iterator
            All iterators are iterables, and generators are iterators that don't need to manully build __iter__ and __next__
        '''
        cur_node = self.head
        while cur_node is not None:
            yield cur_node
            cur_node = cur_node.next

    def __getitem__(self, idx:int) -> Node|None:
        '''index-based, allow to use ll[2]. Alone also allows  for ii in ll:...'''
        if self.head is None:
            return None
        
        if (idx < 0) or idx >= len(self):
            raise IndexError("index out of bound")

        if idx <= (len(self) // 2):
            cur_node = self.head
            count = 0
            while count < idx:
                cur_node = cur_node.next
                count += 1
            return cur_node
        else:
            cur_node = self.tail
            count = len(self)
            while count > idx+1:
                cur_node = cur_node.prev
                count -= 1
            return cur_node

    def is_empty(self) -> bool:
        return self.head is None
    
    def append(self, node:Node):
        self.__idx_counter += 1
        if self.head is None:
            self.head = node
            self.tail = node
            return
        
        cur_node = self.tail
        cur_node.next = node
        node.prev = cur_node
        self.tail = node
    
    def prepend(self, node:Node):
        self.__idx_counter += 1
        if self.head is None:
            self.head = node
            self.tail = node
            return 
        
        cur_node = self.head
        cur_node.prev = node
        node.next = cur_node
        self.head = node

    def insert_at_position(self, idx:int, node:Node):
        if self.head is None:
            print("empty linked list, new node will be head/tail")
            self.append(node)
            return
        
        if idx <= 0:
            self.prepend(node)
            return 
        elif idx >= len(self):
            self.append(node)
            return 
        
        if idx <= (len(self) // 2):
            cur_node = self.head
            count = 0
            while count < idx:
                cur_node = cur_node.next
                count += 1

            nxt_node = cur_node.next
            node.prev = cur_node
            node.next = nxt_node
            cur_node.next = node
            nxt_node.prev = node
            self.__idx_counter += 1
        else:
            cur_node = self.tail
            count = len(self)
            while count > idx:
                cur_node = cur_node.prev
                count -= 1

            prv_node = cur_node.prev
            node.prev = prv_node
            node.next = cur_node
            cur_node.prev = node
            prv_node.next = node
            self.__idx_counter += 1

    def delete_at_position(self, idx:int):
        if self.head is None:
            raise EmptyLinkedListError()
        
        if len(self) == 1:
            self.head = None
            self.tail = None
            self.__idx_counter = 0
            return 

        if idx <= 0:
            self.head = self.head.next
            self.head.prev = None
            self.__idx_counter -= 1
            return 
        elif idx >= len(self):
            self.tail = self.tail.prev
            self.tail.next = None
            self.__idx_counter -= 1
            return 
        
        if idx <= (len(self) // 2):
            cur_node = self.head
            count = 0
            while count < idx:
                cur_node = cur_node.next
                count += 1
            
            nxt_node = cur_node.next
            prv_node = cur_node.prev
            nxt_node.prev = prv_node
            prv_node.next = nxt_node
            self.__idx_counter -= 1
        else:
            cur_node = self.tail
            count = len(self)
            while count > idx+1:
                cur_node = cur_node.prev
                count -= 1

            prv_node = cur_node.prev
            nxt_node = cur_node.next
            prv_node.next = nxt_node
            nxt_node.prev = prv_node
            self.__idx_counter -= 1

    def delete_by_value(self, val:int|float):
        if self.head is None:
            raise EmptyLinkedListError()

        if self.head.data == val:
            nxt_node = self.head.next
            nxt_node.prev = None
            self.head = nxt_node
            self.__idx_counter -= 1
            return

        cur_node = self.head
        while cur_node.next is not None:
            if cur_node.data == val:
                prv_node = cur_node.prev
                nxt_node = cur_node.next
                prv_node.next = nxt_node
                nxt_node.prev = prv_node
                self.__idx_counter -= 1
                return 
            cur_node = cur_node.next

        if self.tail.data == val:
            prv_node = self.tail.prev
            prv_node.next = None
            self.tail = prv_node
            self.__idx_counter -= 1
            return 
        
        print("node with target value not found, linked list returned as it is")
        return 

    def search_by_value(self, val:int|float) -> Node|None:
        if self.head is None:
            print("linked list is empty, return None")
            return None
        
        if self.head.data == val:
            return self.head
        
        cur_node = self.head
        while cur_node.next is not None:
            if cur_node.data == val:
                return cur_node
            cur_node = cur_node.next

        if cur_node.data == val:
            return cur_node
        else:
            print("not found")
            return None

    def display(self) -> None:
        if self.is_empty():
            print("empty list")
            return 
        
        print("None <- ", end='')
        cur_node = self.head
        while cur_node is not None:
            if cur_node.next is None:
                print(f"{cur_node.data} -> None")
            else:
                print(f"{cur_node.data} <-> ", end='')
            cur_node = cur_node.next


if __name__ == "__main__":
        pass






