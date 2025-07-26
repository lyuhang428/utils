# -*- encoding: utf-8 -*-
from node import Node, LinkedList, EmptyLinkedListError

def main():
    ll = LinkedList()
    for i in range(10):
        ll.append(Node(i))

    # print(f"length of link {len(ll)}")
    # ll.display()

    # node = ll.search_by_value(5)
    # print(node.prev.data, node.data, node.next.data)

    # ll.insert_at_position(5, Node(42))
    # print(f"length of inserted link {len(ll)}")
    # ll.display()

    for ii in ll:
        print(ii.data)

    print(list(ll))
    print(ll[2])

if __name__ == "__main__":
    main()