# -*- encoding: utf-8 -*-
from typing import Any
from functools import lru_cache

def timeit(n:int=1) -> callable:
    def inner(func:callable) -> callable:
        def wrapper(*args, **kwargs) -> Any:
            import time
            begin = time.time()
            for _ in range(n):
                res = func(*args, **kwargs)
            end = time.time()
            print(f"{func.__name__} repeated {n} times, took {end - begin:.6f} sec")
            return res
        return wrapper
    return inner

class Cache:
    '''efficiency far worse than lru_cache
        lru_cache can be 40 times faster
    '''
    def __init__(self, func:callable):
        self.func = func
        self.cache = {}

    def __call__(self, *args, **kwargs):
        if args not in self.cache:
            self.cache[args] = self.func(*args, **kwargs)
        return self.cache[args]
    
# @Cache
# @lru_cache
def _fib(n:int) -> int:
    '''true worker'''
    if n <= 1:
        return n    
    return _fib(n - 1) + _fib(n - 2)

@timeit(10)
def fib(n:int) -> int:
    '''if decorate _fib, infinite recursion happens'''
    return _fib(n)




if __name__ == "__main__":
    print(fib(35))

    import time
    from wrapper import wrapper
    begin = time.time()
    for _ in range(10):
        res = wrapper.fib(35)
    end = time.time()
    print(f"{end-begin:.6f}, {res}")