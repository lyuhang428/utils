## 限制性闭壳层Hartree-Fock量子化学计算

工作代码包括：
\- utils.py

\- molecule.py

\- mmd2.py

\- hf.py

\- core.f90


该代码目前只支持闭壳层体系，分子积分（重叠，动能，点电荷吸引，电子排斥）计算使用`McMurchie-Davidson`算法（理论山支持任意角动量量子数，但是只测试到f轨道），电子排斥积分同时支持`Rys quadrature`方法（最高仅支持d轨道）。分子积分部分使用`Fortran`进行计算，通过`f2py`模块在`Python`中调用。分子积分目前仅支持笛卡尔基函数。

计算效率方面比`Psi4`慢近20倍，但是代码易读易修改。

`Fotran`部分依赖`gsl`库进行阶乘以及合流超几何函数（Confluent hypergeometric function）计算。若该库未安装在系统路径，使用前先进行链接：`export LD_LIBRARY_PATH=/path/to/gsl_lib/`. 

使用`multiprocessing`模块以及电子排斥积分八重对称性进行并行计算。
