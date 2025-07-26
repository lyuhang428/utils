#import "@preview/unequivocal-ams:0.1.2": ams-article, theorem, proof
// #set document(author: "lyh",title: "Harmonics Oscillator") // pdf metadata
#set heading(numbering: "1.",outlined: true,bookmarked: true,)
#set par(leading: 0.55em, spacing: 1.55em, first-line-indent: 0.0em, justify: true)
#show heading: set block(above: 1.4em, below: 1em)
#set page(
  margin: 1.0in,
  numbering: "1",
  paper: "a4",
  number-align: center,
  footer: context[#set align(center)
                  #counter(page).display("1 of 1",both:true)])
#set math.equation(numbering: "(1)")
#set math.equation(block: true)
#set math.mat(delim:"[")
#show link: underline
#set text(font: ("Times New Roman", "SimSun"))
#let hbar = math.planck.reduce
#let int = math.integral
#let indent=h(2em)


// content
#align(center,text(21pt)[笔记])
#align(center,text(12pt)[lyh])
#outline(title: "Contents",indent: auto)




_bold letters stand for vectors unless otherwise stated_
#set par(first-line-indent: 2em)

= QC
== 一些概念
可以用平面波描写自由粒子, 平面波的频率与波矢与自由粒子的能量和动量相联系
$ f(bold(r), bold(k))=e^(i bold(k) dot.c bold(r)) e^(-i omega t) $

平面波的频率与波矢不随时间改变, 对应于自由粒子的能量与动量不随时间与位置改变. 一般情况下用一个函数描写粒子的波, 称该函数为波函数. 

由波函数可以得到体系的各种性质, 因此波函数(也称概率幅)描写体系的量子状态(简称状态或态)

波函数在归一化后并不完全确定, 可以用一个常数 $e^{i delta}$ 去乘以波函数. 该常数称为相因子, 归一化波函数可以含有任意相因子. 




== 表象变换
Assume there're two complete orthonormal basis $sum_(n) bold(e_n) bold(e_n)^T=1$ and $sum_n bold(u_n) bold(u_n)^T = 1$, or in Dirac notation $sum_(n) lr(| bold(e_n) lr(angle.r angle.l) bold(e_n) |) = 1$ and $sum_n lr(| bold(u_n) lr(angle.r angle.l) bold(u_n) |)=1$, where $1$ is identity matrix. A vector $bold(f)$ can then be represented as 
$ & bold(f)=sum_(n) a_n bold(e_n) \
  & bold(f)=sum_n bold(u_n) $
向量在某基下的坐标, 指的是在该基下的系数. 

Vector $bold(f)$ is invariant no matter what basis is used, but the coordinate can change. Then
$ a_n = bold(e_n) dot.c bold(f) = bold(e_n) dot.c sum_m b_m bold(u_m) = bold(e_n) dot.c bold(u_1) b_1 + bold(e_n) dot.c bold(u_2) b_2 + dots.c + bold(e_n) dot.c bold(u_n) b_n $

$ mat(
bold(e_1) dot.c bold(u_1), bold(e_1) dot.c bold(u_2), dots.c, bold(e_1) dot.c bold(u_n) ;
bold(e_2) dot.c bold(u_1), bold(e_2) dot.c bold(u_2), dots.c, bold(e_2) dot.c bold(u_n) ;
dots.v, dots.v, dots.down, dots.v ;
bold(e_n) dot.c bold(u_1), bold(e_n) dot.c bold(u_2), dots.c, bold(e_n) dot.c bold(u_n)) dot.c 
mat(b_1 ;
    b_2 ;
    dots.v ;
    b_n ) = 
mat(a_1 ;
    a_2 ;
    dots.v ;
    a_n ) $

and 
$ b_n= bold(u_n) dot.c  bold(f)= bold(u_n) dot.c  sum_(m)a_m bold(e_m)= bold(u_n) dot.c  bold(e_1)a_1+ bold(u_n) dot.c  bold(e_2)a_2+ dots.c + bold(u_n) dot.c  bold(e_n)a_n $


$ mat(
  bold(u_1) dot.c bold(e_1), bold(u_1) dot.c bold(e_2), dots.c, bold(u_1) dot.c bold(e_n) ;
  bold(u_2) dot.c bold(e_1), bold(u_2) dot.c bold(e_2), dots.c, bold(u_2) dot.c bold(e_n) ;
  dots.v, dots.v, dots.down, dots.v ;
  bold(u_n) dot.c bold(e_1), bold(u_n) dot.c bold(e_2), dots.c, bold(u_n) dot.c bold(e_n)
)
dot.c mat(
  a_1 ;
  a_2 ;
  dots.v ;
  a_n ;
)
= mat(
  b_1 ;
  b_2 ;
  dots.v ;
  b_n ;
) $

\
Vector $bold(f)$ under basis ${bold(e_n)}$ has coordinate ${a_n}$, under basis ${bold(u_n)}$ has coordinate ${b_n}$. They are connected by transformation matrix
$ & S = mat(
    bold(e_1^T) ;
    bold(e_2^T) ;
    dots.v  ;
    bold(e_n^T) ;
) dot.c mat(
   bold(u_1), bold(u_2), dots.c, bold(u_n) ;
) \
  & S^T = mat(
  bold(u_1)^T ;
  bold(u_2)^T ;
  dots.v;
  bold(u_n)^T;
 ) 
 dot.c 
 mat(
  bold(e_1), bold(e_2), dots.c, bold(e_n)
) = S^(-1) $

\
Therefore transformation matrix $S$ is orthogonal (or unitary if it's complex). In Dirac notation, making use of completeness condition gives 

$ a_n = lr(angle.l e_n mid(|)f angle.r) 
= lr(angle.l e_n mid(|) sum_m mid(|)u_m mid(angle.r) mid(angle.l) u_m mid(|) f mid(angle.r)) = sum_m lr(angle.l e_n mid(|) u_m angle.r) b_m $

$ & mat(
  angle.l bold(e_1)mid(|)bold(u_1) angle.r, angle.l bold(e_1)mid(|)bold(u_2) angle.r, dots.c, angle.l bold(e_1)mid(|)bold(u_n) angle.r ;

  angle.l bold(e_2)mid(|)bold(u_1) angle.r, angle.l bold(e_2)mid(|)bold(u_2) angle.r, dots.c, angle.l bold(e_2)mid(|)bold(u_n) angle.r ;

  dots.v, dots.v, dots.down, dots.v ;

  angle.l bold(e_n)mid(|)bold(u_1) angle.r, angle.l bold(e_n)mid(|)bold(u_2) angle.r, dots.c, angle.l bold(e_n)mid(|)bold(u_n) angle.r ;
) dot.c mat(
  bold(b_1); 
  bold(b_2); 
  dots.v; 
  bold(b_n); 
) = mat(
  bold(a_1); 
  bold(a_2); 
  dots.v; 
  bold(a_n);
) $

$
  S = mat(
    angle.l bold(e_1)|;
    angle.l bold(e_2)|;
    dots.v; 
    angle.l bold(e_n)|;
  ) dot.c mat(
    |bold(u_1)angle.r, |bold(u_2)angle.r, dots.c, |bold(u_n)angle.r; 
  )
$
Which furnishes the result that transformation matrix $S$ is orthogonal (or unitary) also.


== 算符的表象变换
Assume two basis $\{bold(e_n)mid(|)=1,2,dots.c\}$ and $\{bold(u_n)mid(|)n=1,2,dots.c\}$, operator $H$ under basis $\{bold(e_n)\}$ is $H_(i\j) = angle.l bold(e_i)mid(|)\Hmid(|)bold(e_j) angle.r = sum_(\mn) angle.l bold(e_i)mid(|)bold(u_m)angle.r angle.l bold(u_m)mid(|)\Hmid(|)bold(u_n) angle.r angle.l bold(u_n)mid(|)bold(e_j)angle.r$, then

$
  & mat(
    angle.l bold(e_1)mid(|)bold(u_1)angle.r, angle.l bold(e_1)mid(|)bold(u_2)angle.r, dots.c, angle.l bold(e_1)mid(|)bold(u_n)angle.r; 

    angle.l bold(e_2)mid(|)bold(u_1)angle.r, angle.l bold(e_2)mid(|)bold(u_2)angle.r, dots.c, angle.l bold(e_2)mid(|)bold(u_n)angle.r; 

    dots.v, dots.v, dots.down, dots.v;

    angle.l bold(e_n)mid(|)bold(u_1)angle.r, angle.l bold(e_n)mid(|)bold(u_2)angle.r, dots.c, angle.l bold(e_n)mid(|)bold(u_n)angle.r; 
  ) dot.c 
  
  mat(
    angle.l bold(u_1)mid(|)\Hmid(|)bold(u_1)angle.r, angle.l bold(u_1)mid(|)\Hmid(|)bold(u_2)angle.r, dots.c, angle.l bold(u_1)mid(|)\Hmid(|)bold(u_n)angle.r; 

    angle.l bold(u_2)mid(|)\Hmid(|)bold(u_1)angle.r, angle.l bold(u_2)mid(|)\Hmid(|)bold(u_2)angle.r, dots.c, angle.l bold(u_2)mid(|)\Hmid(|)bold(u_n)angle.r; 

    dots.v, dots.v, dots.down, dots.v;

    angle.l bold(u_n)mid(|)\Hmid(|)bold(u_1)angle.r, angle.l bold(u_n)mid(|)\Hmid(|)bold(u_2)angle.r, dots.c, angle.l bold(u_n)mid(|)\Hmid(|)bold(u_n)angle.r; 
  ) dot.c

  mat(
    angle.l bold(u_1)mid(|)bold(e_1)angle.r, angle.l bold(u_1)mid(|)bold(e_2)angle.r, dots.c, angle.l bold(u_1)mid(|)bold(e_n)angle.r; 

    angle.l bold(u_2)mid(|)bold(e_1)angle.r, angle.l bold(u_2)mid(|)bold(e_2)angle.r, dots.c, angle.l bold(u_2)mid(|)bold(e_n)angle.r; 

    dots.v, dots.v, dots.down, dots.v;

    angle.l bold(u_n)mid(|)bold(e_1)angle.r, angle.l bold(u_n)mid(|)bold(e_2)angle.r, dots.c, angle.l bold(u_n)mid(|)bold(e_n)angle.r; 
  ) \ & = 

  mat(
    angle.l bold(e_1)mid(|)\Hmid(|)bold(e_1)angle.r, angle.l bold(e_1)mid(|)\Hmid(|)bold(e_2)angle.r, dots.c, angle.l bold(e_1)mid(|)\Hmid(|)bold(e_n)angle.r; 

    angle.l bold(e_2)mid(|)\Hmid(|)bold(e_1)angle.r, angle.l bold(e_2)mid(|)\Hmid(|)bold(e_2)angle.r, dots.c, angle.l bold(e_2)mid(|)\Hmid(|)bold(e_n)angle.r; 

    dots.v, dots.v, dots.down, dots.v;

    angle.l bold(e_n)mid(|)\Hmid(|)bold(e_1)angle.r, angle.l bold(e_n)mid(|)\Hmid(|)bold(e_2)angle.r, dots.c, angle.l bold(e_n)mid(|)\Hmid(|)bold(e_n)angle.r; 
  )
$

Let 
$ 
  S = mat(
    angle.l bold(e_1)mid(|)bold(u_1)angle.r, angle.l bold(e_1)mid(|)bold(u_2)angle.r, dots.c, angle.l bold(e_1)mid(|)bold(u_n)angle.r; 

    angle.l bold(e_2)mid(|)bold(u_1)angle.r, angle.l bold(e_2)mid(|)bold(u_2)angle.r, dots.c, angle.l bold(e_2)mid(|)bold(u_n)angle.r; 

    dots.v, dots.v, dots.down, dots.v;

    angle.l bold(e_n)mid(|)bold(u_1)angle.r, angle.l bold(e_n)mid(|)bold(u_2)angle.r, dots.c, angle.l bold(e_n)mid(|)bold(u_n)angle.r; 
  )
$
then $S H_u S^T=H_e$. Here $S$ is also orthogonal or unitary matrix.

== 3D $delta$ function in real- and $k$-space
This content can also be found in "Mathematical physics integral transformation" note. 

$
  & delta(x) = 1/(2pi) int e^(i\kx) \dx \
  
  & (1/(2pi) int e^(i k_x x) \dx) (1/(2pi) int e^(i k_y y) \dy) (1/(2pi) int e^(i k_z z) \dz) = (1/(2pi))^3 int e^(i bold(k) dot.c bold(r)) d^3 bold(r) = delta(bold(k)) \

  & delta(bold(r)) = delta(x)delta(y)delta(z) \
  & delta(bold(k)) = delta(k_x)delta(k_y)delta(k_z) \

  & delta(bold(k)) = delta(k_x)delta(k_y)delta(k_z) = (1/(2pi))^3 int e^(i bold(k) dot.c bold(r)) d^3 bold(r) \

  & delta(bold(k)) = (1/(2pi))^3 int e^(i k_x x)\dx int e^(i k_y y)\dy int e^(i k_z z)\dz \
$

$
    & delta(hbar bold(k)) = (1/(2pi))^3 int e^(i hbar k_x x)\dx int e^(i hbar k_y y)\dy int e^(i hbar k_z z)\dz \

    & = (1/hbar)^3 (1/(2pi))^3 int e^(\ik_x hbar x) \dhbar\x int e^(i k_y hbar y) d hbar y int e^(i k_z hbar z) d hbar z \

    & = (1/hbar)^3 (1/(2pi))^3 int e^(i k_x X) \dX int e^(i k_y Y) \dY int e^(i k_Z Z) \dZ \
    & = (1/hbar)^3 delta(bold(k)) \
$



== Derivative of $delta$ function
$ 
  & int f(x) delta^prime (-x) \dx = -int f(x) \ddelta(-x) = -f(x) delta(-x)|_(-infinity)^(infinity) + int delta(-x) f^prime (x) \dx = f^prime (0) \

  & int f(x)delta^prime (x) \dx = -f^prime (0)  \
  & arrow.double -delta^prime (x) = delta^prime (-x) \
$

$
  & int f(x)delta^((2)) (-x) \dx = -int f(x) d/(\dx) delta^prime (-x) \dx \
  & = int f(x) d/(\dx)delta^prime (x) \dx = int f(x) delta^((2)) (x) \dx \
  &  = f^((2))(0) \
  \
  & int f(x) delta^((2)) (x) \dx = f^((2))(0) \
  & arrow.double delta^((2))(-x) = delta^((2)) (x) \
  \
  & int f(x) delta^((3)) (-x) \dx = -int f(x) d/(\dx) delta^((2)) (-x) \dx \
  & = -int f(x) d/(\dx) delta^((2)) (x) \dx \
  & = -int f(x) delta^((3)) (x) \dx \
  & = f^((3))(0) \
  \
  & int f(x) delta^((3)) (x) \dx = -f^((3)) (0) arrow.double delta^((3)) (-x) = -delta^((3)) (x) \
  & dots.v \
  & delta^((n)) (-x) = (-1)^n delta^((n)) (x) \
$

An example that makes use of the derivative of $delta$ function
$
  1/(2pi) int k^2 e^(i k x) \dk = 1/(2pi) (-partial^2/(partial x^2) int e^(i k x)\dk) = -partial^2/(partial x^2) delta(x) = -delta^((2)) (x)
$



== 坐标表象和动量表象
The position operator is just position itself $hat(r)=r$, the momentum operator is $hat(p)=-i hbar nabla$, for a single component is $hat(p_x)=-i hbar d/(\dx)$. The eigenfunction for momentum operator is $f_(k) (x) = e^(i k_x x)$, and $f_(k)(bold(r)) = A e^(i bold(k) dot.c bold(r))$. $f_k (bold(r))$是动量算符的特征函数, 但是其本身是位矢的函数. $hat(bold(p)) f_k (bold(r)) = hbar bold(k) f_k (bold(r))$, $ hat(bold(p)) f_(k^prime)(bold(r)) = hbar k^prime f_(k^prime) (bold(r))$. $bold(k)$ and $bold(k^prime)$ indicate they are eigenfunctions of different eigenvalues ($bold(k)$ and $bold(k^prime)$) of momentum operator. 其本身当然满足正交性 (考察的是动量特征函数(对应不同的特征值)的正交性, 积掉的是空间坐标!)

$
  & f_bold(k) (bold(r)) = A e^(i bold(k) dot.c bold(r)) \
  & f_bold(k^prime) (bold(r)) = A e^(i bold(k^prime) dot.c bold(r)) \
  & angle.l f_bold(k) (bold(r))mid(|)f_bold(k^prime) (bold(r)) angle.r = |A|^2 int_(bb(R^3)) e^(i (bold(k^prime) - bold(k))dot.c bold(r)) d^3 bold(r)  = (2pi)^3 |A|^2 delta(bold(k^prime) - bold(k)) \
$
这里的积分范围为全空间. 

已知$delta(hbar bold(k)) = delta(bold(k))/hbar^3$, 且$hbar^3 delta(bold(p)) = delta(bold(k)) = delta(bold(p)/hbar)$

$
  & f_bold(k) (bold(r)) = A e^(i bold(k) dot.c bold(r)) \
  & f_bold(k^prime) (bold(r)) = A e^(i bold(k^prime) dot.c bold(r)) \
  & int f^ast_bold(k^prime) (bold(r)) f_bold(k) (bold(r)) d^3bold(r) = int |A|^2 e^(i (bold(k) - bold(k^prime)) dot.c bold(r)) d^3 bold(r) \
  & = int |A|^2 e^(i (bold(p) - bold(p^prime))/(hbar) dot.c bold(r)) d^3 bold(r) \
  & = |A|^2 (2pi)^3 delta((bold(p) - bold(p^prime))/(hbar)) \
  & = |A|^2 (2pi hbar)^3 delta(bold(p) - bold(p^prime)) \
$

特征函数的正教性表现为Dirac $delta$函数，原因在于这里的积分范围是全空间($-infinity$\~$infinity$). 对这里的特征函数要求其归一化为$delta$函数，为非归一化为1
$ |A|^2 (2pi hbar)^3 = 1 arrow.double A = (1/(2pi hbar))^(3/2) $


== 箱归一化
如果考虑周期边界条件. 注意这里的自由粒子并非被无限深势井限制在边长为$L$的盒子中，这里只是要求周期边界条件
$
  cases(
    f_bold(k) (x+L,y,z) = f_bold(k) (x,y,z) \
    f_bold(k) (x,y+L,z) = f_bold(k) (x,y,z) \
    f_bold(k) (x,y,z+L) = f_bold(k) (x,y,z)
  )
$
考虑$x$-分量 $e^(i k_x x) = e^(i k_x (x+L))$, 那么$k_x L = 2pi n_x arrow.double k_x = (2pi)/(L) n_x$, 所以$bold(k) = (2pi)/(L) bold(n)$, 其中$bold(n) = mat(n_x; n_y; n_z)$, $n_i in bb(Z)$. 特征值$bold(k)$不再是连续的，而是离散的. 

在周期边界条件下，动量的特征函数要求归一到1
$
  & f_bold(k) (bold(r)) = A e^(i bold(k) dot.c bold(r)) \
  & f_bold(k^prime) (bold(r)) = A e^(i bold(k^prime) dot.c bold(r)) \
  & k_x = (2pi)/(L) n_x \
  & k_x - k_(x^prime) = (2pi)/(L) n \
  & angle.l f_bold(k^prime) (bold(r))mid(|)f_bold(k) (bold(r))angle.r = |A|^2 int e^(i (bold(k) - bold(k^prime)) dot.c bold(r)) d^3 bold(r) \ 
  & = |A|^2 (int_0^L e^(i (k_x - k_(x^prime))x) \dx) (int_0^L e^(i (k_y - k_(y^prime))y) \dy) (int_0^L e^(i (k_z - k_(z^prime))z) \dz) \
$
注意指数上$k_i - k_(i^prime) = (2pi)/(L) n_i$, $(2pi)/(omega) = L/n$, 盒子的周期$L$是平面波周期$L/n$的整数倍. 又有$1/(2pi) int_0^(2pi) e^(i (n-m) theta) d theta = delta_(\nm)$, 所以$int_0^L e^(i(k_x - k_(x^prime))x) \dx = L delta_(k,k^prime)$. 带入上式可以得到$|A|^2 L^3 delta_(k_x,k_(x^prime)) delta_(k_y,k_(y^prime)) delta_(k_z,k_(z^prime))$, 归一化后得到$A = (1/L)^(3/2)$. 最后得到归一化的满足周期边界条件的动量算符特征函数
$ f_bold(k) (bold(r)) = (1/L)^(3/2) e^(i bold(k) dot.c bold(r)) $



== 几点评论
几点注意事项, 积分上下限只需要差值为 $L$ 即可. 由于周期边界条件的限制, 动量特征值不再构成连续谱, 而只能取分立的值, 因此这里出现的是 Kronecker $delta$, 而动量特征值取连续谱时出现的是 Dirac $delta$. Kronecker $delta$ is discrete version of Dirac $delta$.

在 $hat(O)$ 的特征值构成分立谱的情况下, 特征函数归一化为 Kronecker $delta$. 
$ angle.l f_k mid(|) f_(k^prime) angle.r = delta_(k,k^prime) $

在 $hat(O)$ 的特征值构成连续谱的情况下, 特征函数归一化为 Dirac $delta$
$ angle.l f_k mid(|) f_(k^prime) angle.r = delta(k-k^prime) $

坐标算符与动量算符在非周期边界条件下的特征值构成连续谱, 他们的特征函数正交归一化为 Dirac $delta$ 函数. 在坐标表象下(意指自变量为坐标)的动量算符的特征函数(为简单起见, 接下来只考虑一维情况, 三维情形只是一维的结果乘起来而已)为
$ f_k (x) = sqrt(1/(2pi)) e^(i k x) $
或者
$ f_p (x) = sqrt(1/(2pi hbar)) e^(i p/hbar x) $
该记号指在有确定波矢 $k$ 或确定动量 $p$ 时的动量算符的特征函数在坐标表象下的表现形式. 

动量特征函数归一化为$delta$函数（积掉坐标）
$ angle.l f_k (x)mid(|)f_(k^prime) (x) angle.r = delta(k - k^prime) $
或者
$ angle.l f_p (x)mid(|)f_(p^prime) (x) angle.r = delta(p-p^prime) $


在动量表象下，动量算符就是动量本身 $hat(p) = p$. 动量算符特征函数的自变量在动量表象下为动量$p$
$ hat(p) f_(p_i) (p) = p_i f_(p_i)(p) $
$f_(p_i)(p)$ 指在有确定动量 $p_i$ 时的动量的特征函数. 注意到有 $x delta(x-x_i) = x_i$，得到在动量表象下的动量算符的特征函数为
$ f_(p_i)(p) = delta(p-p_i) $
其归一化到 Dirac $delta$ 函数（在动量表象下自变量为动量，积掉动量）
$ angle.l f_(p_i)(p)mid(|)f_(p_j)(p) angle.r = int delta(p-p_i) delta(p-p_j) \dp = delta(p_i - p_j) $

同理, 在坐标表象下, 坐标算符就是坐标本身, 坐标算符的特征函数在坐标表象下同样也是 Dirac $delta$ 函数
$ x delta(x-x^prime) = x^prime delta(x-x^prime) $

在自身表象下, 动量与坐标的特征函数都是对角的. 


== 动量表象下的位置算符
在纯粹数学领域, 傅里叶变换的系数的选择是人为的, 只要正向与反向的系数的乘积得到 $1/(2pi)$ 即可 (1D)
$
  & F(k) = A int f(x) e^(-i k x) \dx \
  & f(x) = B int F(k) e^(i k x) \dk \
  & A dot.c B = 1/(2pi) \
$

在量子力学领域, 由于波函数归一化的限制, 傅里叶变换的系数的选择不是随意的. 位置算符 $hat(x)$ 的特征函数在位置表象 (实空间表象), 与在动量表象 ($k$ 空间表象) 构成傅里叶变换对
$
  & cases(
    f(k) = sqrt(1/(2pi)) int f(x) e^(-i k x) \dx \
    f(x) = sqrt(1/(2pi)) int f(k) e^(i k x) \dk 
  ) \
  & cases(
    f(p) = sqrt(1/(2pi hbar)) int f(x) e^(-i p/hbar x) \dx \
    f(x) = sqrt(1/(2pi hbar)) int f(k) e^(i p/hbar x) \dk 
  )
$

相应的归一化要求
$
  & angle.l f(k)mid(|)f(k) angle.r = int 1/(2pi) int f^ast (x) e^(i k x) \dx dot.c int f(x^prime) e^(-i k x^prime) d x^prime \dk \
  & = int (1/(2pi) int e^(i k (x-x^prime)) \dk) f^ast (x) f(x^prime) d x^prime \dx \
  & = int delta(x-x^prime) f^ast (x) f(x^prime) d x^prime \dx \
  & = int |f(x)|^2 \dx = 1 \
  \
  & angle.l f(x)mid(|)f(x) angle.r = 1 \
$


位置算符在位置表象下就是坐标本身, 位置算符期望值在位置表象下为
$ angle.l x angle.r = int f^ast (x) hat(x) f(x) \dx = int f^ast (x) x f(x) \dx $

同理, 动量算符在动量表象下就是动量本身. 但是位置算符在动量表象, 以及动量算符在位置表象下具有不同地形式. 在位置表象下动量算符的表达式为 $hat(p) = i hbar d/(\dx)$.

由于期望值不因表象的改变而改变, 故而位置算符的期望值在动量表象下应与在位置表象下等同, 而位置算符在动量表象下并不等于位置本身 (动量表象下的自变量为动量)
$
  & angle.l x angle.r = int f^ast (k) hat(x) f(k) \dk \
  & = int (f^ast (x) e^(i k x)) hat(x) (f(x^prime) e^(-i k x^prime)) \dx \dx^prime \dk \
  & = int f^ast (x) dot.c (int x^prime f(x^prime) delta(x-x^prime) \dx^prime) \dx \
  & = int f^ast (x) x f(x) \dx \
  & arrow.double hat(x) e^(-i k x^prime) = x^prime e^(-i k x^prime) \
$
注意这里 $hat(x)$ 并不作用在 $f(x^prime)$, 因为 $f(x^prime)$ 不是 $k$ 的函数！最后有
$
  & hat(x) = i d/(\dk) \
  & arrow.l.r.double hat(x) = i hbar d/(\dp) \
$

$
  & angle.l x angle.r = 1/(2pi hbar) int f^ast (p) e^(i p/hbar x) hat(x) f(x^prime) e^(-i p/hbar x^prime) \dx^prime \dx \
  & int f^ast (x) x^prime f(x^prime) delta(x-x^prime) \dx^prime \dx \
  & arrow.double hat(x) = i hbar d/(\dp) \
$

For completeness, the momentum operator under position representation is $hat(p) = -i hbar d/(\dx)$.

== A little comments on Dirac $delta$ function
Full description of Dirac $delta$ function can be referred to another notes.

$
  & 1/(2pi) int e^(i (p-p^prime)/hbar x) \dx = hbar delta(p-p^prime) #h(1em) text("scaling property") \
  & 1/(2pi hbar) int e^(i (p-p^prime)/hbar x) \dx = delta(p-p^prime) \
  & 1/(2pi hbar) int e^(i p/hbar x) \dx = delta(p) \
  & 1/(2pi hbar) int e^(i p/hbar x) \dp = delta(x) \
$


= 自由粒子正则系综的密度矩阵
For brevity, here only 1D situation is considered. For a single particle (therefore no particle indistingusable issue) in cubic PBC box with length $L$, the Hamiltonian and corresponding eigenfunctions are
$
  & hat(H) = - hbar^2/(2m) d^2/(\dx^2) \
  & f_k (x) = sqrt(1/L) e^(i k x) \
  & epsilon_k = (hbar^2 k^2)/(2m) \
  & k = (2pi)/L n, #h(1em) n in bb(Z) \
$
$f_k (x)$ 指动量为 $k$ 的，位置表象下的，动量算符的特征函数. 动量算符的特征函数在位置算符特征函数上的投影（即系数）就是位置表象下动量算符的特征函数
$
  angle.l \rmid(|)k angle.r = f_k (r) = sqrt(1/L) e^(i k x)
$

正则系综下的密度算符表达式
$
  & hat(rho) = (e^(-beta hat(H)))/(text("Tr")\{e^(-beta hat(H))\}) \
$
接下来考察密度算符在位置表象下的密度算符矩阵元
$
  & angle.l \rmid(|)e^(-beta hat(H))mid(|)r^prime angle.r = sum_k angle.l\rmid(|)e^(-beta hat(H))mid(|)k angle.r angle.l\kmid(|)r^prime angle.r \
  
  & = sum_k e^(-(beta hbar^2 k^2)/(2m)) angle.l\rmid(|)k angle.r angle.l\kmid(|)r^prime angle.r \

  & = sum_k e^(-(beta hbar^2 k^2)/(2m)) f_k (r) f_(k)^ast (r^prime) \
  & = sum_k e^(-(beta hbar^2 k^2)/(2m)) 1/L e^(i k (x-x^prime))
$

考虑态密度的计算（箱归一化，波矢只能取分立值）
$
  & int_xi^(xi+L) (\dx \dp)/h = L/h hbar \dk = L/(2pi) \dk \
  & sum_k approx L/(2pi) int \dk \
$
对动量的求和过渡到对动量的积分有系数上的差异！将上式带入
$
  & sum_k e^(-(beta hbar^2 k^2)/(2m)) 1/L e^(i k (x-x^prime)) approx L/(2pi) int e^(-(beta hbar^2 k^2)/(2m)) 1/L e^(i k (x-x^prime)) \dk \
  & = 1/(2pi) int e^(-(beta hbar^2 k^2)/(2m)) e^(i k (x-x^prime)) \dk \
  & = sqrt((m)/(2pi beta hbar^2)) e^(-(m)/(2 beta hbar^2) (x-x^prime)^2) \
  & angle.l\rmid(|)e^(-beta hat(H))mid(|)r^prime angle.r = sqrt((m)/(2pi beta hbar^2)) e^(-(m)/(2 beta hbar^2) (x-x^prime)^2) \
$

三维情形下为
$
  angle.l bold(r)mid(|)e^(-beta hat(H))mid(|)bold(r^prime) angle.r = ((m)/(2pi beta hbar^2))^(3/2) e^(-(m)/(2 beta hbar^2) |bold(r)-bold(r)^prime|^2) \
$

The trace of it is then 
$
  & text("Tr")\{e^(-beta hat(H))\} = int_xi^(xi+L) angle.l\rmid(|)e^(-beta hat(H))mid(|)r angle.r \dx \
  & = sqrt(m/(2pi beta hbar^2)) int_xi^(xi+L) \dx \
  & = L sqrt(m/(2pi beta hbar^2)) \
$
因此密度矩阵的矩阵元为
$
  angle.l\rmid(|)rho mid(|)r^prime angle.r = (angle.l\rmid(|)e^(-beta hat(H))mid(|)r^prime angle.r)/(text("Tr")\{e^(-beta hat(H))\}) = 1/L e^(-(m)/(2beta hbar^2)(x-x^prime)^2)
$


Finally let's evaluate the expectation value of Hamiltonian
$
  & angle.l\Hangle.r = text("Tr")\{hat(H)rho\} = sum_(r_i) angle.l\r_i mid(|)hat(H)rho mid(|)r_i angle.r \
  
  & angle.l\r_i mid(|)hat(H)rho mid(|)r_i angle.r = sum_(r_j)angle.l\r_i mid(|)hat(H)mid(|)r_j angle.r angle.l\r_j mid(|)rho mid(|)r_i angle.r \

  & |r_j angle.r = delta(x-x_j) \

  & angle.l\r_i mid(|)hat(H)mid(|)r_j angle.r = int delta(x-x_i) (-(hbar^2)/(2m)) d^2/(\dx^2) delta(x-x_j) \dx \
  & = -(hbar^2)/(2m) int delta(x-x_i) delta^((2))(x-x_j) \dx  = -(hbar^2)/(2m) delta^((2))(x_j - x_i) \

  & angle.l\r_j mid(|)rho mid(|)r_i angle.r = 1/L text("exp")[-(m)/(2beta hbar^2) (x_j - x_i)^2] \

  & sum_(r_j) angle.l\r_i mid(|)hat(H)mid(|)r_j angle.r angle.l\r_j mid(|)rho mid(|)r_i angle.r = int -(hbar^2)/(2m)delta^((2))(x_j - x_i) 1/L text("exp")[-(m)/(2beta hbar^2)(x_j - x_i)^2] \dx_j \

  & = -(hbar^2)/(2m L) int delta^((2))(x_j - x_i) text("exp")[-(m)/(2beta hbar^2)(x_j - x_i)^2] \dx_j \
  & = -(hbar^2)/(2m L) d^2/(\dx_j^2) text("exp")[-(m)/(2beta hbar^2)(x_j - x_i)^2] |_(x_j = x_i) \
  & = -(hbar^2)/(2m L) { (m^2(x_j - x_i)^2)/(beta^2 hbar^4) text("exp")[-(m)/(2beta hbar^2)(x_j - x_i)^2] - m/(beta hbar^2)text("exp")[-(m)/(2beta hbar^2)(x_j - x_i)^2]}|_(x_j = x_i) \
  & = 1/(2beta L) \

  & angle.l hat(H) angle.r = sum_(r_i)sum_(r_j)angle.l\r_i mid(|)hat(H)mid(|)r_j angle.r  angle.l\r_j mid(|)rho mid(|)r_i angle.r \
  & = int_xi^(xi+L) 1/(2beta L) \dx_i = (k_B T) / 2 \
$

上述结果正是能量均分定理. 自由粒子在一维情况只有一个动能二次项, 对总能量的贡献为 $(k_B T)/2$. 

通过内插动量算符在位置表象下的特征函数完备集 (${|\kangle.r}$), 可以从另一个角度计算哈密顿量的期望值
$
  & angle.l x_\imid(|)rho mid(|)x_j angle.r = 1/L text("exp")[-(m)/(2beta hbar^2)(x_i - x_j)^2] \

  & angle.l x_\imid(|)hat(H)rho mid(|)x_i angle.r = sum_(x_j) angle.l x_\imid(|)\Hmid(|)x_j angle.r angle.l x_\jmid(|)rho mid(|)x_i angle.r \

  & angle.l x_\imid(|)hat(H)mid(|)x_j angle.r = sum_k angle.l x_\imid(|)hat(H)mid(|)\kangle.r angle.l\kmid(|)x_j angle.r \

  & = sum_k (hbar^2 k^2)/(2m) angle.l x_\imid(|)k angle.r angle.l\kmid(|)x_\jangle.r \

  & = (hbar^2 k^2)/(2m) sum_k 1/L text("exp")[i k (x_i - x_j)] \

  & approx (hbar^2)/(2m) L/(2pi) 1/L int k^2 text("exp")[i k (x_i - x_j)] \dk \
  & = -(hbar^2)/(2m) 1/(2pi) int d^2/(d x_i^2) text("exp")[i k (x_i - x_j)] \dk \
  & = -(hbar^2)/(2m) d^2/(d x_i^2) 1/(2pi) int text("exp")[i k (x_i - x_j)] \dk \
  & = -(hbar^2)/(2m) d^2/(d x_i^2) delta(x_i - x_j) \
$
剩余步骤同上. 

如果在动量表象下计算则会大大简化, 因为在动量表象下密度算符和哈密顿算符都是对角的
$
  & angle.l\kmid(|)e^(-beta hat(H))mid(|)k^prime angle.r = text("exp")[-(beta hbar^2 k^2)/(2m)] delta_(k,k^prime) = text("exp")[-(beta hbar^2 k^2)/(2m)] \

  & sum_k text("exp")[-(beta hbar^2 k^2)/(2m)] approx L/(2pi) int text("exp")[-(beta hbar^2 k^2)/(2m)] \dk \
  
  & text("Tr"){e^(-beta hat(H))} = L/(2pi) sqrt((2pi m)/(beta hbar^2)) \
  & angle.l\H angle.r = sum_k angle.l\kmid(|)hat(H) rho mid(|)k angle.r \
  & = sum_k sum_(k^prime) angle.l\kmid(|)hat(H)mid(|)k^prime angle.r angle.l\k^prime mid(|)rho mid(|)k angle.r \
  & = sum_k sum_(k^prime) (hbar^2 k^prime^2)/(2m) delta_(k,k^prime) text("exp")[-(beta hbar^2 k^prime^2)/(2m)] delta_(k,k^prime) (2pi)/L sqrt((beta hbar^2)/(2pi m)) \

  & = sum_k (hbar^2  k^2)/(2m) text("exp")[-(beta hbar^2 k^2)/(2m)] 1/L sqrt((2pi beta hbar^2)/(m)) \

  & approx L/(2pi) int (hbar^2 k^2)/(2m) text("exp")[-(beta hbar^2 k^2)/(2m)] 1/L sqrt((2pi beta hbar^2)/(m)) \

  & = 1/(2beta) \
$




= TIPT
== 非简并情况
总哈密顿量可以分为两部分
$ H = H_0 + H^prime $

其中 $H_0$ 是无微扰的哈密顿量，其特征函数和特征值可解，$H^prime$ 是微扰项. 引入一个实参数 $lambda in [0,1]$，通过改变该参数的值可以控制微扰的强度
$ H = H_0 + lambda H^prime $

$lambda$ 最终会被设置成1,即微扰完全打开. 受扰动的特征函数和特征值可以表示成参数 $lambda$ 的幂指数展开
$
  & f = f_0 + lambda f^((1)) + lambda f^((2)) + dots.c \
  & epsilon = epsilon_0 + lambda epsilon^((1)) + lambda epsilon^((2)) + dots.c \
$
其中 $f^((n))$ 与 $epsilon^((n))$ 分别是第 $n$ 级近似特征函数和特征值. $f_0$ 和 $epsilon_0$ 是无微扰情况下定态薛定谔方程的解：$H_0 f_0 = epsilon_0 f_0$. 将微扰展开式带入定态薛定谔方程
$
  (H_0 + H^prime) (f^((0)) + lambda f^((1)) + lambda f^((2)) + dots.c) \
  = (epsilon^((0)) + lambda epsilon^((1)) + lambda epsilon^((2)) + dots.c) (f^((0)) + lambda f^((1)) + lambda f^((2)) + dots.c)
$

注意到左右两边 $lambda$ 的幂次应相等 (这也是引入 $lambda$ 的原因), 有
$
  & (H_0 - epsilon_0) f_0 = 0 #h(1em)text("TISE") \
  & (H_0 - epsilon_0) f^((1)) = (epsilon^((1)) - H^prime)f_0 \
  & (H_0 - epsilon_0) f^((2)) = (epsilon^((1)) - H^prime)f^((1)) + epsilon^((2)) f_0 \
  & (H_0 - epsilon_0) f^((3)) = (epsilon^((1)) - H^prime)f^((2)) + epsilon^((2)) f^((1)) + epsilon^((3))f_0 \

  & (H_0 - epsilon_0) f^((4)) =  (epsilon^((1)) - H^prime)f^((3)) + epsilon^((2)) f^((2)) + epsilon^((3))f^((1)) + epsilon^((4))f_0 \

  & dots.v \
  & (H_0 - epsilon_0) f^((n)) = (epsilon^((1)) - H^prime) f^((n-1)) + sum_(i=2)^n epsilon^((i)) f^((n-i)) \
$

左乘 $angle.l f_0|$ 可以得到能量的各级修正（注意有 $angle.l f_0mid(|)H_0mid(|)f^((n))angle.r = angle.l f^((n))mid(|)H_0f_0angle.r^dagger = epsilon_0 angle.l f^((n))mid(|)f_0angle.r = 0$, $f_0$ 与 $f^((n))$ 彼此正交)
$
  & angle.l f_0mid(|)H_0 - epsilon_0mid(|)f^((n))angle.r = angle.l f_0mid(|)epsilon^((1)) - H^prime mid(|)f^((n-1)) angle.r + sum_(i=2)^n angle.l f_0mid(|)epsilon^((i))mid(|)f^((n-i)) angle.r \
  & 0 = -angle.l f_0mid(|)H^prime mid(|)f^((n-1))angle.r + epsilon^((n)) \
  & arrow.double epsilon^((n)) = angle.l f_0mid(|)H^prime mid(|)f^((n-1))angle.r \
$
可以直接得到能量的一阶修正
$
  epsilon^((1)) = angle.l f_0mid(|)H^prime mid(|)f_0angle.r \
$

无微扰下的哈密顿算符特征函数构成完备集
$
  & H_0 phi_n = epsilon_n phi_n,#h(1em)n=1,2, dots.c \
  & sum_n |phi_n angle.r angle.l phi_n| = 1 \
$

假设无微扰时系统处在第 $k$ 个能级, 即 $f_0 = phi_k$, $epsilon_0 = epsilon_k$. 能量的一阶修正为 $epsilon^((1)) = angle.l phi_k mid(|)H^prime mid(|)phi_k angle.r = H^prime_(\kk)$. 同时一阶微扰近似波函数可以由其展开
$ f^((1))= sum_n a_n phi_n $

带入到 $(H_0 - epsilon_0)f^((1)) = (epsilon^((1)) - H^prime)f_0$ 有
$ (H_0 - epsilon_k) sum_n a_n phi_n = (epsilon^((1)) - H^prime)phi_k $

左乘 $angle.l phi_m|$ 有
$
 & sum_n {angle.l phi_m mid(|) H_0 mid(|) a_n phi_n angle.r - epsilon_k angle.l phi_m mid(|) a_n phi_n angle.r} = epsilon^((1)) angle.l phi_m mid(|)phi_k angle.r - angle.l phi_m mid(|) H^prime mid(|) phi_k angle.r \

 & sum_n {a_n epsilon_n delta_(\mn) - epsilon_k a_n delta_(\mn)} = H^prime_(\kk) delta_(\mk) - H^prime_(\mk) \

  & a_m epsilon_m - epsilon_k a_m = H^prime_(\kk) delta_(\mk) - H^prime_(\mk) \
$

当 $m eq.not k$, 有
$ a_m = (H^prime_(\mk))/(epsilon_k - epsilon_m) $

波函数的一阶修正为
$ f^((1)) =  sum_(n eq.not k) (H^prime_(\nk))/(epsilon_k - epsilon_n)phi_n $


Since the second order correction for energy is $epsilon^((2)) = angle.l\f_0mid(|)H^prime mid(|)f^((1)) angle.r$, where $f_0 = phi_k$. Plugging formula for first order correction wavefunction $f^((1))$ into it furnishes
$
  & epsilon^((2)) = angle.l phi_\kmid(|)H^prime mid(|)sum_(n eq.not k) (H^prime_(\nk))/(epsilon_k - epsilon_n)phi_n angle.r \

  & = sum_(n eq.not k) (H^prime_(\nk))/(epsilon_k - epsilon_n) angle.l phi_\kmid(|) H^prime mid(|)phi_n angle.r \

  & = sum_(n eq.not k) (|H^prime_(\nk)|^2)/(epsilon_k - epsilon_n) \
$

The third order energy correction is $epsilon^((3)) = angle.l f_0mid(|) H^(prime)mid(|) f^((2))angle.r$. Making use of the relation above
$
  & (H_0 - epsilon_0)f^((1)) = (epsilon^((1)) - H^prime)f_0 \

  & (H_0  - epsilon_0)f^((2)) = (epsilon^((1)) - H^prime)f^((1)) + epsilon^((2))f_0 \

  & angle.l f^((2))| text("for first") text("and") angle.l f^((1))| text("for the second") \

  & angle.l f^((2)) mid(|) H_0 - epsilon_0 mid(|) f^((1)) angle.r = -angle.l f^((2)) mid(|) H^prime mid(|) f_0 angle.r \

  & angle.l f^((1)) mid(|) H_0 - epsilon_0 mid(|) f^((2)) angle.r = angle.l f^((1)) mid(|) epsilon^((1)) - H^prime mid(|) f^((1)) angle.r \

  & arrow.double epsilon^((3)) = angle.l f^((1)) mid(|) H^prime - epsilon^((1)) mid(|) f^((1)) angle.r \
$

Which means the third order energy correction can be obtained from first order correction wavefunction, instead of the second order. Then we have 
$
  & epsilon^((3)) = angle.l f^((1)) mid(|) H^prime - epsilon^((1)) mid(|) f^((1)) angle.r \

  & = sum_(n eq.not k) sum_(m eq.not k) (H^prime_(\nk))/(epsilon_k - epsilon_n) (H^prime_(\mk))/(epsilon_k - epsilon_m) angle.l phi_n mid(|) H^prime - epsilon^((1)) mid(|) phi_m angle.r \

  & = sum_(n eq.not k \ m eq.not k) (H^prime_(\nk) H^prime_(\mk))/((epsilon_k - epsilon_n)(epsilon_k - epsilon_m)) {angle.l phi_n mid(|) H^prime mid(|) phi_m angle.r - H^prime_(\kk) delta_(\nm) } \

  & = sum_(n eq.not k \ m eq.not k) (H^prime_(\nk) H^prime_(\mk) H^prime_(\nm))/((epsilon_k - epsilon_n)(epsilon_k - epsilon_m)) - sum_(n eq.not k) (H^prime_(\kk) mid(|)H^prime_(\nk)mid(|)^2)/((epsilon_k - epsilon_n)^2) \
$

做一个小结: 无微扰时系统处于第 $k$ 个能级上, 无微扰时的波函数为 (零阶修正) $f_0 = phi_k$, 能量为 $epsilon_0 = epsilon_k$. 定态微扰的一阶能量修正, 一阶波函数修正, 二阶能量与三阶能量修正为
$
  & epsilon^((1)) = H^prime_(\kk) \

  & epsilon^((2)) = sum_(n eq.not k) (mid(|)H^prime_(\nk)mid(|))/(epsilon_n - epsilon_k) \

  & epsilon^((3)) = sum_(n eq.not k \ m eq.not k) (H^prime_(\nk) H^prime_(\mk)H^prime_(\nm))/((epsilon_k - epsilon_n)(epsilon_k - epsilon_m)) - sum_(n eq.not k) (H^prime_(\kk)mid(|) H^prime_(\nk)mid(|)^2)/((epsilon_k - epsilon_n)^2) \

  & f^((1)) = sum_(n eq.not k) (H^prime_(\nk))/(epsilon_k - epsilon_n) phi_n \
$

又有 $(H_0 - epsilon_0)f^((2)) = (epsilon^((1)) - H^prime)f^((1)) + epsilon^((2))f_0$, 可以得到波函数的二阶修正
$
  & (H_0 - epsilon_k)f^((2)) = (H^prime_(\kk) - H^prime) sum_(n eq.not k) (H^prime_(\nk))/(epsilon_k - epsilon_n)phi_n + sum_(n eq.not k) (mid(|)H^prime_(\nk)mid(|)^2)/(epsilon_k  - epsilon_n)phi_k \

  & sum_l (H_0 - epsilon_k)b_l phi_l = (H^prime_(\kk) - H^prime) sum_(n eq.not k) (H^prime_(\nk))/(epsilon_k - epsilon_n)phi_n + sum_(n eq.not k)(mid(|)H^prime_(\nk)mid(|)^2)/(epsilon_k - epsilon_n)phi_k \

  & sum_l (b_l epsilon_l phi_l - b_l epsilon_k phi_l) = (H^ prime_(\kk) - H^prime) sum_(n eq.not k)(H^prime_(\nk))/(epsilon_k - epsilon_n)phi_n + sum_(n eq.not k)(mid(|)H^prime_(\nk)mid(|)^2)/(epsilon_k - epsilon_n)phi_k \

  & sum_l angle.l phi_m mid(|) b_l epsilon_l - b_l epsilon_k mid(|) phi_l angle.r = sum_(n eq.not k) (H^prime_(\nk))/(epsilon_k - epsilon_n) angle.l phi_m mid(|) H^prime_(\kk) - H^prime  mid(|) phi_n angle.r + sum_(n eq.not k) (mid(|) H^prime_(\nk)mid(|)^2)/(epsilon_k  - epsilon_n)delta_(\mk) \

  & sum_l (b_l epsilon_l delta_(\ml) -b_l epsilon_k delta_(\ml)) = sum_(n eq.not k) (H^prime_(\nk))/(epsilon_k - epsilon_n) \(H^prime_(\kk) delta_(\mn) - H^prime_(\mn)) + sum_(n eq.not k)(mid(|) H^prime_(\nk)mid(|)^2)/(epsilon_k - epsilon_n) delta_(\mk) \

  & b_m(epsilon_m - epsilon_k) = sum_(n eq.not k) (H^prime_(\nk)H^prime_(\kk))/(epsilon_k - epsilon_n)delta_(\mn) - sum_(n eq.not k) (H^prime_(\nk) H^prime_(\mn))/(epsilon_k - epsilon_n) + sum_(n eq.not k) (mid(|) H^prime_(\nk)mid(|)^2)/(epsilon_k - epsilon_n)delta_(\mk) \
$


When $m=k$, we only get trivial result. When $m eq.not k$
$
  & b_m (epsilon_m - epsilon_k) = sum_(n eq.not k) (H^prime_(\nk)H^prime_(\kk))/(epsilon_k - epsilon_n)delta_(\mn) - sum_(n eq.not k) (H^prime_(\nk)H^prime_(\mn))/(epsilon_k - epsilon_n) + sum_(n eq.not k) (mid(|) H^prime_(\nk)mid(|)^2)/(epsilon_k - epsilon_n) delta_(\mk) \

  & b_m (epsilon_m - epsilon_k) = (H^prime_(\mk)H^prime_(\kk))/(epsilon_k - epsilon_m) - sum_(n eq.not k) (H^prime_(\nk)H^prime_(\mn))/(epsilon_k - epsilon_n) +0 \

  & b_m = sum_(n eq.not k) (H^prime_(\nk)H^prime_(\mn))/((epsilon_k - epsilon_n)(epsilon_k - epsilon_m)) - (H^prime_(\mk)H^prime_(\kk))/(epsilon_k - epsilon_m)^2 \
$

and finally
$ f^((2)) =  sum_(m eq.not k){ sum_(n eq.not k) (H^prime_(\nk)H^prime_(\mn))/((epsilon_k - epsilon_n)(epsilon_k - epsilon_m)) - (H^prime_(\mk)H^prime_(\kk))/(epsilon_k - epsilon_m)^2}phi_m $

Using similar trick, the third order wavefunction correction is (maybe wrong)
$ 
  & f^((3)) = sum_(p eq.not k){
   sum_(m eq.not k \ n eq.not k) (H^prime_(\nk)H^prime_(\mn)H^prime_(\pm))/((epsilon_k - epsilon_n)(epsilon_k - epsilon_m)(epsilon_k - epsilon_p)) -  sum_(n eq.not k) (H^prime_(\kk)H^prime_(\nk)H^prime_(\pn))/((epsilon_k - epsilon_n)(epsilon_k - epsilon_p)^2) \

  & - sum_(n eq.not k) (H^prime_(\pk)|H^prime_(\nk)|^2)/((epsilon_k - epsilon_n) (epsilon_k - epsilon_p)^2) - sum_(m eq.not k) (H^prime_(\kk) H^prime_(\mk) H^prime_(\pm))/((epsilon_k - epsilon_m)^2 (epsilon_k - epsilon_p)) + (|H^prime_(\kk)|^2 H^prime_(\pk))/((epsilon_k - epsilon_p)^3) }phi_p \ 
 $


= TDPT
$ H(t) = H_0 + H^prime (t) $
where $i hbar (partial)/(partial t) Phi_n (t) = H_0 Phi_n (t)$, $Phi_n (t) = phi_n (bold(r))text("exp")[-i (epsilon_n)/(hbar) t]$ and $sum_n |Phi_n (t) angle.r angle.l Phi_n (t)| = I$, $angle.l phi_n|phi_m angle.r = delta_(\nm)$, $H_0 phi_n = epsilon_n phi_n$.

Wavefunction $|f(t)angle.r$ satisfies $i hbar (partial)/(partial t)f(t) = (H_0 + H^prime (t))f(t)$. Projecting $|f(t)angle.r$ on basis ${Phi_n (t)| n=1,2, dots.c}$ gives
$
  & sum_n i hbar (partial)/(partial t) (a_n (t) Phi_n (t)) = sum_n (H_0 + H^prime)a_n (t) Phi_n (t) \

  & sum_n {i hbar (partial a_n (t))/(partial t)Phi_n (t) + i hbar a_n (t) (partial)/(partial t)Phi_n (t)} = sum_n (a_n (t)H_0 Phi_n (t) + a_n (t)H^prime Phi_n (t)) \
$

Note that $i hbar a_n (t) partial/(partial t)Phi_n (t) = a_n (t)H_0 Phi_n (0)$, thus
$
 & sum_n {i hbar (partial a_n(t))/(partial t)Phi_n(t)} = sum_n a_n(t)H^prime Phi_n(t) \

 & sum_n i hbar (partial a_n (t))/(partial t) angle.l phi_\mmid(|)phi_n angle.r e^(-i epsilon_n/hbar t) = sum_n a_n (t) angle.l phi_\mmid(|)H^prime mid(|)phi_n angle.r e^(-i epsilon_n/hbar t) \

& i hbar (partial a_m (t))/(partial t) e^(-i epsilon_m/hbar t) =  sum_n a_n(t) H^prime_(\mn) e^(-i (epsilon_n)/(hbar)t) \

& i hbar (partial)/(partial t)a_m (t) =  sum_n a_n (t) H^prime_(\mn)(t) e^(i omega_(\mn) t) \
$

where $omega_(\mn) = (epsilon_m - epsilon_n)/(hbar)$. 上式是含时薛定谔方程的另一种表现形式, 隐含使用相互作用表象. 


同定态微扰, 引入参数 $lambda$ 控制微扰打开的程度 $H^prime (t) = lambda H^prime (t)$, 那么同理在微扰下系数 $a_n (t)$ 可以展开成 $lambda$ 的幂级数
$ a_n (t) = a^((0))_n + lambda a_n^((1)) (t) + lambda^2 a_n^((2))(t) + dots.c $
带入 $i hbar (partial)/(partial t)a_m (t) =  sum_n a_n (t) H^prime_(\mn) (t) e^(i omega_(\mn) t)$ 有
$
  & i hbar (partial)/(partial t) [a_m^((0)) + lambda a_m^((1)) (t) + lambda^2 a_m^((2)) (t) +  dots.c ] \
  & = sum_n [a^((0))_n + lambda a_n^((1)) (t) + lambda^2 a_n^((2)) (t) + dots.c ]lambda H^prime_(\mn) (t) e^(i omega_(\mn)t)
$

两边对应幂次相等, 因此
$ & cases(
    i hbar (partial)/(partial t) a_m^((0)) = 0 \

  i hbar (partial)/(partial t) a_m^((1)) (t) = sum_n a_n^((0)) H^prime_(\mn) (t) e^(i omega_(\mn) t) \

  i hbar (partial)/(partial t) a_m^((2)) (t) = sum_n a_n^((1)) (t) H^prime_(\mn) (t) e^(i omega_(\mn) t) \

  dots.v  \

  i hbar (partial)/(partial t) a_m^((alpha)) (t) = sum_n a_n^((alpha-1)) (t) H^prime_(\mn) (t) e^(i omega_(\mn) t) \
)
$

说明无微扰时的展开系数 $a_m^((0))$ 不含时. 假设微扰在 $t=0$ 时开启, 此时系统处于 $|Phi_k (t)angle.r,r$, $H_0 Phi_k (t) = i hbar partial_t Phi_k (t)$, 所以有 $a_k^((0)) (0) = 1,a_n^((0)) (0) = delta_(\nk)$, 所以有
$
  i hbar (partial)/(partial t) a_m^((1)) (t) = a_k^((0)) H^prime_(\mk) (t) e^(i omega_(\mk) t)= H^prime_(\mk) (t) e^(i omega_(\mk) t)
$

解之可得一级近似
$
  a_m^((1)) (t) = 1/(i hbar) int_0^t H^prime_(\mk) (t^prime) e^(i omega_(\mk) t^prime) \dt^prime \
$

在时刻 $t$ 时系统处于 $|Phi_m (t)angle.r$ 的概率为 $|a_m (t)|^2$. 对微扰后的展开系数取一级近似 $a_m (t) approx a_m^((0)) + a_m^((1)) = a_m^((1))$. 即在微扰作用下系统由初态 $|Phi_k (t) angle.r$ 跃迁到 $|Phi_m (t)angle.r$ 的概率为
$ W_(k arrow.r m) = |a_m (t)|^2 approx |a_m^((1)) (t)|^2 $




//////////////////////////////////////////////////////////////////////////////////////////////////////
= HO
1D harmonic oscillator is second order ordinary homogeneous partial differential equation.
$ frac(d^2y,\dx^2) + a^2y = 0 $

\
The solution to this ODE is superposition of trigonometric functions 
$ & y(x) = A sin a x + B cos a x \
 & y(x) = A e^(i a x) + B e^(-i a x) $

\
In fact the two expressions above are equivalent and connected via Euler equation. And beware the coefficients in the trigonometric form may be complex.

== Series solution of harmonic oscillator
Assume $y(x)$ can be expanded as power series of $x$
#math.equation(block: true,numbering: none)[
  $ y(x) = sum_(n=0)^(infinity) a_n x^n  $
]<label>


== Series solution of quantum harmonic oscillator

== Ladder operator
$hat(a)$ is annihilation operator, $hat(a^dagger)$ is creation operator, $hat(N)$ is number operator (for brevity, '$hat$' is omitted below)
#math.equation(block: true)[
  #block(stroke: red,inset: 1em)[
    #math.equation(block: true,numbering: none)[
      $ & H = -h^2 / (2m) frac(d^2,\dx^2) + (omega^2 m x^2) / 2 \
        & [x,p] = i planck.reduce \
        & cases(a = sqrt(frac(m omega,2 planck.reduce)) (x + frac(i p , m omega)) \
                a^dagger = sqrt(frac(m omega,2 planck.reduce)) (x - frac(i p,m omega))) \ 
        & cases(x = sqrt(frac(planck.reduce,2 m omega)) (a + a^dagger) \
            p = frac(1,i) sqrt(frac(m omega planck.reduce,2)) (a - a^dagger)) \
        & [a,a^dagger] = 1 \
        & N = a^dagger a \
        & [N,a] = -a \
        & [N,a^dagger] = a^dagger \
        & cases(H = planck.reduce omega (N + 1 / 2) \
                N = frac(H,planck.reduce omega) - 1 / 2) \ $
    ]
  ]
]


\
$H$ and $N$ commute, and can be spontaneously diagonalized. Let ${|n angle.r|n=0,1,dots.c}$ be eigenvectors of $H$ and $N$, then 
#math.equation(block: true,numbering: none)[
  $ & N|n angle.r = n |n angle.r \
    & H|n angle.r = (planck.reduce omega N + frac(planck.reduce omega,2)) |n angle.r = (planck.reduce omega n + (planck.reduce omega) / 2)|n angle.r = E|n angle.r \
    & arrow.r.double E = planck.reduce omega (n + 1 / 2), n = 0,1,2,dots.c$
]

\
and we get the discrete energy level
#math.equation(block: true,supplement: "Eq.")[
  #block(stroke: red,inset: 1em)[
    #math.equation(numbering: none)[
      $ E = planck.reduce omega (n + 1 / 2), n = 0,1,2,dots.c $
    ]
  ]
]<energy>

\
$|n angle.r$ is eigenvector of $N$, so are $a|n angle.r$ and $a^dagger|n angle.r$
$ & N a^dagger|n angle.r = ([N,a^dagger] + a^dagger N)|n angle.r = a^dagger|n angle.r + n a^dagger|n angle.r = (n + 1) a^dagger|n angle.r \
  & N a|n angle.r = ([N,a] + a N)|n angle.r = -a|n angle.r + n a|n angle.r = (n - 1)a|n angle.r $

\
let $a^dagger|n angle.r = |p angle.r$, then $N|p angle.r = (n + 1)|p angle.r$, which means $|p angle.r = c|n+1 angle.r$, so does $a^dagger|n angle.r = c|n+1 angle.r$, where $c$ is normalization factor. 

Assume the eigenvectors ${|n angle.r|n=0,1,2,dots.c}$ are orthonormal, then 
#math.equation(block: true,numbering: none)[
  $ & abs(c)^2 angle.l n+1|n+1 angle.r = abs(c)^2 = angle.l n|a a^dagger|n angle.r = angle.l n|N + 1|n angle.r = n+1 \
  & arrow.r.double c = sqrt(n+1) $
]

\
and let $a|n angle.r = |p angle.r$, $N|p angle.r = (n-1)|p angle.r$, $|p angle.r = c|n-1 angle.r$, so as $a|n angle.r = c|n-1 angle.r$
#math.equation(block: true,numbering: none)[
  $ & abs(c)^2 angle.l n-1|n-1 angle.r = abs(c)^2 = angle.l n|a^dagger a|n angle.r = angle.l n|N|n angle.r = n \
  & arrow.r.double c = sqrt(n) $
]

\
therefore we have
#math.equation(block: true, supplement: "Eq.")[
  $ cases(a^dagger|n angle.r = sqrt(n+1) |n+1 angle.r \
        a|n angle.r = sqrt(n) |n-1 angle.r) $
]<ladder>
  
\
$a$ is #text(blue)[lowering operator(annihilation operator)], and $a^dagger$ is #text(blue)[raising operator(creation operator)]. Applying $a^dagger$ on $|0angle.r$ repeatedly gives 
$|n angle.r$
#math.equation(numbering: none,block: true)[
  $ & a^dagger |0angle.r = |1angle.r \
    & a^dagger a^dagger |0angle.r = a^dagger |1angle.r = sqrt(2)|2angle.r \
    & (a^dagger)^3 |0angle.r = sqrt(2) sqrt(3) |3angle.r \
    & (a^dagger)^4 |0angle.r = sqrt(2) sqrt(3) sqrt(4) |4angle.r \
    & dots.v \
    & (a^dagger)^n |0angle.r = sqrt(n!) |n angle.r $
]

\
therefore the nth eigenvector $|n angle.r$ can be evaluated from the ground state
#math.equation(block: true)[
  #block(stroke: red, inset: 1em)[
    #math.equation(block: true,numbering: none)[
      $ |n angle.r = 1 / sqrt(n!) (a^dagger)^n |0angle.r $
    ]
  ]
]

== Eigenfunctions of harmonic oscillator using ladder operator
${|n angle.r|n=0,1,2,dots.c}$ are the eigenvectors of number operator and Hamiltonian operator. The Hamiltonian here is the harmonic oscillator's Hamiltonian, the eigenvectors are the eigenfunction of 
harmonic oscillator, whose energy is given by @energy. However, now the eigenvectors are in number representation, it needs to be converted to the position representation. Under the number representation, 
both number operator and Hamiltonian operator are diagonal due to the orthogonality of ${|n angle.r}$
$ H_(\mn) = angle.l m|H|n angle.r  = planck.reduce omega angle.l m|N+1/2|n angle.r = planck.reduce omega N_(\mn) + planck.reduce omega 1/2 delta_(\mn) =
 planck.reduce omega delta_(\mn) + planck.reduce omega 1/2 delta_(\mn) $

\
since $[H,N]=0$, the number representation is also called energy representation. Under the position representation, the position operator is just position itself, while momentum operator is differential
$ cases(hat(x) = x \
        hat(p)_x = -i planck.reduce d / (\dx) \
        [x,p] = i planck.reduce) $
        
\
and the eigenfunction of position operator under position representation is Dirac $delta$ function
$ cases(x delta(x-x^prime) = x^prime delta(x-x^prime) \
        angle.l delta(x-x')|delta(x-x'') angle.r = integral_infinity^infinity delta(x-x') delta(x-x'') \dx = delta(x''-x')) $

\
the eigenfunction is normalized to Dirac $delta$ function. The eigenfunction of momentum operator under position representation is exponential,
$ cases(f_k (x) = sqrt(1 / (2 pi)) e^(i k x) \
        integral_infinity^infinity 1 / (2 pi) e^(i (k' - k'')x) \dx = delta(k' - k'')) $

\
which is normalized to Dirac $delta$ function. 

\
Position and momentum operator are now connected to ladder operator
#math.equation(block: true,numbering: none)[
  $ & cases(a = sqrt(frac(m omega,2 planck.reduce)) (x + frac(i p , m omega)) \
                a^dagger = sqrt(frac(m omega,2 planck.reduce)) (x - frac(i p,m omega))) \ 
        & cases(x = sqrt(frac(planck.reduce,2 m omega)) (a + a^dagger) \
            p = frac(1,i) sqrt(frac(m omega planck.reduce,2)) (a - a^dagger)) \ $
]

\
and their matrix elements under number representation are then
#math.equation(block: true, supplement: "Eq.")[
  $ & x_(\mn) = angle.l m |sqrt(frac(planck.reduce,2 m omega))(a + a^dagger)|n angle.r \
  & = sqrt(frac(#hbar,2 m omega)) (sqrt(n) angle.l m|n-1 angle.r + sqrt(n+1) angle.l m|n+1 angle.r) \
  & = sqrt(frac(#hbar,2 m omega)) (sqrt(n) delta_(m,n-1) + sqrt(n+1) delta_(m,n+1)) \
  & p_(\mn) = -i sqrt(frac(m omega #hbar,2)) (sqrt(n) delta_(m,n-1) - sqrt(n+1) delta_(m,n+1)) $
]<mat-elem>

\
They are not diagonal under number representation. 
Since $|0 angle.r$ is the ground state of number operator $N = a^dagger a$, it makes sense that applying lowering operator on the ground state gives zero
#math.equation(numbering: none,block: true)[$ a |0angle.r = 0 $]

\
Now let $|x' angle.r$ be one of the eigenvectors of position operator under position representation, and make use of the completeness condition $integral |x'' angle.r angle.l x''| \dx'' = 1$, 
we have 
#math.equation(block: true)[
  $ & angle.l x'|a|0angle.r = 0 \
    & arrow.r.double integral angle.l x'|a|x'' angle.r angle.l x''|0 angle.r \dx'' = 0 \
    & arrow.r.l.double integral angle.l x'|a|x'' angle.r f_0(x'') \dx'' = 0 \
    & = integral (integral delta(x-x') sqrt(frac(m omega,2 #hbar))(x + #hbar / (m omega) d / (\dx))delta(x-x'') \dx) f_0(x'') \dx'' \
    & = integral sqrt(frac(m omega, 2 #hbar)) [integral delta(x-x') x delta(x-x'') \dx + integral delta(x-x') frac(#hbar,m omega) delta'(x-x'') \dx] f_0(x'') \dx'' \ 
    & = sqrt(frac(m omega,2 #hbar)) integral [x'' delta(x''-x') + frac(#hbar,m omega) integral delta(x-x') delta'(x-x'')\dx] f_0(x'') \dx'' \
    & = sqrt(frac(m omega,2 #hbar)) integral x'' delta(x''-x') f_0(x'') \dx'' + sqrt(frac(m omega,2 #hbar)) frac(#hbar,m omega) integral (integral delta(x-x') delta'(x-x'')\dx) f_0(x'') \dx'' \
    & = sqrt(frac(m omega,2 #hbar)) x' f_0(x') - sqrt(frac(m omega,2 #hbar)) frac(#hbar,m omega) integral delta'(x-x')|_(x=x'') f_0(x'') \dx'' \
    & = sqrt(frac(m omega,2 #hbar)) x' f_0(x') - sqrt(frac(m omega,2 #hbar)) frac(#hbar,m omega) integral frac(d,\dx'') delta(x''-x') f_0(x'') \dx'' \
    & = sqrt(frac(m omega,2 #hbar)) x' f_0(x') + sqrt(frac(m omega,2 #hbar)) frac(#hbar,m omega) f'_0(x') = 0$
]

\
Finally, we get 
$ (x + frac(#hbar,m omega) frac(d,\dx)) f_0(x) = 0 $

\
which is the differential equation of the ground state of harmonic oscillator. Solving it gives $f_0(x) = c_1 e^(-frac(m omega,2 #hbar) x^2)$. Normalization gives 
#math.equation(block: true)[
  #block(stroke: red,inset: 1em)[
    #math.equation(numbering: none)[
      $ f_0(x) = (frac(m omega, pi #hbar))^(frac(1,4)) e^(-frac(m omega,2 #hbar)x^2) $
    ]
  ]
]

\
and for excited state eigenfunction
#math.equation(numbering: none,block: true)[
  $ & f_n (x') = angle.l x'|n angle.r \
  & = lr(angle.l x' lr(|frac(1,sqrt(n!)) (a^dagger)^n|)0 angle.r) \
  & = lr(angle.l x' lr(|frac(1,sqrt(n!)) (frac(m omega,2 #hbar))^frac(n,2) (x - frac(i p,m omega))^n|) 0 angle.r) \
  & = lr(angle.l x' mid(|) frac(1,sqrt(n!)) (frac(m omega,2 #hbar))^frac(n,2) (x - frac(i p,m omega))^n (frac(m omega,pi #hbar))^frac(1,4) e^(-frac(m omega,2 #hbar)x^2) angle.r) \
  & = frac(1,sqrt(n!)) (frac(m omega,2 #hbar))^frac(n,2) (frac(m omega,pi #hbar))^frac(1,4) lr(angle.l x' mid(|) (x - frac(i p,m omega))^n mid(|) e^(-frac(m omega,2 #hbar)x^2) angle.r) \
  & = frac(1,sqrt(n!)) (frac(m omega,2 #hbar))^frac(n,2) (frac(m omega,pi #hbar))^frac(1,4) integral delta(x-x') (x - frac(#hbar,m omega) frac(d,\dx))^n e^(-frac(m omega,2 #hbar)x^2) \dx \ 
  & = frac(1,sqrt(n!)) (frac(m omega,2 #hbar))^frac(n,2) (frac(m omega,pi #hbar))^frac(1,4) (x' - frac(#hbar,m omega) frac(d,\dx))^n e^(-frac(m omega,2 #hbar) x'^2) $
]

\
#math.equation(block: true)[
  #block(stroke: red,inset: 1em)[
    #math.equation(numbering: none)[
      $ f_n (x) = frac(1,sqrt(n!)) (frac(m omega,2 #hbar))^frac(n,2) (frac(m omega,pi #hbar))^frac(1,4) lr((x - frac(#hbar,m omega) frac(d,\dx)))^n e^(-frac(m omega,2 #hbar)x^2) $
    ]
  ]
]


== Expectation values
Matrix elements of position and momentum operator under number operator are given by @mat-elem and @ladder
#math.equation(block: true,numbering: none)[
  $ & cases(
    x_(\mn) = sqrt(frac(#hbar,2 m omega)) (sqrt(n) delta_(m,n-1) + sqrt(n+1) delta_(m,n+1)) \
    p_(\mn) = -i sqrt(frac(m omega #hbar,2)) (sqrt(n) delta_(m,n-1) - sqrt(n+1) delta_(m,n+1)) 
  ) \
  & cases(a^dagger|n angle.r = sqrt(n+1) |n+1 angle.r \
        a|n angle.r = sqrt(n) |n-1 angle.r)$
]

\
clearly the expectation values of position and momentum operator at state $|n angle.r$ are zero
#math.equation(numbering: none,block: true)[
  $ & angle.l x angle.r = lr(angle.l n lr(|x|) n angle.r) = 0 \
  & angle.l p angle.r = lr(angle.l n lr(|p|) n angle.r) = 0 $
]

\
but for $x^2$ and $p^2$ we have
$ & lr(angle.l n lr(|x^2|) n angle.r) = frac(#hbar,2 m omega) [ lr(angle.l n lr(|\aa|) n angle.r) + lr(angle.l n lr(|a a^dagger|) n angle.r) + lr(angle.l n lr(|a^dagger a|) n angle.r) + lr(angle.l n lr(|a^dagger a^dagger|) n angle.r) ] \
& = frac(#hbar,2 m omega) (n+1+n) = (n + 1 / 2) frac(#hbar,m omega) \
\
& lr(angle.l n lr(|p^2|) n angle.r) = -frac(m omega #hbar,2) [ lr(angle.l n lr(|\aa|) n angle.r) - lr(angle.l n lr(|a a^dagger|) n angle.r) - lr(angle.l n lr(|a^dagger a|) n angle.r) + lr(angle.l n lr(|a^dagger a^dagger|) n angle.r) ] \
& = -frac(m omega #hbar,2) (-(n+1)-1) = (n + 1 / 2) m omega #hbar \
$

\
For ground state we have $angle.l x^2 angle.r_0 = frac(#hbar,2 m omega)$, and $angle.l p^2 angle.r_0 = frac(m omega #hbar,2)$. Expectation values of ground state kinetic energy and potential energy are then
$ & angle.l E_k angle.r_0 = angle.l frac(p^2,2 m) angle.r = frac(#hbar omega,4) \
  & angle.l E_p angle.r_0 = angle.l frac(m omega^2 x^2,2) angle.r = frac(#hbar omega,4) \
  & angle.l H angle.r_0 = (n + 1 /2) #hbar omega = frac(#hbar omega,2) $

\
and uncertainties
$ & lr(angle.l (Delta x)^2 angle.r) = lr(angle.l lr((x - overline(x))^2) angle.r) = angle.l x^2 angle.r - angle.l x angle.r^2 = (n + 1 / 2) #hbar / (m omega) \
  & angle.l (Delta p)^2 angle.r = (n + 1 / 2) m omega #hbar \
  & angle.l (Delta x)^2 angle.r angle.l (Delta p)^2 angle.r = (n + 1 / 2)^2 #hbar^2 $


\
== Matrix representation of creation and annihilation operator
Assume ${|n angle.r|n=0,1,2,dots.c}$ to be the orthonormal eigenvectors under number representation. The matrix elements of $a^dagger$ are(see also #text(rgb("#85144b"))[#link("https://physics.stackexchange.com/questions/98049/example-of-application-of-creation-annihilation-operators-in-matrix-form")[this link]]) 
#math.equation()[
  #block(stroke: red,inset:1em)[
    #math.equation(supplement: "Eq.",numbering:none)[
  $ & a^dagger_(\mn) = lr(angle.l m lr(|a^dagger|) n angle.r) = sqrt(n+1) lr(angle.l m mid(|) n+1 angle.r) = sqrt(n+1) delta_(\m,n+1)\ 
    & a^dagger = mat(delim:"[", 0,0,0,0,0,0,dots.c; 
                    1,0,0,0,0,0,dots.c;
                    0,sqrt(2),0,0,0,0,dots.c;
                    0,0,sqrt(3),0,0,0,dots.c;
                    0,0,0,sqrt(4),0,0,dots.c;
                    0,0,0,0,sqrt(5),0,dots.c; 
                    dots.v,dots.v,dots.v,dots.v,dots.v,dots.v,dots.down) $
]
  ]
]

\
and ${|n angle.r}$ under number representation are
#math.equation(numbering: none)[
  $ |0angle.r = mat(delim:"[",1;0;0;0;0;0;dots.v), |1angle.r = mat(delim:"[",0;1;0;0;0;0;dots.v), |2angle.r = mat(delim:"[",0;0;1;0;0;0;dots.v), |3angle.r = mat(delim:"[",0;0;0;1;0;0;dots.v), dots.c, |n angle.r = mat(delim:"[",0;dots.v;0;1;0;dots.v;0;) $
]
the $n+1$'th element is one for eigenvector $|n angle.r$. Obviously, they form orthonormal basis vectors in N-dim space(N $arrow.r oo$).

\
Similarly, for annihilation operator $a$ 
#math.equation()[
  #block(stroke: red,inset:1em)[
    #math.equation(numbering: none)[
  $ & a|n angle.r = sqrt(n) |n-1 angle.r \
    & a = mat(delim:"[", 0,1,0,0,0,0,dots.c;
                          0,0,sqrt(2),0,0,0,dots.c;
                        0,0,0,sqrt(3),0,0,dots.c;
                      0,0,0,0,sqrt(4),0,dots.c;
                    0,0,0,0,0,sqrt(5),dots.c;
                  0,0,0,0,0,0,dots.c;
                dots.v,dots.v,dots.v,dots.v,dots.v,dots.v,dots.down)$
]
  ]
]

\
Note that $a$ and $a^dagger$ *#text(red)[are NOT Hermitian]*, which can be seen from their matrices. But $a^dagger a = N$ is Hermitian, and under number representation it's diagonal with diagonal elements being the number of particle
#math.equation()[
  #block(stroke: red,inset:1em)[
    #math.equation(numbering: none)[
  $ N = a^dagger a = mat(delim:"[", 0,,,,,,; 
                        ,1,,,,,; 
                        ,,2,,,,; 
                        ,,,3,,,; 
                        ,,,,4,,;
                        ,,,,,5; 
                        ,,,,,,dots.down,) $
]
  ]
]
The number operator is Hermitian(real symmetric). 

= Numeric solution to HO
== Numeric solution to classic harmonic oscillator

== Numeric solution to quantum harmonic oscillator
Below the Mathematica code 

\
```Mathematica
\[HBar] = 1.;
m = 1.;
\[Omega] = 0.2;
npoints = 100;

x = DiagonalMatrix[Subdivide[-20., 20., npoints - 1]];
dx = x[[2, 2]] - x[[1, 1]]
potential = m*\[Omega]^2/2.*(x . x);

kinetic = -\[HBar]^2/2./
    m*(DiagonalMatrix[Table[-2, npoints]] + 
      DiagonalMatrix[Table[1, npoints - 1], 1] + 
      DiagonalMatrix[Table[1, npoints - 1], -1])/dx^2;

hamiltonian = kinetic + potential;


{val, vec} = Eigensystem[hamiltonian];
ListPlot[val] (* energy difference \hbar\omega *)

norm = Sqrt@Total[vec[[-5]]^2*dx]
ListPlot[vec[[-5]]/norm, PlotRange -> All, Joined -> True] (* numeric eigenvectors *)

(* analytic solution *)
x = Subdivide[-20., 20., 99];
\[Alpha] = Sqrt[m*\[Omega]/\[HBar]];
f[x_, n_] := 
  Sqrt[1/2^n/n!]*(m*\[Omega]/Pi/\[HBar])^(0.25)*
   Exp[-\[Alpha]^2*x^2/2.]*HermiteH[n, \[Alpha]*x];
ListPlot[{f[x, 4]}, PlotRange -> All, Joined -> True]
```

\

We want to numerically solve the quantum harmonic oscillator Schrodinger equation under the position representation. "Position representation" means the eigenfunction is a function of position $x$, and under this representation the position operator $hat(x)$ is diagonal. Analytically, in the position representation $hat(x) = x$; numerically, the position has to be discretized and restricted in a bounded range $x in [x_(min), x_(max)]$, $x = {x_0, x_1, x_2, dots.c, x_i, dots.c, x_N}$, $Delta x = x_2 - x_1$. Then the numeric analogue to analytic $hat(x) = x$ gives the position operator matrix, 
#math.equation()[
  #block(stroke: red,inset:1em)[
    #math.equation(numbering: none)[
      $ hat(x) = mat(delim:"[", x_0,,,,,; 
                    ,x_1,,,,; 
                    ,,x_2,,,;
                    ,,,dots.down,,; 
                    ,,,,x_N;) $
    ]
  ]
]<potential>
which is $N+1 times N+1$ dim matrix. 

\
In the position representation, the eigenfunction of position operator is Dirac $delta$ function(#underline[normalize to $delta$ function])
#math.equation(numbering: none)[
  $ hat(x) delta(x-x') = x delta(x-x') = x' delta(x-x') $ 
]
and the numeric discrete analogue is Kronecker $delta$(#underline[normalize to this $delta$]), thus
#math.equation(numbering: none)[
  $ mat(|x_0 angle.r, |x_1 angle.r, |x_2 angle.r, dots.c, |x_N angle.r) = mat(1,0,0,dots.c,0; 
0,1,0,dots.c,0; 
                                                                              0,0,1,dots.c,0; 
                                                                              dots.v, dots.v, dots.v, dots.down, dots.v; 
                                                                              0,0,0,dots.c, 1;) $
]
Each eigenvector is $N+1$ dim vector. 

\
Since $hat(x)$ and $hat(p)$ do not commute $[hat(x),hat(p)] = i #hbar$. The eigenvectors for $hat(x)$ are not the eigenvectors for $hat(p)$. Analytically, in the position representation the eigenfunction for momentum operator $hat(p)$ is 
$ & f_k (x) = sqrt(1 / (2 pi)) e^(i k x) \ 
  & lr(angle.l f_k' mid(|) f_k angle.r) = 1 / (2 pi) integral_(-oo)^oo e^(i x (k-k')) \dx = delta(k - k') \
  & -i #hbar d / (\dx) f_k (x) = #hbar k f_k (x) = p f_k (x) $

Similarly, the eigenfunction of momentum operator in the position representation is also normalized to Dirac $delta$ function. 

\
To numerically express the momentum operator($hat(p)$) and related kinetic energy operator($(hat(p) dot.c hat(p)) / (2m)$), here finite difference method(FDM) is adopted. For first order derivative, finite difference gives $f'_i = (f_(i+1) - f_(i-1)) / (2Delta)$, while for second order derivative it gives 
$ f''_i = (f_(i+1) + f_(i-1) - 2f_i) / (Delta x^2) $
in matrix form is a tri-diagonal matrix
#math.equation(numbering: none)[
  $ mat(;
      0,dots.c,-2,1,0,dots.c,0;
      0,dots.c,1,-2,1,dots.c,0; 
      0,dots.c,0,1,-2,1,dots.c;
      ;) $
]
We expect the eigenfunction for harmonic oscillator to vanish at infinity(natural boundary condition) $lim_(x arrow.r oo) f(x) = 0$. The position has been discretized $x = {x_0, x_1, x_2, dots.c, x_N}$, therefore we have $f(x_(-1)) = f(x_(N+1)) = 0$, which leads to 
#math.equation()[
  #block(stroke: red, inset: 1em)[
    #math.equation(numbering: none)[
      $ hat(p)^2 = -#hbar^2 d^2 / (\dx^2) =  - #hbar^2 mat(-2,1,,,,,; 
      1,-2,1,,,,; 
      ,1,-2,1,,; 
      ,,1,-2,1,;
      ,,,, dots.down; 
      ,,,,1,-2,1;
      ,,,,,1,-2) $
    ]
  ]
]<kinetic>

Clearly, kinetic operator matrix is real symmetric, and therefore also Hermitian. And the Hamiltonian for harmonic oscillator is then
$ H = - #hbar^2 / (2m) d^2 / (\dx^2) + (m omega^2) / 2 x_i dot.c x_i $
Plugging @potential and @kinetic into it gives the matrix form of Hamiltonian, whose eigenvectors are the numeric solutions to the quantum harmonic oscillator. 

\
贴张图

= 相干态
== 平移算符
to be continued...
== 谐振子相干态
谐振子的特征向量构成正交完备基 ${|n angle.r|n=0,1,2,dots.c}$, 相干态可由其展开
$ |alpha angle.r = e^(-(|alpha|^2)/(2)) sum_(n=0)^oo (alpha^n)/(sqrt(n!)) |n angle.r $
证明待补...

\

湮灭算符的特征函数称为相干态. 由于湮灭算符不是厄米的, 其特征值不一定是实的. 关于相干态在坐标表象下的函数形式, 既可以从特征方程出发通过求解微分方程得到

\
$ & hat(a) |alpha angle.r = alpha |alpha angle.r \
  & hat(a) = sqrt((m omega) / (2 #hbar)) lr((x + (i p) / (m omega))) \
  & sqrt((m omega)/(2 #hbar)) (x + (#hbar)/(m omega)(d)/(\dx)) f(x) = alpha f(x) $

解之并归一化后可得
$ & alpha = Re[alpha] + i Im[alpha]  = |alpha| e^(i theta) = |alpha| (cos theta + i sin theta) \ 
  & |alpha angle.r = ((m omega) / (pi #hbar))^(1/4) exp[-(m omega)/(2 #hbar)(x-sqrt((2 #hbar)/(m omega)) |alpha|cos theta)^2] dot.c [cos(sqrt((2 m omega)/(#hbar))|alpha|sin theta x) + i sin(sqrt((2 m omega)/(#hbar)|alpha|sin theta x))] $


相干态坐标表象下的函数形式还可以通过另一种方式得到
$ & |alpha angle.r = e^(-(|alpha|^2)/(2)) sum_(n=0)^oo (alpha^n)/(sqrt(n!)) |n angle.r \
  & lr(angle.l x mid(|)alpha angle.r) = e^(-(|alpha|^2)/(2)) sum_(n=0)^oo (alpha^n)/(sqrt(n!)) lr(angle.l x mid(|)n angle.r) $

其中 $f_n (x) = lr(angle.l x mid(|)n angle.r)$ 是谐振子第n个特征函数. 
$ & f_n (x) = sqrt((1)/(2^n n!)) ((m omega)/(pi #hbar))^(1/4) e^(-(m omega)/(2 #hbar) x^2) H_n (sqrt((m omega)/(#hbar))x) \
  & e^(2\xt - t^2) = sum_(n=0)^oo (t^n)/(n!) H_n (x) \
  $

\

海森堡运动方程
$ (d)/(\dt) hat(A)(t) = (1)/(i #hbar) [hat(A)(t),hat(H)] $
  
// \
// \
// \
// \
// $ mat(delim:"[",1,2;3,4;)^2 $

// //#theorem[this is a theorem]

// a new paragraph. A new list
// + ff see @ref1
//   - ff
// + gg
// + hh let's cite @shamasundar2002higher


// #figure(
//   image("test.svg",width: 80%),
//   caption: [this is a figure]
// )<ref1>

// #"1 2 3".split()

// #str(123,base:2).split()

// #str("123").len()

// some math $a^2 + b^2 = c^2$, $cos alpha + sin phi = tan psi$, $1/2 * x^2$

// // reference
// #pagebreak()
// #bibliography("test.bib",style: "american-physics-society")
