#import "@preview/unequivocal-ams:0.1.2": ams-article, theorem, proof
// #set document(author: "lyh",title: "Harmonics Oscillator") // pdf metadata
#set heading(outlined: true,bookmarked: true,)
#set par(leading: 0.55em, spacing: 1.55em, first-line-indent: 2.0em, justify: true)
#show heading: set block(above: 1.4em, below: 1em)
#set page(
  margin: 1.0in,
  numbering: "1",
  paper: "a4",
  number-align: center,
  footer: context[#set align(center)
                  #counter(page).display("1 of 1",both:true)],
)
#set math.equation(numbering: "(1)")
#set math.equation(block: true)
#set math.mat(delim:"[")
#show link: underline
#set text(font: ("Times New Roman", "SimSun"))
#let hbar = math.planck.reduce
#let otimes = math.times.circle
#let indent=h(2em)

// https://qiita.com/gomazarashi/items/a7e3d17b13598c1ba143 math fonts
// content
#align(center,text(21pt)[A hopefully gentle (but not
 rigorous) derivation of Lindblad equation])
\
#align(center,text(12pt)[lyh])
\
\
// _bold letters stand for vectors unless otherwise stated_
\
\

/////////////
// main texts
/////////////
= Math preparations
An important identity:
$ lr((A times.circle B)) lr((C times.circle D)) = A C times.circle B D $
this will be used multiple times in derivation. 


eom of density matrix reads
$ & dot(rho(t)) = -frac(i,#hbar)[H,rho] \
  & arrow.l.r.double dot(rho(t)) = cal(L) rho, $
where $cal(L)=-frac(i,hbar)[H,(dot.c)]$. 

In the so-called Liouville space, density matrix $rho$ is 
vectorized: $rho arrow.r |rho angle.r angle.r$, and $cal(L)$ is 
a matrix that acts on a matrix $cal(L) arrow.r hat(hat(cal(L)))$, 
which is called superoperator. 

Liouville space and all related algebra do not introduce new physics, it's merely a 
mathematical trick, to make use of linear algebra techniques. 

Let ${|i angle.r| i=1.2.dots.c}$ be orthonormal basis, then operator
 $hat(rho)$ can be decomposed as 
 $ hat(rho) = sum_(i,j) lr(|i angle.r angle.l i|)hat(rho)lr(|j angle.r angle.l j|) = sum_(i,j) rho_(\ij) lr(|i angle.r angle.l j|), $

which is equivalent to a vectorized form
$ sum_(\ij) rho_(\ij)lr(|i angle.r angle.l j|) <=> sum_(\ij) rho_(\ij) |i angle.r times.circle |j angle.r $

and matrix acting on it can be rewritten as 
$ & H |rho angle.r angle.r = sum rho_(\ij) H|i angle.r times.circle |j angle.r equiv sum rho_(\ij) H times.circle I |i angle.r times.circle |j angle.r \
  & => H |rho angle.r angle.r = H times.circle I |rho angle.r angle.r \
  \
  \
  \
  & rho H = sum rho_(\ij)|i angle.r angle.l j|H = sum_(\ij)rho_(\ij)|i angle.r lr((H^dagger |j angle.r)^dagger) equiv sum rho_(\ij)|i angle.r times.circle H^dagger|j angle.r  equiv sum rho_(\ij) I times.circle H^dagger|i angle.r times.circle |j angle.r  \
  & => I times.circle H^dagger |rho angle.r angle.r \
  \
  \
  \
  & H rho - rho H = [H,rho]_- equiv lr((H times.circle I - I times.circle H^dagger))|rho angle.r angle.r  $
Therefore the Liouville superoperator is
$ cal(L) = - i / hbar [H,(dot.c)] = -i/hbar lr((H times.circle I - I times.circle H)) $
and $cal(L)|rho angle.r angle.r = -i/hbar lr((H times.circle I - I times.circle H)) |rho angle.r angle.r$. 

This is how Liouville superoperator is constructed in practice. 

The construction of anti-commutator superoperator is similar: ${s^dagger s,(dot.c)}_+ = s^dagger s times.circle I + I times.circle s^dagger s$. Another type is $s rho s^dagger$, which can be constructed by
$ & s rho s^dagger = sum rho_(\ij) s|i angle.r angle.l j|s^dagger \
  & = sum rho_(\ij) s|i angle.r lr((s|j angle.r))^dagger \
  & equiv sum rho_(\ij) s|i angle.r times.circle s|j angle.r \
  & = s times.circle s |rho angle.r angle.r \ $

In summary
#set math.equation(numbering: none)
#table(  
  columns: 3,
  row-gutter: 1em,
  column-gutter: 1em,
  stroke: none,
  table.header[$[H,(dot.c)]_-$][${s^dagger s,(dot.c)}_+$][$s(dot.c)s^dagger$],
  [$H times.circle I - I times.circle H$], [$s^dagger s times.circle I + I times.circle s^dagger s$], [$s times.circle s$]
)


#set math.equation(numbering: "(1)")
= Lindbladian
Hamiltonian is factored into two parts $H = H^0 + H^prime (t)$, where $H^0$ is easily solvable and all difficult parts are absorbed into $H^prime (t)$. The eom of density matrix under interaction picture is given by
$ & dot(rho^(text(\t\ot))_I) = partial/(partial t) lr((e^(i (H^0 t)/hbar) rho_S^text(t\ot) e^(-i (H^0 t)/hbar))) \
  & = (i H^0) / hbar e^((i H^0 t)/hbar) rho e^(-i (H^0 t)/hbar) + e^((i H^0 t)/hbar) dot(rho) e^(-i (H^0 t)/hbar) - e^((i H^0 t)/hbar) rho (i H^0)/(hbar) e^(-i (H^0 t)/hbar) \
  & dot(rho_I^text(t\ot) (t)) = i/hbar e^((i H^0 t)/hbar) H^0 rho e^(-i (H^0 t)/hbar) - i/hbar e^((i H^0 t)/hbar) [H^0 + H^prime,rho] e^(-i (H^0 t)/hbar) - e^((i H^0 t)/hbar) rho e^(-i (H^0 t)/hbar) (i H^0)/hbar \
  & = i/hbar H^0_I rho_I^(t\ot) - i/hbar e^((i H^0 t)/hbar) H^0 rho e^(-i (H^0 t)/hbar) - i/hbar e^((i H^0 t)/hbar) H^prime rho e^(-i (H^0 t)/hbar) + i/hbar e^((i H^0 t)/hbar) rho H^0 e^(-i (H^0 t)/hbar) + e^((i H^0 t)/hbar) rho H^prime e^(-i (H^0 t)/hbar) - rho_I^(t\ot) H^0_I i/hbar \
  & = i/hbar (H_I^0 rho_I^(t\ot) - H_I^0 rho_I^(t\ot) - H^prime_I rho_I^(t\ot) + rho_I^(t\ot)H_I^0 + rho_I^(t\ot) H_I^prime - rho_I^(t\ot) H^0_I) \
  & = -i/hbar [H_I^prime (t), rho_I^(t\ot)(t)] \ $
Under interation picture, the eom of density matrix has the same form as that under Schrodinger picture. 

Since density matrix is time-dependent, the formal solution to its eom involves time-ordering
$ & rho_I^(t\ot)(t) - rho_I^(t\ot)(0) = -i/hbar integral_0^t [H_I^prime (t^prime), rho_I^(t\ot)(t^prime)] \dt^prime \
  & rho_I^(t\ot)(t) = rho_I^(t\ot)(0) - i/hbar integral_0^t [H_I^prime (t^prime), rho_I^(t\ot)(0) - i/hbar integral_0^(t^prime)[H_I^prime (t^(prime prime)),rho_I^(t\ot)(t^(prime prime))] \dt^(prime prime)] \dt^(prime) \
  & dots.v \ $
The expansion stops at the second order:
$ rho_I^(t\ot)(t) = rho_I^(t\ot)(0) - i/hbar integral_0^t [H_i^prime (t^prime), rho_I^(t\ot)(0)]\dt^prime - 1/hbar^2 integral_0^t \dt^prime [H_I^prime (t^prime), integral_0^(t^prime)[H_I^prime (t^(prime prime)),rho_I^(t\ot)(t^(prime prime))]\dt^(prime prime)] $
Taking time derivative to get rid of one integral
$ dot(rho_I^(t\ot))(t) = -i/hbar [H_I^prime (t),rho_I^(t\ot)(0)] - 1/hbar^2 integral_0^t [H_I^prime (t),[H_I^prime (t^prime),rho_I^(t\ot) (t^prime)]] \dt^prime $

So far, no approximations are made. We assume the system and bath are weakly coupled, thus the total density matrix is approximately a direct product of system density matrix and bath density matrix: $rho_I^("tot")(t) approx rho_I^S (t) times.circle rho_I^B (t)$; and we assume the bath is too large to be affected: $rho_I^B (t) = rho_I^B (0) = R_0$, therefore $rho_I^("tot")(t) = rho_I^S (t) times.circle R_0$. 

The total Hamiltonian reads $H = H_S + H_B + H^prime$, where $H_S$ and $H_B$ act solely on system and bath, while the coupling part $H^prime$ acts on both. This writing is misleading as the Hamiltonian is actually direct product:
$ & H_S = H_S times.circle I \
  & H_S rho^"tot" = H_S times.circle I rho_S times.circle rho_B = H_S rho_S times.circle rho_B \
  \
  \
  & H_B = I times.circle H_B \
  & H_B rho^"tot" = I times.circle H_B rho_S times.circle rho_B = rho_S times.circle H_B rho_B \
  \
  \
  & H^prime = s times.circle Gamma \
  & H^prime rho^"tot" = s rho_S times.circle Gamma rho_B, \ $

where $s$ and $Gamma$ are system and bath operator, respectively. Under interaction picture, $s$ and $Gamma$ are transformed accordingly
$ & s -> tilde(s) \
  & Gamma --> tilde(Gamma) \
  & H_I^prime = tilde(s) times.circle tilde(Gamma) \ $

To get the information of system, we trace out the bath part:
$ & dot(rho_(S,I) (t)) = "Tr"_"B" dot(rho_I^("tot")(t)) \
  & = -i/hbar "Tr"_"B" {H_I^prime (t) rho_(I)^S (0) times.circle R_0 - rho_I^S (0) times.circle R_0 H_I^prime (t)} - 1/hbar^2 "Tr"_"B" integral_0^t [H_I^prime (t),[H_I^prime (t^prime),rho_I^(t\ot) (t^prime)]] \dt^prime \
  & = -i/hbar (tilde(s)(t) rho_I^S (0) "Tr"_"B" {tilde(Gamma)(t) R_0 } - rho_I^S (0) tilde(s)(t) "Tr"_"B" {R_0 tilde(Gamma)(t)}) - dots.c \ $

Here $"Tr"_"B" {R_0 tilde(Gamma)(t)}$ is the expectation value of bath, we assume it to be zero. Then
$ & dot(rho_(S,I) (t)) = "Tr"_"B" dot(rho_I^("tot")(t)) \
  & = - 1/hbar^2 "Tr"_"B" integral_0^t [H_I^prime (t),[H_I^prime (t^prime),rho_I^(t\ot) (t^prime)]] \dt^prime, \ $
which states the state of system at time $t$ is related to the state of system at time $t^prime$ that is prior to $t$. If the thermal-equilibrium of bath is much faster than the evolution of system(since bath is huge), the state of system is not determined by its history, so we get
$ dot(rho_(S,I) (t)) = - 1/hbar^2 "Tr"_"B" integral_0^t [H_I^prime (t),[H_I^prime (t^prime),rho_I^(t\ot) (t)]] \dt^prime $<tag>

// In Schrodinger picture, operators are time-independent, but the eom of density matrix under Schrodinger picture is given by $dot(rho) = - i/hbar [H,rho]$, which seems to be contradict. The traditional operators are functions of momentum and/or position, which is not the case for density matrix. Density matrix is not traditional operator. 

Interaction Haimltonian is direct product of system operator and bath operator $H^prime = hbar sum_i s_i otimes Gamma_i$, and under interation picture it has the same form $H_I^prime = hbar sum_i tilde(s_i) otimes tilde(Gamma_i)$. So #ref(<tag>) can be further simplified:
$ & dot(rho_(S,I))(t) = -1/hbar^2 integral_0^t "Tr"_"B" [hbar sum_i tilde(s_i) (t) otimes tilde(Gamma_i) (t), [hbar sum_j tilde(s_j) (t^prime) otimes tilde(Gamma_j) (t^prime), rho_i^"tot" (t)]] \dt^prime \
  & = -sum_(\ij) integral_0^t "Tr"_"B" [tilde(s_i)(t) otimes tilde(Gamma_i)(t), tilde(s_j)(t^prime) otimes tilde(Gamma_j)(t^prime) rho_I^"tot" (t) - rho_I^"tot" (t) tilde(s_j) (t^prime) otimes tilde(Gamma_j) (t^prime)] \dt^prime \ $
Toal density matrix is approximately direct product of system density matrix and bath density matrix: $rho_I^"tot" (t) approx rho_(S,I) (t) otimes rho_(B,I) (t) = rho_(S,I) (t) otimes R_0$. So
$ dot(rho_(S,I))(t) = & -sum_(\ij) integral_0^t  "Tr"_"B" { tilde(s_i) (t) otimes tilde(Gamma_i) (t) tilde(s_j) (t^prime) otimes tilde(Gamma_j) (t^prime) rho_(S,I) (t) otimes R_0 \
  & - tilde(s_i) (t) otimes tilde(Gamma_i) (t) rho_(S,I) (t) otimes R_0 tilde(s_j) (t^prime) otimes tilde(Gamma_j) (t^prime) \
  & -tilde(s_j) (t^prime) otimes tilde(Gamma_j) (t^prime) rho_(S,I) (t) otimes R_0 tilde(s_i) (t) otimes tilde(Gamma_i) (t) \
  & + rho_(S,I) (t) otimes R_0 tilde(s_j) (t^prime) otimes tilde(Gamma_j) (t^prime) tilde(s_i) (t) otimes tilde(Gamma_i) (t)} \dt^prime \ $

$ dot(rho_(S,I)) (t) = & -sum_(\ij) integral_0^t "Tr"_"B" { \
  &  (tilde(s_i) (t) tilde(s_j) (t^prime) otimes tilde(Gamma_i) (t) tilde(Gamma_j) (t^prime)) rho_(S,I) (t) otimes R_0 \
  & - (tilde(s_i) (t) rho_(S,I) (t) otimes tilde(Gamma_i) (t) R_0) tilde(s_j) (t^prime) otimes tilde(Gamma_j) (t^prime) \
  & - (tilde(s_j) (t^prime) rho_(S,I) (t) otimes tilde(Gamma_j) (t^prime) R_0) tilde(s_i) (t) otimes tilde(Gamma_i) (t) \
  & + (rho_(S,I) (t) tilde(s_j) (t^prime) otimes R_0 tilde(Gamma_j) (t^prime)) tilde(s_i) (t) otimes tilde(Gamma_i) (t)} \dt^prime \
  & = -sum_(\ij) integral_0^t "Tr"_"B" { \
  &  tilde(s_i) (t) tilde(s_j) (t^prime) rho_(S,I) (t) otimes tilde(Gamma_i) (t) tilde(Gamma_j) (t^prime) R_0 \
  & - tilde(s_i) (t) rho_(S,I) (t) tilde(s_j) (t^prime) otimes tilde(Gamma_i) (t) R_0 tilde(Gamma_j) (t^prime) \
  & - tilde(s_j) (t^prime) rho_(S,I) (t) tilde(s_i) (t) otimes tilde(Gamma_j) (t^prime) R_0 tilde(Gamma_i) (t) \
  & + rho_(S,I) (t) tilde(s_j) (t^prime) tilde(s_i) (t) otimes R_0 tilde(Gamma_j) (t^prime) tilde(Gamma_i) (t) } \dt^prime \
  & = -sum_(\ij) integral_0^t angle.l tilde(Gamma_i) (t) tilde(Gamma_j) (t^prime) angle.r (tilde(s_i) (t) tilde(s_j) (t^prime) rho_(S,I) (t) - tilde(s_j) (t^prime) rho_(S,I) (t) tilde(s_i) (t)) \
  & - angle.l tilde(Gamma_j) (t^prime) tilde(Gamma_i) (t) angle.r (tilde(s_i) (t) rho_(S,I) (t) tilde(s_j) (t^prime) - rho_(S,I) (t) tilde(s_j) (t^prime) tilde(s_i) (t)) \dt^prime, \ $
where $angle.l tilde(Gamma_i) (t) tilde(Gamma_j) (t^prime) angle.r$ is the correlation function of bath operator. Correlation between different bath operators are zero; correlation within the same bath operator behaviours like delta function:
$ angle.l tilde(Gamma_i) (t) tilde(Gamma_j) (t^prime) angle.r prop delta_(\ij) delta(t-t^prime) $<corr>

System operators are hermitian. Replacing $s_i$ in the first term with $s_i^dagger$, and $s_j$ with $s_j^dagger$ in the second term renders
$ dot(rho_(S,I)) (t) = & -sum_(\ij) integral_0^t angle.l tilde(Gamma_i) (t) tilde(Gamma_j) (t^prime) angle.r (tilde(s_i^dagger) (t) tilde(s_j) (t^prime) rho_(S,I) (t) - tilde(s_j) (t^prime) rho_(S,I) (t) tilde(s_i^dagger) (t)) \
  & - angle.l tilde(Gamma_j) (t^prime) tilde(Gamma_i) (t) angle.r (tilde(s_i) (t) rho_(S,I) (t) tilde(s_j^dagger) (t^prime) - rho_(S,I) (t) tilde(s_j^dagger) (t^prime) tilde(s_i) (t)) \ $<derivation>

If the system operator satisfies
$ & [H_"sys", s_i] = -hbar omega_i s_i \
  & [H_"sys", s_i^dagger] = hbar omega_i s_i^dagger, \ $ 
for example the system Hamiltonian is qho without zpe and system operator is creation/annihilation operator, then under interaction picture system operator reads
$ & tilde(s_i) = e^(i H_"sys"/hbar t) s_i e^(-i H_"sys"/hbar t) \
  & tilde(s_i^dagger) = e^(i H_"sys"/hbar t) s_i^dagger e^(-i H_"sys"/hbar t) \ $
Here comes operator exponential $e^A B e^(-A)$, and Baker-Hausdorff formula states
$ e^A B e^(-A) = B + [A,B] + 1/2! [A, [A,B]] + 1/3! [A, [A, [A,B]]] + dots.c $ 
#proof[
  #set math.equation(numbering: none)
  Let $F(lambda) = e^(lambda A) B e^(-lambda A)$, its Taylor expansion at $lambda = 0$ is:
  $ & F^prime (lambda)|_(lambda=0) = (A e^(lambda A) B e^(-lambda A) - e^(lambda A) B e^(-lambda A) A)|_(lambda=0) = (\AF - \FA)|_(lambda=0) = [A,B] \
    & F^(prime prime) (lambda)|_(lambda=0) = (A F^prime - F^prime A)|_(lambda=0) = [A, [A, F]]|_(lambda=0) = [A, [A,B]] \
    & F^(prime prime prime) (lambda)|_(lambda=0) = (A F^(prime prime) - F^(prime prime) A)|_(lambda=0) = (A[A, [A, F]] - [A, [A, F]]A)|_(lambda=0) = [A, [A, [A, F]]] = [A, [A, [A, B]]] \
    & dots.v \ $
    Finally let $lambda=1$. 
]

Then the system operators $s_i$ and $s_i^dagger$ under interaction picture look like
$ & tilde(s_i) = s_i + (-hbar omega_i (\it)/hbar) s_i + ((-hbar omega_i (\it)/hbar))^2/2! s_i + ((-hbar omega_i (\it)/hbar))^3/3! s_i + dots.c  = e^(-i omega_i t) s_i \
  & tilde(s_i^dagger) = e^(i omega_i t) s_i^dagger \
  & \ $<bath>

Combine #ref(<derivation>), #ref(<corr>), and #ref(<bath>) together, and let $gamma_i$ depict the coupling strength, we get
$ & dot(rho_(S,I)) (t) = -sum_i integral_0^t gamma_i delta(t-t^prime) (tilde(s_i^dagger) (t) tilde(s_i) (t^prime) rho_(S,I) (t) - tilde(s_i) (t^prime) rho_(S,I) (t) tilde(s_i^dagger) (t) \
  & - tilde(s_i) (t) rho_(S,I) (t) tilde(s_i^dagger) (t^prime) + rho_(S,I) (t) tilde(s_i^dagger) (t^prime) tilde(s_i) (t)) \dt^prime \
  & = -sum_i integral_0^t gamma_i delta(t-t^prime) (s_i^dagger (t) s_i (t^prime) rho_(S,I) (t) e^(i omega_i (t-t^prime)) - s_i (t^prime) rho_(S,I) (t) s_i^dagger (t) e^(i omega_i (t-t^prime)) \
  & - s_i (t) rho_(S,I) (t) s_i^dagger (t^prime) e^(i omega_i (t^prime - t)) + rho_(S,I) (t) s_i^dagger (t^prime) s_i (t)) e^(i omega_i (t^prime - t)) \dt^prime \ $
Note that $integral_0^t e^(i omega (t-tau)) delta(t-tau) d tau = 1/2$, then
$ & dot(rho_(S,I)) (t) = -sum_i gamma_i/2 (s_i^dagger (t) s_i (t) rho_(S,I) (t) - s_i (t) rho_(S,I) (t) s_i^dagger (t) - s_i (t) rho_(S,I) (t) s_i^dagger (t) + rho_(S,I) (t) s_i^dagger (t) s_i (t)) \
  & = -sum_i gamma_i/2 ({s_i^dagger (t) s_i (t), rho_(S,I) (t)} - 2 s_i (t) rho_(S,I) (t) s_i^dagger (t)) \ $

This is how system-bath coupling affects the dynamics of system density matrix. Going back to Schrodinger picture, and adding on unitary part finally gives
$ dot(rho_(S,S)) (t) = -i/hbar [H_"sys", rho_(S,S) (t)] $






