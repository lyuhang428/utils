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
#align(center,text(21pt)[Some Notes])
\
#align(center,text(12pt)[lyh])
\
\
_bold letters stand for vectors unless otherwise stated_
\
\

/////////////
// main texts
/////////////
= Math preparations
- MathWorld generalized Pauli algebra: #link("https://mathworld.wolfram.com/GeneralizedGell-MannMatrix.html"). 

- Direct product and partial trace: #link("https://zhuanlan.zhihu.com/p/653816083")
- Math fonts: #link("https://typst.app/docs/reference/math/variants/")

\
If ${mid(|)phi_i^A angle.r mid(|) i=1,2,dots.c,N}$ and ${|phi_i^B angle.r mid(|)i=1,2,dots.c,N}$ are basis sets of space $A in bb(C^N)$ and $B in bb(C^N)$, respectively. Let new space $C = A otimes B in bb(C^{N times N})$, and basis set in this space will be ${|phi_i^A angle.r otimes |phi_j^B angle.r mid(|) i=1,2,dots.c,N; j=1,2,dots.c,N}$. Vectors on $C$ can be expanded on this basis set $|f^(\AB) angle.r in C = sum_(\ij) c_(\ij) |phi_i^A angle.r otimes |phi_j^B angle.r$. When we write down $|f^(\AB) angle.r$ as state vector, we implicitly specify a pure state, becasue mixed state cannot be expressed as state vector in Hilbert space. The density matrix is 

$ & rho^(\AB) = |f^(\AB) angle.r angle.l f^(\AB)| \
  & = (sum_(\ij)c_(\ij)|phi_i^A angle.r otimes|phi_j^B angle.r) (sum_(\kl)c_(\kl)^ast angle.l phi_k^A|otimes angle.l phi_l^B|) \
  & = sum_(\i\j\kl)c_(\i\j\kl) |phi_i^A angle.r angle.l phi_k^A otimes |phi_j^B angle.r angle.l phi_l^B| \ $
Sometimes in a confusing notation $rho_(\AB) = sum_(\i\j\kl) |i angle.r_A angle.l j| otimes |k angle.r_B  angle.l l|$. 

// Assume $C = A\otimes B$, tracing out $A$ will give us $B$. Let $\{|\phi_i\rangle|i=1,2,\cdots\}$ be arbitrary complete orthonormal basis set, and $I$ the identity matrix.
// $$\begin{split}
//     & \sum_i \left(\langle\phi_i|\otimes I\right) C \left(|\phi_i\rangle\otimes I\right) \\
//     & = \sum_i \left(\langle\phi_i|\otimes I\right) \left(A\otimes B\right) \left(|\phi_i\rangle\otimes I\right) \\
//     & = \sum_i \langle\phi_i|A\otimes IB |\phi_i\rangle \otimes I \\
//     & = \sum_i \langle \phi_i | A | \phi_i \rangle \otimes IBI \\
//     & = \sum_i \langle \phi_i | A | \phi_i \rangle \otimes B \\
//     & = \text{Tr}A \otimes B = \text{Tr}A\cdot B \\
// \end{split}$$

// Similarily,
// $$\begin{split}
//     & \sum_i \left( I \otimes \langle\phi_i| \right) A\otimes B \left( I \otimes |\phi_i\rangle \right) \\
//     & = \sum_i I A \otimes \langle\phi_i| B I\otimes |\phi_i\rangle \\
//     & = \sum_i IAI \otimes \langle \phi_i| B |\phi_i \rangle \\
//     & = \text{Tr}B\cdot A \\
// \end{split}$$

// This is how partial trace is performed. Partial trace is an inverse operation of tensor product. If $A$ and $B$ are known, and $A = B\otimes C$, then $C = \frac{\text{Tr}_\text{B}A}{\text{Tr}B}$. 
