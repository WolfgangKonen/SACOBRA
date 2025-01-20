# SACOBRA
Self-adjusting constrained optimization in R.

SACOBRA is a package for numerical **constrained optimization** of expensive black-box functions under **severely limited budgets**. The problem to solve is: 
```math 
\mbox{Minimize}\quad  f(\vec{x}) , \vec{x} \in [\vec{a},\vec{b}] \subset \mathbf{R}^d 
```
$$ \mbox{subject to}\quad g_i(\vec{x}) \le 0, i=1,\ldots,m    $$
$$ \mbox{~~~~~~~~~~}\quad\quad h_j(\vec{x}) = 0, j=1,\ldots,r.    $$

SACOBRA performs optimization with a minimum of true function evaluations. It has proven to work well on problems with high dimensions (e.g. $d=124$) and many constraints (e.g. $m+r=60$). It is usable for all kind of **numerical** optimization of continuous functions, but not for combinatorial optimization or discrete optimization.

SACOBRA can be used for unconstrained optimization problems as well, but its focus is on constrained optimization. Other packages might be more suitable for unconstrained optimization.

< SACOBRA logo >

## Installation

To install from GitHub:
```
install.packages("devtools")    # if not yet installed

require(devtools)
install_github("WolfgangKonen/SACOBRA", dependencies = NA)
```

## Usage
< how to set up function and constraints >

### Examples
< how to run examples and where the examples are >
< an example with a nice image>



