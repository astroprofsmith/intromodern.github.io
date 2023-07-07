---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{code-cell}
:tags: ["remove-cell"]

from IPython import display
import numpy as np
import matplotlib.pyplot as plt

```

# Chapter 4: Lorentz Transformations

We have already seen that requiring the size of a four-vector
to be invariant across reference frames leads to interesting
implications (time dilation, for example).  However, there are
more insights that can be gained by figuring out a method to
predict what the components of a four vector will be in a second
reference frame, given measurements of that four vector's components
in a first reference frame.  The change of a four vector from
one reference frame to another with a constant relative velocity
is called a **Lorentz Transformation**.


## Setting Up the Problem

The Lorentz transformation converts the components of a 4-vector
in one inertial reference frame to the components as measured in
a second frame moving at a constant velocity relative to the first.
This method works for **any** four-vector, but we will define it
using the displacement four-vector, as that is easy to visualize.

We start by recreating Figure 3.2 as Figure 4.1 -- to make sure
you are applying the transformation properly, it is useful to always
start by recreating this diagram and making sure that you are assigning
the correct numbers to the elements, as defined in this diagram.



```{code-cell}
:tags: ["remove-input"]
# 3D plot of a spacetime diagram with x, ct, and y
plt.figure(figsize=(5,5))
plt.arrow(0,0,1,0,head_width=0.1)
plt.arrow(0,0,0,1,head_width=0.1)

plt.arrow(-0.5,-0.5,1,0,head_width=0.1)
plt.arrow(-0.5,-0.5,0,1,head_width=0.1)

plt.arrow(0.2,-0.25,0.75,0,head_width=0.05)

plt.plot([0.2],[0.1],'ro')
plt.plot([0.9],[1.1],'bo')
plt.arrow(0.2,0.1,.7,1.0,head_width=0.15)

ax = plt.gca()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.axis([-1,2,-1,2])
ax.text(0.15, 1, "S")
ax.text(-0.05, 1.2, "ct")
ax.text(1.2, 0, "x")
ax.text(1.1, -0.25, "v_R")
ax.text(-0.35, .5, "S'")
ax.text(-.55, .7, "ct'")
ax.text(.7, -0.5, "x'")
ax.text(1.0,1.0,"dx_4")
plt.show()
```
```{note}
Figure 4.1 -- The spacetime diagram that defines a Lorentz transformation.
There are two events, indicated by red and blue dots, and we define a
displacement four vector between them (indicated by an arrow).  These
components are defined in the reference frame $S$, and we want to know
the components in a second frame $S^\prime$, which is moving with a
constant relative velocity $\vec{v}_R$.  We **define** $x$ such that it lies
parallel to $\vec{v}_R$.
```


The 4-vector components are determined by the two events shown as
colored dots in Figure 4.1. Since we have the choice of which is the
x-axis, we will always pick the relative velocity of the two observers
to lie along the $x$-axis, to make the math simpler. Putting the
diagram of Figure 4.1 into 4-vector notation, we need a tool that
will take us from $[dx_4]$ to $[dx_4]^\prime$, like so:
```{math}
:label: eqlort
[dx_4]^\prime = {\cal L}_x(\beta_R)[dx_4]
```

The Lorentz transformation can be mathematically represented as a
matrix, because it transforms one vector into another. The Lorentz
transformation, therefore, is a $4 \times 4$ matrix because there are
4 elements in the vector. If $L_{\alpha\beta}$ is a an arbitrary
element of this matrix (row $\alpha$ column $\beta$, where the
Greek letters stand for 0, 1, 2, and 3), the transformation
matrix has the form:
```{math}
:label: lormat
{\cal L}_{\alpha\beta} = 
\begin{bmatrix}
L_{00} & L_{01} & L_{02} & L_{03}\\
L_{10} & L_{11} & L_{12} & L_{13}\\
L_{20} & L_{21} & L_{22} & L_{23}\\
L_{30} & L_{31} & L_{32} & L_{33}
\end{bmatrix}
```
The problem then is to figure out what the 16 numbers in ${\cal L}$
are that will achieve our goal.

```{warning} 
Note that the $\beta$ used as the index counter of the matrices has
nothing to do with the relative velocity $\beta_R$.  Sometimes other
Greek letters are used for counting indices, $\mu$ and $\nu$ in
particular.  You must notice how the letter is being used to be
able to tell whether $\beta$ is a speed or a matrix index.
```

## Determining the Elements

Is there a way that we can determine some of the elements of the of
the Lorentz transformation? We must figure out how the transformation
works in specific cases, and these cases will then point to general
values that could be used in any situation.  All that we really know
about the Lorentz transformation is that the 4-vector that is produced
by its operation on the original 4-vector must have the same size as
the first 4-vector, or the Michelson-Morley experimental results would
be contradicted.

To determine some of the elements of the Lorentz transformation
matrix, consider what happens when the displacement 4-vector
```{math}
:label: eqdisp4
[dx_4] =
\begin{bmatrix}
icdt\\
dx\\
dy\\
dz
\end{bmatrix}
```
is transformed from one inertial reference frame to another. If the
relative velocity is in the $x$-direction, then displacements in the
$y$ and $z$ directions won't change.  This derives from the independence
of motion along perpendicular axes: if I am 10 m south of you, and I
take one step East, I will still be 10 m South of you (assuming we are
far from the poles, of course).

```{note}
Although the displacements perpendicular to the direction of relative
motion are not affected by a Lorentz transformation, this independence
does not apply to **any** object that you might perform a Lorentz transformation
on.  In electromagnetism, for example, the electric and mangetic fields in
the $y$ and $z$ direction are different in a frame moving in the $x$ direction.
We use the dimensional independence of $dy$ and $dz$ to derive values for
${\cal L}$, but once we do that we don't have to insist that $y$ and
$z$ are unaffected in all circumstances at all times.  Be cautious.
```

This means that $dy^\prime = dy$.  However, if we write out the third
row of Equation {eq}`eqlort`, we get $$dy^\prime = L_{20}icdt+L_{21}dx
+L_{22}dy + L_{32}dz.$$
The only way this equality can hold for all displacement four vectors
is if $L_{22}=1$ and the other elements are zero.  The same logic hold
for the fourth row, so we know the matrix must look like
```{math}
:label: lormat2
{\cal L} = 
\begin{bmatrix}
L_{00} & L_{01} & L_{02} & L_{03}\\
L_{10} & L_{11} & L_{12} & L_{13}\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
```
That takes care of half the numbers in one go!

To get the next set, we have to look at a particular case.
Consider the examples at the end of Chapter 3.  They all
involved switching from a reference frame where two events
were at rest into a reference frame in which there would have
to be motion for something to get from one event to the other.
We know therefore that the original four displacement is
```{math}
:label: eqrestdx4
[dx_4] =
\begin{bmatrix}
icdt_0\\
0\\
0\\
0
\end{bmatrix}
```
We can also figure out what the transformed components have to be.  We
know that the time interval will be dilated in the primed frame: $dt =
\gamma_R dt_0$.  Furthermore, if the events are at rest in the
original frame, and the primed frame is moving to the right at speed
$\beta_R$, then from the point of view of the primed frame, the
location of the events would have to be moving left at the same speed,
so $dx^\prime = -v_Rdt = -\beta_R cdt = -\beta_R\gamma_R cdt_0$, where
we have used time dilation to turn the $dt$ into $dt_0$.

That gives us two further relations that will enable us to get four of the
numbers.  We look at the first two rows of Equation {eq}`eqlort`:
$$icdt^\prime = L_{00}icdt_0+L_{01}dx_0$$
and
$$dx^\prime = L_{10}icdt_0+L_{11}dx_0.$$
First of all, notice that $dy$ and $dz$ do not enter into these
equations.  This means that the elements that would multiply
$dy$ and $dz$ ($L_{02}$, $L_{03}$, $L_{12}$, $L_{13}$) must
all be zero!  Second, we use the fact that $dx_0=0$ (events at rest)
to drop those terms out of the equations.  Divide what's left
by each other to get
```{math}
:label: eqLrat
\frac{dx^\prime}{}
```