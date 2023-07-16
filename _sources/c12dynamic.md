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
from myst_nb import glue

```

# Chapter 12: Relativistic Electrodynamics

Now that we have an expression for the electromagnetic
field tensor, which treats $\vec{E}$ and $\vec{B}$ as
different parts of a single object, and we know how that
tensor transforms under Lorentz transformations, we can
start to see how the rest of electrodynamics plays out
within an SR context.  We will start with examining
charge and current.

## Charge/Current Four Vector and Continuity

We spent much of the first ten chapters figuring out how to cast the
three-vectors of Classical Mechanics into four-vector form.  We cn
play the same game with the vectors of Classical Electrodynamics.  One
of the most important is the current density $\vec{J}$, which
describes how much charge flows through a given area in a given amount
of time.  The dimensions are therefore charge per area per time, or
current per area, which is a kind of density.  Mathematically, we
write this as
```{math}
:label: currdenJ
\vec{J} = \rho \vec{v}
```
which is a volume charge density (charge/volume) multiplied by a
velocity (length/time).  One factor of length cancels out,
leaving charge per area per time, as required.

Knowing what we know now, we can recognize that the lengths between
these moving charges has been contracted along the direction of
motion, so we can go ahead and write this in terms of what the density
would be in a reference frame where the charges are at rest:
```{math}
:label: currdenJ0
\vec{J} = \gamma \rho_0 \vec{v}
```
where $\gamma =1/\sqrt{1-v^2/c^2}$, of course.  Although the density
is charge over volume, only the dimension along the direction of
motion gets contracted, so the density only changes by one factor of
$\gamma$.

This vector field can make up the spacial part of the current density
four-vector, but what do we do with the time part?  Well, we have
a scalar $\rho_0$ times $\gamma \vec{v}$, and $\gamma \vec{v}$ is
the spatial part of the four-velocity.  Maybe we can make the
four-current by multiplying $[v_4]$ by $\rho_0$?  And indeed we can:
```{math}
:label: fourcurr
\boxed{
[J_4] =
\rho_0 \gamma 
\begin{bmatrix}
ic\\
v_x\\
v_y\\
v_z
\end{bmatrix}
}
```

And its size is therefore $\rho_0^2[v_4]^2 = -\rho_0^2c^2$.  which is
certainly Lorentz invarient.  Since we know a four-vector times a
scalar is a four-vector, this must be *a* four-vector, but how do we
know it's *the* four vector that we want?  Two reasons: it slots very
smoothly into Maxwell's Equations, which we will explore in the next
section.  But first, we will see that it gives a nice way to express
charge conservation in the form of the "continuity equation", which
says that the rate at which charge is flowing into or out of a
location is equal and opposite to the rate at which the amount of
charge at that location is changing.

In classical notation, we would say
```{math}
:label: conteqclass
\vec{\nabla}\cdot\vec{J} = \frac{\partial \rho v_x}{\partial x}+ \frac{\partial \rho v_y}{\partial y}+ \frac{\partial \rho v_z}{\partial z} = -\frac{\partial \rho}{\partial t}
```
If we move the right side over the to left and multiply top and bottom by $c$,
it should start to suggest something:
```{math}
:label: movetimeover
\frac{\partial (\gamma \rho_0 c)}{\partial (ct)}+\frac{\partial \gamma \rho_0 v_x}{\partial x}+ \frac{\partial \gamma \rho_0 v_y}{\partial y}+ \frac{\partial \gamma \rho_0 v_z}{\partial z} = 0
```
You will note that each successive term being differentiated is indeed one of the
terms of our current four-vector.  Furthermore, each successive derivative is being
done with respect to successive variables of the four-displacement!  So this suggests
a kind of dot product, which we could write this way:
```{math}
:label: fourdivergence1
\left(\frac{\partial}{\partial (ct)},\frac{\partial}{\partial x},\frac{\partial}{\partial y},\frac{\partial}{\partial z}\right)
\begin{bmatrix}
\gamma \rho_0 c\\
\gamma \rho_0 v_x\\
\gamma \rho_0 v_y\\
\gamma \rho_0 v_z
\end{bmatrix}
= 0
```
If we define a **four-divergence** as the row vector above and write it as
```{math}
:label: fourdivdef
\partial_\alpha = \frac{\partial}{\partial x^\alpha} = \left(\frac{\partial}{\partial (ct)},\frac{\partial}{\partial x},\frac{\partial}{\partial y},\frac{\partial}{\partial z}\right)
```
then the continuity equation can be written as simply as
```{math}
:label: continuitySR
\boxed{
\partial_\alpha J^\alpha = 0
}
```
Here you can see how the Einstein notation starts to save you a lot of writing.  Equations
{eq}`continuitySR` and {eq}`movetimeover` are the exact same equation, with the same
information content, but Equation {eq}`continuitySR` requires *so* much less writing,
it's much more popular than other methods.

```{note}
You will note that the four-gradient $\partial_\alpha$ is written in
covariant form, but it is using the contravariant components for the
derivatives, so it is written with the subscript, but when you write
it out in terms of components, the terms in the denominator get the
superscript.  But there's no minus sign in this case, because the
components are the contravariant form.
```

So, this four vector makes intuitive sense (charge time velocity gives
current), it has the properties of a four-vector (because it's a
scalar $\rho_0$ times a four vector $v^\alpha$), it maintains charge
conservation through the continuity equation, and it sets us up to
write Maxwell's Equations in a more compact form, too (you'll recall
Maxwell's Equations depend on the charge and current densities on the
right side and the $\vec{E}$ and $\vec{B}$ fields on the left -- see
the next section of this chapter).  We therefore use this as the
current four-vector $J^\alpha$, or sometimes the "charge-current
four-vector".

## Maxwell's Equations as Four Vectors

## EM Forces

