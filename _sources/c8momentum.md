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

# Chapter 8: The Momentum Four-Vector

You have been introduced to the displacement 4-vector, the Lorentz
transformation and to the velocity 4-vector. These quantities can be
used to locate objects and determine their motion. Borrowing an idea
introduced by Newton to predict the behavior of dynamical systems, we
introduce the momentum 4-vector. The 4-velocity is a proper 4-vector
in that its compo- nents in two different reference frames are related
by the Lorentz transfor- mation, and its size is a Lorentz scalar as
demanded by the Michelson- Morley experiment.


Knowing that if I multiply a proper 4-vector by a Lorentz scalar, I will get another proper 4-vector, I introduce:
```{math}
:label: p4def
[p_4]=m_0[v_4]
```
where $m_0$ is the mass as measured by an observer in the rest frame
of the object. As is true with the time interval as measured in the
rest frame, the 'rest mass' , the mass measured in the rest frame of
the object, is a Lorentz invariant.

## Momentum

The components of the momentum 4-vector introduced in equation
{eq}`p4def` are
```{math}
:label: p4comp
[p_4] = m_0[v_4] = m_0 \gamma \frac{[dx_4]}{dt}
\begin{bmatrix}
icm_0\gamma\\
m_0\gamma v_x\\
m_0\gamma v_y\\
m_0\gamma v_z
\end{bmatrix}
```

```{margin}
Example 8.1:

If a particle of rest mass $m_0 = \\1.00\times10^{-27}$ kg
is travelling at a speed of $\beta = 0.866$ in the $x$
direction, what is its momentum?

At that speed, $\gamma=2$ and therefore
$\vec{p} = 2.0(1.00\times10^{-27}~{\rm kg})(0.866)$
$(3\times 10^8~{\rm m/s})\hat{x}$
$=5.20\times 10^{-19}~{\rm kg~m/s}~\hat{x}$.

This is twice what
Newton would have expected the momentum of the particle to be.
```

The three spatial components are almost the classical components of the
momentum of a particle, which is usually defined as $\vec{p}\equiv
m\vec{v}$.  What Newton missed was the factor of $\gamma$. It is not
surprising that this happened as the value of $\gamma$ for the things
Newton could observe
is different from the value 1 by just a few parts in $10^{12}$,
hardly worth taking into account.

So, let us re-define the physical momentum of a particle as:
```{math}
:label: pdef
\vec{p}\equiv m_0 \gamma \vec{v}
```
which means Eqaution {eq}`p4def` becomes
```{math}
:label: p4full
[p_4] =
\begin{bmatrix}
icm_0\gamma\\
p_x\\
p_y\\
p_z
\end{bmatrix}
```
To check if this 4-vector does follow the dictates of the
Michelson-Morley experiment, we use the form for the
4-momentum found in Equation {eq}`p4comp` and take the dot
product of the 4-vector with itself:
```{math}
:label: pmag
[p_4]^2 = m_0^2\gamma^2(-c^2+v_x^2+v_y^2+v_z^2)
= -m_0^2\gamma^2(1-v^2/c^2) = -m_0^2c^2
```
Since $m_0$ is an invariant and so is $c$, the length of the
momentum 4-vector is a Lorentz invariant. Once again, it's
negative.

f a system is made up of more than one particle, the total
4-momentum for the system is the sum of the momentum 4-vectors
of each of the pieces. It is the size of the net (sum of
momentum of each of the pieces) 4-momentum that is a Lorentz
scalar.

## Momentum and Lorentz Transformations


The components of the momentum 4-vector transform from those
measured by one inertial observer to those measured by a second
inertial observer moving with respec to the first with $\beta_R$
in the $+x$ direction in the same manner as any other 4-vector:
```{math}
:label: p4lort
[p_4]' = {\cal L}_x(\beta_R)[p_4]
=
\begin{bmatrix}
\gamma_R & -i\beta_R\gamma_R & 0 & 0\\
i\beta_R\gamma_R & \gamma_R & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
icm_0\gamma\\
m_0\gamma v_x\\
m_0\gamma v_y\\
m_0\gamma v_z
\end{bmatrix}
=
\gamma_Rm_0\gamma
\begin{bmatrix}
i(c-\beta_Rv_x)\\
-\beta_Rc+v_x\\
v_y\\
v_z
\end{bmatrix}
```
or
```{math}
:label: p4p
[p_4]' =
\gamma_R\gamma m_0c
\begin{bmatrix}
i(1-\beta_R\beta_x)\\
\beta_x-\beta_R\\
\beta_y\\
\beta_z
\end{bmatrix}
=
\begin{bmatrix}
icm_0\gamma'\\
p'_x\\
p'_y\\
p'_z
\end{bmatrix}
```
From this we can deduce predictions like
```{math}
:label: gamp
\gamma' = \gamma_R\gamma(1-\beta_R\beta_x)
```
or
```{math}
:label: pxp
p'_x = \gamma_R\gamma m_0c(\beta_x-\beta_R)=\gamma'm_0v'_x
```
Divide Equation {eq}`pxp` by Equation {eq}`gamp` to get the
velocity addition formula again!

## Example of Transforming Momentum

Suppose that two observers in different reference frames were
measuring the momentum of the particle in Example 5.1 above.
Say the primed observer is moving with $\beta_R=0.866$ in the
$+x$ direction with respect to the unprimed observer.  The
unprimed observer would measure a momentum 4-vector of:
```{math}
[p_4] =
\begin{bmatrix}
i6.00\times10^{-19}~{\rm kg~m/s}\\
5.20\times10^{-19}~{\rm kg~m/s}
\end{bmatrix}
```
The components of the momentum 4-vector in the primed co-ordinate
system are found by using the Lorentz transformation:
```{math}
[p_4]' =
\begin{bmatrix}
2 & -i1.73\\
i1.73 & 2
\end{bmatrix}
\begin{bmatrix}
i6.00\times10^{-19}~{\rm kg~m/s}\\
5.20\times10^{-19}~{\rm kg~m/s}
\end{bmatrix}
```
```{math}
[p_4]' =
\begin{bmatrix}
i3.00\times10^{-19}~{\rm kg~m/s}\\
0.0
\end{bmatrix}
```

The primed frame of reference is moving at the same speed as the
particle, so the primed observer is in the rest frame of the
particle, hence the physical momentum should be zero. Which it is!
The time component in the rest frame should be the rest mass times
$c$, which you can easily check is the case here.


## Momentum and Energy

