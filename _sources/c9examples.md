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

# Chapter 9: Examples of Relativity in Action

In this chapter we will explore how to use the energy-momentum four
vector to understand three example situations.  First, in the
collision of two particles, how much energy is available to
potentially create more particles?  Second, how to apply the Lorentz
transformation to the momentum four vector when the particle in the
system is a photon?  And finally, what happens to the energy of a
photon when it collides with (also called "scatters off of") a free
electron?

## The Center of Momentum Frame

The traditional scattering experiment that physicists use to probe the
nature of the interaction between pieces of a system is typically done
in a laboratory where one particle is accelerated and smashed into a
stationary target particle.  The energy, momentum, etc. of the debris
from the interaction are then measured.  Other scattering experiments
use colliding beams, where two beams of particles traveling in
opposite directions are allowed to interact. While this gives more
energy to do interesting things, the probability of interactions
occuring are smaller.

In the simplest case, there is one incident particle (with a rest
energy of $E_{oi}$, but moving with a momentum $p_i$ and a kinetic
energy of $KE_i$) and one target particle (at rest, with a rest energy
of $E_{0t}$) that interact in a particular place at a particular time.
Take the $x$ axis to be parallel to the momentum vector of the
incident particle.  Then we can write the four-momentum of the
incident particle as
```{math}
:label: p4inc
[p_4]_i =
\begin{bmatrix}
i(E_{0i} + KE_i)/c\\
p_i\\
0\\
0
\end{bmatrix}
```
and the four-momentum of the target particle is
```{math}
:label: p4tar
[p_4]_t =
\begin{bmatrix}
iE_{0t}/c\\
0\\
0\\
0
\end{bmatrix}
```

The total 4-momentum for the scattering as determined by an observer
in the laboratory reference frame is (writing the incident momentum in
terms of the KE using Equation {eq}`pfrome`):
```{math}
:label: p4tot9
[p_4]_{\rm lab} =
\begin{bmatrix}
i\frac{1}{c}(E_{0i} + KE_i + E_{0t})\\
\frac{1}{c}\sqrt{KE_i^2 + 2 E_{0i}KE_i}\\
0\\
0
\end{bmatrix}
```

When the scattering interaction happens, new particles can be formed,
but the momentum 4-vector describing the system remains constant (is
conserved) as long as you stay in the same inertial reference
frame. All of the appropriate quantum numbers must also be conserved,
e.g. charge, baryon number, etc.

Suppose that there is a second observer, watching this same scattering
process take place. This observer is moving in the same direction as
the incident particle, but at a slower speed. This observer sees the
incident particle moving slower, but now the target particle has some
physical momentum that is pointing in the opposite direction that the
observer is traveling.


At some relative speed $\beta_{\rm com}$ the second observer will
measure the incident particle and the 'target' particle to have the
same size physical momentum, but pointing in opposite directions. The
total physical momentum in this case is 0. This reference frame is
known as the zero-momentum or center-of-momentum (com) reference
frame. The 4-momentum for this frame is:
```{math}
:label: p4com
[p_4]_{\rm com} =
\begin{bmatrix}
iE_{\rm com}/c\\
0\\
0\\
0
\end{bmatrix}
```

The center-of-momentum frame is a very special case.  After the
interaction takes place, the physical momentum of all of the particles
still must equal 0. If this scattering process were to make different
particles than the two incident particles, the threshold (the lowest
energy at which it could occur) would have the new particle produced
each with 0 momentum (and KE).  This means that $E_{\rm com}$ is
the maximum energy available to make new particles in the interaction.
We say maximum because the KE of the daughter particles won't actually be zero
in a real situation.  If the two parent particles survive the interaction,
their rest mass must be accounted for out of the total energy, too.
The value $E_{\rm com}$ is a budget of energy that must be allocated to
the rest energies and kinetic energies of everything that comes out of
the interaction.

The time component of the 4-momentum in the com reference frame (the
total energy in the center-of-momentum) tells which new particles can
be produced in a scattering reaction. If there is any energy left over
after creating the rest energies of the new particles, it goes into
the kinetic energy of the new particles. The kinetic energy is
distributed in such a way that the total physical momentum of the
particles is 0.

If the 4-momentum in the lab is known, then the
4-momentum in the com frame can be found using the Lorentz
transformation:
```{math}
:label: p4comlort0
[p_4]_{\rm com} = {\cal L}_x (\beta_{\rm com})[p_4]_{\rm lab}
```
Plugging in Equation {eq}`p4tot9` (and ignoring $y$ and $z$) we get
```{math}
:label: p4comlort1
[p_4]_{\rm com} =
\begin{bmatrix}
\gamma_{\rm com} & -i\beta_{\rm com}\gamma_{\rm com} \\
i\beta_{\rm com}\gamma_{\rm com} & \gamma_{\rm com} \\
\end{bmatrix}
\begin{bmatrix}
i\frac{E_{\rm lab}}{c}\\
p_i\\
\end{bmatrix}
```
Matrix multiplication yields:
```{math}
:label: p4com2
\begin{bmatrix}
iE_{\rm com}/c\\
0
\end{bmatrix}
=
\gamma_{\rm com}
\begin{bmatrix}
i(E_{\rm lab}/c-\beta_{\rm com}p_i)\\
p_i-\beta_{\rm com}E_{\rm lab}/c\\
\end{bmatrix}
```

The second row of Equation {eq}`p4com2` tells how to find the relative
speed of the com reference frame if you know the total energy and momentum in
the lab reference frame.  The two terms must cancel, if the total momentum
is to be zero, and that can only happen if
```{math}
:label: betacom
\beta_{\rm com} = \frac{cp_i}{E_{\rm lab}} =
\frac{\sqrt{KE_i^2 + 2 E_{0i}KE_i}}{E_{0i} + KE_i + E_{0t}}
```
There really isn't any further simplification one can do to
Equation {eq}`betacom`.  Every quantity on the right would be
known as inputs for a given experiment, so you could calculate
the necessary speed of the com for any target particle being
hit by any incident particle with any kinetic energy you choose.

We can, however, apply a limit check and consider what happens as the
KE of the incident particle gets very, very large.  In such a limit,
the rest energies become negligible and $\beta_{\rm com} \rightarrow
1$.  This makes sense because the faster the incident particle is
going, the closer its speed gets to one, the more it dominates the sum
that makes up the total momentum, so effectively the target plays no
role, and the com is going right along with the incident particle.

If the target and the incident particle are the same kind of particle,
such that $E_{0i}=E_{0t}\equiv E_0$, then we can do one further
simplification -- pull out a $\sqrt{KE_i}$ in the numerator:
```{math}
:label: betacom2same
\beta_{\rm com} = \sqrt{KE_i}
\frac{\sqrt{KE_i + 2 E_0}}{KE_i + 2E_0}
=  \sqrt{\frac{KE_i}{KE_i + 2E_0}}
= \left(1+\frac{2E_0}{KE_i}\right)^{-1/2}
```
In this case you can see, again, that as $KE_i$ gets large,
$\beta_{\rm com}\rightarrow 1$, but it's easier to see.

The first row of Equation {eq}`p4com2`, on the other hand,
tells us the total energy available in the com frame:
```{math}
:label: Etotcom
E_{\rm com} = \gamma_{\rm com}(E_{\rm lab} - \beta_{\rm com} c p_i)
```
The fact that we are subtracting something from $E_{\rm lab}$ on the
right side of Equation {eq}`Etotcom` makes it look like the amount of
energy in the com frame is less than the amount of energy in the lab
frame, but you have to be cautious with such generalizations, because
$\gamma_{\rm com}$ is also getting larger as $p_i$ increases, so it's
not so trivial to see which effect will dominate.

Need to make a graph of Ecom vs. KEi for two protons.


## Photons and the Doppler Shift

## Compton Scattering
