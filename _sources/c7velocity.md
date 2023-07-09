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

# Chapter 7: The Velocity Four-Vector


## Background

In the slow moving world of dump trucks and physics professors, the
physical velocity of an object as measured by an observer is given by
the vector quantity:
```{math}
:label: oldvel
\vec{v} = \frac{d\vec{r}}{dt} = (dx/dt,dy/dt,dz/dt) = v_x\hat{x}
+v_y\hat{y}+v_z\hat{z}
```
A second observer, moving at velocity $\vec{v}_R$ in some direction,
measures this same object moving with velocity:
```{math}
:label: olddvel
\vec{v}' = \vec{v}-\vec{v}_R = \frac{d\vec{r}'}{dt} = (dx'/dt,dy'/dt,dz'/dt) 
```

The value of the time interval dt is the same for both observers.  Prior
to Einstein, it would have never occured to anyone to question whether
$dt$ might be different for different observers.

In Equation {eq}`oldvel`, the derivative can be understood in two
ways. It is operator that transforms the displacement vector, that can
be a function of time, into the velocity vector (that also can be a
function of time). A second way of understanding the derivative is
that it is a ratio of two elements, a very very small vector change in
the displacement divided by a very very small change in scalar time.

Physical velocity is a vector because the displacement, $d\vec{R}$ is a
3-vector and $dt$ is a scalar. A vector divided by a scalar is a vector.

## The Four Velocity

We cannot define a four velocity by simply taking the displacement
four vector $[dx_4]$ and dividing by $dt$, because $dt$ is not a
scalar.  It changes under a Lorentz transformation, and therefore
$[dx_4]/dt$ will not have the properties of a four vector as defined
in Chapter 2.  However, the **proper** time interval, $dt_0$, **is** a
scalar!  So we can define a four velocity as
```{math}
:label: eq4vdef
[v_4] = \frac{[dx_4]}{dt_0}
```
where $dt_0$ is the Lorentz scalar as described in Chapter 6. Since
$dt_0$ is a scalar, $1/dt_0$ is also a scalar.  The
4-velocity is a proper 4-vector because the displacement 4-vector is a
proper 4-vector and $1/dt_0$ is a Lorentz scalar. A Lorentz scalar times a
proper 4-vector is a proper 4-vector.

Wtien you multiply a vector by a scalar, what really happens is that
you multiply each of the components of the vector by that scalar
quantity. Putting in the components of the displacement 4-vector into
equation {eq}`eq4vdef` and multiplying each by $1/dt_0$ gives:
```{math}
:label: eq4v0
[v_4] = \frac{[dx_4]}{dt_0} =
\begin{bmatrix}
ic\frac{dt}{dt_0}\\
\frac{dx}{dt_0}\\
\frac{dy}{dt_0}\\
\frac{dz}{dt_0}
\end{bmatrix}
```

The derivatives that are the three spatial components of the velocity
4-vector are not the velocity of the object, as $dx$ is in measured in
a different reference frame than $dt_0$.  Within a given reference
frame the motion of an object would be written as $dx/dt$, where $t$
is measured **in that reference frame** (an observer in one reference
frame does not "just know" what is going on in another reference
frame, after all).  This would be the velocity as Isaac Newton would
have understood the term.

However, we can use the time-dilation formula from Equation {eq}`timedilation`
to convert $dt_0$ to $dt$.  Equation {eq}`eq4v0` then can be
written as:
```{math}
:label: eq4vgam
[v_4] = \frac{[dx_4]}{dt_0} = \gamma \frac{[dx_4]}{dt}
\begin{bmatrix}
ic\gamma\\
\gamma v_x\\
\gamma v_y\\
\gamma v_z
\end{bmatrix}
```

According to the correspondence principle, the theory we call special
re(ativity has to reduce to Newton's mechanics when $\beta\ll 1$. As
$\beta$ approaches 0, $\gamma$ approaches 1. Notice that the three
space terms of the 4-velocity are just the components of the physical
velocity as Newton would know it. The time component is just the speed
of light. At this point, the 4-velocity meets the requirements of the
correspondence principle. If Newton had been able to observe things
that were moving very fast with respect to him, he probably would have
noticed that the factor of $\gamma$ would have to be accounted for in
the definition of the physical velocity of the object.


```{note}
From here on, there are going to be several, perhaps even many,
different instances of the factor $\gamma$.  In all cases, $\gamma =
1/\sqrt{1-\beta^2}$, but $\beta$ is a speed, and therefore must be
relative to something else.  In the case of Equation {eq}`eq4vgam`,
the speed in the $\gamma$ is the speed of whatever it is that has the
velocity this four-vector is describing.  It relates to the time
interval between two events, compared with those two events at rest.
This is **not** necessarily the same thing as the relative speed
between two reference frames, which I will call $\beta_R$ and
$\gamma_R$.  You will soon see cases where there are three reference
frames of interest: the rest frame of an object, a frame in which the
object is moving, and then a third frame in which the object will be
moving at a different speed.  There are therefore five velocities that
could be of interest (the velocity of the object in the second frame,
the velocity of the object in the third frame, the relative velocity
between frame two and three, and the two relative velocities between
frames two and three, each, to the rest frame).  Each of these
velocities could have a useful $\gamma$ associated with it.  From here
on, you **must** pay attention to **which** $\beta$ speed goes with
which $\gamma$!
```