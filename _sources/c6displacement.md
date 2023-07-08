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

# Chapter 6: The Displacement Four Vector

In this chapter, we will take a deep dive into the properties of
the displacement four vector, using all the tools developed in the
previous chapters.  Some of the material will repeat what has been
presented already, but in more depth.

## The Difference Between Two Events

The model we will be using to predict the behavior of systems moving
with speeds significant compared to the speed of light require that
the physical process be one such that there are two events occurring
that can be viewed by different observers. Some examples of physical
processes that can be understood this way are:


1. A strange particle being produced (first event), and at some later
time, decaying into other particles (second event).

2. A rocket ship leaving the earth (first event), and arriving at
  another star (second event).

3. A high energy electron striking a proton (first event) and then
  being scattered and leaving the interaction area at a different
  energy and direction of motion from the initial electron (second
  event).

4. A photon with energy in the radar range sent out by a cop (first
  event) that is reflected by a moving car, and then returns to the
  cop at a different energy (second event)

5. A star cruiser traveling at 0.900 the speed of light (first event)
  and then shooting of a laser pulse at an angle of $30^\circ$ with
  respect to the direction of motion of the star cruiser (second
  event).

6. A high energy proton entering the magnetic field of a bubble chamber
  (first event) and this same proton leaving the magnetic field
  traveling in a different direction (second event).

This model (which of course is special relativity) suggests how an
inertial observer, such as one that is moving at a constant speed (not
accelerating), should analyze what is different about the two
events. It also provides a way (called the Lorentz transformation) of
calculating what other inertial observers, moving with respect to the
first observer, would measure as the difference between these two
events.

## A Summory of The Story So Far

Before diving into the analysis, a quick summary of what we have
establishd so far.  Figure 6.1 reproduces Figure 4.1, showing two
events in spacetime diagram with a displadcement four vector between
them.



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
plt.arrow(0.2,0.1,.7,1.0,head_width=0.15,length_includes_head=True)

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
Figure 6.1 -- Two events in a spacetime diagram with a displacement
four vector between them.
```

The displacement between these two events is to be analyzed by two
inertial observers traveling at a relative velocity of $\beta_R$ with
respect to each other. By convention, the $x$-axis is chosen such that
the relative velocity points in the $+x$ direction. Each of the
observers measures the physical displacement $(dx,dy,dz)$ and the time
interval $(dt)$ between the two events. Using the measured values, each
observer then constructs the displacement 4-vector that describes
these two events:
```{math}
:label: eqtwo4disp
[dx_4] =
\begin{bmatrix}
icdt\\
dx\\
dy\\
dz
\end{bmatrix}
\hspace{1cm}
[dx_4]' =
\begin{bmatrix}
icdt'\\
dx'\\
dy'\\
dz'
\end{bmatrix}
```

Now, as an example, let us consider the displacement 4-vector as
measured by two different observers. For simplicitly, we assume the
two events as measured by the unprimed observer (shown in Figure 6.1)
have real displacements $dy$ and $dz$ equal O. We can therefore ignore
the last two components and write the displacement 4-vector in the
unprimed co-ordinate system as:
```{math}
:label: equnprime
[dx_4] =
\begin{bmatrix}
icdt\\
dx
\end{bmatrix}
```
The elements of the displacement 4-vector as measured by the primed
observer are found by operating on the unprimed displacement 4-vector
(equation {eq}`equnprime`) with the Lorentz transformation (equation
{eq}`lormat`) for relative velocity of $\beta_R$ in the $+x$ direction.
$$[dx_4]^\prime = {\cal L}_x(\beta_R) [dx_4]$$
```{math}
:label: lortrans
\begin{bmatrix}
icdt^\prime\\
dx'
\end{bmatrix}
=
\begin{bmatrix}
\gamma_R & -i\beta_R\gamma_R\\
i\beta_R\gamma_R & \gamma_R
\end{bmatrix}
\begin{bmatrix}
icdt\\
dx
\end{bmatrix}
```


```{margin}
Why didn't Newton notice this?

In Newton's time,a relative velocity of $1 \times 10^4$ m/s (the
orbital speed of the moon) was about as large as a relative velocity
he could see. This implies $\beta_R = 3.3 \times 10^{-5}$ and
$\gamma_R = (1 + 5\times10^{-11})$. The value for $\gamma_R$ is
really very close to 1, so it is not surprising that Newton missed
this small effect!  The time interval in any frame would have
been less than a trillionth of a fractional change away from any
other frame.  It's no wonder Newton assumed time intervals were
the same in any frame!
```

Multiplying out the matrix equation {eq}`lortrans` to find the
components of the displacement 4-vector as measured by the primed
observer in terms of the values of the components measured by the
unprimed observer gives:
```{math}
:label: eqdtp
cdt' = \gamma_R(cdt - \beta_R dx)
```
```{math}
:label: eqdxp
dx' = \gamma_R(dx - \beta_R c dt)
```

In particular, we note that if $dx=0$ and therefore $dt=dt_0$, then
equations {eq}`eqdtp` and {eq}`eqdxp` simplify to $dt' = \gamma_R
dt_0$ and $dx'=-\gamma_R\beta_R c dt_0 = -v_R dt' \rightarrow
dx'/dt'=-v_R$.  The first equation is **Time Dilation**, or that the
time interval between two events in their rest frame is the shortest
it can be, compared with any other frame in relative motion.  The
latter equation is simply indicating that if the primed frame is
moving right, the events that were at rest will be measured to have
moved to the left at the same speed.


The size of this four vector, also called the "interval"
between the events, is $$[dx_4]^2 = dx^2-c^2dt^2,$$ which
in the rest frame reduces to $-c^2dt_0^2$.  The interval in
the primed frame is therefore $$dx'^2-c^2dt'^2,$$
or
$$ \gamma_R^2(dx - \beta_R c dt)^2- \gamma_R^2(cdt - \beta_R dx)^2$$
Plugging in the definition for $\gamma_R$ and expaning the
binomial squares:
```{math}
:label: eqinvar
\frac{dx^2-2dx\beta_Rcdt+\beta_R^2c^2dt^2-c^2dt^2+2cdt\beta_Rdx-\beta_R^2dx^2}{1-\beta_R^2}
```
The cross terms cancel:
```{math}
:label: eqinvar2
\frac{dx^2+\beta_R^2c^2dt^2-c^2dt^2-\beta_R^2dx^2}{1-\beta_R^2}
```
regroup and pull out a factor of $-\beta_R^2$:
```{math}
:label: eqinvar3
\frac{dx^2-c^2dt^2-\beta_R^2(dx^2-c^2dt^2)}{1-\beta_R^2}
```
Pull out the common factor and the $1-\beta_R^2$ will cancel,
leaving $dx^2-c^2dt^2$, which is the same interval we started
with in the unprimed frame!  The interval does not change under
the Lorentz transformation.  This should of course be no surprise,
because that is the condition under which we **derived** the
elements of the matrix.  In particular, the interval in the
rest frame is therefore $-c^2dt_0^2$, which means whatever the
interval works out to be in any other frame, it must also
end up being equal to $-c^2dt_0^2$.

Hence, the size of the displacement 4-vector is an invariant under a
Lorentz Transformation.  Therefore, we say that it is a Lorentz
**scalar**.  It is interesting to note that the size of this invariant
quantity is negative! This is quite different from the traditional 3
vectors used in Newtonian mechanics.

## Muons in the Atmosphere

## Length Contraction

