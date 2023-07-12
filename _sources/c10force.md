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

# Chapter 10: Relativistic Force and Acceleration

## Background


Dynamics is the study of how motion of an object changes with time. In
the late 1600's Issac Newton proposed a model for predicting the
changes in motion of a slow moving (in 21st century terms)
object. From observed changes in motion, a physicist could deduce the
interactions that the object experienced. Each interaction is
characterized by a vector, called a force, that is determined by
where the object is, how fast it is traveling, and values of
quantities such as mass, electrical charge, etc.

Newton suggested that if one calculates the vector sum of
all the forces acting on the object, then this net force
```{math}
:label: Newton2nd
\vec{F}_{\rm net} = \Sigma_n \vec{F}_n =
\frac{d\vec{p}}{dt} = \frac{d(m\vec{v})}{dt}
```

where $m$ is the mass of the slow moving object, $\vec{v}$ is the
velocity relative to the observer of the object and $dt$ is the
(infinitesimal) time interval. While this model worked well for the
speeds that were known in the 1700s, its application to objects moving
near the speed of light has some problems. One problem is with the
$dt$ term. As we know now, $dt$ depends on the reference frame of the
observer and is therefore not actally a scalar, as Newton would have
understood it to be. If the velocity is a vector, and $dt$ is not a
scalar, then the net force is not a vector.

Second, there is nothing inherent in Equation {eq}`Newton2nd` to
inhibit a net force from changing the objects speed to an arbitrarily
fast value.  This is no problem for Newton, who never had to deal with
very large speeds, but we now know that nothing can move faster than
$c$, the speed limit of the universe.  Equation {eq}`Newton2nd` has
no built in speed limit.

To fix these problems, we need to cast Equation {eq}`Newton2nd` in
4vector form.  Then we can predict the acceleration of objects near the
speed of light and characterize their motion.

## Newton's Second Law with Four Vectors

To fix Equation {eq}`Newton2nd`, we keep the underlying meaning
of the equation but recast it in relativistic terms, as we did to
get four velocity from four displacement.  Instead of $d\vec{p}/dt$,
we write
```{math}
:label: minkforce
[K_4] = \frac{d}{dt_0}[p_4]
```
The four force is the derivative (with respect to the proper time, a
scalar) of the four momentum.  Since $[p_4]$ is a proper 4- vector and
$dt_O$ is a scalar, then $[K4]$ is also a proper 4-vector. The letter
K is in the 4-force to make it easier to distinguish between it and
the vector force used by Newton. This form of the force is often
called the Minkowski Force in honor of H. Minkowski who introduced the
4-vector notation in the early 1900s.

Using the definition of the momentum 4-vector in terms of the
Newtonian velocity $\vec{v}$, with components $v_x$, $v_y$, and $v_z$
gives:
```{math}
:label: fourforce
[K_4] = \frac{d(m_0[v_4])}{dt_0} =
\begin{bmatrix}
i \frac{1}{c} dE/dt_0\\
d\vec{p}/dt_0
\end{bmatrix}
=
\begin{bmatrix}
i d(\gamma m_0c)/dt_0\\
d(\gamma m_0\vec{v})/dt_0
\end{bmatrix}
```
where we again re-define the velocity and momentum vectors to include
the relativistic $\gamma$ term that Newton missed because he never
measured objects moving near the speed of light.

## The Minkowski Force

The Newtonian force vector is related to the momentum by Equation
{eq}`Newton2nd` In this equation, $dt$ is the time as measured by an
observer who sees the object moving with velocity $\vec{v}$, not the
observer in the rest frame of the object, If we convert the times in
Equation {eq}`fourforce` from the rest frame to the moving reference
frame using the time dilation Equation {eq}`timedilation`, Equation
{eq}`fourforce` becomes:
```{math}
:label: M4force
[K_4] = 
\begin{bmatrix}
i \frac{\gamma}{c} dE/dt\\
\gamma d\vec{p}/dt
\end{bmatrix}
=
\begin{bmatrix}
K_0\\
\gamma \vec{F}
\end{bmatrix}
```
The three spatial terms of the Minkowski force are the
components of the Newtonian force, multiplied by a factor
of $\gamma$.

The time component needs a little more careful attention to interpret
it properly.  Clearly, it will be related to the rate of change of the
energy of the object experiencing the force.  This rate of change is
called the power delivered to the object.  How is this related to the
force?  In Newtonian terms, the work done on an object by a force as
the object displaces by some amount $d\vec{x}$ is given by $dW =
\vec{F}\cdot d\vec{x}$.  The power is therefore $P = dW/dt =
\vec{F}\cdot \vec{v}$.  Let us therefore see what happens when we take
the dot product of the four velocity with the four force:
```{math}
:label: v4dotF41
[v_4]\cdot[K_4] = [v_4]\cdot \frac{d(m_0[v_4])}{dt_0}
```
However, we can take advantage of the fact that $d(x^2)=2xdx$
to write this as
```{math}
:label: v4dotF42
[v_4]\cdot[K_4] = \frac{d}{dt_0}\left(\frac{m_0}{2}[v_4]\cdot[v_4]\right)
```
But we know that $[v_4]^2$ is the Lorentz invariant $-c^2$,
which is also a constant in time.  So
```{math}
:label: v4dotF43
[v_4]\cdot[K_4] = \frac{d}{dt_0}\left(-\frac{m_0c^2}{2}\right) = 0
```
The result of the four dimensional dot product must therefore be
zero, so let's write it out in terms of components:
```{math}
:label: v4dotF4comp
[v_4]\cdot[K_4] = \gamma^2(v_xF_x + v_yF_y +v_zF_z) + i\gamma c K_0 = 0
```
Therefore
```{math}
:label: findingK0
\gamma^2(\vec{v}\cdot\vec{F}) = -i\gamma c K_0
```
and
```{math}
:label: Ktime
\boxed{
K_0 = i\gamma(\vec{F}\cdot\vec{\beta})
}
```
So the Minkowski force becomes
```{math}
:label: M4forcefin
[K_4] = 
\begin{bmatrix}
i\gamma (\vec{F}\cdot\vec{\beta})\\
\gamma \vec{F}
\end{bmatrix}
```
Again, we know from Newtonian mechanics that the rate of
change of the energy of an object is $\vec{F}\cdot\vec{v}$,
and if we compare the forms of $K_0$ in Equations {eq}`fourforce`
and {eq}`M4forcefin`, we see that
```{math}
:label: Kpower
K_0 = i\gamma(\vec{F}\cdot\vec{\beta}) = i\frac{\gamma}{c}\frac{dE}{dt}
= i\gamma\frac{d}{dt}(\gamma m_0 c)
```
so $\vec{F}\cdot\vec{\beta} = d(\gamma m_0 c)/dt$, as you might
expect, and the time component of the Minkowski force is telling
us the power being delivered to the object experiencing the force.

```{note}
It's interesting to note that although $\vec{F}\cdot\vec{v}$ is
power, and has dimensions of energy per time, $\vec{F}\cdot\vec{\beta}$
has dimensions of force, and yet it is **telling** us the power,
because all you need to do is multiply $K_0$ by $c/\gamma$ and you
get the power.
```



### Example 7.1

Calculate the Minkowski force exerted on an electron traveling
vertically upward near the surface of the earth with $\beta = 0.866$.

Solution: Let the $z$ axis be vertical with the positive value
pointing up. Near the surface of the earth, the electron experiences a
gravitational force $\vec{F}_g = m\vec{g} = -m\times9.81~{\rm
N/kg}~\hat{z}$.  The mass of the electron has a value of
$(9.10938188\pm0.00000072)\times 10^{-31}$ kg.  The force exerted by
gravity on the electron is therefore:
```{math}
:label: Fgemin
\vec{F}_g = (9.11\times10^{-31}~{\rm kg})\times(9.81~{\rm N/kg})(-\hat{z})
= 8.93\times10^{-30}~{\rm N}(-\hat{z})
```
```{margin}
Always good to remember that for $\beta=0.866$, $\gamma=2$
```
Since the force is pointing only in the $-\hat{z}$ direction, the
$F_x$ and $F_y$ comonents are zero. Using Equation {eq}`M4forcefin` we
now can find the spatial components of $[K_4]$ by calculating the
value of $\gamma$ when $\beta=0.866$.  $K_x$ and $K_y$ are therefore
zero, but $K_z = \gamma m g$, or $-1.78\times10^{-29}~{\rm N}$.  Since
the motion is in the $+\hat{z}$ direction, the time component becomes:
```{math}
:label: Fgtime
K_0 = i\gamma F_g(-\hat{z}) \cdot (\beta \hat{z}) =
-i (2.00 \times 8.93\times10^{-30}~{\rm N} \times 0.866) =
- i 1.54\times10^{-29}~{\rm N}
```
Combining these results gives the Minkowski force:
```{math}
:label: Fgmink
[K_4] =
\begin{bmatrix}
-i 1.54\times10^{-29}~{\rm N}\\
0\\
0\\
-1.78\times10^{-29}~{\rm N}
\end{bmatrix}
```
Remember that the actual force experienced by the electron will be
$K_3/\gamma$, or half the number given in Equation {eq}`Fgmink`,
because $\gamma=2$ at this speed.  The time component still has units
of newtons, because $K_0$ has a factor of $1/c$ in it.  If you want
to know the power delivered to the electron, you have to take the
$K_0$ term as given in Equation {eq}`Fgmink` (without the $i$, of
course) and multiply it by $c/\gamma$.  This would imply the electron
is losing energy at a rate of $2.31\times10^{-21}$ watts.

## The Acceleration Four Vector

Although Special Relativity demands that any inertial reference frames
we consider be moving at constant relative velocities, there is nothing
to prevent objects within a particular reference from from accelerating.
You do have to be careful about switching frames, though, because you
can't use the SR rules to switch into a rest frame for an accelerating
object!

The relativistic acceleration, a well defined but strange quantity, can
be modeled in two ways. First, it would be the proper time derivative
of the 4-velocity:
```{math}
:label: acc4dv
[a_4] = \frac{d([v_4])}{dt_0} = \gamma \frac{d([v_4])}{dt}
```
or one could use the relativistic analog of Newton's law for a
particle of constant rest energy (mass):
```{math}
:label: acc4ma
[a_4] = [K_4]/m_0
```


Consider an observer in the lab who sees a particle of mass $m_0$,
traveling in the $+x$ direction with some speed $v$.  A force is then
applied to the particle in some arbitrary direction, causing the
observer to see a change in the particle's velocity in the $x$, $y$
and $z$ directions ($dv_x$, $dv_y$, and $dv_z$). This observer will
measure an acceleration $\vec{a} = d\vec{v}/dt$.  We can write down
the four-acceleration as observed in the lab as:
```{math}
:label: acc4lab
[a_4] =
\gamma
\begin{bmatrix}
i\frac{d(\gamma c)}{dt}\\
\frac{d(\gamma v_x)}{dt}\\
\frac{d(\gamma v_y)}{dt}\\
\frac{d(\gamma v_z)}{dt}
\end{bmatrix}
```
Remember that $\gamma\equiv1/\sqrt{1-\beta^2}$ has all three
components of velocity in it: $\beta^2\equiv \beta_x^2 + \beta_y^2
+\beta_z^2$.  Since both $\gamma$ and the velocity are changing
with time, we must use the product rule of calculus to evaluate
Equation {eq}`acc4lab`:
```{math}
:label: acc4prodrule
[a_4] =
\gamma
\begin{bmatrix}
ic\frac{d\gamma}{dt}\\
v_x\frac{d\gamma}{dt}\\
v_y\frac{d\gamma}{dt}\\
v_z\frac{d\gamma}{dt}
\end{bmatrix}
+
\gamma^2
\begin{bmatrix}
i\frac{dc}{dt}\\
\frac{dv_x}{dt}\\
\frac{dv_y}{dt}\\
\frac{dv_z}{dt}
\end{bmatrix}
```
Taking the derivative of $\gamma$ works like this:
```{math}
:label: dgamma
d\gamma = d\left((1-\beta^2)^{-1/2}\right) = -\frac{1}{2}\left(1-\beta^2\right)^{-3/2}(-2\beta)d\beta = \gamma^3\beta d\beta
```
This is a complicated relationship that depends on the three components
of the velocity as well as $dv_x$, $dv_y$, and $dv_z$.

Placing this result back into Equation {eq}`acc4prodrule` and using
the fact that $dc/dt=0$ gives us:
```{math}
:label: acc4dgam
[a_4] =
\gamma^2
\begin{bmatrix}
i\gamma^2\beta\frac{d\beta}{dt}c\\
\gamma^2v_x\beta\frac{d\beta}{dt}+\frac{dv_x}{dt}\\
\gamma^2v_y\beta\frac{d\beta}{dt}+\frac{dv_y}{dt}\\
\gamma^2v_z\beta\frac{d\beta}{dt}+\frac{dv_z}{dt}
\end{bmatrix}
```
Equation {eq}`acc4dgam` is a very complex equation, not least because
each component of $[a_4]$ depends on all three components of both the
velocity (through $\beta$) and the acceleration measured by the
observer in the lab.  This is very new to us at this point: that the
acceleration in $x$ should depend on both the velocity and
acceleration in $y$ and $z$.  We are used to the spatial dimensions
being independent of each other.  That is apparently no longer the
case.  Of course, the correspondance principle still applies, and you
can see that if $\beta\rightarrow 0$, the first terms all drop out,
$\gamma\rightarrow 1$ and we get back acceleration is the derivative
of the velocity, with independent dimensions, just as we expect in
Newtonian mechanics.

It is hard, perhaps impossible, to get any kind of intuition out of
Equation {eq}`acc4dgam`, but we can look at a special case to see what
happens.  Consider the case where a constant force is being exerted on
a particle in the $x$ direction only ($\vec{F} = F_x\hat{x}$), and the
particle starts from rest.  Then the velocity will only have an $x$
component, so we can just use $\vec{v} = v\hat{x}$ and drop the last
two components of the four vector, since there will never be any
motion in the $y$ or $z$ directions.  Then Equation {eq}`acc4dgam`
simplifies a bit:
```{math}
:label: acc4Fx
[a_4] =
\gamma^2
\begin{bmatrix}
i\gamma^2\frac{v}{c}\frac{dv}{dt}\\
\gamma^2\left(\frac{v}{c}\right)^2\frac{dv}{dt}+\frac{dv}{dt}\\
\end{bmatrix}
```
The relativistic form for Newton's Second law (Equation {eq}`acc4ma`)
then becomes
```{math}
:label: newton2x
m_0[a_4] = [K_4] \rightarrow
m_0\frac{dv}{dt}\gamma^4
\begin{bmatrix}
i\frac{v}{c}\\
1
\end{bmatrix}
=
\gamma F_x
\begin{bmatrix}
i\frac{v}{c}\\
1
\end{bmatrix}
```
Since the four-vectors are the same, we can compare the factors out
front to find that
```{math}
:label: relacc
\frac{dv}{dt} = \frac{1}{\gamma^3}\frac{F_x}{m_0}
```
The corrspondance principle check is satisfied because for
$\gamma\approx 1$, this is just $a=F/m$, as Newton would expect.
However, as $v\rightarrow c$, $\gamma \rightarrow \infty$, which means
the acceleration is going to plummet to zero as the speed gets close
to $c$, regardless of the fact that the force has not changed.  Yet,
again, the universe conspires to make sure that nothing can go
faster than $c$.  The closer you get to $c$, the harder it is to
get any further acceleration.  This limit is displayed in Figure 10.1 --
two initial speeds are shown to indicate that it doesn't matter how
fast the object is moving to begin with.  Even under a constant force,
it will never move faster than the universal speed limit.

```{code-cell}
:tags: ["remove-input"]
Foverm = 9.81
c = 2.995e8
v0b = 0.8*c
v0 = 0.0

tmax = 200000000
dt = 10000.0
t = 0.0
nn = int(tmax/dt)
tarr = np.zeros((nn,1))
varr0 = np.zeros((nn,1))
varr1 = np.zeros((nn,1))
varr0[0] = v0
varr1[0] = v0b
i = 0

while i<nn-1:
  gamma = 1.0/np.sqrt(1.0-varr0[i]**2/c**2)
  a = Foverm/gamma**3
  varr0[i+1] = varr0[i]+a*dt
  gamma = 1.0/np.sqrt(1.0-varr1[i]**2/c**2)
  a = Foverm/gamma**3
  varr1[i+1] = varr1[i]+a*dt
  tarr[i+1] = tarr[i]+dt
  i=i+1

tarr = tarr/3.15e7
plt.figure(figsize=(14,6))
plt.plot(tarr,varr0/c,'b-',label='Start from 0.0')
plt.plot(tarr,varr1/c,'m-',label='Start from 0.8')
plt.plot([tarr[0],tarr[-1]],[1.0,1.0],'r:')
plt.xlabel('Time (years)')
plt.ylabel('Beta')
plt.legend()
plt.show()
```
```{note}
Figure 10.1 -- A graph of $\beta$ as a function of time (in years) in
the lab frame for an object of mass 1 kg, experiencing a constant
force of $9.81$ N in the $x$ direction.  The blue line shows the
object starting from rest, while the magneta line shows the object
starting from $\beta=0.8$.  In neither case does the speed increase
beyond 1, which is indicated by the red horizontal dotted line.
```

Figure 10.1 also shows that if you could push an object with a force
equal to its weight on the surface of the Earth, within a few years,
the speed of that object would approach the speed of light, even if it
started out from rest.  However, in practical terms, maintaining even
such a modest thrust for years without fair is hard to do without
running out of fuel.  But that's another topic for another day.

The general case of applying the relativistic Newton's law to a
particle having laboratory accelerations and velocities in all three
directions, as well as the force having all three components is very
complicated, as $\beta$ and $\gamma$ have all three components of the
velocity in them. In general, the force and the acceleration will not
be in the same direction, something that does not happen in classical
physics.

