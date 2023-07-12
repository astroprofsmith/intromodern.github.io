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
[p_4] = m_0[v_4] = m_0 \gamma \frac{[dx_4]}{dt}=
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
is different from the value 1 by just a few parts in $10^{12},$
hardly worth taking into account.

So, let us re-define the physical momentum of a particle as:
```{math}
:label: pdef
\vec{p}\equiv m_0 \gamma \vec{v}
```
which means Equation {eq}`p4def` becomes
```{math}
:label: p4full
\boxed{
[p_4] =
\begin{bmatrix}
icm_0\gamma\\
p_x\\
p_y\\
p_z
\end{bmatrix}
}
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

If a system is made up of more than one particle, the total
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

Example 8.2: Suppose that two observers in different reference frames
were measuring the momentum of the particle in Example 5.1 above.  Say
the primed observer is moving with $\beta_R=0.866$ in the $+x$
direction with respect to the unprimed observer.  The unprimed
observer would measure a momentum 4-vector of:
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

The size of the momentum 4-vector can be found by squaring each of the
components in Equation {eq}`p4comp` and adding them up.
```{math}
:label: p4size
[p_4]^2 = p_x^2+p_y^2+p_z^2-(\gamma m_0 c)^2 = p^2 - (\gamma m_0 c)^2
```
where we take $p$ to be the physical momentum (including the $\gamma$!).
Multiplying both sides by $c^2$ and doing a bit of algebra gives
```{math}
:label: eqcp1
\boxed{
(cp)^2 +(m_0 c^2)^2 = (\gamma m_0 c^2)^2}
```
The term on the right of Equation {eq}`eqcp1` can be identified
by making use of the correspondence principle. Expanding the $\gamma$
on the right using a Taylor expansion about $\beta = 0$ gives:
```{math}
:label: gammataylor
\gamma = (1-\beta^2)^{-1/2} \rightarrow 1 +\frac{1}{2}
\beta^2 + \frac{3}{8} \beta^4  + \frac{15}{48} \beta^6 + ...
```
If we assume that $\beta$ is small enough that we can ignore
all the higher order powers, then $\gamma \approx \beta^2/2$
and $\gamma m_0 c^2$ becomes approximately
```{math}
:label: totalEexp
\gamma m_0 c^2 \approx \left(1+\frac{1}{2}\beta^2\right) m_0 c^2
= m_0c^2 + \frac{1}{2}m_0 v^2
```
The last term on the right is recognizable as Newton's definition
of the kinetic energy of the particle. The dimensions of all the terms
in Equation {eq}`totalEexp` are energy.  As the velocity becomes exactly
equal to 0, the term on the right (and all the terms in the expansion,
of course) disappear, but there still is a term left in this 'energy'
equation.  The term $m_0 c^2$ can therefore be interpreted as an energy
associated with the particle when it is at rest.  This is of course the
rest energy of the particle, giving rise to that most famous of equations,
```{math}
:label: emc2
\boxed{
E_0 = m_0 c^2}
```
To be completely accurate, then, the kinetic energy should really
include all those higher order terms.  Not just $1/2mv^2$, but
a $v^4$ term and a $v^6$ term and a $v^8$ term and so on.  This
expansion means that the closer $v$ gets to $c$, the more the KE
must increase for small increases in $v$.  Or to put it another way,
the same increase in KE will yield a much smaller increase in $v$
if $v$ is already close to $c$.  That energy goes more and more into
the higher order terms rather than into increasing $v$.  The closer
you get to the speed of light, the more energy it takes to go any
faster.  This is another way of saying that nothing can move faster
than light.  It would take an infinite amount of energy.


```{note}
How did Newton miss the rest energy of a particle? The rest energy of
a person with a (rest) mass of 100 kg is $E =100~{\rm
kg}(3\times10^8~{\rm m/s})^2 = 9\times10^{18}$ joules.  This is an
unbelievably large number.  The total yearly electrical demand of the
entire US is on the order of $10^{19}$ joules (according to the
[EIA](https://www.eia.gov/electricity/annual/html/epa_01_02.html)).
This is about the same order of magnitude!!!  Two reasons Newton
didn't notice this massive amount of energy.  Most importantly, only a
change in energy is connected with anything happening.  Energy that
doesn't change has no measurable implications, so you would never know
it was there.  Secondly, a 100 kg person does not just vanish, in
Newton's time or now, so that $10^{19}$ joules of energy isn't
available to be associated with anything else you could measure.
Until people discovered that the rest energy of particles can change
(through processes like fusion or radioactivity), there was no reason
to notice an energy associated with rest mass, because it never
changed.  I even hear people today talk about a "Law of Conservation
of Mass," even though there is no such law.  Mass can and does change,
freeing up that rest energy to be associate with motion.  This energy
is the cause of making the Sun shine, nuclear power plants, or nuclear
bombs.  Mass is most definitely not conserved, but Isaac Newton had no
idea.
```

If the total energy $\gamma m_0 c^2$ is just the rest energy when the
particle is not moving,
the fact that it increases with the particle's speed ($\gamma$ gets bigger
than one) suggests that the remaining energy above and beyond the rest
energy would be the kinetic energy.  This conclusion is supported by
the fact that when $\gamma$ is small but not zero, $\gamma m_0 c^2$ is
approximately the rest energy plus the classical kinetic energy.
The $1/2 mv^2$ rule for kinetic energy that we all learned in introductory
physics is therefore only an approximation.  A full equation would be
```{math}
:label: etotal
\boxed{
E_{\rm total} =  \gamma m_0 c^2 = E_0 + KE }
```
We can use Equation {eq}`eqcp1` to say that
```{math}
:label: etotal2
\boxed{
E_{\rm total}^2 =  (\gamma m_0 c^2)^2 = E_0^2 + (cp)^2 }
```
or we could turn it around and solve for the momentum
of a particle given its energies:
```{math}
:label: pfrome1
cp = \sqrt{E_{\rm total}^2 - E_0^2}
```
Expanding the square using Equation {eq}`etotal`:
```{math}
:label: pfrome
cp = \sqrt{E_0^2+2E_0KE+KE^2 - E_0^2} = \sqrt{KE^2+2E_0KE}
```
You could pull the KE out of the square root, to get
```{math}
:label: pfrome2
\boxed{
p = \frac{KE}{c}\sqrt{1+2\frac{E_0}{KE}}}
```


In the world of fast moving objects, the momentum and energy of the
system are connected as shown in equation {eq}`etotal2`.  

```{note}
One of the subjects where the rules of relativity cannot be ignored
is in the realm of high energy particle physics.  In that regime, it
is often useful to use units of eV for energy (or keV or even MeV).
An eV is the amount the energy of an electron changes when moving
through a potential difference of one volt.  One volt is a joule per
coulomb, so 
because the charge on an electron is $1.602\times10^{-19}$ C, the
value of one eV in SI units is $1.602\times10^{-19}$ J.  Because
of Equation {eq}`etotal2`, we know that $cp$ has dimensions of energy,
so it is often more convenient when you know the energy in eV to
express the momentum in eV/c.  Then you don't have to worry about
converting everything to kg m/s.  See the example 8.3 in the margin,
below, for how easy eV units can be.
```


This means
that the momentum 4-vector can be written in terms of the physical
momentum and total energy of the system as:
```{math}
:label: p4fullE
\boxed{
[p_4] =
\begin{bmatrix}
iE_{\rm tot}/c\\
p_x\\
p_y\\
p_z
\end{bmatrix}
}
```
which means the time component of the momentum four vector is telling
you the total energy of the system.  Many people therefore call the
momentum four vector "the energy-momentum four vector".

```{warning}
Please make a careful note that Equations {eq}`etotal` and {eq}`etotal2`
are NOT the same equation!  Do NOT think that $cp$ is the kinetic energy.
Each term in Equation {eq}`etotal2` is **squared**, and the sum of squares
is **not** the same thing as the square of a sum!  The term $cp$ is clearly
related to the kinetic energy -- they have the same dimensions and they
both go up and down together -- but they are not the same thing.  You
can see that they are related but not equal in Equation {eq}`pfrome2`.
```

Equation {eq}`etotal` can be used to find the speed of a particle of
known rest and kinetic energy.  If the total energy is rest plus kinetic,
and the total energy is $\gamma$ times the rest energy, then
```{math}
:label: gamE
E_{\rm tot} = \gamma E_0 = E_0 + KE \rightarrow \boxed{\gamma = 1 + \frac{KE}{E_0}}
```
This suggests an easy test for deciding if you have to use the
relativistic formulae to describe the motion of an object.  If you
know the KE and the $E_0$, use Equation {eq}`gamE` to get $\gamma$.
If this number is significantly greater than one, then you need to use
the relativistic equations, but if it's very close to one, you can get
away with using Newtonian mechanics.

```{margin}
Example 8.3:

An electron gun in the lab can accelerate electrons through a potential
difference of up to 500 V.  In energy terms, this means the electron will
gain 500 eV of kinetic energy.  The rest energy of an electron is about
500 keV.  Therefore, the Lorentz factor for the electrons in this device
will never be greater than $\gamma = 1 + 5\times10^2/5\times10^5 = 1.001$,
and it is not important to use relativistic equations to describe the
motion of the electrons in this apparatus.

At $\gamma=1.001$, we can figure out the speed, and we get $\beta
= \sqrt{1-1/\gamma^2} = 0.0447$, which works out to $1.3\times10^7$ m/s.
That's still something like 25 million miles an hour, which shows you
just how fast "close to the speed of light" really is.  Even 25 million
miles an hour isn't fast enough to really need to use relativity theory.
```

## The Rest Energy

The rest energy is the energy measured in the rest frame of a particle.
It can be determined from the rest mass using Equation {eq}`emc2`.
The elementary particles have fixed rest energies.  For the proton --
938.26 MeV, the neutron -- 939.55 MeV, the electron -- 0.511 Mev, etc.
However, not all particles have rest energy. A photon (gamma ray) does
not have any rest energy, and a neutrino's is so small it can barely be
measured.

Equation {eq}`gamE` would therefore imply that the value for the $\gamma$
of these particles is infinite. If $\gamma$ is infinite, then $\beta= 1$,
and these particles must be traveling at the speed of light.  Always.
This is not surprising for the gamma ray, as it is light. What is
surprising is that other particles such as the neutrino, which are
not photons, also travel at the same speed.

Also surprising is that all of these particles have momentum even
though they do not have any rest mass.  You might think, if you think
of momentum as mass times velocity, that no mass means no momentum.
However, if you look at Equation {eq}`etotal2`, the only way this
equation can work for particles that have energy but no mass is if
$p=E/c$.  Particles with energy but no mass **still have momentum**!!!

```{note}
You can also see this in Equation {eq}`pfrome2`: if the rest energy
goes to zero, $p\rightarrow KE/c$.  For a photon, the KE is $hf$, Planck's
constant times the frequency, so the momentum of a photon is $hf/c$,
or $h/\lambda$, where $\lambda$ is the wavelength.  See also the
treatment in the section of Chapter 9 on the Doppler Shift.
```

On the other hand, for a particle with $E_0 \neq 0$ to travel **at**
the speed of light, Equation {eq}`gamE` implies it must have
infinite kinetic energy 
because $\gamma$ is infinite for an object traveling the
speed of light. No matter how much energy is given to this massive particle,
it will always be traveling slower than the speed of light.
Consequently, objects with $E_0 \neq 0$ travel slower than the speed
of light and objects with $E\approx 0$ travel at the speed of light.
Nothing that travels at or lower than the speed of light can ever travel
faster than the speed of light because it would
require an infinite amount of energy, which is not readily available
to us, to make the transition from slightly less than the speed of
light to slightly more than the speed of light.  As the saying goes...


**The speed of light. It's not only a good idea, it's the law!**

## Conservation Laws and the Momentum Four-Vector

If an interaction happens that is isolated from its surroundings, the
total energy $E_{\rm tot}$ (as measured by a particular inertial
observer) of that system will be constant before, during, and after
the interaction.  This is the same $E$ that appears in the time
component of the total momentum 4-vector of this system. If the
interaction experienced by each of the pieces of the system is only
due to interactions with other pieces of the system, then the net
physical momentum (as measured by a particular inertial observer) for
this system is constant. The physical momentum is a 3-vector, so both
its size and its direction must remain constant. These laws of physics
hold for any particular inertial observer.


When a certain physical parameter remains constant while time changes,
then that parameter is said to be conserved. Under the conditions
stated in the previous paragraph, the sum of spatial components of
the momentum 4-vector, when added as 3-vectors, is conserved. In a
similar manner, the time component of the momentum 4-vector is also
conserved. For a particular inertial observer, the momentum 4-vector
is conserved.

If a particular inertial observer measures the total energy of the
system, this value will remain constant as long as the system observed
is isolated. A second inertial observer, traveling with some relative
velocity $\beta_R$ with respect to the first, will also measure a
constant total energy, but the value that that observer measures will,
in general, be different from the value measured by the first
observer. The Lorentz transformation shows how the two different total
energy values are related. What remains constant between two two
reference frames is the size of the momentum 4-vector, as it is a
Lorentz invariant.

We therefore try to use the words consistently: to be "conserved" is
to remain constant within a particular reference frame, during a
time interval, but to be "invariant" is to be the same in different
reference frames moving with respect to each other.  Shifting from
one reference frame to another is not a process that conserves energy.
An object at rest in one frame will have only its rest energy,
but in another frame it will have rest and kinetic energy.  Shifting
between these frames does **not** do work upon the object.  If I
change my walking speed relative to you, your speed changes relative
to me, but I have not exerted a force on you to change your speed.
The kinds of analysis that depend on conservation of energy or the
work-energy theorem are meant to be carried out in a single frame
of reference.  If you switch to another reference frame, you must
recalculate the energy number somehow -- it will not in general be
the same.

## Example of Momentum Conservation

Example 8.4: A neutral pion (rest energy $E_0=135$ MeV) decays into
two photons. Show that in the rest frame of the pion, the two photons
move in opposite directions. Find the energy of these two photons.

In the rest frame of the pion, the momentum 4-vector for the pion is:
```{math}
:label: p4pi
[p_4]_{\rm before} =
\begin{bmatrix}
i135~{\rm MeV}/c\\
0.0\\
0.0\\
0.0
\end{bmatrix}
```
The observer in the rest frame of the pion watches the decay take place.
The $x$-direction is defined as the direction that one of the dautgher photons
takes. The second photon travels in some unknown direction. The observer in
the rest frame of the pion measures a momentum 4-vector of the two photons as:
```{math}
:label: p42pho
[p_4]_{\rm after} =
\begin{bmatrix}
i(E_1+E_2)/c\\
p_{x1}+p_{x2}\\
p_{y2}\\
p_{z2}
\end{bmatrix}
```
For two 4-vectors to be equal, each of their equivalent components have to
be equal.  Therefore the second photon cannot have any momentum in the
$y$ and $z$ directions, and the momentum of the second photon in the
$x$ direction must be equal and opposite to that of the first photon,
so that the momentum adds up to equal the zero it was before the collision.
This shows that the photons, whichever direction they go, must go opposite
to each other.

```{note}
The ease with which we can set things to zero and cancel them makes this
particular frame of reference, in which the momentum was and remains
zero, particularly useful.  It is referred to as the **center of momentum
frame**, and we will use it in a detailed example in the next chapter.
```

We can examine the energy of the two photons by looking at the time
component.  Since the four momentum is conserved, the time component
must also be conserved, and therefore the sum of the two photon
energies must be 135 MeV/c.  Since the photons have no mass, their
momenta and their energies must be proportional (by Equation
{eq}`eqcp1`), so the energy must be split evenly between them.  Each
photon has an energy of 67.5 MeV, and therefore a momentum of 67.5
MeV/c, but each in an opposite direction to the other.