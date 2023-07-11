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

However, Figure 9.1 shows how much energy is available in the COM
frame, according to Equation {eq}`Etotcom` (after subtracting off the
rest energies of the original two particles) as a function of how much
kinetic energy the incident particle has.  You can see that the amount
of available energy is significantly less than the amount of KE put in
(the black dotted line represent equality).  It appears that some of
the energy in the lab frame goes into moving the center of momentum of
the system.

```{code-cell}
:tags: ["remove-input"]
E0 = 938
KE = np.arange(0,100,1)*1000.0
beta = np.sqrt(KE/(KE+2.0*E0))
gam = 1.0/np.sqrt(1-beta**2)
Elab = 2.0*E0+KE
pi = np.sqrt(KE**2+2.0*E0*KE)
Ecom = gam*(Elab-beta*pi)
plt.figure(figsize=(9,5))
plt.plot(KE,Ecom-2.0*E0,'b-',label='E in COM')
plt.plot(KE,KE,'k:',label = 'Eout = Ein')
plt.xlabel('Kinetic Energy of Incident Particle (MeV)')
plt.ylabel('Available Energy in COM Frame (MeV)')
plt.axis([0,1.0e5,0,1.25e4])
plt.legend()
plt.show()
```
```{note}
Figure 9.1 -- Graph of energy available in the COM frame
(after subtracting the rest energy of the two particles involved)
as a function of the KE of the incident particle in the lab
frame.  There is significantly less energy available in the
COM frame to make particles than in the lab frame, if you
put all the KE into accelerating one particle at the stationary
target.
```

Need to add Examples 1 and 2.

## Photons and the Doppler Shift


The Doppler shift tells what happens to the measured energy of a
photon when it is viewed by different inertial observers. The velocity
of the photon is the same for both inertial observers, but the energy
changes. The momentum 4-vector for a photon as seen by the first
observer is
```{math}
:label: p4phot
[p_4]_1 =
\begin{bmatrix}
iE_1/c\\
E_1/c\\
0\\
0
\end{bmatrix}
=
\begin{bmatrix}
ih/\lambda_1\\
h/\lambda_1\\
0\\
0
\end{bmatrix}
```
where we use $E=hf$ and $c=\lambda f$ to get the latter form.

Notice that the size of the momentum 4-vector for a photon (or any
zero rest energy particle) is 0. That is because the size of the
momentum 4-vector = $-m_0c^2$, and the particle does not have any rest
energy.

A second inertial observer moving with $\beta_R$ in the $+x$ direction
with respect to the first inertial observer, watching this same photon,
would measure a momentum 4-vector given by (again ignoring $y$ and
$z$):
```{math}
:label: p4photlort
[p_4]_2 =
\begin{bmatrix}
iE_2/c\\
E_2/c
\end{bmatrix}
=
\begin{bmatrix}
\gamma_R & -i\beta_R\gamma_R \\
i\beta_R\gamma_R & \gamma_R \\
\end{bmatrix}
\begin{bmatrix}
iE_1/c\\
E_1/c
\end{bmatrix}
=
\gamma_R
\begin{bmatrix}
i(E_1/c-\beta_RE_1/c)\\
E_1/c-\beta_RE_1/c
\end{bmatrix}
```
So
```{math}
:label: p4phot2
[p_4]_2 = 
\begin{bmatrix}
iE_2/c\\
E_2/c
\end{bmatrix}
=
E_2/c
\begin{bmatrix}
i\\
1
\end{bmatrix}
=
\gamma_RE_1(1-\beta_R)/c
\begin{bmatrix}
i\\
1
\end{bmatrix}
```
The part in brackets is clearly the same, so $E_2=\gamma_R(1-\beta_R)E_1$.
Plug in the definition of $\gamma_R$ to get
```{math}
:label: Edopp
E_2 = E_1 \frac{1-\beta_R}{\sqrt{1-\beta^2}}
=
 E_1 \frac{1-\beta_R}{\sqrt{(1-\beta)(1+\beta)}}
= E_1\sqrt{\frac{1-\beta}{1+\beta}}
```
Since $\beta\leq 1$, the numerator in Equation {eq}`Edopp` will
always be less than the denominator, so as long as the second observer
is moving to the right (receding from the source of the light),
she will measure a smaller energy than the first observer would.
If the second observer switches direction and heads toward the source
of light, the $\beta$ will change sign, which effectively flips the
ratio under the square root.  In that case, $E_2\geq E_1$ -- the
energy goes up.

This is qualitatively not different from the case of a massive particle:
if you throw a ball at me, and I run toward you, the ball will have more
energy in my reference frame.  If I run away, the ball will have less
energy.  However, for a massive object like a ball or a proton, that
energy is associated with speed.  For a massless particle like a photon,
the speed does not change -- it always moves at $c$.  The energy in that
case is associated with the frequency (or equivalently, the wavelenth)
of the light.  So, for light, we often write Equation {eq}`Edopp` in
terms of frequency ($f$) or wavelength ($\lambda$):
```{math}
:label: Fdopp
f_2 = f_1\sqrt{\frac{1-\beta}{1+\beta}}
```
and
```{math}
:label: Ldopp
\lambda_2 = \lambda_1\sqrt{\frac{1+\beta}{1-\beta}}
```
Since $\lambda = c/f$, these two equations are reciprocals.

Since Equation {eq}`Fdopp` says the photons will be shifted toward
smaller frequencies if the relative velocity of the source and
observer is positive (receding), and the smallest frequencies of
visible light are red, we refer to this kind of change as a "redshift",
and if the source and receiver are approaching (negative relative
velocity), we call this change a "blueshift".  This can be confusing
if you are considering photons outside the visible range of light.
An x-ray, for example, subject to a blueshift, changes to an even higher
frequency, **away** from blue light on the electromagnetic spectrum.
So if the "blue" and "red" confuse you, just think of going to
higher or lower frequency (which means going the opposite direction
in wavelength).

```{warning}
It doesn't make sense to ask whether the source or the receiver is
"really" moving, because all motion is relative.  For a sound wave,
the medium it moves through provides a context in which it might
make sense to say one or the other is "really" moving (relative to the
medium), but in Chapter 1 we showed that the Michaelson Morely experiment
disproves the hypothesis that light has a medium.  So any motion
of source or receiver is equivalent and indistinguishable.
```




## Compton Scattering
