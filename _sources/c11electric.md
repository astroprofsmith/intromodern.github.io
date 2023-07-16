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

# Chapter 11: The Electromagnetic Tensor

Many introductory relativity textbooks stop with dynamics, but I think
it's important to include some of the implications of realtivity
theory for the understanding of electricity and magnetism.  There are
two reasons why it's important to see how SR affects E&M: first, for
mechanics, the effects of relativity are completely negligible in our
world of physics professors and dump trucks.  The weird effects only show
up at speeds near the speed of light, and for most day-to-day processes
you can blissfully ignore them.

Within the context of electricity, however, relativistic effects
occur, and are important, at *any* speed.  Even for the literal
snail's pace of the average drift speed of a 1 A current in a copper
wire, the differing observations in reference frames in relative
motion are obvious and unavoidable.

Second, the theory of electrodynamics that culminates in Maxwell's
Equations is automatically consistent with SR.  Newtonian mechanics
needs massive modifications to work with SR, but the bulk of
relativistic E&M is learning how to write it in four-vector notation
and remarking, "Wow!  It works!"  For a full treatment of this topic,
I recommend an advanced electrodynamics book such as Griffiths.  In
this eBook, I only plan to introduce the subject, not explore all the
implications.

## Magnetic Force as a Relativistic Effect

### Reminder of Classical Situation

```{figure} images/classicIBF.png
:alt: classicalEMBforce
:class: bg-primary mb-1
:width: 800px
:align: center
:name: classicfig

Schematic of the classical understanding of the force on a charge
moving antiparallel to a current-carrying wire.  A description of
all the pieces of this figure is in the text.

```

Before we dive into tackling how relativity implies we must rethink
E&M, let's first walk through a classical understanding of a simple
problem.  Then we will recast that problem through a relativistic lens
and see how it invites a different understanding of the phenomenon.
We will start with the force on an electric charge moving parallel to
a current-carrying wire.  The situation is shown in {numref}`classicfig`
as a schematic.  A wire carries a current $I$ to the right.  Meanwhile
a charge $q$ is moving to the right with a velocity $\vec{u}$ at a
distance $r$ from the center of the wire.  We observe, empirically, that
the charge experiences an upward force (toward the wire) with a
magnitude of
```{math}
:label: Fclassic
\vec{F} = -\hat{r} \frac{qu\mu_0 I}{2\pi r}
```
and this force remains perpendicular to the motion of the particle,
no matter how we rotate the system around the wire.

We therefore invent a field we  call $\vec{B}$ and we say there's
a magnetic force $\vec{F}_m = q\vec{u}\times\vec{B}$.  To make that
force line up with the observed force, we say the field $\vec{B}$
points tangent to a circle around the wire and has a magnitude of
$\mu_0 I/2\pi r$.


### The Relativistic Model

As with most relativity problems, this one comes down to setting up
the situation in different reference frames and then switching between
them with Lorentz transformations.  There are three frames we will
need to compare, because there are two different kinds of motion (the
current and the charge), and we need to compare them both to a common
reference frame.  Or to put it another way, we have the initial frame
in which the charge has a velocity $\vec{u}$, we have the rest frame
of the charge, and then we will need the rest frame of the *current*,
as well.

Let us begin with the "lab frame", which I will call "Frame A", which
corresponds to the situation described in the previous section.  There
is a current $I$ to the right, but I am going to break this down into
two parts: a positive current going right and a negative current going
left.  Each current comprises half the charge density, but they are
equal and opposite densities, so the total charge density is zero, but
the total current is $I$ to the right.  This way of conceptualizing the
setup is shown in {numref}`frameAfig`, with the positive current in
red and the negative current in blue.

```{figure} images/frameAfig.png
:alt: frameApicture
:class: bg-primary mb-1
:width: 800px
:align: center
:name: frameAfig

Schematic of a particular way of visualizing the current $I$ in a
reference frame where the total charge density of the wire is zero.
We define the current as $2\lambda v$, but there is a $+\lambda$ moving
to the right (red) and a $-\lambda$ moving to the left (blue).  The sum
of the two $\lambda$ is zero, but since the $v$ values are also opposite,
$\lambda v$ adds up instead of cancelling.

```

We then need to transform into a frame moving with $v_R=u$ to the
right.  In this frame, which I will call Frame B, the charge $+q$ is
at rest.  According to the classic magnetic force model (Equation
{eq}`Fclassic`), the magnetic force on the charge should be zero.  It
is not moving, so it cannot experience a magnetic force.  However,
since an observer in Frame A observes the charge being attracted to
the wire, an observer in Frame B must also observe the charge being
attracted to the wire.  But how, if there is no force?  No *magnetic*
force!  There could be a *different* force.

```{note}
This is actually the exact sort of dilemma that led Einstein to develop
SR in the first place.  How could you have a magnetic force in one
reference frame, and an electric force in another reference frame?
Shouldn't the laws of Physics be the same in all reference frames?
His 1905 paper that introduced SR to the world was not called "Special
Relativity" after all.  The name of that paper was "On the Electrodynamics
of Moving Bodies".
```

In Frame B, lengths along the direction of relative motion will be
contracted.  When we talk about a charge density $\lambda$, we are
talking about a certain number of charges *per* *length*.  So the
charge densities in the wire will change, and since these sets of
charges have different velocities, relative to $\vec{v}_R$, they will
contract by different amounts.  Given our formula for length
contraction, we know
```{math}
:label: lambdapm
\lambda_\pm =\pm\frac{dq}{dz} =
\pm\gamma_\pm\frac{dq}{dz_0} = \pm\gamma_\pm\lambda_0
```
In this frame, the two charge densities will **not**
cancel out, since they will no longer be equal.  The wire will **not**
be electrically neutral.
```{math}
:label: labdatot
\lambda_{\rm totB} = \lambda_+-\lambda_- = \lambda_0(\gamma_+-\gamma_-)
```


To be more specific, instead of talking about just one $\lambda$, we
must now talk about a $\lambda_+$ and a $\lambda_-$, since they will
be different.  (Note that $\lambda_\pm$ is explicitly defined in Frame
B).  To know how much they are contracted, we have to know their
speeds.  We have a formula for that:
```{math}
:label: vplusminus
v_\pm = \frac{v\mp u}{1\mp uv}
```
where $v$ (the speed of the current) and $u$ (the relative velocity
between A and B) are measured in Frame A.  On the other hand,
```{math}
:label: gammapm
\gamma_\pm = \frac{1}{\sqrt{1-v_{\pm}^2/c^2}}
```
where $v_\pm$ is measured in Frame B.

But what are $v_\pm$ measured with respect *to*?  These velocities
represent the velocities of the charges in the wire, the $+$ and $-$
charges separately, but velocity relative to what?

We need a Frame C, in which the positive charges in the wire are at
rest.  In Frame C, the charge density of the positive charges will be
$\lambda_0$ (we don't actually need to know what the charge density of
the negative charges in Frame C is -- we just need an anchor point
with which to compare the others.  If the positive charges are at rest
in Frame C, that means that Frame A is moving with speed $v$ with
respect to Frame C!  The original speed of the charges in Frame A was
$v$.  So this means that $\lambda = \gamma\lambda_0$, where
$\gamma=1/\sqrt{1-v^2/c^2}$.

We now have all the pieces we need to pull it all together.  To make
the math easier to write, I will take $c=1$ for a while, and then put
the $c$ back later.  From comparing Frames A and C, we know $\lambda =
\gamma\lambda_0$, where $\gamma = 1/\sqrt{1-v^2}$.  From comparing Frames
B and C, we know $\lambda_\pm = \pm\gamma_\pm\lambda_0$, where $\gamma_\pm
= 1/\sqrt{1-v_\pm^2}$.  We also know that $v_\pm = (v\mp u)/(1\mp vu)$.
If we plug in the velocity formula for $v_\pm$, we get
```{math}
:label: gampmuv
\gamma_\pm = \frac{1}{\sqrt{1-\left(\frac{v\mp u}{1\mp uv}\right)^2}}
```
find a common denominator and flip it upside down to get
```{math}
:label: gampmflip
\gamma_\pm = \sqrt{\frac{(1\mp uv)^2}{(1\mp vu)^2 - (v\mp u)^2}}
```
Expand out the squared terms in the denominator
```{math}
:label: gampmexpand
\gamma_\pm = \frac{1\mp uv}{\sqrt{(1\mp 2vu+v^2u^2 - v^2\pm 2vu-u^2)}}
```
Note the factors of $2vu$ will cancel, as whatever sign they have, the
other will be opposite.  That will leave $1+v^2u^2-v^2-u^2$ in the
denominator, but that's just the same as $(1-v^2)(1-u^2)$
and $1/\sqrt{1-v^2}$ is just $\gamma$!  So this becomes
```{math}
:label: gampmcollapse
\gamma_\pm = \gamma\frac{1\mp uv}{\sqrt{(1-u^2)}}
```
Now we have equations for $\gamma_\pm$, we can figure out the total
charge density in Frame B:
```{math}
:label: lamtotB
\lambda_{\rm totB} = \lambda_++\lambda_- = \lambda_0(\gamma_+-\gamma_-)=
\gamma\lambda_0\left[\frac{1-vu}{\sqrt{1-u^2}}-\frac{1+vu}{\sqrt{1-u^2}}\right] =
\lambda\left[\frac{-2vu}{\sqrt{1-u^2}}\right]
```
Note that the total density in Frame B is negative!  There are some
important pieces in this equation we can pull out.  Notice the
$2\lambda v$.  This is what we would call $I$ in the original Frame A.
Furthermore $1/\sqrt{1-u^2}$ is just the $\gamma_R$ for the
transformation between Frames A and B.

We know that a line charge density in Frame B will have an *electric*
field given by
```{math}
:label: Eline
\vec{E} = \hat{r}\frac{\lambda_{\rm totB}}{2\pi \epsilon_0 r}
```
But we have an equation for $\lambda_{\rm totB}$, it's Equation {eq}`lamtotB`!
```{math}
:label: Elinelam
\vec{E} = \hat{r}\frac{-Iu}{2\pi \epsilon_0 rc^2}\gamma_R
```
where I put the $c$ back in, now.  Since $\vec{F}_e = q\vec{E}$, and using
$c^2=1/\mu_0\epsilon_0$, we can
write this as
```{math}
:label: EforceB
\vec{F}_E = q\vec{E} = -\hat{r}\frac{qu\mu_0I}{2\pi r}\gamma_R
```
for the *electric* force experienced in Frame B.  If we want to switch
back into Frame A, we use the formula we found for the Minkowski force
(Equation {eq}`M4forcefin`) and drop the $\gamma_R$ to get
```{math}
:label: EfofrB
\vec{F}_{EA} = \vec{F}_{EB}/\gamma_R = -\hat{r}\frac{qu\mu_0I}{2\pi r}
```
This is the exact same force that, classically, we would have said was a
*magnetic* force!  We had the same equation back at the start of the
chapter, as Equation {eq}`Fclassic`!

So which is it?  Magnetic or Electric?  The answer is both and
neither!  We are heading toward an understanding that there is just
one "thing" that we currently call an electromagnetic field, and
whether we *name* it "electric" or "magnetic" depends on which
reference frame you're observing in.  If SR had been developed before
magnetism was discovered, we probably wouldn't even have a word for
"magnetic field" -- there would just be result of transforming an
electric field into another reference frame.  Just like we now talk
about "spacetime", and the space and time parts get mixed up when you
change frames, and just like how $x$ and $y$ get mixed up when you
rotate your coordinates around the $z$ axis, it seems that $vec{E}$
and $\vec{B}$ get mixed up when you switch reference frames, too.

But we can't just make a four vector out of them!  They have six
components, not four.  To understand how to fold $E$ and $B$ into the
SR formalism, we must first figure out how they transform when you
switch frames.  Then we can see if and how the Lorentz transformation
is relevant to them.


## How the Fields Transform

## The Electromagnetic Field Tensor

