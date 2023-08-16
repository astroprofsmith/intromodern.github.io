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

** Are there any figures I can make for this stuff?**

# Chapter 12: Relativistic Electrodynamics

Now that we have an expression for the electromagnetic field tensor,
which treats $\vec{E}$ and $\vec{B}$ as different parts of a single
object, and we know how that tensor transforms under Lorentz
transformations, we can start to see how the rest of electrodynamics
plays out within an SR context.  This chapter is almost entirely just
showing how you can write classical electrodynamical models to show
that they are intrinsically consistent with the rules of special
relativity already.  We will start with examining how you represent
the Minkowski Force for a charge moving with velocity $\vec{u}$ in a
region where there are both $\vec{E}$ and $\vec{B}$ fields.

## The Lorentz Force

You learned in introductory E&M that a charge
moving with velocity $\vec{u}$ in a region where
there are both $\vec{E}$ and $\vec{B}$ fields experiences
a force given by
```{math}
:label: lorforce
\vec{F} = q\vec{E} + q\vec{u}\times\vec{B}
```
Given what we now know about how $\vec{E}$ and
$\vec{B}$ are parts of a larger tensor, how does
that affect the way we conceptualize this force?

We start with our definition of the Minkowski Force
from Equation {eq}`M4forcefin`:
```{math}
:label: M4forceEM
K^\alpha = 
\begin{bmatrix}
\frac{\gamma}{c} (\vec{F}\cdot\vec{u})\\
\gamma \vec{F}
\end{bmatrix}
```
and if we plug in Equation {eq}`lorforce`, we get
```{math}
:label: M4EM
K^\alpha = 
\begin{bmatrix}
\frac{\gamma}{c}q (\vec{E}\cdot\vec{u})\\
\gamma q (\vec{E} + \vec{u}\times\vec{B})
\end{bmatrix}
```
Note that $\vec{B}$ has disappeared from the time component because
the cross product of $\vec{u}\times\vec{B}$ always must be orthogonal
to $\vec{u}$, so if you take the dot product of $\vec{u}\cdot(\vec{u}
\times\vec{B})$, you will always get zero.

If we write out all four components of Equation {eq}`M4EM`, we get
```{math}
:label: M4EMall4
\begin{bmatrix}
K^0\\
K^1\\
K^2\\
K^2
\end{bmatrix}
=
\gamma q
\begin{bmatrix}
(E_xu_x + E_yu_y + E_zu_z)/c\\
E_x +u_yB_z - u_zB_y\\
E_y +u_zB_x - u_xB_z\\
E_z +u_xB_y - u_yB_x
\end{bmatrix}
```

Let's examine these terms carefully, and compare them to the EM field
tensor (Equation {eq}`EMtensorfin`).  Note that the time component is just
what you would get if you multiplied each term across the top row of
the EM tensor with the corresponding term of the velocity four vector
(the time component of $u^\alpha$ is eliminated by the zero on the
diagonal of the tensor).  That suggests this equation might include
the product of the EM tensor with the four-velocity.

What happens if you multiply the terms across the second row of the EM
tensor with the corresponding terms of the four velocity?  The $c$
will cancel when multiplying $-E_x/c$ by $\gamma c$, then the $u_x$
will be eliminated by the zeros along the diagonal.  That leaves $B_z$
times $u_y$ and $-B_y$ times $u_z$, which is exactly what we see in
Equation {eq}`M4EMall4`!  You can check the other two components, but
this is nothing more or less than writing:
```{math}
:label: fullLorF
\boxed{
K^\alpha = qu_\beta F^{\alpha\beta}
}
```
Amazing that all those terms can be packed into such a small equation!

The Einstein notation needs a bit of unpacking, here.  The fact that
$\beta$ appears twice, once up and once down, means that we are
summing over all the possible values of $\beta$, which is the index
for the column in the EM tensor.  The $\alpha$ index refers to the
row.  You can tell this because you have the contravariant version of
$K^\alpha$, which means it will be a column vector, and the $\alpha$
counts which row each component is in.  If $\alpha$ is the row on the
left side of the equation, it must also be the row on the right side.
So $F^{13}$ is the second row, fourth column ($-B_y$).

If we work out the matrix multiplication, the first row of this matrix
equation would look like
```{math}
:label: LFrow0
K^0 = q(u_0F^{00} + u_1F^{01} + u_2F^{02} + u_3F^{03})
= q\gamma (0.0 + u_xE_x + u_yE_y + u_zE_z)
```
which matches Equation {eq}`M4EMall4`, and the second row is
```{math}
:label: LFrow1
K^1 = q(u_0F^{10} + u_1F^{11} + u_2F^{12} + u_3F^{13})=
q\gamma ((-c)(-E_x/c) + 0.0 + u_yB_z - u_zB_y)
```
```{math}
:label: LFrow1b
K^1 = \gamma q(\vec{E} + \vec{u}\times\vec{B})\cdot\hat{x}
```
which is the $x$ component of the Lorentz Force! (times $\gamma$)
You can work out $K^2$ and $K^3$ and show they are the $y$ and
$z$ components of the Lorentz force yourself.

Note that the index $\beta$ for the four velocity in Equation
{eq}`fullLorF` is a subscript.  When you lower the index, you also
make the time component negative.  Hence the $-c$ in Equation
{eq}`LFrow1`.  The first column of the EM tensor ($F^{\alpha 0}$, see
Equation {eq}`EMtensorfin`) has negative signs on all the $E$ components, so you
need the $\gamma c$ component of the four-velocity to be negative,
which makes the final answer positive, as it needs to be for Equation
{eq}`M4EMall4`.  If you use it consistently, Einstein notation will
always get you the right answer, but you have to be careful about
using it right.

```{warning}
You have to watch the greek letters carefully in Einstein notation.
If you just looked at Equation {eq}`fullLorF` and thought "oh,
that's a row vector times a matrix, I'll multiply across the row
vector and down the columns of the matrix," you would get the wrong
answer.  Even though the $u$ is first, you are still multiplying
each term in $u$ by the terms across each row (counting the columns
with $\beta$) in $F$, because you have to write it out like Equations
{eq}`LFrow0` and {eq}`LFrow1` to see which terms are really multiplying
which other terms.
```

The notation of SR therefore gives us a way to write the electric and
magnetic forces on a particle as a *single* Minkowski force, that is
the result of interactions between a moving charge and a single
electromagnetic field tensor.  There is no need to conceptualize
electric and magnetic fields and forces as separate entities -- they
are aspects of a single phenomenon.

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
four-vector, but what do we do with the time part?  Well, we have a
scalar $\rho_0$ times $\gamma \vec{v}$, and $\gamma \vec{v}$ is the
spatial part of the four-velocity.  Maybe we can make the four-current
by multiplying the four-velocity $v^\alpha$ by $\rho_0$?  And indeed
we can:
```{math}
:label: fourcurr
\boxed{
J^\alpha =
\rho_0 \gamma 
\begin{bmatrix}
c\\
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
The right side is the rate of change of the charge density at a particular
location, and the left side is the divergence of the vector field.  So if
the amount of charge is decreasing, $\partial \rho/\partial t <0$ and the
divergence is positive, which means the current is flowing outward.  If the
current is flowing inward, then the divergence is negative, and the time
derivative is positive, which means the amount of charge at that location
is increasing, which makes sense.

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
You will note that the four-divergence $\partial_\alpha$ is written in
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

We have combined $\rho$ and $\vec{J}$ into a single four-vector.
We have also combined $\vec{E}$ and $\vec{B}$ into a single tensor.
Does this mean we can combine the four Maxwell's Equations into a
single equation, too?  Well, almost.  In this section, I will show you
how you can write them as two (shorter) equations.  To collapse them
down to a single equation requires learning about the four-vector form
of the potential, and that is slightly beyond the scope of this book.
See more advanced texts like Griffiths that continue to four-potentials,
Gauge transformations, and the D'Alembert Equation.

For now, let's remind ourselves of the four classic Maxwell's Equations:
```{math}
:label: gauss
\vec{\nabla}\cdot\vec{E} = \frac{\rho}{\epsilon_0}
```
```{math}
:label: nogaussB
\vec{\nabla}\cdot\vec{B} = 0
```
```{math}
:label: faraday
\vec{\nabla}\times\vec{E} + \frac{\partial \vec{B}}{\partial t} = 0
```
```{math}
:label: ampere
\vec{\nabla}\times\vec{B} - \mu_0\epsilon_0 \frac{\partial \vec{E}}{\partial t} = \mu_0\vec{J}
```

Now, I have already introduced the four-divergence in the continuity
equation (Equation {eq}`fourdivdef`).  What happens if we take the
four-divergence of the EM field tensor?  Using Equation
{eq}`EMtensorfin`, that works out to
```{math}
:label: div4EM
\partial_\nu F^{\mu\nu} = \left(\frac{\partial}{\partial (ct)},\frac{\partial}{\partial x},\frac{\partial}{\partial y},\frac{\partial}{\partial z}\right)
\begin{pmatrix}
0 & E_x/c & E_y/c & E_z/c\\
-E_x/c & 0 & B_z & -B_y\\
-E_y/c & -B_z & 0 & B_x\\
-E_z/c & B_y & -B_z & 0
\end{pmatrix}
```

Let's start with the first row.
```{math}
:label: divrow0
\partial_\nu F^{0\nu} = \partial_0 F^{00} + \partial_1 F^{01}
+ \partial_2 F^{02} + \partial_3 F^{03}
```
The zero in the upper left corner of the EM tensor ($F^{00}$) kills
the time derivative ($\partial_0$).  For the rest, we have the
successive derivatives of each component of $\vec{E}$ with respect to
each successive space dimension.  That's a three-divergence.  So that
first row works out to $\vec{\nabla}\cdot\vec{E}/c$.  But Equation
{eq}`gauss` tells us what the divergence of $\vec{E}$ is!  It's
$\rho/\epsilon_0$!  So we get
```{math}
:label: zerothterm1
\partial_\nu F^{0\nu} = \vec{\nabla} \cdot \vec{E}/c = \frac{\rho}{c\epsilon_0}
```
But $\rho$ is just the zeroth term of the four-current divided by $c$!
So this becomes
```{math}
:label: zerothterm2
\partial_\nu F^{0\nu} = \frac{J^0}{c^2\epsilon_0}
```
and you should also remember that $c^2=1/\mu_0\epsilon_0$, so if we take
on the idea that
```{math}
:label: SRmaxwell1
\boxed{
\partial_\nu F^{\mu\nu} = \mu_0 J^\mu
}
```
We have shown that the time component of this equation is none other than
Gauss's Law (Equation {eq}`gauss`)!
```{margin}
Note that $\mu_0$ is the magnetic
permeability of free space, and do not confuse it with the $\mu$ we are using
to indicate the four-vector components!
```

What about the space components?  Lets look at $\mu=1$:
```{math}
:label: div4EM_1
\partial_\nu F^{1\nu} = \frac{\partial}{\partial x^0} F^{10}
+  \frac{\partial}{\partial x^1} F^{11}
+  \frac{\partial}{\partial x^2} F^{12}
+  \frac{\partial}{\partial x^3} F^{13}
```

```{math}
:label: div4EM_2
\partial_\nu F^{1\nu} = \frac{\partial}{\partial (ct)} F^{10}
+  \frac{\partial}{\partial x} F^{11}
+  \frac{\partial}{\partial y} F^{12}
+  \frac{\partial}{\partial z} F^{13}
= -\frac{1}{c^2}\frac{\partial E_x}{\partial t} + 0.0 +
\frac{\partial B_z}{\partial y} - \frac{\partial B_y}{\partial z} 
```

The first term is just the $x$ component of $-\mu_0\epsilon_0
\frac{\partial \vec{E}}{\partial t}$, while the last two terms are
just the $x$ component of the curl of $\vec{B}$.  The $x$ component of
$\mu_0J^\nu$ is just $\mu_0 J_x$ (the $x$ component of the current
density), so all of this is just the $x$ component of Ampere's Law,
Equation {eq}`ampere`!  You can verify that doing this for the third
and fourth rows will get you the other two components of Ampere's Law,
so the four-divergence of the EM field tensor is no more and no less
than Gauss's Law and Ampere's Law, together, and by writing Equation
{eq}`SRmaxwell1`, we have collapsed those two of Maxwell's Equations
down into one.

What about the other two?  Our first clue is that both Equation
{eq}`nogaussB` and Equation {eq}`faraday` are equal to zero.  So our
final equation will also be equal to zero.  The second clue is that if
you compare the left sides of Equations {eq}`nogaussB` and
{eq}`faraday` to {eq}`gauss` and {eq}`ampere`, you can see that they
look the same, just swapping $\vec{E}$ and $\vec{B}$ around.  At the
end of Chapter 11, we constructed *two* versions of the EM tensor, one
of which had the $E$ and $B$ components switched.  So, if we bring
back the dual tensor $G^{\mu\nu}$ at this point (Equation
{eq}`Dualtensorfin`), and take its four-divergence, we will see that
```{math}
:label: dual4div
\boxed{
\partial_\nu G^{\mu\nu}=0
}
```
will give us the other two Maxwell's Equations.

Again, let's start with the time component: the first row.  The
zero in position $G^{00}$ will kill the time derivative, and then
we have $\vec{\nabla}\cdot\vec{B}$ from the other three components,
just like we got the divergence of $\vec{E}$ before.  We know that
$\vec{B}$ does not have a divergence, so setting this to zero will
take care of Equation {eq}`nogaussB`.

For the spatial components, we take $\mu=1$ as an example, again,
and we get the same result as Equation {eq}`div4EM_1`, just with
$E$ and $B$ swapped:
```{math}
:label: div4dual_1
\partial_\nu G^{1\nu} = \frac{\partial}{\partial (ct)} G^{10}
+  \frac{\partial}{\partial x} G^{11}
+  \frac{\partial}{\partial y} G^{12}
+  \frac{\partial}{\partial z} G^{13}
= -\frac{1}{c}\frac{\partial B_x}{\partial t} + 0.0 -
\frac{\partial E_z}{c\partial y} + \frac{\partial E_y}{c\partial z} 
```
when you include the $y$ and $z$ components, this becomes
```{math}
:label: spacedual
-\frac{1}{c}\frac{\partial\vec{B}}{\partial t} +\frac{1}{c}\vec{\nabla}\times\vec{E}
```
set this to zero and multiply through by $c$ to get Faraday's Law,
Equation {eq}`faraday`.

So Equations {eq}`SRmaxwell1` and {eq}`dual4div` represent all four
of Maxwell's Equations in only two compact versions.  While this could
arguably be described as simply a notational trick, since these two
equations are really representing the same eight separate equations
that the traditional format writes in four (the two divergence equations
are just one equation each, while the two curl equations are three
component equations each).  There is still something aesthetically
pleasing about seeing so much of nature, the entirety of what we
call electromagnetic phenomena, collapsed into two equations (and
eventually into one).

Aside from aesthetics, why might we want to do this?  The main reason
is that we know four-vectors obey the rules of special relativity.  By
expressing the laws of electricity and magnetism in these terms, we
know that they are already consistent with SR.  We don't need to make
major modifications to accomodate relativity, as we did with newtonian
mechanics.  There is nothing equivalent to the speed of light limit
that pops up in relativistic E&M that wasn't there in classical E&M.
Maxwell's Equations automatically incorporate the insights of SR.
What changes when you understand this is the way you envision the
world, not the way the world itself works.

## Summary?

## Problems

(rough ideas)

1) What is the Minkowski force on a moving charge near a current-carrying
wire?  Use a Lorentz transformation to convert this four vector into
the rest frame of the charge.  Compare with derivation in last chapter
that didn't use Lorentz transformations and four-vectors.

2) The plane-wave solution to Maxwell's Equations is a travelling E-M wave
that would take the following form:
```{math}
\vec{E}(x,t) = E_0 \cos{(kx-\omega t)}\hat{y}
```
```{math}
\vec{B}(x,t) = \frac{E_0}{c} \cos{(kx-\omega t)}\hat{z}
```
The intensity of the light wave is proportional to $E_0^2$, and the frequency
is $\omega/2\pi$. The speed is $c=\omega/k$.

a) Write out the EM field tensor for this situation.

b) Carry out the double Lorentz transformation to express this solution in a
reference frame moving at a relative speed of $\beta$ in the $x$ direction.
Find $E$, $k$, $\omega$, and $c$ in the new frame.

c) Verify that $c'=c$ and that the wavelength and frequency change according
to the Doppler formula.

d) What happens to the light intensity?


3) Write a Python program to apply the momentum principle and
iteratively update the velocity and position of a particle with rest
mass $m_0 = 0.5$, that starts at the origin, moving with initial
velocity $\vec{\beta}=0.001\hat{x}$, when $\vec{E} = 10.0~\hat{y}$ and
$\vec{B}=10.0~\hat{z}$ (using units with $c=1$).  Verify that the
trajectory looks like {numref}`circularpathfig`.

4) Derive $\partial_\mu J^\mu = 0$ from $\partial_\mu F^{\mu\nu} = \mu_0 J^\mu$.

5) Show the other components work for Maxwell's Equations.


