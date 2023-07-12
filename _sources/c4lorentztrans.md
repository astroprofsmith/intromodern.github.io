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

# Chapter 4: Lorentz Transformations

We have already seen that requiring the size of a four-vector
to be invariant across reference frames leads to interesting
implications (time dilation, for example).  However, there are
more insights that can be gained by figuring out a method to
predict what the components of a four vector will be in a second
reference frame, given measurements of that four vector's components
in a first reference frame.  The change of a four vector from
one reference frame to another with a constant relative velocity
is called a **Lorentz Transformation**.


## Setting Up the Problem

The Lorentz transformation converts the components of a 4-vector
in one inertial reference frame to the components as measured in
a second frame moving at a constant velocity relative to the first.
This method works for **any** four-vector, but we will define it
using the displacement four-vector, as that is easy to visualize.

We start by recreating Figure 3.2 as Figure 4.1 -- to make sure
you are applying the transformation properly, it is useful to always
start by recreating this diagram and making sure that you are assigning
the correct numbers to the elements, as defined in this diagram.



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
Figure 4.1 -- The spacetime diagram that defines a Lorentz transformation.
There are two events, indicated by red and blue dots, and we define a
displacement four vector between them (indicated by an arrow).  These
components are defined in the reference frame $S$, and we want to know
the components in a second frame $S^\prime$, which is moving with a
constant relative velocity $\vec{v}_R$.  We **define** $x$ such that it lies
parallel to $\vec{v}_R$.
```


The 4-vector components are determined by the two events shown as
colored dots in Figure 4.1. Since we have the choice of which is the
x-axis, we will always pick the relative velocity of the two observers
to lie along the $x$-axis, to make the math simpler. Putting the
diagram of Figure 4.1 into 4-vector notation, we need a tool that
will take us from $[dx_4]$ to $[dx_4]^\prime$, like so:
```{math}
:label: eqlort
[dx_4]^\prime = {\cal L}_x(\beta_R)[dx_4]
```

The Lorentz transformation can be mathematically represented as a
matrix, because it transforms one vector into another. The Lorentz
transformation, therefore, is a $4 \times 4$ matrix because there are
4 elements in the vector. If $L_{\alpha\beta}$ is a an arbitrary
element of this matrix (row $\alpha$ column $\beta$, where the
Greek letters stand for 0, 1, 2, and 3), the transformation
matrix has the form:
```{math}
:label: lormat1
{\cal L}_{\alpha\beta} = 
\begin{bmatrix}
L_{00} & L_{01} & L_{02} & L_{03}\\
L_{10} & L_{11} & L_{12} & L_{13}\\
L_{20} & L_{21} & L_{22} & L_{23}\\
L_{30} & L_{31} & L_{32} & L_{33}
\end{bmatrix}
```
The problem then is to figure out what the 16 numbers in ${\cal L}$
are that will achieve our goal.

```{warning} 
Note that the $\beta$ used as the index counter of the matrices has
nothing to do with the relative velocity $\beta_R$.  Sometimes other
Greek letters are used for counting indices, $\mu$ and $\nu$ in
particular.  You must notice how the letter is being used to be
able to tell whether $\beta$ is a speed or a matrix index.
```

## Determining the Elements

Is there a way that we can determine some of the elements of the of
the Lorentz transformation? We must figure out how the transformation
works in specific cases, and these cases will then point to general
values that could be used in any situation.  All that we really know
about the Lorentz transformation is that the 4-vector that is produced
by its operation on the original 4-vector must have the same size as
the first 4-vector, or the Michelson-Morley experimental results would
be contradicted.

To determine some of the elements of the Lorentz transformation
matrix, consider what happens when the displacement 4-vector
```{math}
:label: eqdisp4
[dx_4] =
\begin{bmatrix}
icdt\\
dx\\
dy\\
dz
\end{bmatrix}
```
is transformed from one inertial reference frame to another. If the
relative velocity is in the $x$-direction, then displacements in the
$y$ and $z$ directions won't change.  This derives from the independence
of motion along perpendicular axes: if I am 10 m south of you, and I
take one step East, I will still be 10 m South of you (assuming we are
far from the poles, of course).

```{note}
Although the displacements perpendicular to the direction of relative
motion are not affected by a Lorentz transformation, this independence
does not apply to **any** object that you might perform a Lorentz transformation
on.  In electromagnetism, for example, the electric and mangetic fields in
the $y$ and $z$ direction are different in a frame moving in the $x$ direction.
We use the dimensional independence of $dy$ and $dz$ to derive values for
${\cal L}$, but once we do that we don't have to insist that the $y$ and
$z$ for all vectors are unaffected in all circumstances at all times.  Be cautious.
```

This means that $dy^\prime = dy$.  However, if we write out the third
row of Equation {eq}`eqlort`, we get $$dy^\prime = L_{20}icdt+L_{21}dx
+L_{22}dy + L_{32}dz.$$
The only way this equality can hold for all displacement four vectors
is if $L_{22}=1$ and the other elements are zero.  The same logic hold
for the fourth row, so we know the matrix must look like
```{math}
:label: lormat2
{\cal L} = 
\begin{bmatrix}
L_{00} & L_{01} & L_{02} & L_{03}\\
L_{10} & L_{11} & L_{12} & L_{13}\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
```
That takes care of half the numbers in one go!

To get the next set, we have to look at a particular case.
Consider the examples at the end of Chapter 3.  They all
involved switching from a reference frame where two events
were at rest into a reference frame in which there would have
to be motion for something to get from one event to the other.
We know therefore that the original four displacement is
```{math}
:label: eqrestdx4
[dx_4] =
\begin{bmatrix}
icdt_0\\
0\\
0\\
0
\end{bmatrix}
```
We can also figure out what the transformed components have to be.  We
know that the time interval will be dilated in the primed frame: $dt =
\gamma_R dt_0$.  Furthermore, if the events are at rest in the
original frame, and the primed frame is moving to the right at speed
$\beta_R$, then from the point of view of the primed frame, the
location of the events would have to be moving left at the same speed,
so $dx^\prime = -v_Rdt = -\beta_R cdt = -\beta_R\gamma_R cdt_0$, where
we have used time dilation to turn the $dt$ into $dt_0$.

That gives us two further relations that will enable us to get four of the
numbers.  We look at the first two rows of Equation {eq}`eqlort`, using
Equation {eq}`eqrestdx4` as the four vector:
$$icdt^\prime = L_{00}icdt_0$$
and
$$dx^\prime = L_{10}icdt_0.$$
The first of these two equations will get you back the time dilation
if you take $L_{00}=\gamma_R$, and then the second equation matches with
the equation for $dx^\prime$ as long as $L_{10}=i\beta_R\gamma_R$.

This leaves us with (so far):
```{math}
:label: lormat3
{\cal L} = 
\begin{bmatrix}
\gamma_R & L_{01} & L_{02} & L_{03}\\
i\beta_R\gamma_R & L_{11} & L_{12} & L_{13}\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
```

The number of unknowns in the Lorentz transformation matrix has
decreased from 16 to 6. The rest can be found using similar tricks.
For example, you could set up the events in Figure 4.1 such that the
primed observer was in the rest frame of the events and work out the
values for the rest of the unknowns. When all such tricks are
exhausted, the Lorentz transformation between two frames of reference
(Fig 4.1) that have a relative speed of $\beta_R$ in the positive x
direction is given by the matrix:
```{math}
:label: lormat
\boxed{
{\cal L}_x(\beta_R) = 
\begin{bmatrix}
\gamma_R & -i\beta_R\gamma_R & 0 & 0\\
i\beta_R\gamma_R & \gamma_R & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
}
```
The subscript $x$ is to remind us that the relative velocity is
in the $x$ direction.  If that were not the case, the matrix
would look a lot more complicated!  Note that this formulation
is based on the diagram in Figure 4.1 -- if the velocity were in
the opposite direction, the $\beta_R$ would become $-\beta_R$,
the **inverse Lorentz transformation**.

## Examples of Lorentz Transformations

### Example 4.1: Transform a Displacement Four-vector

An observer measures the 4-displacement between two events:
```{math}
:label: eqdx41
[dR_4] =
\begin{bmatrix}
i3.00~{\rm m}\\
1.00~{\rm m}\\
-0.250~{\rm m}\\
0.00~{\rm m}
\end{bmatrix}
```
A second observer, traveling with speed $1.31\times10^8$~m/s in the
$x$ direction with respect to the first observer measures the
4-displacement for these same two events. What values of the elements
of the 4-displaceinent does this second observer measure?

To solve this problem, you use the form of the Lorentz transformation
in equation {eq}`lormat`, plugging in the numbers.  We know that
$\beta_R = v_R/c = 1.31/3.00 = 0.438$, and therefore $\gamma_R=1.11$.
Multiply the two to get $\beta_R\gamma_R=0.486$.

Now, we are ready to find the elements of $[dR_4]'$
```{math}
:label: ex41
[dR_4]' = {\cal L}_x(\beta_R)[dR_4] = 
\begin{bmatrix}
1.11 & -i0.486 & 0 & 0\\
i0.486 & 1.11 & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
i3.00~{\rm m}\\
1.00~{\rm m}\\
-0.250~{\rm m}\\
0.00~{\rm m}
\end{bmatrix}
```
multiplying out the matrix gives
```{math}
[dR_4]' =
\begin{bmatrix}
1.11(i3.00~{\rm m})-i0.486(1.00~{\rm m})+0+0\\
-0.486(3.00~{\rm m})+1.11(1.00~{\rm m})+0+0\\
0+0+-0.250~{\rm m}+0\\
0+0+0+0.00~{\rm m}
\end{bmatrix}
=
\begin{bmatrix}
i2.79~{\rm m}\\
-0.351~{\rm m}\\
-0.250~{\rm m}\\
0.00~{\rm m}
\end{bmatrix}
```
As a check, you can note that the size of the two four-vectors is the
same:
$$[dR_4]^2 = (-3.00^2+1.00^2+(-0.250)^2+0.0^2)~{\rm m}^2 = -8.009~{\rm m}^2$$
and
$$[dR_4]'^2 = (-2.79^2+(-0.351)^2+(-0.250)^2+0.00^2)~{\rm m}^2 = -8.009~{\rm m}^2$$
Note which minus signs are squared (and therefore go away) and which
ones are not.  Also note that all components have units of length.

### Example 4.2: Angles in Space

Find the physical angle that the displacement given in example 4.1
makes with respect to the $x$ axis. Determine the angle as measured by
the primed observer.

The angle that tlie unprimed observer makes is:
$$\tan{\theta} =\frac{dy}{dx} = \frac{1.00~{\rm m}}{-0.250~{\rm m}} = -4.00
\rightarrow \theta = -76.0^\circ$$
The angle is in the 4th quadrant.  
This makes sense since dy is negative and dx is positive.

The second observer moving with respect to the first measures an angle:
$$\tan{\theta^\prime} =\frac{dy'}{dx'} = \frac{-0.351~{\rm m}}{-0.250~{\rm m}} = 1.40
\rightarrow \theta^\prime = 234.0^\circ$$
This angle is in the third quadrant as both $dx'$ and $dy'$ are negative.
Essentially, the second observer is overtaking the first, so the
horizontal displacement flips around and points the other way.

### Example 4.3: Using an Inverse Lorentz Transformation

The primed observer in example 4.1 believes that the unprimed observer
is moving with speed $v_R = 1.3\times10^8$ m/s in the $-x$ direction
with respect to him. Show that this is a reasonable assumption by
finding the elements of the inverse Lorentz transformation that
transform the coinponents she measured for $[dR4]'$ into those
measured by the observer in the unprimed reference frame.  In
other words, imagine you didn't know $\beta_R$ and find it from
the given four vectors.

Going back to the transformation in example 4.1, but setting it
up as an inverse transformation, we get:
```{math}
:label: eqinvlort
[dx_4] = {\cal L}_x(-\beta_R)[dx_4]'
```
```{math}
\begin{bmatrix}
i3.00~{\rm m}\\
1.00~{\rm m}\\
-0.250~{\rm m}\\
0.00~{\rm m}
\end{bmatrix}
=
\begin{bmatrix}
\gamma_R & i\beta_R\gamma_R & 0 & 0\\
-i\beta_R\gamma_R & \gamma_R & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
\begin{bmatrix}
i2.79~{\rm m}\\
-0.351~{\rm m}\\
-0.250~{\rm m}\\
0.00~{\rm m}
\end{bmatrix}
```
Note the minus sign has been switched.
The $y$ and $z$ components don't change, but if you work out the
first two components and then take their ratio, you get
$$\frac{dx}{cdt}=\frac{1.00~{\rm m}}{3.00~{\rm m}} =
\frac{\gamma_R[-0.351~{\rm m} + \beta_R(2.79~{\rm m})]}{\gamma_R
[2.79~{\rm m}-\beta_R(0.351~{\rm m})]}$$
The factors of $\gamma_R$ cancel, and then you can sove for $\beta_R$
to get $\beta_R = 0.44$ which means $v_R = 1.3\times10^8$ m/s,
as expected.

```{note}
Programming a Lorentz Transformation

It would be a very good idea for you to write your own code that
takes a $\beta_R$ and a four vector and returns the transformed
four vector.  You could also write it to return the size of the
four vector before and after the transformation, so you could verify
that it hasn't changed.  You could write the code to set up and
multiply matrices, or you could write the four equations out
separately.  Which method you use would be irrelevant to the output.
```