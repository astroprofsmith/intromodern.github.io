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

# Chapter 3: Four Vectors

## Time as a Dimension

At the end of the last chapter, we argued that if the speed of light
was to be invariant -- to be measured the same in any reference frame
-- then the quantity $$dx^2 + dy^2 + dz^2 - c^2dt^2$$ must also be
invariant across reference frames.  This formulation looks a lot like
a kind of 4 dimensional Pythagorean theorem, except for the minus sign
in the last term.  If it were not for that minus sign, we could well
imagine we were looking at the magnitude of a four-dimensional
displacment vector, just adding time as a new component.  If a
spacial displacement vector is $d\vec{r} = (dx,dy,dz)$, and its
magnitude is $dx^2+dy^2+dz^2$, perhaps it would make sense to have
a four-dimensional vector $dx_4 = (cdt,dx,dy,dz)$, except that
if we were to square each term and add them up, we would get
$+c^2dt^2$, not minus.  We need a way to remember to include that
minus sign.

There are two ways of keeping track of the minus sign.  One is simpler
for the beginner, and the other makes the more complicated material
later easier to manage.  The simpler way is to introduce complex
numbers as a bookkeeping trick.  If $i^2 = -1$, then we can write
the time component of the vector as $icdt$, and its square will be
$-c^2dt^2$.  This method ensures that the minus sign always goes in
the right place.  There is nothing physical or mysterious about the
$i$; it is no more and no less than a handy way to keep track of where
to put the minus sign.  It always shows up in the time component of the
four-vector, and it never shows up in the space components.

```{warning}
If you ever see a time component without an $i$, or find an $i$ in
a space component, then you have made a mistake in your math somewhere.
```

With this tool, we can define a displacement four-vector as
\begin{equation}
[dx_4] =
\begin{bmatrix}
icdt\\
dx\\
dy\\
dz
\end{bmatrix}
\end{equation}
The "size" of this four-dimensional object would then be $$[dx_4]^2 =
dx^2 + dy^2 + dz^2 - c^2dt^2.$$ We need a better word than "size", but
"magnitude" is not a good choice, because we are so used to how
magnitudes work with three-dimensional vectors.  The magnitude of a
three-vector is always positive (unless all three of its components
are zero), but $[dx_4]^2$ could be positive, negative, or zero!
For a displacement four-vector, we call the size the "interval", but
there isn't a better word than "size" for other four vectors.

What a person measures when light propagates is a set of four things:
the three components of the physical displacement and the time
interval it took to make this displacement. These values, when
arranged as above, can be thought of as 4 components of a
vector.  The vector is called the displacement 4-vector for the events.

```{warning}
There is no consensus on what order to put the four components.
Some people put the three space components first and time in the
fourth position.  My own preference is to count from zero and put
time in the zeroth (first) position.  Then the space components
are still 1, 2, and 3, but time is a bit different, just like zero
is a bit different from the positive integers.  Having time first
also saves writing.  We often don't need to use $y$ or $z$,
but if you put time last, you always need to write all four.  If you
put time first, you can sometimes only write the first two.  Finally,
if you find yourself in the unusual circumstance of needing to have
more than three spatial dimensions, it is easy to add 4, 5, 6, and so
on, but if time is in space four, then you have space dimensions
1, 2, 3, 5, 6, 7, etc., which is confusing and harder to program in
a computer.
```

In general, a four vector is an object with four components,
one of which gets a minus sign when you square the components,
and the other three have the properties of normal three vectors.
The size of any four vector is the same in any inertial reference
frame, and to convert a four vector from one frame to another,
there is a specific procedure one must follow (see Chapter 4).

The second way involves being aware of a distinction between covariant
and contravariant four-vectors.  The mathematical details are not
necessary at this stage (see Appendix?), but suffice to say that the
contravariant version of a four vector should be written as a column,
as in Equation 6, only without the $i$.  The standard notion is to use
a superscript greek letter to refer to the components, where the Greek
letter could stand for 0, 1, 2, or 3.  So
\begin{equation}
dx^\alpha =
\begin{bmatrix}
cdt\\
dx\\
dy\\
dz
\end{bmatrix}
\end{equation}
```{warning}
The superscript here does NOT mean "raised to the power of".  In this
notation, $dy$ would be written $dx^2$, but the 2 does not mean squared,
it means the third component of the contravariant four-vector.  If
people are using this notation, you have to be aware from the context
whether they mean a contravariant component or an exponent.
```
There is also a covariant version of the four-vector, which is written
as a row with a subscript Greek letter, but also includes the minus sign: 
$dx_\alpha = [-cdt,dx,dydz]$.  For the four vectors we are considering,
$dx_0 =-dx^0$, but $dx_1=dx^1$ and so forth.  So if you multiply the
two together, according to the rules of multiplying matrices,
$$dx_\alpha dx^\alpha = dx_0dx^0 + dx_1dx^1+dx_2dx^2+dx_3dx^3$$
$$dx_\alpha dx^\alpha = -dx^0dx^0 + dx^1dx^1+dx^2dx^2+dx^3dx^3$$
$$dx_\alpha dx^\alpha = -c^2dt^2 + dx^2+dy^2+dz^2$$
where in that last equation, the 2 does mean squared.

In this notation, whenever there is a Greek letter that appears twice,
once up and once down, that is a shorthand for "multiply across the
row and down the column and add them up".  It plays the same role as a
dummy variable in integration, and therefore does not appear in the
result.

Once you have the concept of a four-vector, mathematically, it is
useful to also have a method of displaying them, graphically.  Such a
graph is called a "spacetime diagram".  I will give a brief
introduction to them here, and then we will explore their properties
in much more depth in the next chapter.

%The number of Greek letters you need to specify a particular
%component is called the "rank".  A four vector only needs one number,
%so it is a first-rank object.  First-rank, but four-dimensional, because
%there are four choices for that one number.  A second rank object could
%be written as a matrix, like
%\begin{equation}
%g_{\alpha\beta} =
%\begin{bmatrix}
%-1 & 0 & 0 &0\\
%0 & 1 & 0 & 0\\
%0 & 0 & 1 & 0\\
%0 & 0 & 0 & 1
%\end{bmatrix}
%\end{equation}
%In this case, if you saw $g_{\alpha\beta}dx^\beta$, the repetition
%of beta would mean multiply across the columns of $g$ and down the
%rows of $dx$, but the $\alpha$ would stay the same, so you would get
%$$dx_\alpha = g_{\alpha\beta}dx^\beta = [-cdt,dx,dy,dz]$$

## Spacetime Diagrams

Superficially, a spacetime diagram hardly seems worthy of its own
distinct name: the convention is to use the vertical axis to represent
time measurements (clock readings), and the horizontal axis to represent
one of the spatial dimensions (ruler readings).  Traditionally, we choose
to call this space axis $x$.  These are usually all the dimensions we
have on a two-dimensional display like a computer screen, although we
can cheat a little by adding a second spatial dimension in projection,
coming out of the screen.  A three-dimensional representation of such a
spacetime diagram is shown in Figure 3.1, along with some common features.

It is also traditional to multiply the time axis by $c$, using $ct$ as
the dimension, rather than just $t$.  This has two benefits: first,
this turns the dimensions of the axis from time into space, which
allows it to have the same scale as the horizontal axis.  This makes
it possible to compare lengths on both axes, which otherwise would
have different dimensions and therefore be incomparable (which is
bigger, one second or one meter?).  Also, the displacement four vector
has a time component of $cdt$, so having the $c$ included in the axis
means that distances on the diagram can accurately represent the
components of the displacement four-vector.

```{code-cell}
:tags: ["remove-input"]
# 3D plot of a spacetime diagram with x, ct, and y
url='https://glowscript.org/#/user/dasmith/folder/Public/program/SRspacetimediag'
display.IFrame(src=url,width=800,height=600)
```
```{note}
Figure 3.1 -- Three-dimensional spacetime diagram, which you can
rotate to observe from different angles.  Events are represented as
two spheres: orange and purple.  The displacement four vector between
them is black.  The projections of this four vector are also drawn
within each of the three planes.  The black arrow in the magenta plane
would be the spatial displacement (in $x$ and $y$ only, not $z$),
while the cyan and yellow planes also include the elapsed time between
the events.  Click and drag to rotate the "cube" so that you are
looking face-on to the side bounded by $x$ and $ct$, with the origin
in the lower left corner.  This is the standard way to draw a 2D
spacetime diagram.  The $y$ axis would then be diagonal in projection.
There is nowhere to put a $z$ axis on the computer screen.
```

Much of special relativity amounts to comparing a set of observations
in one reference frame to that of a second reference frame moving
relative to the first.  By convention, we usually define the $x$ axis
to point along the direction of the relative motion, so that the
second frame is moving to the right.  This is the default assumption,
although it is certainly not **necessary**.  If you must choose $x$
to lie in a different direction than the relative velocity, that
will make the math more complicated, so always make sure to draw
a picture before you start doing math, to make sure the equations
actually match the setup.  Figure 3.2 shows the standard way to
represent two spacetime diagrams in relative motion.

Events, being infinitesimal in duration and extent, are represented on
the spacetime diagram as dots.  A collection of events that represent
the motion of a single object through space over time are called a
worldline.  The tip of my nose at any given moment could be said to be
an event, and if I am not moving, then the worldline of my nose would
be a vertical line on a spacetime diagram.  If I were to walk at a
constant speed, my nose's worldline would be a tilted straight line,
since for every time interval $cdt$ I would have moved my nose a
distance $dx=vdt=\beta cdt$.  The slope of the line, rise over run,
would therefore be $1/\beta$, since we put the time axis on the
vertical.  As my speed goes to zero, the slope becomes infinite: a
vertical line.  The lowest slope the worldline of a real object can
have is 1, since that would imply $1/\beta=1$ or $v=c$, the fastest a
real object can move.

We'll explore the implications of these properties more in the next chapter.
For now, it is sufficient to note that we can define a displacement
four-vector between two dots on a spacetime diagram, and that the
components of this four vector represent horizontal and vertical sides
of that triangle.  Be aware, however, that although this looks exactly
like a triangle, that pesky minus sign is still there, which means the
"size" of the hypoteneuse is not constrained to be positive, as it would
with a normal triangle!  If the sides are equal, the length is zero, not
$dx\sqrt{2}$.

```{code-cell}
:tags: ["remove-input"]
# 3D plot of a spacetime diagram with x, ct, and y
plt.figure(figsize=(5,5))
plt.arrow(0,0,1,0,head_width=0.1)
plt.arrow(0,0,0,1,head_width=0.1)

plt.arrow(-0.5,-0.5,1,0,head_width=0.1)
plt.arrow(-0.5,-0.5,0,1,head_width=0.1)

plt.arrow(0.2,-0.25,0.75,0,head_width=0.05)
ax = plt.gca()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.axis([-1,2,-2,2])
ax.text(0.15, 1, "S")
ax.text(-0.05, 1.2, "ct")
ax.text(1.2, 0, "x")
ax.text(1.1, -0.25, "v_R")
ax.text(-0.35, .5, "S'")
ax.text(-.55, .7, "ct'")
ax.text(.7, -0.5, "x'")
plt.show()
```
```{note}
Figure 3.2 -- Spacetime diagram of two reference frames in relative motion.
The primed frame S' (drawn to the left and down a little) is moving to the
right at speed v_R, relative to the unprimed frame S, which is above and to
the right.  The x and x' axes are defined to be parallel to the vector direction
of v_R.
```


## Examples of Four Vectors

