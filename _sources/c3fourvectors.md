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
```{math}
:label: eqdx4
[dx_4] =
\begin{bmatrix}
icdt\\
dx\\
dy\\
dz
\end{bmatrix}
```

```{margin} Example 3.1

An observer sees two events occuring. She measures a physical
displacement of $dx = 1.20$ m, $dy = 2.00$ m and $dz = 0.00$ m. The
time interval between these events is measured to be $0.0100$
microseconds. Write down the eleinents of the displacement 4 vector
measured by this observer and find the size of this displacement 4
vector.

The three spatial elements are given in the problem. What I need to do
is to calculate the size of the time component: $cdt = 3\times10^8$
m/s $(1.00\times 10^{-8}$ s) = 3.00 m.  The displacement four vector
is therefore
\begin{equation*}
[dx_4] =
\begin{bmatrix}
i3.00~{\rm m}\\
1.20~{\rm m}\\
2.00~{\rm m}\\
0.00~{\rm m}
\end{bmatrix}
\end{equation*}
The size of this displacement 4 vector is:
$$(1.20^2+2.00^2)~{\rm m}^2 = -3.56~{\rm m}^2$$
```

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
$dx_\alpha = [-cdt,dx,dy,dz]$.  For the four vectors we are considering,
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
result.  This is called "Einstein summation notation" or just
"Einstein notation", and it won't come up again in this book until we
apply SR to Electromagnetism in Chapter XX.  I mention it here to stress
that using the $i$ is not the **only** way to keep track of the minus
sign.

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
bigger, one second or one meter?).  Second, the displacement four vector
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
them is black.  The components of this four vector are also drawn
as arrows parallel to the three axes.  Note that there could well be
enormous displacement in $z$ that would not show up on this diagram.
Click and drag to rotate the "cube" so that you are
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


## Examples of Displacement Four Vectors

The key to using this model of special relativity to describe
physically real happenings is to be able to write down the
displacement 4-vector for that happening. To do this, there must be
two events occurring, events which both observers can detect. Then,
each observer measures the physical displacements andthe time interval
between the two events and converts them into the components of the
displacement 4-vector as shown in Equation {eq}`eqdx4`.

### Muon Decay

For example, consider the following interesting happening. A pion is
traveling along in a bubble chamber, leaving a track of bubbles. At
some time and place, it decays into a muon. I know that something
happens because the track of bubbles changes direction. The muon moves
on for a while and then decays into an electron. Again, I know this
happens because the bubble track changes direction. The two events are
the change in direction of the tracks as shown in Figure 3.3.  Note
that the figure on the left is a diagram of what happens in space, while
the figure on the right is a spacetime diagram!

```{code-cell}
:tags: ["remove-input"]
p = np.array([[0,0,0],
[1,1,5],
[2,1.5,7],
[3,0,15]])
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig.suptitle('Creation and Destruction of a Muon')
ax1.plot(p[:,0], p[:,1], 'k:')
ax1.plot([p[1,0],p[2,0]], [p[1,1],p[2,1]], 'r.-')
ax1.set_xlim([0, 3])
ax1.set_ylim([-1,3])
ax1.set_title('Space Trajectory')
ax1.set_xlabel('x')
ax1.set_ylabel('y')

ax2.plot(p[:,0], p[:,2], 'k:')
ax2.plot([p[1,0],p[2,0]], [p[1,2],p[2,2]], 'r.-')
ax2.set_xlim([0, 3])
ax2.set_ylim([-1,15])
ax2.set_title('Spacetime Diagram')
ax2.set_xlabel('x')
ax2.set_ylabel('ct')

plt.show()
```
```{note}
Figure 3.3 -- Two diagrams to display the motion of particles.  On the left,
we have a plot of $y$ vs. $x$, which shows the trajectories of the particles.
A pion enters from the lower left going up and to the right, turns into a
muon and goes stright to the right, and then turns into an electron that goes
down and to the right.  This diagram tells you nothing about how fast the
particles are going, or how much time passes between the events.  On the right,
we have a spacetime diagram, which tells us nothing about the motion in the $y$
direction on the left plot (the tracks on this diagram do not turn around),
but does convey the velocity in the $x$ direction (not the $y$ or $z$!).
Can you describe how this velocity component changes throughout
this sequence of events?  (Note: this is not realistic.)
```

If we consider the two events of the creation and disappearence of the muon,
observers in the reference frame shown in Figure 3.3 could measure the displacement
and time duration between these events (and we can read them off the graph,
albiet with no units), and then put them into the four-vector format:
\begin{equation}
[dx_4] =
\begin{bmatrix}
icdt\\
dx\\
dy\\
dz
\end{bmatrix}
 =
\begin{bmatrix}
i2\\
1\\
0.5\\
0
\end{bmatrix}
\end{equation}
I put zero as the $z$ displacement although technically we simply have no information
about whether there was any displacement in the $z$ direction.

### Space ship

Often in SR problems, it is important to make a graph showing the
events under consideration, and to indicate the relative motion of
the reference frames, as in Figure 3.2.  This will help you visualize
what's happening and ensure that you set up your equations properly.

Consider a space ship that leaves the Earth and travels 4.2 LY to the
star Alpha Centauri at a constant speed of $\beta=0.75$.  For this
scenario, there are two reference frames that are of interest: the one
that is at rest with respect to the space ship and another that is at
rest with respect to the earth. The two events will be the departure
from earth and the arrival at alpha centauri. The relative speed
between these two observers is just the speed of the space ship in the
direction of Alpha Centauri. The events can be diagrammed as shown in
Figure 3.3.

An observer that is at rest with respect to the earth will measure
(evenually -- once all the rulers and clocks are returned and
collated) a physical displacement in the $x$-direction of the distance
between earth and Alpha Centauri: 4.2 light years. There will be no
physical displacement in the $y$ or $z$ directions. There will be some
time interval $dt_\oplus$ that the earth based observer will measure
for the trip. Since we know the velocity that the space is traveling
with respect to the earth, the time it takes to get to the star in the
Earth's reference frame is $$dt_\oplus = \frac{dx_\oplus}{v_{\rm
rocket}} = \frac{4.2~{\rm LY}}{0.75c} \rightarrow cdt_\oplus =
5.6~{\rm LY}$$ Note that $cdt_\oplus$ has units of length.


```{code-cell}
:tags: ["remove-input"]
dx=4.2
cdt = 5.6
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig.suptitle('Journey to Alpha Centauri')
ax1.arrow(0,0,5,0,head_width=0.2)
ax1.arrow(0,0,0,6,head_width=0.2)

ax1.arrow(-1.5,-1.5,2,0,head_width=0.2)
ax1.arrow(-1.5,-1.5,0,2,head_width=0.2)

ax1.arrow(2,-0.75,1.75,0,head_width=0.1)

ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)
ax1.axis([-2,6,-2,8])
ax1.text(1.0, 6, "S")
ax1.text(-0.05, 6.5, "ct")
ax1.text(5.5, 0, "x")
ax1.text(4.5, -0.8, "v_R")
ax1.text(-1.1, .1, "S'")
ax1.text(-1.55, 1.0, "ct'")
ax1.text(1.0, -1.5, "x'")
ax1.set_title('Rest Frame of Earth')

ax1.plot([0.5],[0.5],'ro')
ax1.plot([0.5+dx],[0.5+cdt],'bo')

ax2.arrow(0,0,5,0,head_width=0.2)
ax2.arrow(0,0,0,6,head_width=0.2)

ax2.arrow(-1.5,-1.5,2,0,head_width=0.2)
ax2.arrow(-1.5,-1.5,0,2,head_width=0.2)

ax2.arrow(4,-0.75,-1.75,0,head_width=0.1)

ax2.get_xaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)
ax2.axis([-2,6,-2,8])
ax2.text(1.0, 6, "S'")
ax2.text(-0.05, 6.5, "ct'")
ax2.text(5.5, 0, "x'")
ax2.text(4.5, -0.8, "v_R")
ax2.text(-1.1, .1, "S")
ax2.text(-1.55, 1.0, "ct")
ax2.text(1.0, -1.5, "x")
ax2.set_title('Rest Frame of Rocket')

beta = 0.75
gam = 1/np.sqrt(1-beta**2)
cdtp = cdt/gam
ax2.plot([3],[0.5],'ro')
ax2.plot([3],[0.5+cdtp],'bo')

plt.show()
```
```{note}
Figure 3.3 -- Spacetime diagrams of two events in two different
reference frames in relative motion.  The left diagram represents the
reference frame at rest with respect to the Earth.  In this reference
frame, the rocket leaves at the red dot on the lower left, travels 4.3
LY and arrives at the blue dot 5.6 years later.  The primed frame in
this case is the rest frame of the rocket.  The right diagram shows
the same two events from the primed frame: the Earth falls away behind
the rocket while Alpha Centauri comes to meet it, while the rocket remains
in the same place.
```

We can therefore write a displacement four vector in the reference
frame at rest with respect to the Earth:
\begin{equation}
[dx_4]_\oplus =
\begin{bmatrix}
icdt_\oplus\\
dx\\
0\\
0
\end{bmatrix}
 =
\begin{bmatrix}
i5.6~{\rm LY}\\
4.2~{\rm LY}\\
0\\
0
\end{bmatrix}
\end{equation}

The observer aboard the space craft experiences these same two events:
leaving the earth and arriving at Alpha Centauri. They both seem to
occur just outside of the window of the space craft. There is no
displacement between these two events! $dx_{\rm rocket} = 0$. The observer
aboard the space craft will measure some time interval, $dt_{\rm rocket}$
between the time he sees the earth outside of his window and sees
alpha centauri outside the window. The observer aboard the space craft
measures a displacement 4-vector:
\begin{equation}
[dx_4]_{\rm rocket} =
\begin{bmatrix}
icdt_{\rm rocket}\\
0\\
0\\
0
\end{bmatrix}
\end{equation}
The Michaelson-Morely experiment demands that the size of these two
four-vectors be the same.  Clearly, the observer on the spacecraft will
measure a different duration for the trip than the clocks in the
reference frame at rest with respect to the Earth will!  Note this has
nothing to do with the first observer being physically **on** the Earth --
the $dt_\oplus$ is measured by the infinite lattice of clocks at
rest with respect to the Earth.

The sizes of the two four vectors are:
\begin{equation}
dx_\oplus^2 - (cdt_\oplus)^2 = -(cdt_{\rm rocket})^2
\end{equation}
Solve for $dt_{\rm rocket}$:
\begin{equation}
dt_{\rm rocket} = dt_\oplus\sqrt{1-\left(\frac{dx_\oplus}{cdt_\oplus}\right)^2}
= dt_\oplus \sqrt{1-\beta^2}
\end{equation}
Note that the quantity under the square root will always be less than one,
so the time interval on the rocket will always be less than the time
interval in the Earth's reference frame!

```{margin} Example 3.2

The half-life for a muon is $2.20 \times 10^{-6}$ seconds as ineasured
in the rest frame of the muon. A particle accelerator, in the
laboratory reference frame, produces a beam of a large number of muons,
traveling at a speed $v = 2.900 \times 10^8$ m/s as observed in the
laboratory. Find the half-life of the muons as measured by the
observer in the laboratory frame of reference.

The value of $dt_0 = 2.20\times 10^{-6}$ s. What I am looking for is $dt$.
Equation {eq}`eqtimed` relates the
two. What I need is to figure out a value for $\beta$.  Since $\beta
= v/c,$ it must be $2.900\times 10^8$ m/s/($3.000\times 10^8$ m/s)=.9666.
The Lorentz Factor is therefore $\gamma = 3.906$, so
$dt = \gamma dt_0 = 8.59\times 10^{-6}$ s
```

This is the famous **Time Dilation**, which is sometimes characterized as
"moving clocks run slow", but I find this phrase confusing, as all clocks
could be in motion, relative to something else.  A more precise formulation
is "a clock in a reference frame at rest with respect to two events will
measure the shortest possible time interval between those two events."
This time interval, being unique, also gets a name, and is called the
**proper time interval**, and is often designated $dt_0$.
If you are not clear on which clock is moving, draw a spacetime diagram and
see whether the two events you are considering are in the same place.
Often, the confusion arises because someone is not carefully considering
what the two events in question actually are, and mistakenly compares
two different time intervals.

The quantity $\sqrt{1-\beta^2}$ comes up so often it has a name: the
Lorentz Factor, although it is usually more covenient to move it to the
other side of the equation, since the clock on the rocket is the one
at rest and therefore $dt_{\rm rocket}$ is the proper time interval.
In this case,
```{math}
:label: eqtimed
dt_\oplus = \frac{dt_{\rm rocket}}{\sqrt{1-\beta^2}} = dt_{\rm rocket}\gamma,
```
where $\gamma \equiv 1/\sqrt{1-\beta^2}$ is the Lorentz Factor.
The Lorentz Factor is always greater than or equal to one (equality
in the case of $\beta=0$ -- no relative motion means the clocks will agree),
so the clocks in the Earth's frame will always measure a longer time
interval, no matter the speed of the rocket.

In the case presented here, $\beta = 0.75$, so $\gamma = 1.512$,
and therefore the clocks on the rocket will measure a time interval
of 3.7 years.

### Muons in the Atmosphere

Muons are created in the upper atmosphere (around 100 km up) of the
Earth when cosmic rays from space collide with atoms in the air.
These muons continue toward the Earth at high speed.  However, the
half-life for muons at rest is only about 2.2 microseconds.  By the
time $22~\mu$s have passed, the vast majority of them will have
decayed.  Even at the speed of light, they could only travel about 6.6
km in $22~\mu$s.  There should not be any muons reaching the ground,
and yet you can see many of them in a cloud chamber at sea level.
Even close to the speed of light, they would need at least 330
microseconds to get down here.

```{margin}
Equation {eq}`eqlorentz` is a very useful relationship to know.
```
Much as in Example 3.2, the halflife of the muons needs to be dilated
in the reference frame at rest with respect to the Earth.  For the
muons to make it to the Earth's surface, the time dilation factor
would have to be $\sim 330/2=115$.  We can solve for $\beta$ from the
Lorentz Factor:
```{math}
:label: eqlorentz
\gamma^2 = \frac{1}{1-\beta^2} \rightarrow \beta =
\sqrt{1-\frac{1}{\gamma^2}}.
```
Plug in 115 for $\gamma$ and get $\beta
\geq 0.999962$ -- this is the minimum speed the muons would have to be
going if they are to reach the ground before they decay.  The
fact that this actually happens is further solid evidence to support
that time dilation is real.

The spacetime diagrams for these events would look just like Figure
3.3, only the $x$ axis would stand for the height above the ground,
pointing down.  In that case, the left diagram would be the Earth's
reference frame, while the right diagram would be the muon's reference
frame.  The red dot would be the creation of the muon in the upper
atmosphere, while the blue dot would be the muon reaching the ground.
The vertical separation between the dots on the right would be $cdt_0
= 6.6$ km, while the horizontal separation in the diagram on the left
would be the 100 km height of the atmosphere.  The vertical separation
in the left diagram would be $cdt = 100$ km (technically 100 km
divided by 0.999962 -- therefore the ratio of $dx$ to $cdt$ would be
$\beta$, as it should be.).

Writing all those nines is a little boring, so usually when the speeds
start to get close to the speed of light, we just write the speed **as
the Lorentz Factor**, rather than converting it back to $\beta$, or
even $v$.  So it would be perfectly acceptable to just say something
like "travelling with a speed of $\gamma = 115$", even though technically
$\gamma$ is not a speed.

## Summary


