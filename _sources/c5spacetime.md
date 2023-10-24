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

# Chapter 5: Properties of Spacetime Diagrams

Now that we have the tools of the displacement four vector, the
spacetime diagram, and the Lorentz transformation, it is worth taking
some time to look at what we can learn about the properties of space
and time by applying these tools to the spacetime diagram.

## Summary so Far

So far in our journey through Special Relativity, we have been
considering the relative displacement in spacetime of two events, as
measured in two different inertial frames of reference that are in
relative motion with some constant velocity $v_R$ (or $\beta_R$).  We
draw these events on a spacetime diagram, such as in
{numref}`fig2events` -- the left side shows two events at rest in a
reference frame S, while the right side shows where those same two
events would land in a frame S' which is moving to the right at speed
$v_R$ relative to S.  Since S' is moving right, the events will be
measured as moving to the left, separated by some displacement $dx'$
and some time interval $dt'$ (which is different from the time
interval $dt_0$ in the rest frame).

```{code-cell}
:tags: ["remove-cell"]
dx=4.2
cdt = 5.6
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
fig.suptitle('Two Events in Two Reference Frames')
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
ax2.set_title('Frame in Relative Motion')

ax2.plot([4.5],[0.5],'ro')
ax2.plot([4.5-dx],[0.5+cdt],'bo')

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
ax1.set_title('Rest Frame of the Events')

beta = 0.75
gam = 1/np.sqrt(1-beta**2)
cdtp = cdt/gam
ax1.plot([3],[0.5],'ro')
ax1.plot([3],[0.5+cdtp],'bo')
glue("twoevents", fig, display=False)

```

```{glue:figure} twoevents
:figwidth: 800px
:name: fig2events

Two events at rest in an unprimed frame and in a primed frame
moving to the right (which means the events displace to the left) with a speed
of $\beta_R = 0.75$.
```

For two events at rest, we call the time interval between them the proper time, $dt_0$,
and the square of the displacement four vector is $-c^2dt_0^2$.  Because the
size of a four vector is invariant through a Lorentz transformation, this must
be the size of the displacement four vector in the S' frame as well, which means
```{math}
-c^2dt_0^2 = dx^{'2} - c^2dt^{'2}
```
Factor out the $cdt'$ and divide through by $c^2$ to get
```{math}
dt_0^2 = dt^{'2} \left(1 - \left(\frac{dx'}{cdt'}\right)^2\right)
```
But we take $dx'/cdt'\equiv-\beta_R$, so
```{math}
dt_0^2 = dt^{'2} (1 - \beta_R^2)
```
and if we take $\gamma_R^2\equiv 1/(1-\beta_R^2)$, then
```{math}
:label: timedilation5
\boxed{
dt' = \gamma_R dt_0
}
```
which we call **Time Dilation**, and that is why the blue dot on the
right is higher than the blue dot on the left.

What I wish to do in this chapter is to apply this kind of analysis to sets
of events, to present patterns in spacetime that will help you improve your
intuition about the implications that SR demands, which are so contrary to
our day to day experience.

## Types of Intervals

In the list of principles that SR demands is the insistance that an
observer can never "just know" what is going on somewhere else.  Since
in many frames of reference two events are separated by some spatial
displacement $dx$, it is worth considering just how the observer is to
know that this displacement is indeed $dx$ and not some other number.
In the primed frame shown in {numref}`fig2events`, an observer cannot
move from the red dot to the blue dot -- such an observer would not in
fact be **in** the primed frame.  Such an observer would be at rest
with respect to the unprimed frame, and the blue dot would be at the
observer's location.  If the observer is at the location of the red
dot in the primed frame they cannot be at the location of the blue
dot.

How, then, is our primed observer to measure the distance to the blue
dot, if they cannot walk over there?  One imaginative method that is
often used is to claim that the observer waits until everything is
over, and then they collapse the whole lattice of rulers and clocks
and piece together measurements of $dx'$ and $cdt'$ after the fact
from the records of these clocks.

While that would certainly work in principle, it's not very
emotionally satisfying.  Breaking down an infinite lattice of infinite
clocks takes a long time.  Is there a way an observer could determine
$dx'$ in a more immediate fashion?  Indeed there is.  The observer
could send a pulse of light from their location (in the same place as
the red dot, but not necessarily at the same time) in such a way that
it reflects off a strategically placed mirror at the same location as
the blue dot (at the exact moment represented by the location of the
blue dot on the $cdt'$ axis), and returns to the location of the
observer at some later time.

Such a setup is illustrated in {numref}`figinterval`.  The worldline
of the (stationary) observer is a vertical red line.  The red dot is
stationary along this line, while the blue dot is off to the side.
Orange lines represent the worldlines of the light that the observer
must send and receive if they are to get information about the blue
event back to their location.  There is a slider at the bottom of the
diagram that allows you to move the blue dot up and down, relative to
the red dot.  Note that "up and down" means "earlier or later" in
time.



```{code-cell}
:tags: ["remove-input"]
# A VPython tool to show types of intervals
url1 = "https://glowscript.org/#/user/dasmith/folder/Public/program/SRintervals"
geroch = display.IFrame(src=url1,width=880,height=700)
glue("intervalfig",geroch, display=False)

```

```{glue:figure} intervalfig
:figwidth: 800px
:name: figinterval

Interactive spacetime diagram.  The red line represents
the worldline of a stationary observer.  The blue dot represents some
event $p$ on that worldline.  The green dot represents some other event
$q$ that most definitely does not reside on the worldline with $p$.
To get information about $q$, therefore, the observer must send a
light ray out to $q$ and get the reflection back.  The worldlines of
these light rays are shown in orange.  The speed of light is taken to
be 1.  The observer can therefore
define two time intervals that represent the elapsed time between
event $p$ and the events when the light was emitted and received.
From $t_1$ and $t_2$ the observer can calculate a $\Delta x$ and
a $\Delta t$, as described in the text.  These four values as
well as the interval you get from the displacement four vector
are shown on the diagram.  Move the slider to move $q$ up and down
relative to $p$, and click the box to turn light cones for $p$ on
and off.  
```

Given such a setup, the observer can define two time intervals, which
we will call $t_1$ and $t_2$.  The first is the time from the red dot
until the light returns, and the second is the time from the moment
the light is sent out until the time of the red dot.  These time
intervals are represented in {numref}`figinterval` by a white and a
magenta arrow, respectively.  The value of $t_1$ is positive if the
red dot happens before the light returns, and the value of $t_2$ is
positive if the light is sent out before the red dot.  These two
values, $t_1$ and $t_2$, are printed out on the diagram, and you can
see by moving the slider that if you shift the order of the red and
blue events, either $t_1$ or $t_2$ (but not both!) will switch to
negative.

Given these two numbers and our knowledge about the speed of light, we
can help our observer calculate values for $dx'$ and $cdt'$, without
ever going over to the blue dot!  If the red and blue dots are
simultaneous, then $dt'=0$, and $t_1$ must equal $t_2$ (this is the
initial setup for {numref}`figinterval`, and you can see the symmetry
yourself).  The later the blue dot shifts, the larger $dt'$ should
get, and the earlier the blue dot shifts, the more negative $dt'$
should get.  Therefore, our observer can conclude that the temporal
displacement between the red and blue dots is
```{math}
:label: cdtp5
dt' = \frac{t_1-t_2}{2}
```

The total time to go out and back is $t_1+t_2$.  At the speed of
light, the distance traveled would be duration times speed, so
```{math}
:label: cdxp5
dx' = c\frac{t_1+t_2}{2}
```
These numbers are also displayed in {numref}`figinterval` (without the primes),
using units where $c=1$, for simplicity.
The numbers in the diagram have been calculated from $t_1$ and $t_2$,
not measured directly from the graph.

```{margin}
To show that $t_1t_2$ is the same as the interval, multiply
Equation {eq}`cdtp5` by $c$ and then add it to Equation {eq}`cdxp5`:
$$cdt'+dx' = ct_1$$ (the $t_2$ factors will cancel).  Then
subtract {eq}`cdxp5` from {eq}`cdtp5` to cancel the $t_1$ terms
and get
$$dx'-cdt' = ct_2$$
Multiply the two to get
$$c^2t_1t_2 = (dx'-cdt')(dx'+cdt')$$
which is
$$c^2t_1t_2 = dx'^2-c^2dt'^2$$
and the latter is just the square of the displacement four
vector that we have been using.  So in units of $c=1$, the
product of $t_1$ and $t_2$ is the interval!
```

The interval of the displacement four vector, $dx'^2-c^2dt'^2$ is also
shown to the right of the red dot.  What is interesting is that we can
take Equations {eq}`cdtp5` and {eq}`cdxp5` and solve them for $t_1$
and $t_2$, and if you multiply $t_1$ by $t_2$, you can show this
equals $dx'^2-c^2dt'^2$!  The product $t_1t_2$ (shown in the diagram
to the right of the blue dot) is just another way of writing the
square of the displacement four vector!

Try sliding the blue dot up and down and verify that these two ways of
writing the interval are always the same (to within possible rounding
errors).  The lesson here is that it **is** possible to measure the
displacement between two events without an observer actually going from
one event to the other (which can only happen in the rest frame of the
observer).

Both these ways of writing the interval suggest three possible ranges
of interest that an interval might fall into.  As you move the slider
back and forth, you should be able to identify these three possibilities.
Either both $t_1$ and $t_2$ are positive, in which case the interval is
positive, or one of the two is positive and the other negative, in
which case the interval is negative, or one of the two is zero, in
which case the interval is zero.  These three options correspond to,
in four-vector notation, to $dx'>cdt'$, $dx'<cdt'$, or $dx'=cdt'$,
respectively.

It is useful to classify a displacement four-vector into one of these
three categories, because there are specific properties that each of
these types of displacement four-vector have.  Move the slider such
that the interval displays zero.  In this case, the displacement from
red to blue will be one of the orange lines -- the worldline of the
light that travels either out from or back to the observer's
worldline.  The displacement between red and blue in this case must be
just like the displacement that light would follow, so this kind of
interval is called a "lightlike interval".  For a lightlike interval,
either $t_1$ or $t_2$ is zero (depending on whether the red or the
blue event happened first), which means that $dx'=cdt'$, or
$dx'/dt'=c$, which means anything moving along that worldline has to
be moving at the speed of light.

If the interval is negative, that means that $dx'<cdt'$.  In the sum
of squares that represent the size of the four vector, the time
component has a minus sign, up to and including the extreme case of
the proper interval $-c^2dt_0^2$.  A negative interval is therefore
more like a time displacement than it is like a space displacement, so
such intervals are called "timelike".  As you might probably guess at
this point, displacements where $dx'>cdt'$ have positive intervals and are
therefore called "spacelike".

As you will see as this chapter develops, it is useful to group
displacement four-vectors in these categories, because all timelike
intervals share certain properties that it is useful to remember, as
do all spacelike intervals.  All lightlike intervals are in some sense
the same, as they all equal zero, and they all have to equal zero in
all reference frames, although of course the time and space components
can change individually.

For example, consider all the possible events that are lightlike
displaced from the red dot.  These are all events that lie along
diagonal lines that cross at the red dot and make $45^\circ$ angles
with the horizontal (if $dx'=cdt'$, then the slope is 1 and it makes a
$45^\circ$ angle).  If the red dot sends out light, that light will go
up and away from the red dot at a $45^\circ$ angle.  Any event that
sends light to the red dot must lie below the red dot at a $45^\circ$
angle.  All the lightlike intervals that connect to the red dot
therefore make an $\times$ across this diagram.  However, this set of
events is actually referred to as a "light cone".  Why a cone?
Because if we do include the $y$ dimension as pointing into the
computer screen, then the $\times$ can be rotated around the vertical
axis, and instead of an $\times$, we get a cone.  Click the button on
{numref}`figinterval` to see the light cones associated with the red
dot.

Of course, it's only a cone if we include two space dimensions, $x$
and $y$.  If we could include $z$, the light would be travelling in a
sphere, either expanding out from the red dot or collapsing to it.  We
can't make a four-dimensional graph, though, so we represent the
sphere as a cone, and the term "light cone" has stuck.

```{warning}
We will always talk about the "light cones" associated with any event,
but please remember that it's really an expanding sphere in 3D space.
```

If any real thing wanted to get from the red dot to a point on the
upper light cone, or from a point on the lower light cone to the red
dot, this thing would have to travel at the speed of light to do so.
However, for any event **inside** these cones, it would in principle
be possible to get to or from the red dot without hitting light speed.
Therefore, the set of events inside the upper light cone are all the
events on which the event at the red dot could **possibly** exert any
kind of influence.  We therefore call the events inside the cone the
"future" of the red dot.  All the events in the lower cone could
**possibly** influence what happens at the red dot, so we call this
set of events the "past" of the red dot.  Every single event has
its own light cones, and therefore its own set of past and future.

Events that lie outside these light cones are neither past nor future,
but some "other" that we don't have a good name for.  Events that are
spacelike separated cannot influence each other, and one cannot be
either cause or effect for the other.  It takes light from the Sun
about eight minutes to get here.  If the Sun were to somehow vanish,
that event would be outside the light cone of the Earth right now.  It
would take eight minutes for the light cone of that event to intersect
the world line of the Earth, and only then would the horrific darkness
and bitter cold ensue.  So enjoy your eight minutes!


## How do the Axes Change?

Consider the vertical and horizontal axes in a spacetime diagram.
Each point on either axis is itself an event, with a coordinate of
$x=0$ or $t=0$.  So an axis on a spacetime diagram is the set of
events that has one of the two coordinates be zero.  The origin, of
course, is the event where both coordinates are zero.  One of the
principles described in Chapter 2 is that all events in spacetime have
to be observed in all frames.  Switching reference frames can't just
make something not happen.

Therefore, we can ask how the coordinates of the set of events that
make up the axes of the unprimed frame change when we shift into the
primed frame.  Or, to turn it around, which events will land on the
axes in the primed frame, and where were they originally in the
unprimed frame?  I can explain this conceptually for the vertical
axis, and then we will see what the math for the Lorentz
transformation tells us.

The vertical axis is all the events that occur at $x=0$.  A vertical
worldline represents an object at rest.  Let's say I am an observer at
rest at the origin, so my worldline is the vertical axis.  If I switch
into a primed reference frame moving to the right at $v_R$, each
successive event will shift further to the left, as I pass further and
futher to the right.  So the original axis will tilt left in the
primed frame.  Furthermore, as I walk to the right, successive events
that were previously to my right will now move right next to me.  If I
draw my worldline on the original spacetime diagram, it will be a
rightward tilting line, and so this line will be the vertical axis in
the new frame.

Therefore, if I want to draw the set of events on the unprimed frame
that will be the vertical axis in a primed frame, I will draw a tilted
axis with a slope of $1/\beta_R$.  My speed is distance over time, but
the spacetime diagram puts time on the vertical axis, so "rise over
run" implies a slope of one over speed.  The faster the relative speed
of the primed frame, the more tilted the new axis line will be, up to
a speed (and therefore slope) of 1, since I can't go faster than
light.

```{warning}
I stress that this axis is only tilted in the original, unprimed,
frame.  If I were to redraw the spacetime diagram in the primed frame,
the axes would be horizontal and vertical, as normal.  The set of
events in the unprimed frame that will be vertical in the primed frame
make up a tilted line in the unprimed frame.  {numref}`fig16events` shows a case
where the axes are implicitly redrawn every time the reference frame
is shifted, while {numref}`figaxestilt` is keeping the same original axes while
showing which events would be on the horizontal and vertical axes, if
you were to redraw them.
```

We can show this intuitive prediction mathematically by performing
a Lorentz transformation.  We can ask, of all the points $(ict,x)$
in the spacetime plane, which ones will end up having $x'=0$ in the
primed frame?  The transformation looks like:
```{math}
:label: axistrans5
\begin{bmatrix}
ict'\\
x'
\end{bmatrix}
=
\begin{bmatrix}
\gamma_R & -i\beta_R\gamma_R\\
i\beta_R\gamma_R & \gamma_R
\end{bmatrix}
\begin{bmatrix}
ict\\
x
\end{bmatrix}
=
\gamma_R
\begin{bmatrix}
i(ct-\beta_R x)\\
x-\beta_R ct
\end{bmatrix}
```
According to Equation {eq}`axistrans5`, the set of points that will be
at $x'=0$ will have the condition that $x=\beta_R ct$.  Since the time
axis is vertical, this equation describes a straight line through the
origin with a slope of $1/\beta_R$, as predicted by intuition, above.

However, Equation {eq}`axistrans5` lets us go further, to ask a
question that is much harder to describe intuitively: which events
will land on the horizontal axis in the primed frame?  These are the
points where $t'=0$, and we can find them by setting the time
component of Equation {eq}`axistrans5` equal to zero.  In that case,
$ct=\beta_R x$, and that is the equation of a straight line through
the origin with a slope of $\beta_R$.  The primed horizontal axis
will be an upward tilting line that shift to a steeper slope as you
increase $\beta_R$, up to 1, since that is as fast as you can go.

The results of these two equations are depicted graphically in
{numref}`figaxestilt`.  You can shift the slider below the figure to
change $\beta_R$.  The yellow arrows that represent which events will
be on the axes in the primed frame will tilt as you change the
relative speed.  Again, the actual unprimed axes remain unchanged --
the tilted axes are telling you which events will end up on the
perpendicular axes in the primed frame, if you were to redraw them.

```{code-cell}
:tags: ["remove-input"]
# Interactive spacetime diagram to allow the user to tilt the axes
url1 = "https://glowscript.org/#/user/dasmith/folder/Public/program/SRlorentzaxes"
tilting = display.IFrame(src=url1,width=880,height=650)
glue("tilttheaxes",tilting, display=False)

```

```{glue:figure} tilttheaxes
:figwidth: 800px
:name: figaxestilt

A spacetime diagram where the events that lie on the
axes are indicated with yellow arrows.  A slider below the diagram
allows you to change the relative speed of a primed reference frame,
and the arrows will tilt to indicate which events will end up on the
axes of the primed frame, were you to redraw the diagram in the new
frame.  Note that the sets of events tilt together as you increase the
relative velocity to the right.  The original axes remain where they
were, indicated by white arrows.
```

## Simultaneity is Relative

The behavior of the vertical axis is easy to understand just by
imagining someone walking.  The faster you walk, the further apart in
space the events that you pass by will be, and therefore your worldline
(which is by definition the vertical axis in the primed frame) will
tilt further and further.

Consider the event represented by the blue dot in
{numref}`figaxestilt`: as the diagram is first drawn, with the
observer at rest in the unprimed frame, the event is to the right of
the vertical axis.  As you move the slider to the right and the axis
tilts, you can choose a speed such that the yellow arrow lies on the
blue dot.  This represents walking just fast enough that you get from
the origin to the location of the blue dot just as it happens.  If you
walk faster, you will pass the blue dot before it happens, and the
event will happen to the left of the axis in the primed frame.

It is therefore possible, by altering the relative speed, to choose a
reference frame in which the blue dot event is either left, right, or
at the same location as the origin.  This amounts to simply choosing a
speed so that you fall short, overtake, or precisely reach the event
as it happens.  Since you can get to the blue dot from the origin by
moving at a speed less than $c$, the blue dot lies inside the light
cone of the origin, and the interval between the origin and the blue
dot is timelike.  This leads to the conclusion that through a careful
choice of reference frame, a later event can be right, left, or at the
same position in space as an earlier one, if the two are timelike
separated.

This hopefully seems rather intuitive.  You have overtaken or reached
events many times in your life, so hopefully it is not hard to imagine
the implications of shifting the vertical axis as shown in
{numref}`figaxestilt`.  However, it is much harder to imagine how the
red dot interacts with the horizontal axis.  The mathematical
description is almost identical to the case of the blue dot, only
rotated $90^\circ$.  I will repeat what I said above, in only slightly
different words: through a careful choice of reference frame, an event
to the right can be before, after, or at the same time as an event to
the left, if the two are spacelike separated.

Play with the slider and convince yourself that this is accurate.  By
changing the relative speed, you can put the red dot above the
horizontal axis (the red event happens after the event at the origin),
below the horizontal axis (the red event happens before the event at
the origin), or on the horizontal axis (the red event happens at the
same time as the event at the origin).  This is very hard to accept,
but this means that simultaneity is relative.  Whether events happen
at the same time, or in which order they happen, depends on your
choice of reference frame, as long as the events are spacelike
separated.  The fastest you can go is $\beta_R=1$, which is a line
with a $45^\circ$ angle slope, so you can never go fast enough to get
the horizontal axis to reach, say, the blue dot.

One of the reasons this may bother you is it may seem like this may
contradict causality.  If I can arbitrarily switch the order of events by
choosing a different reference frame, can I make effects happen before
causes?  Thankfully, no.  You can only change the temporal sequence of
events that are spacelike separated, and events that are spacelike
separated can never be cause or effect for each other, because
something would have to move faster than light to carry the influence
from one event to the other.  Nature preserves causality, but what
you think of as "simultaneous" will depend on what reference frame
you are in.

## Conclusions from Sets of Events

To sum up what we have learned about spacetime from exploring how
these diagrams change under Lorentz transformations, I have created
one more interactive diagram, shown in  {numref}`fig16events`.  This diagram
also has a slider that lets you change the relative speed of a primed
reference frame, but every time you change it, the program redraws the
diagram with the new axes.

```{warning}
The axes in {numref}`fig16events` are always labeled $x$ and $ct$,
even though as soon as you change $\beta_R$ away from zero, they
strictly should change to $x'$ and $ct'$.
```

On this diagram are shown sixteen spheres to represent a set of events.
In the original, unprimed frame, the events make up a square, so you
could think of them as four people standing still, snapping their fingers
four times in unison.  Each vertical column of four events are happening
at the same place, and each horizontal row of events are happening at
the same time.

```{code-cell}
:tags: ["remove-input"]
# Spacetime diagram to show how event shift under transformation
url1 = "https://glowscript.org/#/user/dasmith/folder/Public/program/SRstgrid"
strch = display.IFrame(src=url1,width=880,height=650)

glue("stretchevents",strch, display=False)

```

```{glue:figure} stretchevents
:figwidth: 800px
:name: fig16events

Interactive spacetime diagram that shows 16 events.  The
default is to have the events in a grid.  You can think of this as
four objects at rest, separated by some distance along the $x$ axis,
with four events along each vertical worldline.  There is a slider to
change a relative velocity $\beta_R$.  However, what happens when you
change the slider is that the program applies a Lorentz transformation
to each event and then plots the 16 events in a new reference frame.
Each time you change the slider, you are changing to a new reference
frame.  Certain points are marked in color, and the implications of
how they shift are discussed in the text.
```

As you change the relative speed of the primed frame, these sixteen
events will shift around according to the Lorentz transformation.
This is the opposite operation to pinching the axes in
{numref}`figaxestilt` -- here, the points between the axes are being
stretched.  Play with sliding the marker back and forth and see how
the events move.  You could imagine this as the "pinched" axes in
{numref}`figaxestilt` being streched back to perpendicular, and
therefore all the events between them being pulled outward, too, as if
they were on a rubber sheet.

```{note}
It is useful
to start thinking about space and time as being able to be stretched
and squeezed, as this is an important feature of General Relativity.
```

I have marked several spheres with color, and made three of them leave
trails behind, to draw your attention to them and gain particular
insights.  The cyan sphere marks the origin, and note that it does not
move.  If all the components of a four vector are zero, the Lorentz
transformation will not change them.  You can use the cyan sphere as
an anchor by which you can compare the locations of the others.

The orange event on the time axis will move left or right (since it is
timelike separated from the origin) as well as up, leaving a
bowl-shaped trail behind.  This is Time Dilation.  The shortest
possible time interval will be measured in the frame where the two
events (cyan and orange) are at rest.  If you zoom in on the orange
sphere, you will see that the trail gets flatter and flatter, the more
you zoom in.  This is a representation that for small $\beta_R$, the
Lorentz factor gets close to one, in which case $dt'\approx dt_0$.  If
the intervals are equal, the trail would be flat, and the closer you
zoom in, the flatter it gets.

The red spheres are all lightlike separated from the origin, and you
can see that when you move the slider, they **stay** along the
$45^\circ$ line through the origin.  Since the speed of light is the
same in all reference frames, events that are lightlike separated will
be lightlike separated in all reference frames.  All the spheres above
the red spheres are timelike separated from the cyan sphere, and all
the spheres below the red spheres are spacelike separated from the
cyan sphere.  You can see that no matter how you move the slider, you
can never change a timelike interval into a spacelike interval, or
vice versa.

The blue sphere represents an event that is simulataneous with the
cyan event at the origin in the original frame.  By moving the slider,
you can see that the temporal order of these events can be switched.
In fact, pick any two events you like that are separated by a spacelike
interval, and you can find ranges of $\beta_R$ where the order is switched.
What counts as "simultaneous" depends on your frame of reference.

Finally, note that the blue event leaves a trail that looks just like
the orange trail, only rotated $90^\circ$.  You might think from this
similarity that there must also be a length dilation to match the
time dilation illustrated by the orange trail.  This is not the case.
Instead, we talk about a length **contraction**.  How can this be?
The difference lies in what we mean by "length."

Imagine that the leftmost and rightmost column of events lie upon
vertical worldlines that represent the left and right ends of an
object at rest.  Then the magenta sphere and the blue sphere are in
the same place, and are separated from the cyan sphere by the same
spatial distance.  However, the cyan and blue events are simultaneous,
so we define the length of the object by the locations in space of
these events.  The magenta event is the same distance away, and therefore
also measures the length of the object.

If you increase the relative speed $\beta_R$ to about $0.67$, you will
see that although the blue event is much further away, it also happens
much earlier than the cyan event.  The spatial displacement between
blue and cyan can no longer be considered to be the length of the
object!  In this reference frame, the object is moving, and if you
consider the distance between the back end of a moving car ten seconds
ago and the front end of the car now, then you are including the
distance driven by the car during those ten seconds as part of the
length of the car, which is not the way we usually would think of the
length of the car.

Instead, we must consider the magneta event, which is simultaneous in
**this** reference frame with the cyan event.  Since it was at the end
of the object in the original frame, it must also be at the end of the
object in this frame, so in **this** frame, we would consider the
"length" of the object to be the spatial displacement between the cyan
and magenta events, which you can see from the trails is smaller than
the spatial displacement between cyan and blue in the original frame
of reference.  Lengths contract.  We will work out examples of length
contraction more thoroughly in the next chapter.

```{note}
It might be useful to look at [this interactive
version](https://alexonscience.com/projects/spacetimeglobe/) of
{numref}`fig16events` -- it lets you add your own events and see how
they shift around under Lorentz transformations.
```

## Hyperbolic Rotation and Rapidity

It might, when you move the slider in {numref}`figaxestilt`, remind
you of rotation.  Certainly there is rotation going on -- the axes are
swinging around the origin.  But it's also clearly not the kind of
rotation you're used to.  When we usually rotate the $x$ and $y$
coordinate axes around the $z$ axis, the two arrows move together.
However, in {numref}`figaxestilt` , the arrows are moving
symmetrically in opposite directions.  In the familiar kind of
rotations, you can keep swinging the axes around and around and
around, increasing the angle as far as you like.  However, in a
spacetime diagram, there is an asymptotic limit to how far the axis
will swing -- they will both move toward a slope of 1.

This is still a rotation, but it's called a hyperbolic rotation.  You
can even mathematically make it look like a rotation.  Recall that when
we rotated the coordinate axes around the $z$ axis, we determined how
that would affect a vector by multiplying the vector by a matrix:
```{math}
:label: normalrot
\begin{pmatrix}
v_x'\\
v_y'\\
v_z'
\end{pmatrix}
=
\begin{pmatrix}
\cos{\theta}&-\sin{\theta}&0\\
\sin{\theta}&\cos{\theta}&0\\
0&0&1
\end{pmatrix}
\begin{pmatrix}
v_x\\
v_y\\
v_z
\end{pmatrix}
```
When you multiply that out, you
get
```{math}
:label: rotang
\begin{pmatrix}
v_x'\\
v_y'\\
v_z'
\end{pmatrix}
=
\begin{pmatrix}
v_x\cos{\theta}-v_y\sin{\theta}\\
v_x\sin{\theta}+v_y\cos{\theta}\\
v_z
\end{pmatrix}
```
Note that the $x$ and $y$ components of the vector get "mixed up" when you
rotate the coordinate system.  By rotation, you are turning part of $x$ into $y$
and vice versa.

The periodic nature of the sinusoidal functions corresponds to the angle
being able to go around and around and around.  To have the angle approach
an asmptote, as in a spacetime diagram, we need the hyperbolic trig functions,
$\tanh$, $\sinh$, and $\cosh$.  In a normal $x-y$ plane, if you had a vector
in that plane, the angle the vector would make with the $x$ axis would be
$\tan{\theta} = v_y/v_x$, opposite over adjacent.  However, for spacetime,
it's the hyperbolic tangent: $\tanh{\phi} = dx/cdt = v/c = \beta$.  Note
that we put $x$ over $y$ to get $\beta$.

There's a trig identity that says 
```{math}
:label: tanh
\tanh^2{\phi} = 1 - \frac{1}{\cosh^2{\phi}} = \beta^2
```
But we already know that $\beta^2 = 1 - 1/\gamma^2$, from the definition
of $\gamma$, so $\cosh{\phi} = \gamma$!

Finally, we know $\cosh^2+\sinh^2 = 1$, so if $\cosh{\phi} = \gamma$, then
$\sinh^2{\phi} = 1-\cosh^2{\phi}$
```{math}
:label: sinh
\sinh^2{\phi} = 1-\gamma^2 = 1-\frac{1}{1-\beta^2} = \frac{\beta^2}{1-\beta^2}
= \beta^2\gamma^2
```
So $\sinh{\phi} = \beta\gamma$.  But these are just the elements of the Lorentz
transformation matrix!  This means we can write the Lorentz matrix as
```{math}
:label: lortrapid
{\cal L}_x(\phi) =
\begin{bmatrix}
\cosh{\phi} & -i\sinh{\phi} & 0& 0\\
i\sinh{\phi} & \cosh{\phi} & 0& 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix}
```
If you compare Equation {eq}`lortrapid` and {eq}`normalrot`, you can see that
except for the hyperbolic part, and the $i$, they look very similar.
This means, in a very real but very odd sense, changing speeds actually rotates
space and time into each other, in much the same way that rotating axes rotates
$x$ and $y$ into each other.

The letter $\phi$ in these equations is called the "rapidity", and it depends on
$\beta$ as
```{math}
:label: rapidity
\phi = \ln{\left(\sqrt{\frac{1+\beta}{1-\beta}}\right)}
```
However, despite
the trig functions, don't get confused and think that $\phi$ is an angle in the
regular sense of how you think of that term.  In Equation {eq}`normalrot`, the letter
$\theta$ corresponds to the angle through which we rotate the axes of the diagram.
However, in the case of Equation {eq}`lortrapid`, if $v\rightarrow c$, then
the axis will pinch to an angle of $45^\circ$, while $\phi\rightarrow\infty$
(the $1-\beta$ in the denominator of Equation {eq}`rapidity` goes to zero).


Why might you want to do this?  Two reasons.  First of all, it's kind
of neat to think about a "boost" (the name of the operation of
changing speed into a new reference frame) as a kind of rotation between
space and time.  Secondly, should you need to apply two Lorentz
transformations in a row, you may recall there are trig identities
that let you write the product of trig functions as a trig function of
the sum of the angles.  In this case, you can work out that two
successive Lorentz transformations, if you write them like Equation
{eq}`lortrapid`, work out to a single Lorentz transformation using the
sum of the rapidities of the two original transformations: ${\cal
L}_x(\phi_1) {\cal L}_x(\phi_2) = {\cal L}_x(\phi_1+\phi_2)$.  This
could save a lot of number crunching.

## Problems

1) According to observers in a reference frame at rest with respect to
both the Earth and a star 3.25 LY away, an astronaut on a spaceship
near the star eats breakfast at 7 am, while their spouse on Earth
eats lunch at five hours later.  How fast would an observer in a rocket
need to travel, and in what direction, to conclude that the spouse ate
lunch before the astronaut ate breakfast?

2) George Gamov, in his delightful but dated book, *Mr. Tompkins in
Wonderland*, posits a world where the speed of light is only a few
km/hr, and all these odd effects of relativity are commonplace for
people in that world.  *Need to look in the book and recreate his
argument about the murder observed from the train*

3) Show that even when you tilt the axes as in {numref}`figaxestilt`,
the properties of the Lorentz transformation are such that the new
axes are still orthogonal, even though they are clearly not
perpendicular in the figure.

4) The following questions are with regard to {numref}`fig16events`: 

a) Can you find relative speeds for which each of the events below the
red dots happen before the cyan event?

b) Can you find a speed at which any of the events above the red dots
happen before the cyan event?

c) At $\beta_R=+0.67$, which of the dots lies on the vertical axis
above the cyan dot?  That is, where was this dot originally when
$\beta_R=0$?

d) At $\beta_R=+0.67$ consider the duration of time between the cyan
event and the event from part c).  Is this time interval longer or
shorter than the duration of time (as measured when $\beta_R=0$)
between the cyan event and the orange event?  Does this make sense?

4) Use the definition of rapidity $\phi(\beta)$ (Equation {eq}`rapidity`)
to prove that $\sinh{\phi} = \gamma$.

5) This one could end up being rather annoying, so probably best to just
do it for the time component.  Show that two successive boosts, when you
use rapidity and the hyperbolic trig functions, are the same as a single
boost with the sum of the rapidities.  In other words, show that
```{math}
{\cal L}_x(\phi_2){\cal L}_x(\phi_1) = {\cal L}_x(\phi_1+\phi_2)
```

6) Note that in considering two successive boosts, you can't just
conclude that the final speed $\beta_{1+2}$ will be the sum of the two
individual boosts ($\beta_{1+2}\neq\beta_1+\beta_2$ -- if you doubt
this, consider the case where $\beta_1=\beta_2=0.75$.  Do you see how
adding those together would be a problem?).  However, you can just add
rapidities, because they can go as high as you like.  Start with
$\phi_{1+2} = \phi_1+\phi_2$ and plug in Equation {eq}`rapidity` to
get a formula for $\beta_{1+2}$ as a function of $\beta_1$ and
$\beta_2$.  You will see this formula again in Chapter 7.

7) An event is at $(ct,x)=(1~{\rm m},2~{\rm m})$ in some reference
frame.  Is the displacement from the origin to this event timelike,
lightlike, or spacelike?  How fast would another reference frame have
to travel relative to this one to observe this event on the horizontal
axis?  Perform a Lorentz transformation at this value for $\beta_R$
and verify that this is so.

8) Consider the case where a bolt of lightning creates a peal of
thunder.  You are 5 km away from this event.  Write down displacement
four vectors for the events of your seeing the lighting and your
hearing the thunder, relative to the event of their creation.  Find
their intervals and determine if they are timelike, lightlike, or
spacelike.  With appropriate scaling, you can use
{numref}`figinterval` to represent these situations.  Do your
classifications make intuitive sense, and how do they compare with the
light cones from the original event?