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

# Chapter 13: Concepts of General Relativity

Although a full treatment of General Relativity (GR) is far beyond the
scope of this book, I would like to leave you with some hints at
the major concepts you must grapple with, if you extend the theory
beyond the special case of reference frames moving with constant
relative velocities.

## The Twin Paradox Revisited

To begin our overview of GR, let us return to our old friend, the twin
paradox, and consider the following situation: we have two events in
the same location in space, one event later than the other in time.
We have two individual clocks that have been syncronized, and they are
both at the same location as the two events.  At the instant of the
first event, one of the clocks begins moving away.  You could consider
this to be the first event, if you like.  The other clock remains
stationary.  At some later time, the clock returns to the location of
the stationary clock.  This return marks the second event.  This
situation is shown in {numref}`twopaths`, where the red line represents
the worldline of the stationary clock and the blue curve represents
a possible path that the moving clock might take.  For our purposes,
the exact path the moving clock takes is unimportant, only that it
goes away and comes back.  It must also, of course, be a physically
possible trajectory: no doubling back along the time axis or moving
faster than $c$!

From what we learned about SR, we can consider the displacement
four-vector between the two events, and it would have a size given
by $-\Delta t_0^2$, the proper time between these two events, as measured
by the red clock.  However, it is not at all clear that the blue clock
will measure the same time interval.  Consider each worldline to be
made up of infinitely many, infinitely small four-diplacements between
events that are infinitesimally far apart along the colored lines.
Each step along the worldline would have its own four-displacement,
which would have its own size $-ds^2=-c^2dt^2+dx^2$.  For the purposes
of this analysis, I am going to multiply through by $-1$ and use
$ds = \sqrt{c^2dt^2-dx^2}$.  It's not strictly necessary, but it makes
the argument easier to follow.

The total interval along each worldline would then be the integral
of all the $ds$ intervals for each of the steps.
```{math}
:label: intint
\Delta S = \int_{\rm path} ds = \int_{\rm path} \sqrt{c^2dt^2-dx^2}
```
Factor out the $dt$ to make it easier to understand:
```{math}
:label: intint2
\Delta S = \int_{\rm path} cdt\sqrt{1-\frac{dx^2}{c^2dt^2}} = \int_{\rm path}
cdt\sqrt{1-\frac{v^2}{c^2}}
```
What we have here is a path integral of some function, call it $f$,
that depends on the speed of the clock at each moment along the path:
$f = \sqrt{1-v^2/c^2}$.  For the red path, $v=0$, so $f=1$ and the
integral is trivially just the proper time.  However, along the blue
path $0\geq v^2<1$, so $0<f\leq 1$, which means the result of this
integral (which is the time interval elapsed on the blue clock, as we
are integrating over $dt$) is always going to be less than or equal to
the proper time.  In particular, for the twin paradox, if the blue path
were to involve traveling to some other star and back to the Earth,
the clock on the spaceship (the blue clock) would measure less time
than the clock on the Earth.  This is the same result we derived before
by only considering a one-way trip at constant speed.

No matter what path the blue clock takes, it will always measure less
time than the stationary (red) clock.  In this context, we can accurately
say that moving clocks run slow, because we are comparing two clocks
within the same reference frame.  When we first encountered time dilation,
we were considering the elapsed time between two events as measured by
infinite lattices of synchronized clocks, one set in motion relative to the
other set.  In this context, within a single reference frame, we can say
that a clock at rest with respect to two events will measure the LONGEST
time interval between those two events, as compared with any other clock
that follows a different path between the events.  This sounds like it
contradicts the previous formulation of time dilation, but you have to
consider carefully what is being compared.


```{code-cell}
:tags: ["remove-cell"]
# ST plot of two paths
fig = plt.figure(figsize=(5,5))
plt.arrow(-0.5,-0.5,1,0,head_width=0.1)
plt.arrow(-0.5,-0.5,0,1,head_width=0.1)

plt.plot([0.2],[-0.1],'ko')
plt.plot([0.2],[1.1],'ko')
plt.plot([0.2,0.2],[-0.1,1.1],'r-')

ax = plt.gca()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.axis([-1,2,-1,2])
ax.text(-0.55, 0.8, "ct")
ax.text(0.8, -0.5, "x")

y1 = -0.1
y2 = 1.1
w1 = 2*np.pi/(y2-y1)
w2 = 1.5*w1
y = np.linspace(y1,y2,100)
x = 0.2+0.3*np.sin(w1*(y-y1))+0.22*np.sin(w2*(y-y1))
plt.plot(x,y,'b-')


glue("twopathfig", fig, display=False)

```

```{glue:figure} twopathfig
:figwidth: 800px
:name: twopaths


Two events at the same place in a spacetime diagram, but
showing two possible paths between them.  The red vertical
line is the worldline of an object at rest, while the blue
line shows the path of a moving object that returns to the
same place in time for the second event.
```

In SR, there is no question of returning to the same location
as the first event.  This would not be possible, as to turn
around and come back would require changing the relative velocity,
and SR only considers unchanging relative velocities.  In this
new context, there is an asymmetry between the moving clock and
the stationary clock, because only the moving clock *changes*
its motion.  The only way to describe these occurances from the
perspective of the blue clock would be to allow for accelerating
reference frames.

## Accelerating Reference Frames

Imagine you are driving along in a car at a constant velocity.
Everything in the car would behave according to Newton's laws.
Objects sitting on the dashboard would remain at rest (from your
point of view -- an observer sitting on the side of the road
would see the objects moving at the same velocity as the car),
and if you threw a ball into the back seat, you would see it
move in a nice parabola.  If you had a helium balloon, it would
float upward to the length of its string.  The fuzzy dice hang
straight down from the rear view mirror.  All is calm.

Then suddenly the driver slams on the brakes.  The objects in the car
fly forward.  You feel thrown forward against the strap of your safety
belt.  The fuzzy dice swing forward toward the windshield.  From your
point of view, thinking in terms of Newton's laws, a mysterious force
has just appeared that gives all the objects in the car an additional
forward velocity.  What happens to the helium balloon?  Take a moment
to think about that and make a prediction.  After you have decided, go
to watch [this video](https://www.youtube.com/watch?v=y8mzDvpKzfY) to
see what the balloon actually does.

From the side of the road, there is no mysterious force.  All objects
continue moving forward, as their velocity would demand, but the car
slows down, as does therefore the safety strap, which pushes backwards
against you and slows you down to match the car.  The balloon is
lighter than air, so while the air in the car slides forward along
with you and the fuzzy dice, that pushes the lighter balloon toward
the back of the car.

If the driver had, instead of braking, turned the wheel hard to
the left, so that the car began moving in a circle, everything
in the car would behave as if there were a mysterious force pushing
to the right.  You would lean right, the fuzzy dice would dangle to
the right, and the balloon would lean to the left.  Although a person
on the side of the road would say everything in the car was trying
to go in a straight line, but the car was pushing to the left,
within the reference frame of the car, everything would be acting
as if a force to the right had suddenly appeared.

Let's consider that circular trajectory in the context of
{numref}`twopaths`.  What can we say about a clock that moves on a
circular path through spacetime?  Such a path is represented in
{numref}`circularSTfig`.  To get from the orange event to the purple
event, a clock at rest would just follow a vertical line, but the
black curve represents a clock that moves around in a circle.  Based
on the analysis of the previous section, such a clock would have to
measure a shorter time interval than the clock at rest (if the system
were sitting on a merry-go-round, we could consider the clock at rest
to be at the center axis of the rotating disc -- such a clock would be
synchronized with the lattice of clocks making up the reference frame
at rest with the orange and purple events, even if the clock at the
center did not pass through those two events).

```{code-cell}
:tags: ["remove-cell"]
# 3D spacetime diagram of circular motion
url1 = "https://glowscript.org/#/user/dasmith/folder/Public/program/SRspacetimecircle"
test = display.IFrame(src=url1,width=800,height=700)
glue("stcircfig",test, display=False)

```

```{glue:figure} stcircfig
:figwidth: 800px
:name: circularSTfig

3D spacetime diagram of an object moving in a circle in the xy plane.
The worldline of the object is shown as a black curve.  The object starts
at the orange event and ends at the purple event.  Both events are in the
same place, but the purple event is later than the orange event.  Rotate
the diagram to see projections of the worldline onto the different planes.
In particular, rotate the diagram to look straight down the $ct$ axis to
see the circular path in space.
```

Within the reference frame of the purple and orange events (the clock
at rest), for a given constant angular velocity $\omega$, the larger
the radius of the circle, the faster the clock will have to move to
complete the circle (according to $v=\omega r$) and therefore the
larger the discrepency in the time intervals, according to Equation
{eq}`intint2`.  A clock further out on the merry-go-round will be
moving faster and therefore measure a shorter time interval than a
clock closer to the center.  Furthermore, we also know that lengths
contract at faster speeds, so if a person at the edge of the
merry-go-round were to lay down meter sticks to measure the
circumference of the circle, they would need more meter sticks than
they would have expected without relativity.  The radius of the
circle, being perpendicular to the motion, is not affected, so $2\pi
r$ is now less than the actual circumference, according to the
(contracted) meter sticks.

The circumference being equal to $2\pi r$ is proven within a flat (or
Euclidian) space.  If space is curved, that equation no longer holds.
Consider the circle formed by a latitude line near the South Pole:
because the Earth is curved, the measured radius of the Earth along
the surface is going to be longer than the radius of the circle you
would get by slicing a flat circle along the latitude line.  Therefore,
the circumference of the latitude circle will be shorter than $2\pi r$,
because the Earth is curved.  Due to the length contraction of the
meter sticks, our scientist on the merry-go-round would also conclude
that she is in a curved space, as the geometry is not behaving according
to the rules of flat space.

The effects on time and geometry are tightly linked to the
acceleration of the reference frame of the merry-go-round.  The larger
the radius, the larger the acceleration needs to be ($a=\omega^2 r$)
to maintain the circular path.  To someone on the merry-go-round,
there would seem to be a mysterious outward (centrifugal) force,
and to move toward the center would be to move "up" against this
force, and to move away from the center would be to move "down"
with this force.  The further down they go, the larger the effects
of time dilation and length contraction would get.  They would
conclude their spacetime was deviating more and more from Euclidean
flatness.  The same chain of logic holds for the acceleration
due to gravity, because the forces are equivalent.

## The Equivalence Principle

Consider a rocket ship in deep space, far from any source of gravity.
If that rocket were moving at a constant speed, everything in the
rocket would drift around inside, because in the reference frame
moving with the rocket, everything would be moving together, and would
therefore be the same as if the rocket were not moving at all.  This
insight was first expressed by Galileo, although he talked about
sailing ships, rather than rockets.  It is only through the rocking of
the waves, he pointed out, the *deviations* from constant speed, that
we can tell the ship is moving.  On our rocket, a laser affixed to one
wall would shoot a horizontal beam straight across to the other side,
regardless of the rocket speed, because the photons depart the wall
with the same upward momentum as the rocket.  Whether that upward
velocity is zero or not, as long as it is constant, the result will be
the same.

However, if the rocket is accelerating, the objects in the rocket will
tend to "fall" to the floor.  To an observer outside the rocket,
watching the rocket accelerate on by, the objects are measured to be
moving at a constant velocity, while the floor is accelerating up to
meet them.  Until they actually touch the floor, they are just
drifting through deep space, unaware that anything else is happening.
From inside the rocket, though, the objects seem to be moving down to
the floor.  If the rocket is accelerating up at 9.8 m/s/s, it will
seem to anyone in the rocket that everything inside is experiencing a
downward force exactly equivalent to the gravity on Earth.

This insight, that free-falling objects are no different from objects
floating in deep space, was once described by Einstein as his happiest
thought.  The story goes that he was strolling through town and
saw painters working on a clock tower.  It occurred to him that if
one of the painters were to fall, the painter would effectively feel
no gravity until he hit the ground (not so happy for the painter).
In the SR part of this book, we insisted on only dealing with reference
frames moving at constant speed, which we considered to be "inertial
reference frames".  However, a reference frame in gravitational free
fall is *also* an inertial frame, because it is no different from
a reference frame that is not moving in the absence of gravity.

There is a deep mystery in introductory Physics that is often glossed
over, and perhaps not even mentioned.  Gravity is introduced through
Newton's Law of Universal Gravitation, that any two masses in the
universe attract each other in proportion to the product of their
masses and inversely proportional to the square of the distance
between them, or
```{math}
:label: newUgrav
\vec{F} = -\frac{GMm}{r^2}\hat{r}.
```
Teachers have also by this point probably introduced Newton's
famous second law, often expressed that the total force on a system is
proportional to the rate of change of the velocity of a system,
with the constant of proportionality being called the mass, or
```{math}
:label: new2law
\vec{F} = m\vec{a}.
```
For a falling object near the Earth's surface, these two forces are
equated and the $m$ term cancels out, indicating that the acceleration
of any object near the surface of the Earth is the same, and we use
$GM/R^2$ with values for the mass and radius of the Earth to get
$g\equiv GM/R^2 = 9.8$ m/s/s.  The velocity of any object falling
near the Earth's surface will change by 9.8 m/s every second, regardless
of that object's mass.

Buried in the chain of logic was a big assumption, so big and obvious
that most students don't even notice, and often are confused when it
is actually pointed out to them: there is nothing in that description
that demands that the two times $m$ is used, they must be referring to
the same thing.  The $m$ in Equation {eq}`newUgrav` is defined in
terms of how strongly masses pull on each other through gravity.  The
$m$ in Equation {eq}`new2law` is defined in terms of how much an
object's velocity changes when you push on it.  There's nothing in
those two definitions that says those two *have* to be the same thing!
It's a huge assumption and a complete mathematical slight of hand
trick to define them with the same letter $m$ in the first place.

Say something about experimental tests.

Einstein's theory of General Relativity insists that they are indeed
the same thing, and provides a reason why that might be the case.
If what we call gravity is indeed a kind of illusion caused by being
in a curved spacetime, then falling objects aren't really falling,
but rather just drifting through spacetime as they normally would,
until the ground, or the floor of the rocket, comes up from below and
hits them.  Everything falls at the same rate because everything
is drifting together through a curved spacetime.  GR suggests that
you will never be able to do any experiment to be able to tell
the difference between these two situations (standing on the Earth
vs. standing on the floor of an accelerating rocket), because they
are equivalent.

Consider again the accelerating car.  However, now imagine the car is
parked.  Everything in the car will behave according to Newton's laws:
the dice will hang down, the safety strap will not press against you,
and the balloon will rise straight up to the top of its string.  If,
now, the Incredible Hulk were to grab the back of the car and lift it
up, such that the car were tilted at a large angle, consider what
would happen inside the car: the dice would shift forward, you would
feel pushed against the safety belt, and the balloon would tilt back
toward the back of the car.  An observer outside would say everything
is simply staying oriented to gravity's "down", but inside the car, it
would seem like a mysterious force was pushing foward, *exactly* as if
the driver had slammed on the brakes.  In fact, if the two kinds of
mass are truly equivalent, there is absolutely no difference between
braking the car and tilting the car.  In one case, we would say a
mysterious "fictional" force was adding a horizontal compontent to the
total force on objects in the car, and in the other case, we would say
that the force of gravity was rotating, so that it gained a horizontal
component that wasn't there before.  In both cases, the total force
vector, from the point of view of everything inside the car, would
swing from vertically down to diagonally forward.

This insight is so important, that this "Equivalence Principle"
occupies a similar place in the foundation of GR that the two
postulates of SR do for that theory.  Within the formulation of GR,
what we consider to be a force that we call Gravity is really an
illusion of objects trying to move through curved spacetime, just like
the balloon and the fuzzy dice aren't really being pushed by invisible
forces in the accelerating car.

## The Pound-Rebkha Experiment

## The Weight of a Box of Photons

## Einstein's Equation

## Metrics

## Geodesics
