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

```

# Chapter 1:  Roots of Relativity

Before we dive into the set of predictions and interpretations that
make up what we call "Special Relativity" (or SR), it's important to
lay some groundwork and establish some definitions for what we will be
talking about.  Much of SR will twist and challenge what you thought
you knew about space and time, so it is important to define our terms
carefully and precisely, or else it can get very confusing as to what
we are actually talking about.  In this chapter, we will describe the
rules by which the game of SR is played, and describe some of the
reasons why those rules have been accepted.  It's important to know,
right from the beginning, that the predictions of relativity, in both
the general and specific forms, have passed every single test to which
they have been subjected.  As far as we can tell, the universe really
does work this way.

## Experimental Foundation

Before the 1900s, a popular theory about the propagation of light
waves was that they traveled in some medium called the 'ether'. The
'ether' plays the analogous role that the air plays in the propagation
of sound. The physical properties of the air determine the speed of
the sound waves. In the case of the 'ether', it supported a speed of
propagation of $3.0 \times 10^8$ m/s. Following the usual practices of
physics, the validity of this theory would be determined by
experiments. As a test of this theory, in 1881, Michelson did an
experiment to measure the relative velocity of the ether with respect
to the earth. The 'best' ether theory of the time suggested that the
ether was at rest with respect to the fixed stars. The earth orbits
the sun once a year at an orbital speed of about $3.0 \times 10^4$
m/s.  If a light beam on earth travels in the same direction as the
orbital motion about the sun, then that beam of light is in a frame of
reference that is moving at $3.0 \times 10^4$ m/s with respect to the
speed of the ether. If the beam of light is heading towards the sun,
it is moving in a frame of reference that is at respect to the ether.

Michelson designed an interferometer (schematically shown in Figure
1.1) that used partially silvered mirrors to split a beam of light
into two legs that were at right angles to each other. By careful
examination of the interference of the light after it traveled along
the two beam lines and was recombined by the partially silvered
mirrors, he could measure a difference in speed of light in the two
arms of about 300 m/s, much smaller than the orbital speed of the
earth about the sun. The interferometer was floated on a pot of
mercury and allowed to rotate slowly. During a rotation of 360
degrees, at one time one of the legs would be pointing in the
direction of motion of the earth (moving with respect to the ether)
and the other at rest with respect to the ether. A quarter of a turn
later, the light in the arm that was at rest with respect to the ether
would be moving and vice versa. This anticipated change in the speed
of light when the either was at rest with respect to the mirrors or
moving with respect to the mirrors in the other arm could be detected
in the changing of the interference pattern of the combined beams.

Much to Michelson's surprise, he measured no difference in speed as
the arm was at rest with respect to the ether or moving with respect
to the ether. Michelson had this terse description of the experiment:
"The interpretation of these results is that there is no displacement
of the interference bands. ... The result of the hypothesis of a
stationary ether is thus shown to be incorrect." (A. A. Michelson,
Am. J. Sci, 122, 120 (1881)).  Michelson was encouraged by some
important English physicists to continue his research. Along with
E.W. Morley, he built a much bigger interferometer that had multiple
mirrors and was 10 times as sensitive to the relative motion. He could
measure a relative velocity of 30 m/s. Again, they found no
difference in the speed in the two arms. (A. A. Michelson and
E. W. Morley, Am. J. Sci., 134, 333 (1887)). The experiment was
repeated often over the next few years, at even greater precision, but
always gave the same null result: there is no stationary ether. As a
consequence of his careful and brilliant work about the nature of the
propagation of light, Michelson was awarded the Nobel prize in 1907.

As with most interesting experiments, Michelson got a result that did
not agree with his intuition about the outcome of the
experiment. After a number of very bright physicists examined these
results, a relatively simple and rather clever model was developed
that agreed with the experimental results.

```{code-cell}
:tags: ["remove-input"]
# Insert VPython simulation of a Michaelson Interferometer
# Allow user to rotate system, relative to ether
# Have radio button to include/remove ether
url1 = "https://glowscript.org/#/user/dasmith/folder/Public/program/SRinterferometer"
display.IFrame(src=url1,width=700,height=700)
```
```{note}
Figure 1.1 -- Animation of a simplified Michaelson interferometer.
```


## Example of Binary Star Systems

Another piece of evidence for the constancy of the speed of light
comes from observations of binary star systems.  If two stars,
orbiting each other, lie in a plane that is shared with the Earth,
then we will observe eclipses every half orbit, as one star blocks at
least some of the light from the other.  Half of each orbit, one star
is moving toward the Earth, and the other is moving away.

If the speed of light depended on the relative speed of the emitter
and/or the receiver, the light from the star approaching the Earth
would be moving faster than the light from the receeding star.  Over
the years of travel between the stars and the Earth, the faster light
would catch up to the slower light, perhaps even overtaking it.  The
observations of the moments in the orbit, such as one star passing in
front of the other, would arrive at the Earth at completely distorted
times.

Of course, no such distortion has ever been observed.  Eclipses arrive
like clockwork.  Arrival times of the eclipses' beginnings and endings
might need to be corrected for the Earth's motion around the Sun
(light will take up to sixteen minutes longer to arrive at the Earth
when it is on the far side of the Sun from the star), but there is no
evidence that some of the light from the stars is getting here faster
than other light.

```{code-cell}
:tags: ["remove-input"]

# Insert simulation of an eclipsing Red-blue binary light curve
# Allow user to change c and see how light curve distorts

```

## Model to describe the experimental results.

One interpretation of Michelson's experimental results is to say that
every inertial (moving at constant speed... not accelerating) observer
measures that light travels at the same speed. (This is quite in
contrast to how sound moves through the air. If I am moving towards
the source, the sound seems to be traveling faster than when I am
moving away from the source of the sound.) If an inertial observer
measures the speed of the photons that come from a flashlight that she
holds in her hand, she would get the value $c = 3.0 \times 10^8$
m/s. If a second observer, moving at a speed of $2.8 \times 10^8$ runs
away from (or towards) the flash light would also measure the speed of
the photons emitted by that same flashlight to be $c$, precisely the
same value as the other observer. This is intuitatively absurd, but it
agrees with the experiments, so it is a 'not wrong' model for the
propagation of light.

The next step in developing the model is to transform this literary
statement into mathematical terms so that quantitative predictions can
be made. Consider a point source emitting light waves in three
dimension. An observer would see the light waves traveling out in the
shape of a sphere whose radius is dependent of time (as seen in
Figure 1.3). As the light is emitted, in some very small
time interval $dt$, the wave front makes a spherical shape of radius
$r = c dt$. The 3-dimensional equation for a spherical shape centered
at the origin is:

```{math}
:label: eqsph
x^2+y^2+z^2= r^2
```
Using the distances displayed Figure 1.3, equation {eq}`eqsph` becomes:
```{math}
:label: eqsph2
(dx)^2+(dy)^2+(dz)^2 = (cdt)^2,
```
where $dx$, $dy$, and $dz$ are how much the observer measured the wave
front to propagate in time $dt$, and $c$ is the speed of light as
measured by the observer.

Consider a second reference frame, which we designate "the primed
frame" (all variables measured in this frame will have a prime on
them, like $x^\prime$), traveling at a constant speed $v_R$ in the
$\hat{x}$ direction with respect to the first observer (the first
observer find the source of the light flash to be at rest) were to
observe the same event (the wave front of the light propagating away
from the point source), what should she see?

```{code-cell}
:tags: ["remove-input"]
# Insert animation of expanding sphere of light.  Show radial
# arrow and overlay x, y, z, and ct
url = "https://glowscript.org/#/user/dasmith/folder/Public/program/SRExpandSphere"
display.IFrame(src=url,width=600,height=600)
```
```{note}
Figure 1.3 -- Animation of an expanding sphere of photons from an
initial flash.  The radius of the sphere, with length $ct$
is shown as a black arrow, with the cartesian components indicated
by red arrows.  The animation will rotate to display the three-dimensional
nature of the diagram.
```

Common sense says that the wave front should not longer seem to be
spherical, but should appear to be oblate (a squashed
sphere). However, what Michelson found was that this spherical wave
had to look the same to this observer as it appeared to the first
observer. Otherwise, the second observer would measure a different
velocity for the light wave. What the observer in the primed reference
frame measures is a spherical wave that propagates outward at the same
speed as shown in Equation {eq}`eqsph`!. The second (primed) observer
also sees a spherical wave, but this time the observer measures:
```{math}
:label: eqsphp
(dx^\prime)^2 + (dy^\prime)^2 + (dz^\prime)^2 = (c dt^\prime)^2
```
where $dx^\prime$, $dy^\prime$, and $dz^\prime$ are how much the
second observer measured the wave front to propagate in time
$dt^\prime$, and $c$ (not $c^\prime$!) is the speed of light.

Equations {eq}`eqsph2` and {eq}`eqsphp` can be combined to be:
\begin{equation}
(dx)^2+(dy)^2+(dz)^2 - (cdt)^2 = (dx^\prime)^2 + (dy^\prime)^2 + (dz^\prime)^2 - (c dt^\prime)^2
\end{equation}(eqn:ds)

It appears that there is a sum of squares of measured displacements
and time intervals for the measurement of the propagation of the wave
front of a light wave that are the same for each of the two observers
(moving with respect to each other) that remains constant. This
statement can be written as a conservation law. The sum:
\begin{equation}
(dx^\prime)^2 + (dy^\prime)^2 + (dz^\prime)^2 - (c dt^\prime)^2 
\end{equation}(eqn:interv)
must remain constant. This conservation law is quite different from
what Newton and company thought was 'conserved' in the measurement of
the motion of an thing. In the Newtonian world, the value of $dt$
would be the same for both observers, implying a $c^\prime \neq c$,
not this strange looking sum. But, if this model agrees with the
experiments, it is not wrong; even if it disagrees with Newton's
model.


## The Postulates of Relativity

All of relativity follows from two surprisingly simple postulates.
Once you accept these postulates, you can use the tools of mathematics
to derive predictions for how space and time will behave, and how
objects will move within them.  

1. The speed of light is the same in all inertial reference frames.

This is the most surprising one, but as explained above, nature seems
to have forced this idea onto us.  No other speed behaves like this.
When someone walks forward on an airplane in flight, a fellow
passenger would perceive them as walking at a meter per second or so,
but someone on the ground would perceive them as moving at hundreds of
meters per second.  This seems obvious.  The premise that light from a
flashlight on an airplane would be measured as having the same speed
on the plane as from the ground seems absurd.  However, as introduced
above, this is what actually does happen in nature.

2. The laws of physics are the same in all inertial reference frames.

Another way to express this postulate is to assert that there is no
universally preferred, or absolute, reference frame.  If you are on a
train moving East, and I am on a train moving West, there is no
physical experiment we can perform that will tell us who is "really"
moving.  You might try to assert that while both trains are moving,
the ground is not, but of course, the Earth is spinning and orbiting a
star that is itself orbiting a galaxy that is drifting through the
universe.  The closest thing there is to a universal reference frame
is the reference frame of the cosmic microwave background radiation,
which permeates the whole universe as a record of that moment when
most of space became transparent about 14 billion years ago.  That's a
very large system to choose as a reference, but the choice is
arbitrary.  There's nothing about the universe that says even that
choice is more natural than any other.  Once you choose a reference
frame, any set of experiments or measurements you carry out will
follow the same set of rules as they would in any other reference
frame.  The answers could be different, but they wouldn't be "more
right" or "more wrong" than measurements carried out in some other
inertial reference frame.

## Inertial Reference Frames and Events

The term "inertial reference frame" has come up several times already.
Before we go any further, we need to be very clear on what this term
means.  In essence, a "reference frame" refers to a single set of
space and time coordinates.  Pick an origin in space and time, and
then clock readings mark locations along a time axis, and ruler
readings mark locations along spacial axes, which usually are the
standard cartesian x, y, and z.  Four numbers can therefore define a
moment location in space at a moment in time, relative to the chosen
origin.

The difference between two ruler readings along an axis is
called a displacement (if the two points in question are not along the
same axis, we must use the Pythagorean Theorem to find the amount of
displacement along the diagonal, or we can express the total
displacement as a vector with a displaceent along each axis).  The
difference between two clock readings is called a time interval or a
duration.  If you consider smaller and smaller durations and
displacements, the region described by these numbers will shrink to be
a mathematical point in space and an instantaneous moment in time.
This is sometimes represented as a finger snap or a lightning flash,
although in reality neither of those actions are infinitesimal in
extent.  An infinitesimal chunk of space and time is called an
"event".  The total collection of all events is considered to be
"reality".

If reality is the collection of all events that happen throughout
space and time, the question that concerns us in this course on
Special Relativity is how these events are measured by different
observers in relative motion.  We need to be a bit more specific about
what it means to establish a system of coordinates in space and time.
Typically by the time you reach the age at which you are reading this
book, you will already take for granted the idea that space and time
are a kind of fixed lattice underlying reality as we know it, and that
locations in space and time have an absolute meaning independent of
who is trying to measure them.  One of the goals of this book is to
convince you this intuitive picture does not actually match reality.

Therefore, we must be more specific about what these four numbers are
and how they might be measured.  The typical image that is called
forth to represent this process is a lattice of rulers and clocks.
Imagine you had an infinite number of clocks and also an infinite
number of rulers.  You painstakingly synchronize all these clocks, and
confirm that the rulers are all not different in length.  Now, you
very slowly move the clocks out through space, putting one at the end
of each ruler, and then extending the rulers out from each placed
clock, so that in the end, you have a lattice of cubes made by the
rulers, with a clock at each vertex.  The location of each clock is
therefore the number of rulers in each direction from the origin, and
each clock has a reading, based on the overall synchronization.  A
single event, therefore, can be recorded as happening at the nearest
clock (as long as these rulers are imaginary, we can imagine them
being as small as we need them to be to achieve the desired spacial
resolution), marked at that clock reading: x, y, z, and t.  Those
numbers can be collected after the events by bringing the clocks back
together and collating their readings.

This concept of an infinite lattice of synchronized clocks represents
our idea of a "reference frame".  As long as we are imagining an infinite
lattice, we can imagine a second infinite lattice (or a third, or as
many as we like) that is passing through the first, moving at a constant
relative velocity.  This velocity could be measured by recording when
certain clocks pass each other, and when the clocks are gathered up again
later, the number of rulers traveled divided by the durations on the
relevant clocks would reveal the relative velocity, which we will
designate as $v_R$.  For special relativity, we require that $v_R$
be constant.

```{code-cell}
:tags: ["remove-input"]
# Insert animation of expanding sphere of light.  Show radial
# arrow and overlay x, y, z, and ct
url = "https://glowscript.org/#/user/dasmith/folder/Public/program/SRrefframes"
display.IFrame(src=url,width=800,height=800)
```
```{note}
Figure 1.4 -- Animation of two reference frames in relative motion.
The spheres represent clocks, and the rods represent rulers.  The viewer
must imagine that the lattice continues on indefinitely an all directions.
The camera is at rest with respect to the reference frame colored red,
while the blue frame is in constant relative motion.  The animation will
loop, to represent more clocks and rods coming in from the side.
```



Ultimately, we are interested in how objects in the universe move and
interact with each other.  To study this motion, we need to make
measurements of location and duration, which needs to happen within a
chosen reference frame.  Measurements recorded in different reference
frame will be compared, but they should never be mixed together.  A
displacement measure in one reference frame divided by a duration
measured in another reference frame would not represent a real
velocity.  An object could be at rest (velocity of zero) in one
reference frame, but to observers recording in a second reference
frame, moving at a relative velocity to the first, the object would be
recorded in different places at different times.  In other words, to
the observers at rest with respect to the second reference frame, the
object and the lattice of rulers and clocks with it would be moving,
not their own lattice of rulers and clocks.

Another important assumption we must make about the reference frames
we consider in special relativity is that they be "inertial frames".
What this means, simply expressed, is that Newton's laws are accurate
when objects are observed within the frame.  Objects' velocities do
not change unless a force acts upon them, momentum and energy are
conserved, and so forth.  An example of a non-inertial reference frame
would be if you put the reference frame lattice on a carousel.  Using
this rotating lattice of rulers and clocks, objects could be observed
to move in complicated, changing patterns, without any observable
forces causing the changes in motion (so we make up fictitious forces
and call them names like "centrifugal" to make Newton's laws keep
working).  The rules of special relativity only work in inertial
reference frames.  Usually, but not always, this means reference
frames that have constant velocity, including zero.

## Important Correlaries to the Postulates

Nothing physical can travel faster than $c$.  We will see why in later
chapters of this book, but it's important to be aware of that from the
beginning.  Once you demand that a speed be the same in all reference
frames, as the first postulate does, the universe will conspire to
make sure nothing moves faster than that.  It's important to note that
we call this speed "the speed of light", but that's because light,
having no inertial mass, moves at the fastest speed possible.  Light
will slow down in materials like water or glass, so if we're being
careful, we should call this speed "the speed of light in a vacuum".
But really, it's the speed light moves in a vacuum simply because
that's the fastest speed the universe allows.  Anything with no
inertial mass will move at that speed (the particles called neutrinos
also have almost no mass, and they move that fast, too).  So it is not
possible for anything to move from one place to another faster than
$c$.

For precisely that reason, you can't "just know" what's happening at
some other location than your own.  It takes time for anything, anything
that might convey information or have an effect, to move.  Even when
you look across the room, you are seeing the room as it was a few
nanoseconds ago.  The Moon you see in the sky is the Moon as it was
about a second ago.  The Sun that shines on us is the Sun from eight
minutes ago.  We often describe physical situations as if we have an
instantaneous overview of everything everywhere all at once.  This is
impossible.  Even to measure the ends of an object with a ruler demands
either moving one's self to the other side of the object, or sending a
light beam to the other side and back.  Many misunderstandings and
confusions in special relativity arise from people assuming they just
know what is going on at another location, without allowing for the
time it would take to get from one location to another.

Most of special relativity is about comparing measurements according
to two different reference frames that are in motion relative to each
other.  *Special* relativity demands the assumption that the relative
velocity between the frames is constant (a "special" case).  If the
motion of the reference frames is not constant, the theory must be
modified to a more general form, which makes it General Relativity.

Once you have two different reference frames, it is important to state
three characteristics that might seem self-evident, but they have
important implications that are worth articulating explicitly.  First,
imagine each frame of reference has a set of observers associated with
it, or acting within it.  Each set of observers has to agree that the
same events happened.  Switching reference frames does not change what
actually happens.  People in different reference frames may well
disagree on when and where events happen, but they should not disagree
on *whether* the events happen.

Second, each set of observers must agree on what the results of
measurements are.  Not only must there be agreement within the set of
observers within a particular reference frame, but each set of
observers must agree on what measurements are made in another
reference frame.  Those measurements might well contradict each other,
but all observers should agree on what they are.  For example, if a
set of observers in reference frame 1 agree that a room is 12 m
across, a set of observers in reference frame 2 should not be able to
truthfully assert that the first set of observers actually measured
the room to be 11 m across.  The second set of observers could
feasibly measure the room to be 11 m across, themselves, but they
should agree that the first set of observers measured 12 m.

All sets of observers should agree that causes come before effects.
We will find that switching reference frames will change many aspects
of space and time that you are used to thinking of as universal.
However, causes must still come before effects.  This is, as some have
suggested, a kind of definition of time.  As the old saying goes,
“Time... is what keeps everything from happening at once.”
(R. Cummings, "The Girl in the Golden Atom", 1919).  See also [this
Science Asylum YouTube
video](https://www.youtube.com/watch?v=7HBKEDyFTv8).  Not an
**operational** definition, mind you, but a useful description.

## Then Let's Begin...

Armed with these postulates, concepts, and principles, we are
ready to begin characterizing the universe of space and time,
and to discover the surprises waiting therein...