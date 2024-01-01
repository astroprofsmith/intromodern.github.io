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

# Chapter 14: Some Applications of General Relativity

## Mercury's Orbit

General Relativity demands a complete reconceptualization of how
gravity works.   Why go to so much effort?  What about Newtonian
gravity is so inadequate as to merit such a radical move?  Well,
as successful as Newton's gravity was, it was not without its
critics.  Leibnitz, for example, was never satisfied that Newton
provided no mechanism by which two objects in space might exert
forces on each other.  How does the Earth "know" where the Sun
is, that it should be pulled in the right direction?  Leibnitz
suggested the intervening space might be filled with some kind of
vortices to mediate the interaction, but such vortices were never
found, and the success of Newtonian gravity in mapping the orbits
of the comets and planets convinced everyone of its utility.

In fact, in the mid-19th century, Newtonian calculations were
successful in predicting the location of a new planet that no one
hitherto had known was there.  Careful measurements of the movement of
Uranus showed that it was not following the predictions of Newtonian
calculations, unless there was an unknown object out there exerting
forces on it that were not being accounted for.  Laborious
calculations from two different people, independently, predicted where
this unknown object should be, and the planet Neptune was discovered
in a location consistent with their predictions.

```{note}
Although German astronomer Johann Gottfried Galle is rightly credited
with discovering Neptune in 1846, he is actually not the first person to
*observe* the planet.  That honor goes to none other than Galileo
Galilei, who in 1612 was tracking the movement of the moons of Jupiter.

{numref}`neptune` shows a scan from Galileo's observing notebook,
which he sketched on the night of 27 Dec 1612.  It shows Jupiter
in the middle, with three moons to the right and a single moon to
the left.  He also included a background star, above and to the
left, or so he thought.  If you take a planetarium application
and set the date to 27 Dec 1612, you can scan through the night
to discover exactly when the moons line up with Galileo's drawing,
and you will see that the nearby "star" is actually Neptune!


```{figure} images/Galileo.png
:alt: galileoneptune
:class: bg-primary mb-1
:width: 700px
:align: center
:name: neptune

Sketch from Galileo's observing notebook on the night of 27 Dec 1612.
The dotted line leads to a nearby star that Galileo was using as a
reference, but back extrapolation shows that this "star" was actually,
unbeknownst to anyone at the time, the planet Neptune!
```



By the mid-19th century, it was becoming clear that Mercury did not
orbit strictly as Newtonian gravity said it should, either.  If you
ignore the other planets and just consider the force from the Sun
(being the largest and closest body to Mercury, this is not
unreasonable), Newtonian gravity would predict that Mercury should
travel in an elliptical orbit.  However, Mercury moves as if its
ellipse is drifting around the Sun: when Mercury returns to the point
of its ellipse that is closest to the Sun, the perihelion, that point
is not in the same place it was at the previous orbit.  The measured
drift of the perihelion is about 566$^{\prime\prime}$ per century.
Mercury takes about three months to complete one orbit, so in a single
orbit, the location of one perihelion only differs by about one and a
half arcseconds from the previous one.  That is a tiny change (one
arcsecond being 1/3600th of a degree, or about a tenth the width of a
sheet of paper held at arm's length).

Some of this shift could be explained by the effects of the other
planets, dominated by Jupiter, tugging on Mercury from other
directions.  Careful calculations (before electronic computers!) added
up to a prediction of 523$^{\prime\prime}$ per century.  This is not a
large discrepancy -- the difference is only 43 seconds of arc every
century -- but that was enough to be measured.  Scientists at the time
thought that since the prediction of Neptune was so successful,
perhaps there was another planet inside Mercury's orbit that was not
being accounted for in the calculations.  This hypothetical planet was
dubbed Vulcan, and for over 50 years, astronomers tried to find it.
{numref}`vulcan` shows the orbit of Vulcan, as factual as all the
other orbits.

```{figure} images/vulcan.jpg
:alt: vulcan1846
:class: bg-primary mb-1
:width: 700px
:align: center
:name: vulcan

Solar System map produced for education in 1846.  The
planet Vulcan is shown as orbiting closer to the Sun
than Mercury.  See [Library of Congress](https://www.loc.gov/resource/g3180.ct003790) for more details.
```
## Gravitational Lensing

If objects, including photons, really follow the paths through
spacetime that maximize the elapsed time, then their paths through
curved space would not be what we would consider to be straight lines.
Airline pilots presumably follow the shortest possible flight paths,
to conserve fuel, but if you drew a flight path from New York City to
Los Angeles on a flat map, the line would "curve up" perhaps even into
Canada, rather than going "straight" across Missouri.  This is because
the Earth is, of course, not flat, and the shortest distance on a
sphere is a great circle.  Draw the same flight path on a globe, and
it will look straight when you look directly down on it.  One of the
best ways we have to determine if we are in a curved spacetime is to
observe the trajectories of photons and see if they curve.

In the previous chapter, we discussed how objects in an accelerating
rocket ship will seem to those inside the ship to be falling to the
floor, while to those outside, it is the floor that is rushing up to
catch the objects.  Let us examine how that argument applies to
photons.  Consider first a stationary rocket, or equivalently, one
moving at constant velocity.  Imagine our astronauts have mounted a
laser pointer on one wall of the rocket.  It fires photons directly
across the room at the opposing wall, and we can mark where they hit
the wall.  To observers both inside and outside the rocket, the
photons will seem to fly in straight lines across the rocket.  Even if
the rocket is moving at a constant (non-zero) velocity, the photons
leaving the laser pointer will gain the same upwards momentum and
therefore will keep up with the rocket, traveling horizontally to
those inside the rocket, and along a straight diagonal line to those
outside.

However, if the rocket accelerates, the photons in flight will not
pick up the additional velocity gained by the rocket, and will therefore
hit the other wall below the mark we placed when the rocket was not
accelerating.  To the external observer, the photons will fly in a
straight line, while the rocket moves upward ever faster, ensuring
that they hit the wall below the mark.  To an observer inside the rocket,
the photon beam will appear to fall downward, just as a fired bullet
would, to strike the wall below the mark.

How far would the photons fall?  At $c$, it would take the photons ten
nanoseconds to cross a 3 m room in the rocket.  At an acceleration of
10 m/s/s, the photons would fall $5\times10^{-16}$ m in ten
nanosecond.  This is on the order of the size of subatomic particle,
like a proton.  That's not an easily measurable amount, in real life.
However, if it were, the Equivalence Principle demands that light moving
across a 3 m room on the surface of the Earth should also fall by
$5\times10^{-16}$ m.  If not, we could determine whether or not we
are actually sitting on the Earth's surface by measuring whether or
not a horizontal beam of photons falls while travelling across the room.
If it hits the mark, we are on the Earth.  If it hits below the mark,
we are in an accelerating rocket.

The Equivalence Principle says no such experiment can distinguish the
two situations, so the photons near the Earth must also fall down.  If
the photons are traveling along a geodesic, that means the spacetime
near the Earth must be curved.

In practice, the Earth is just not massive enough to provide strong
enough spacetime curvature to measure the falling of a beam of
photons.  If we want to see whether light is deflected by gravity,
we need to consider much, much more massive objects, like stars, or
even better, entire galaxies.  If we consider the trajectories of
photons traversing the vast distances of interstellar or intergalactic
space, we can ask whether those trajectories are distorted away from
straight lines by the gravity of massive objects.  If mass distorts the
spacetime around itself, one would expect a photon path near a massive
object to bend, and not continue straight as one would expect the path
in flat spacetime to do.

In 1915, when Einstein was developing GR, he had one ready-made
real-world phenomenon that he could explain: the precession of the
perihelion of Mercury.  GR could explain those 43 arcseconds per
century without the need of a new planet.  However, in science, it's
never completely convincing to explain a phenomenon which is already
known -- it is too easy to rig the game to reach a conclusion you have
already decided in advance.  It is much more convincing if your theory
can predict an observation which no one has made yet, particularly if
it is an observation no one is expecting.  Buoyed by the agreement of
his calculations with the Mercury observations, Einstein was looking
for a new phenomenon he could predict.  He turned to the bending of
light by gravity.

Since the Earth was too small to generate measurable bending, Einstein
suggested the light passing by a much bigger object, such as the Sun,
be measured.  If the light from distant stars showed a deflection in
direction, his theory could be tested through quantitative prediction
and measurement.  The problem of course is that it is very difficult
to observe starlight passing very close to the Sun, since we cannot
usually see stars during the day time.  A major exception to this
restriction is during a total solar eclipse.  Then the Moon blocks
enough of the sunlight that the position of stars could be measured.
As luck would have it, not only was a total solar eclipse coming up in
1919, but that eclipse would occur near the Hyades star cluster, so
there would be many background stars available for comparison.
Careful measurements of the background stars' positions, when compared
with images of the same cluster taken months earlier, at night, when
the Sun was in the other direction, could reveal whether the stars'
positions shifted, and if so by how much.

There are actually two alternate hypotheses to the GR prediction.
If spacetime is not curved at all, and if massless photons experience
no force, one would expect no shift in their positions.  It is also
possible to approach this problem from an almost Newtonian paradigm
and hypothesize that photons, while massless, might still experience
an acceleration $a=GM/r^2$, just like a massive object would, regardless
of the fact that you can't divide both sides of an equation by $m=0$.
This Newtonian prediction is different from the Einsteinian prediction
because in the Newtonian paradigm, time remains a universal constant,
but in the GR paradigm, the "clocks" of the photons will slow down
as they approach the Sun and then speed up again as they receed,
increasing the deflection to twice what a Newtonian deflection
would predict.

This slowing down of time near the massive object causes the
trajectory of the light wave to shift, much like the slowing down of
light in glass or water causes the trajectory of light to shift in
these media.  Controlling this deflection by shaping the amount light
slows down is the basis for lens design and the whole science of
optics, and it is in this sense that the GR phenomenon is called
Gravitational Lensing.  Unlike, say, a glass lens, the light
trajectories around a star are not bent so as to get them to form
an image.  The bending is greatest near the center and least far from
the object, which is the opposite of a glass optical lens.  To make
a glass lens behave like a gravitational lens, its curvature would
have to be greatest near the center, much like the base of a stemmed
wine glass.  {numref}`glasslens` shows an example of a glass lens that
was constructed deliberately to behave like a gravitational lens.  The
glass in the photograph was about 45 cm in diameter, and the curved
middle extended roughly 10 cm toward the viewer.

```{figure} images/glasslens.jpg
:alt: glasslens
:class: bg-primary mb-1
:width: 700px
:align: center
:name: glasslens

A glass lens distorts the view of a rope light, suspended horizontally
behind it.  The lens is shaped like that of the base of a wine glass,
curving forward at the center from a flat edge, with the stem in the
middle shaved off.  The flat rope light appears curved above the center
of the lens, and a second image of the light appears below the center,
just as a point gravitational lens would create two images of a
background star, if you considered the lensing to be happening along
a vertical plane slice through the middle of this picture.
```

In 1919, multiple expeditions from the Royal Society of London set out
to South America and Africa to wait along the path of the Moon's
shadow and image the sky during the eclipse.  Frank Dyson and Arthur
Eddington led these expeditions, and when they returned, they
published their measured shifts in the stars' positions.  The graph
from their 1920 report is shown in {numref}`Edd19`.  Stars nearer to
the Sun are clearly deflected more, and although the uncertainties are
fairly large, they clearly rule out the Newtonian prediction and are
consistent with the Einsteinian calculations.  This confirmation of a
German scientist's predictions by English observers, overthrowing the
long-established theories of that premeire English scientist, Isaac
Newton himself, so shortly after the two countries were at war,
enthralled the world and catapulted Einstein into the international
spotlight, making him the pop culture figure he remains today.

```{figure} images/Eddington.png
:alt: edd1919
:class: bg-primary mb-1
:width: 700px
:align: center
:name: Edd19

Dyson and Eddington's published report of star image deflection during
the 1919 solar eclipse.  (Dyson, F. W.; Eddington, A. S., Davidson
C. (1920). "A determination of the deflection of light by the Sun's
gravitational field, from observations made at the total eclipse of 29
May 1919". Philosophical Transactions of the Royal Society 220A:
291–333) The angular distance from the Sun is printed along the
horizontal axis, but note that the distance increases to the left.
The y-axis shows the total amount of deflection, decreasing for
increasing distance from the Sun.  The lighter solid line is the
best-fit linear relationship, while the darker solid line in the
middle is the prediction from GR.  The prediction from Newtonian
physics is shown by the dotted line below, a factor of two less than
what Einstein predicted, and clearly inconsistent with the measured
results.
```

How can we understand this deflection in more depth?  The easiest
approximation is to use the thin lens approximation.  With the
distances between the observer, the lens, and the source being
literally astronomical, we can assume that the bending happens within
a space along the path that is insignificantly short compared with the
other distances in the problem.  The light travels in straight lines
far from the lensing object, and then bends as it passes through a
plane at the lens location.  This scenario is presented in
{numref}`thingravlens`, as viewed from a location to the side.  The
source is far away to the left, a distance SL from the lens.  The
observer is on the right, a distance LO from the lens.  The distance
from the source to the observer is SO=SL+LO.  The light path reaches a
closest distance to the lens of $b$, and then is deflected through an
angle $\alpha$.  The angles $\beta$ and $\theta$ are measured with
respect to the baseline that connects the observer to the lens:
$\beta$ is the angular position of the original source (the magenta
dotted line), and $\theta$ is the angular position of the shifted
image (the cyan solid line) where the light rays appear to be coming
from.  The angles in {numref}`thingravlens` are greatly exaggerated so
they can be seen by eye.  In real life, all these angles are very
small, so we use the small angle approximation for all of them.




```{code-cell}
:tags: ["remove-cell"]
# The thin lens diagram

fig = plt.figure(figsize=(14,8))

plt.plot([-0.9],[1.0],'ro')
plt.plot([1.5],[-0.5],'bo')
plt.plot([0.0],[-0.5],'ko')

xdir = np.linspace(-0.9,1.5,100)
ydir = 1.0 - (1.5/2.4)*(xdir+0.9)
plt.plot(xdir,ydir,'m:')

ybent = np.ones(np.shape(xdir))
plt.plot(xdir,ybent,'c:')
ybent[xdir>0.0] = 1.0 - xdir[xdir>0.0]
plt.plot(xdir,ybent,'c-')
ybent = 1.0 - xdir
plt.plot(xdir,ybent,'c:')

plt.plot([-0.9],ybent[0],'co')

plt.plot([0.0,0.0],[-0.5,1.0],'k:')
plt.plot([-0.9,1.5],[-0.5,-0.5],'k:')

plt.text(-0.1,0.1,'b',fontsize=18)
plt.text(0.2,0.85,r'$\alpha$',fontsize=18)
plt.text(-0.2,1.2,r'$\alpha$',fontsize=18)
plt.text(0.8,-0.3,r'$\theta$',fontsize=18)
plt.text(1.2,-0.45,r'$\beta$',fontsize=18)
plt.text(0.7,-0.7,'LO',fontsize=18)
plt.text(-0.4,-0.7,'SL',fontsize=18)

ax = plt.gca()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.axis([-1,2,-1,2])

glue("thinlensfig", fig, display=False)
```
```{glue:figure} thinlensfig
:figwidth: 800px
:name: thingravlens


The thin gravitational lens diagram.  The lens is the black dot
at the bottom center.  The source is the red dot to the left, and
the observer is the blue dot to the right.  The cyan line shows
the path that photons take to actually reach from the source to the
observer, deflected through an angle $\alpha$ from their original
trajectory.  Photons that started out along the magenta
dotted line would have reached the observer if the lens weren't there,
but in this case would be deflected far below the observer (this
deflection is not shown).  With the lens, the photons appear to the
observer to be coming from the location of the cyan dot, and therefore
the image of the source is shifted by an angle $\theta-\beta$.
```

We split the vertical height of the right side into two pieces
and use the small angle approximation to add them in terms
of the triangles involved: $\theta (SL+LO) = \alpha SL + \beta (SL+LO)$



```{code-cell}
:tags: ["remove-cell"]
# Newton vs. Einstein: angular deflection in passing point mass
def simp(yy,dx):
  w = np.zeros(np.shape(yy))
  w[0] = 1
  w[1::2] = 4
  w[2::2] = 2
  w[-1] = 1
  return(np.matmul(w,yy)*dx/3)
  
def einstein(x,m,r):
  return(m*r*(4*x*x+r*r)/(x*x+r*r)**2.5)
  
def newton(x,m,r):
  return(m*r/(x*x+r*r)**1.5)
  
nn = 1.0e7
zmax = 30
dz = zmax/nn
z = np.linspace(0,zmax,int(nn))
mass = 1.0
radius = 1.0

ez = einstein(z,mass,radius)
nz = newton(z,mass,radius)

edef = simp(ez,dz)*2.0
ndef = simp(nz,dz)*2.0

fig = plt.figure(figsize=(14,8))
plt.plot(z[z<5]+1,nz[z<5],'r-',label='Newton')
plt.plot(z[z<5]+1,ez[z<5],'b-',label='Einstein')
plt.title('Contribution to deflection')
plt.xlabel('Distance from Star ($R_\odot$)')
plt.ylabel('Angular deflection')
plt.legend()
plt.text(2.5,1,'Newton total: {:4.2f}'.format(ndef),fontsize=18)
plt.text(4.5,1,'Einstein total: {:4.2f}'.format(edef),fontsize=18)

glue("angdeffig", fig, display=False)
```
```{glue:figure} angdeffig
:figwidth: 800px
:name: angdef


Calculation of angular deflection of each step for a photon passing by
an object of mass $1~M_\odot$.  Only the receding half of the
trajectory from the point of closest approach ($1~R_\odot$) is shown.
The approaching trajectory, not shown, is symmetric.  The curves, when
numerically integrated and then doubled to account for the approaching
half of the trajectory, show the total amount of angular deflection.
These answers are shown superimposed on the graph.  Note that the GR
calculation yields exactly double the Newtonian calculation for the
total deflection, although the blue curve is not twice the red curve
at every point. *NOTE*  Need to check units of deflection.
```



## Acceleration While Standing Still?

(The Veritasium Video)

https://www.youtube.com/watch?v=XRr1kaXKBsU

Not sure this is worth including.  Maybe.  It's about 10 minutes in.

## Black Holes

## The Big Bang

