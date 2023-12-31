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
(being the largest and closest body to Mercury, this is not unreasonable),
Newtonian gravity would predict that Mercury should travel in an
elliptical orbit.  However, Mercury moves as if its ellipse is drifting
around the Sun: when Mercury returns to the point of its ellipse
that is closest to the Sun, that point is not in the same place it
was at the previous orbit.  This is not a large effect -- the shift
is only 43 seconds of arc every century -- but that was enough to be
measured.  Scientists at the time thought that since the prediction
of Neptune was so successful, perhaps there was another planet inside
Mercury's orbit that was not being accounted for in the calculations.
This hypothetical planet was dubbed Vulcan, and for over 50 years,
astronomers tried to find it.  {numref}`vulcan` shows the orbit of
Vulcan, as factual as all the other orbits.

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




```{code-cell}
:tags: ["remove-cell"]
# Newton vs. Einstein: angular deflection in passing point mass
#def simp(yy,dx):
#  w = np.zeros(np.shape(yy))
#  w[0] = 1
#  w[1::2] = 4
#  w[2::2] = 2
#  w[-1] = 1
#  return(np.matmul(w,yy)*dx/3)
  
#def einstein(x,m,r):
#  return(m*r*(4*x*x+r*r)/(x*x+r*r)**2.5)
  
#def newton(x,m,r):
#  return(m*r/(x*x+r*r)**1.5)
  
nn = 1.0e7
zmax = 30
dz = zmax/nn
#z = np.linspace(0,zmax,int(nn))
mass = 1.0
radius = 1.0

#ez = einstein(z,mass,radius)
#nz = newton(z,mass,radius)

#edef = simp(ez,dz)*2.0
#ndef = simp(nz,dz)*2.0

fig = plt.figure(figsize=(14,8))
#plt.plot(z[z<5],nz[z<5],'r-',label='Newton')
#plt.plot(z[z<5],ez[z<5],'b-',label='Einstein')
#plt.title('Contribution to deflection')
#plt.xlabel('Distance from closest approach')
#plt.ylabel('Angular deflection')
#plt.legend()
#plt.text(2,1,'Newton total: {:4.2f}'.format(ndef))
#plt.text(4,1,'Einstein total: {:4.2f}'.format(edef))
plt.plot([1.1],[1.1],'b.')

glue("angdeffig", fig, display=False)

```

```{glue:figure} angdeffig
:figwidth: 800px
:name: angdef


Calculation of angular deflection of each step for a photon passing
by an object of mass 1 M$_\odot$.
```

## Acceleration While Standing Still?

(The Veritasium Video)

## Black Holes

## The Big Bang

