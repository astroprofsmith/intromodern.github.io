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
:label: intint
\Delta S = \int_{\rm path} cdt\sqrt(1-\frac{dx^2}{c^2dt^2}} = \int_{\rm path}
cdt\sqrt{1-\frac{v^2}{c^2}}
```
What we have here is a path integral of some function, call it $f$,
that depends on the speed of the clock at each moment along the path:
$f = \sqrt{1-v^2/c^2}$.  For the red path, $v=0$, so $f=1$ and the
integral is trivially just the proper time.  However, along the blue
path $0\geq v^2<1$, so $0<f\leq 1$, which means the result of this
integral (which is the time interval elapsed on the blue clock, as we
are integrating over $dt$) is always going to be less than or equal to
the proper time.

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



## A Rotating Reference Frame

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

## The Equivalence Principle

## The Pound-Rebkha Experiment

## The Weight of a Box of Photons

## Einstein's Equation

## Metrics

