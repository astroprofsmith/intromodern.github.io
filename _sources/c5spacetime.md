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

# Chapter 5: Properties of Spacetime Diagrams

Now that we have the tools of the displacement four vector, the
spacetime diagram, and the Lorentz transformation, it is worth taking
some time to look at what we can learn about the properties of space
and time by applying these tools to the spacetime diagram.

## Types of Intervals

```{code-cell}
:tags: ["remove-input"]
# Insert VPython simulation of a Michaelson Interferometer
# Allow user to rotate system, relative to ether
# Have radio button to include/remove ether
url1 = "https://glowscript.org/#/user/dasmith/folder/Public/program/SRintervals"
display.IFrame(src=url1,width=880,height=700)
```
```{note}
Figure 5.1 -- Interactive spacetime diagram.  The red line represents
the worldline of a stationary observer.  The blue dot represents some
event $p$ on that worldline.  The green dot represents some other event
$q$ that most definitely does not reside on the worldline with $p$.
To get information about $q$, therefore, the observer must send a
light ray out to $q$ and get the reflection back.  The worldlines of
these light rays are shown in orange.  The observer can therefore
define two time intervals that represent the elapsed time between
event $p$ and the events when the light was emitted and received.
From $t_1$ and $t_2$ the observer can calculate a $\Delta x$ and
a $\Delta t$, as described in the text.  These four values as
well as the interval you get from the displacement four vector
are shown on the diagram.  Move the slider to move $q$ up and down
relative to $p$, and click the box to turn light cones for $p$ on
and off.
```

## How do the Axes Change?

