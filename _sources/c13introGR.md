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

Although a full treatment of General Relativity is far beyond the
scope of this book, I would like to leave you with some hints at
the major concepts you must grapple with, if you extend the theory
beyond the special case of reference frames moving with constant
relative velocities.

## The Twin Paradox Revisited



```{code-cell}
:tags: ["remove-cell"]
# ST plot of two paths
fig = plt.figure(figsize=(5,5))
plt.arrow(-0.5,-0.5,1,0,head_width=0.1)
plt.arrow(-0.5,-0.5,0,1,head_width=0.1)

plt.plot([0.2],[-0.1],'ko')
plt.plot([0.2],[1.1],'ko')

ax = plt.gca()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.axis([-1,2,-1,2])
ax.text(-0.05, 1.2, "ct")
ax.text(1.2, 0, "x")
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

## The Equivalence Principle

## The Pound-Rebkha Experiment

## The Weight of a Box of Photons

## Einstein's Equation

## Metrics

