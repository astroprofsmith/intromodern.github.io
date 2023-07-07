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

# Chapter 6: The Displacement Four Vector

In this chapter, we will take a deep dive into the properties of
the displacement four vector, using all the tools developed in the
previous chapters.  Some of the material will repeat what has been
presented already, but in more depth.

## The Difference Between Two Events

The model we will be using to predict the behavior of systems moving
with speeds significant compared to the speed of light require that
the physical process be one such that there are two events occurring
that can be viewed by different observers. Some examples of physical
processes that can be understood this way are:


1. A strange particle being produced (first event), and at some later
time, decaying into other particles (second event).

2. A rocket ship leaving the earth (first event), and arriving at
  another star (second event).

3. A high energy electron striking a proton (first event) and then
  being scattered and leaving the interaction area at a different
  energy and direction of motion from the initial electron (second
  event).

4. A photon with energy in the radar range sent out by a cop (first
  event) that is reflected by a moving car, and then returns to the
  cop at a different energy (second event)

5. A star cruiser traveling at 0.900 the speed of light (first event)
  and then shooting of a laser pulse at an angle of $30^\circ$ with
  respect to the direction of motion of the star cruiser (second
  event).

6. A high energy proton entering the magnetic field of a bubble chamber
  (first event) and this same proton leaving the magnetic field
  traveling in a different direction (second event).

This model (which of course is special relativity) suggests how an
inertial observer, such as one that is moving at a constant speed (not
accelerating), should analyze what is different about the two
events. It also provides a way (called the Lorentz transformation) of
calculating what other inertial observers, moving with respect to the
first observer, would measure as the difference between these two
events.



