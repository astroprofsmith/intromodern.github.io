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
from myst_nb import glue
```

(introduction)=
# Particle Physics


This is how to have the image change with mouse over.  The images
have to be in the "_static" directory.

<img src="_static/frameAfig.png" onmouseover="this.src='_static/frameBfig.png'"
onmouseout="this.src='_static/frameAfig.png'"/>



## Chickenhouse Experiments


Example of how to include VPython simulation

```{code-cell}
:tags: ["remove-cell"]
# Insert VPython simulation of a Michelson Interferometer
# Allow user to rotate system, relative to ether
# Have radio button to include/remove ether
url1 = "https://glowscript.org/#/user/dasmith/folder/Public/program/SRinterferometer"
test = display.IFrame(src=url1,width=800,height=700)
glue("interfig",test, display=False)

```

```{glue:figure} interfig
:figwidth: 800px
:name: michelfig

Animation of a simplified schematic of a Michelson
interferometer.  A laser, represented by the red cylinder to the left,
shines a beam to the right.  The light is split by a diagonal
half-silvered mirror.  Half the beam continues to the right, while
half the beam goes up.  Each of these half-beams is reflected from a
flat mirror back along the incoming path.  In the animation, the
reflected beam is shifted slightly toward the viewer -- rotate the
image to see the difference clearly.  The reflected beams then
recombine at the splitter, resulting in an outgoing downward beam that
is the sum of the two beams.  The animation begins in a world with no
ether.  If you click on the checkbox below the animation, you can turn
on an ether, and the slider will rotate the direction of the relative
motion between the ether and the apparatus.  Waves going with the
ether are stretched out, while waves heading upstream are squished.  By
rotating the direction, you can see that the outgoing beam varies
dramatically based on the angle.  No such dependence was ever
observed.
```



## Book Overview

Example of how to include an image in a figure.


```{image} images/SR_LC.png
:alt: lightcurveSR
:class: bg-primary mb-1
:width: 700px
:align: center
:name: srlc
```
