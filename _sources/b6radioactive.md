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
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
```

(halflives)=
# Radioactivity

## Half lives

```{code-cell}
:tags: ["remove-input"]

mybutt = widgets.Button(
    description='Click me',
    disabled=False,
    button_style='', # 'success', 'info', 'warning', 'danger' or ''
    tooltip='Click me',
    icon='check' # (FontAwesome names without the `fa-` prefix)
)

myslide = widgets.IntSlider(
    value=7,
    min=0,
    max=10,
    step=1,
    description='Test:',
    disabled=False,
    continuous_update=False,
    orientation='horizontal',
    readout=True,
    readout_format='d'
)


fig = plt.figure(figsize=(9,5))
plt.plot([0.0],[0.0],'ro',label='Test')
plt.xlabel('Test X')
plt.ylabel('Test Y')
plt.legend()
plt.show()
mybutt
```


## Radiometric Dating




