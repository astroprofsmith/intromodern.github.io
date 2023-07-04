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

# Chapter 3: Four Vectors

## Time as a Dimension

At the end of the last chapter, we argued that if the speed of light
was to be invarient -- to be measured the same in any reference frame
-- then the quantity $$dx^2 + dy^2 + dz^2 - c^2dt^2$$ must also be
invariant across reference frames.  This formulation looks a lot like
a kind of 4 dimensional Pythagorean theorem, except for the minus sign
in the last term.  If it were not for that minus sign, we could well
imagine we were looking at the magnitude of a four-dimensional
displacment vector, just adding time as a new component.  If a
spacial displacement vector is $d\vec{r} = (dx,dy,dz)$, and its
magnitude is $dx^2+dy^2+dz^2$, perhaps it would make sense to have
a four-dimensional vector $dx_4 = (cdt,dx,dy,dz)$, except that
if we were to square each term and add them up, we would get
$+c^2dt^2$, not minus.  We need a way to remember to include that
minus sign.

There are two ways of keeping track of the minus sign.  One is simpler
for the beginner, and the other makes the more complicated material
later easier to manage.  The simpler way is to introduce complex
numbers as a bookkeeping trick.  If $i^2 = -1$, then we can write
the time component of the vector as $icdt$, and its square will be
$-c^2dt^2$.  This method ensures that the minus sign always goes in
the right place.  There is nothing physical or mysterious about the
$i$; it is no more and no less than a handy way to keep track of where
to put the minus sign.  It always shows up in the time component of the
four-vector, and it never shows up in the space components.

```{warning}
If you ever see a time component without an $i$, or find an $i$ in
a space component, then you have made a mistake in your math somewhere.
```

With this tool, we can define a displacement four-vector as
\begin{equation}
[dx_4] =
\begin{bmatrix}
icdt\\
dx\\
dy\\
dz
\end{bmatrix}
\end{equation}
The "size" of this four-dimensional object would then be $$[dx_4]^2 =
dx^2 + dy^2 + dz^2 - c^2dt^2.$$ We need a better word than "size", but
"magnitude" is not a good choice, because we are so used to how
magnitudes work with three-dimensional vectors.  The magnitude of a
three-vector is always positive (unless all three of its components
are zero), but $[dx_4]^2$ could be positive, negative, or zero!
For a displacement four-vector, we call the size the "interval", but
there isn't a better word than "size" for other four vectors.

What a person measures when light propagates is a set of four things:
the three components of the physical displacement and the time
interval it took to make this displacement. These values, when
arranged as above, can be thought of as 4 components of a
vector.  The vector is called the displacement 4-vector for the events.

```{warning}
There is no consensus on what order to put the four components.
Some people put the three space components first and time in the
fourth position.  My own preference is to count from zero and put
time in the zeroth (first) position.  Then the space components
are still 1, 2, and 3, but time is a bit different, just like zero
is a bit different from the positive integers.  Having time first
also makes writing easier.  We often don't need to use $y$ or $z$,
but if you put time last, you always need to write all four.  If you
put time first, you can sometimes only write the first two.  Finally,
if you find yourself in the unusual circumstance of needing to have
more than three spacial dimensions, it is easy to add 4, 5, 6, and so
on, but if time is in space four, then you have space dimensions
1, 2, 3, 5, 6, 7, etc., which is confusing and harder to program in
a computer.
```

In general, a four vector is an object with four components,
one of which gets a minus sign when you square the components,
and the other three have the properties of normal three vectors.
The size of any four vector is the same in any inertial reference
frame, and to convert a four vector from one frame to another,
there is a specific procedure one must follow (see Chapter 4).

The second way involves being aware of a distinction between covariant
and contravariant four-vectors.  The mathematical details are not
necessary at this stage, but suffice to say that the contravariant
version of a four vector should be written as a column, as in Equation
6, only without the $i$.  The standard notion is to use a superscript
greek letter to refer to the components, where the Greek letter could
stand for 0, 1, 2, or 3.  So
\begin{equation}
dx^\alpha =
\begin{bmatrix}
cdt\\
dx\\
dy\\
dz
\end{bmatrix}
\end{equation}
```{warning}
The superscript here does NOT mean "raised to the power of".  In this
notation, $dy$ would be written $dx^2$, but the 2 does not mean squared,
it means the third component of the contravariant four-vector.  If
people are using this notation, you have to be aware from the context
whether they mean a contravariant component or an exponent.
```
There is also a covariant version of the four-vector, which is written
as a row with a subscript Greek letter, but also includes the minus sign: 
$dx_\alpha = [-cdt,dx,dydz]$.  For the four vectors we are considering,
$dx_0 =-dx^0$, but $dx_1=dx^1$ and so forth.  So if you multiply the
two together, according to the rules of multiplying matrices,
$$dx_\alpha dx^\alpha = dx_0dx^0 + dx_1dx^1+dx_2dx^2+dx_3dx^3$$
$$dx_\alpha dx^\alpha = -dx^0dx^0 + dx^1dx^1+dx^2dx^2+dx^3dx^3$$
$$dx_\alpha dx^\alpha = -c^2dt^2 + dx^2+dy^2+dz^2$$
where in that last equation, the 2 does mean squared.

In this notation, whenever there is a Greek letter that appears twice,
once up and once down, that is a shorthand for "multiply across the
row and down the column and add them up".  It plays the same role as a
dummy variable in integration, and therefore does not appear in the
result.

%The number of Greek letters you need to specify a particular
%component is called the "rank".  A four vector only needs one number,
%so it is a first-rank object.  First-rank, but four-dimensional, because
%there are four choices for that one number.  A second rank object could
%be written as a matrix, like
%\begin{equation}
%g_{\alpha\beta} =
%\begin{bmatrix}
%-1 & 0 & 0 &0\\
%0 & 1 & 0 & 0\\
%0 & 0 & 1 & 0\\
%0 & 0 & 0 & 1
%\end{bmatrix}
%\end{equation}
%In this case, if you saw $g_{\alpha\beta}dx^\beta$, the repetition
%of beta would mean multiply across the columns of $g$ and down the
%rows of $dx$, but the $\alpha$ would stay the same, so you would get
%$$dx_\alpha = g_{\alpha\beta}dx^\beta = [-cdt,dx,dy,dz]$$