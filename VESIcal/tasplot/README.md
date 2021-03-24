# README #

tasplot is a Python module that adds the names of volcanic rock types to a plot of total alkali versus silica data

### Installation ###

* Clone or [Download](https://bitbucket.org/jsteven5/tasplot/downloads) the git repository
* Install with `> pip install -e .` in the repository directory (your
  local changes will be reflected when you import the module

Alternatively, you can install the master version straight from
Bitbucket: simply install using `pip` with the https repository URL:

`> pip install git+https://bitbucket.org/jsteven5/tasplot.git`

This will automatically link the module onto your path, allowing you to
use it from any python code within your environment.

### Running ###

*  The script assumes that you are making a plot with Matplotlib.
*  It is used as follows:


```
#!python

import tasplot
import matplotlib.pyplot as plt

silica = [50, 60, 70]
total_alkalis = [4, 5, 6]

fig = plt.figure()
ax1 = fig.subplot(111)
tasplot.add_LeMaitre_fields(ax1)
plt.plot(silica, total_alkalis)
```



### Contribution guidelines ###

* Please feel free to suggest or contribute improvements.  You can contact me via Bitbucket, or [@volcan01010 on Twitter](https://www.twitter.com/volcan01010).
