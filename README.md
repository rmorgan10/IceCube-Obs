# IceCube-Obs
Scripts to determine DECam observability of IceCube alerts

## Installation

- Step 0: Download an installation of [miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Step 1: Build a virtual environment with the following command:

```python
conda create --name des-icecube python numpy astropy pyephem matplotlib pandas scipy healpy basemap=1.0.7
```

- Step 2: Activate the virtual environment

```python
conda activate des-icecube
```

- Step 3: Fix errors in the `basemap` source code

The module `basemap` is great for making detailed plots, but recently has been superseded by `cartopy`. After `basemap` stopped being maintained, `matplotlib` depreciated some of the function calls that `basemap` uses. There are two locations where a minor change is required.

Open the file `/Users/<USERNAME>/miniconda/envs/des-icecube/lib/python2.7/site-packages/mpl_toolkits/basemap/__init__.py` in your favorite text editor.

We need to replace all calls of `ax.get_axis_bgcolor()` with `ax.get_fc()`, as described in this [Stack Overflow post](https://stackoverflow.com/questions/50691151/axessubplot-object-has-no-attribute-get-axis-bgcolor). The two lines are 1623 and 1767.

**Line 1623** inside the function `drawmapboundary()`: `fill_color = ax.get_axis_bgcolor()`

Change to `fill_color = ax.get_fc()`

**Line 1767** inside the function `fillcontinents()`: `axisbgc = ax.get_axis_bgcolor()`

Change to `axisbgc = ax.get_fc()`

- Step 4: You should be good to go! Try `python test_icecube_observability.py`