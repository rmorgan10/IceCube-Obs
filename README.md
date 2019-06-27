# IceCube-Obs
Scripts to determine DECam observability of IceCube alerts

## Installation

- Step 0: Download an installation of [miniconda](https://docs.conda.io/en/latest/miniconda.html)
- Step 1: Clone this repository

```
cd <path/to/where/you/want/to/work>
git clone https://github.com/rmorgan10/IceCube-Obs.git
```

- Step 2: Build a virtual environment with the following command:

```
conda create --name des-icecube python=2.7 numpy astropy pyephem matplotlib pandas scipy healpy basemap=1.0.7
```

- Step 3: Activate the virtual environment

```
conda activate des-icecube
```

- Step 4: Fix errors in the `basemap` source code

The module `basemap` is great for making detailed plots, but recently has been superseded by `cartopy`. After `basemap` stopped being maintained, `matplotlib` depreciated some of the function calls that `basemap` uses. There are two locations where a minor change is required.

Open the file `/Users/<USERNAME>/miniconda/envs/des-icecube/lib/python2.7/site-packages/mpl_toolkits/basemap/__init__.py` in your favorite text editor.

We need to replace all calls of `ax.get_axis_bgcolor()` with `ax.get_fc()`, as described in this [Stack Overflow post](https://stackoverflow.com/questions/50691151/axessubplot-object-has-no-attribute-get-axis-bgcolor). The two lines are 1623 and 1767.

**Line 1623** inside the function `drawmapboundary()`: `fill_color = ax.get_axis_bgcolor()`

Change to `fill_color = ax.get_fc()`

**Line 1767** inside the function `fillcontinents()`: `axisbgc = ax.get_axis_bgcolor()`

Change to `axisbgc = ax.get_fc()`

- Step 5: You should be good to go! Try `python test_icecube_observability.py`

## Usage

These scripts are designed to determine the observability of a single IceCube alert from CTIO. The file `test_icecube_observability.py` is where the properties of an alert can be entered.

To test the observability of an alert:

- Step 1: Open `test_icecube_observability.py` in your favorite text editor
- Step 2: Comment out existing `event_dict`s
- Step 3: Create a new `event_dict` with the properties of the alert you are interested in. It should look like this:

```python
event_dict = {'766165': {'ra': 65.7866,
                         'dec': -37.4431,
                         'eventid': 766165,
                         'time': '2019/5/25 00:27:07'}}
```

- Step 4: Save your changes to the file
- Step 5: Run the scripts with the command `python test_icecube_observability.py`

## Output

The scripts will both print diagnostic information to the terminal and produce plots for assessing the observability of the IceCube alert.

**Terminal Output Example**

```
Event
  Event ID = 132707
  (ra, dec) = (342.7772, 10.0547)
Date
  Now = 2019/6/27 17:07:12 (UTC)
  Search time = 2019/6/27 17:07:09 (UTC)
  Optimal time = 2019/6/28 09:10:36 (UTC)
  Airmass at optimal time = 1.31
Sun
  Angular separation = 107.57 (deg)
  Next rising = 2019/6/28 11:39:49 (UTC)
  Next setting = 2019/6/27 21:52:54 (UTC)
Moon
  Illumination = 0.22
  Angular separation = 55.76 (deg)
  Next rising = 2019/6/28 07:16:45 (UTC)
  Next setting = 2019/6/27 18:11:09 (UTC)
  Next new moon = 2019/7/2 19:16:11 (UTC)
  Next full moon = 2019/7/16 21:38:11 (UTC)
Galactic
  (l, b) = (80.7101, -42.7402)
  E(B-V) = 0.02
```

The main function of the observing scripts is to determine the optimal time to perform observations based on the Sun, Moon, and rotation of the Earth. That optimal time is shown near the top of the output. At this optimal time, details on the position of the Sun, position and illumination of the Moon, and extinction casused by dust in the Milky Way (Galactic) are displayed.

**All-Sky Map Example**

![](./output/skymap.example)

The black band is areas of the sky where the plane of the Milky Way causes a lot of light extinction. The blue outline shaped like a tank is the DES footprint, a large part of the sky that DES has spent about 6 years covering. The neutino is represented by the red dot, the Sun is represented by the yellow dot, and the Moon is shown as a grey/purple dot.

**Spherical Rendering of the Skymap Example**

![](./output/ortho.example)

The center of the green circle is the zenith, or point on the sky that is directly overhead, from CTIO in Chile. The green circle itself is a limiting angle where we do not want to observe events outside of it since we would have to look through a lot of the Earth's atmosphere. For reference, the amount of atmosphere you look through is known as the airmass. Airmass = 1.0 is directly overhead at the zenith, and the green circle is at airmass = 2.0. Mathematically, the airmass is the secant of the angle between the line of sight and the zenith.

**Airmass Plot Example**

![](./output/airmass.example)

The black curve in the plot is the airmass as a function of time. This plot is useful for determining whether the airmass is increasing or decreasing with time, and how rapidly it is changing, at the point of starting observations.

