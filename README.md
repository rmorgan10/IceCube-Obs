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

(Add information on how to interpret plots)