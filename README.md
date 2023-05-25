# Jump detection using feature-Based time series clustering for NEMS

This jump detection library is directed towards an application in nanomechanical mass measurments as 
described in `Nanomechanical mass measurements through feature-based time series clustering`. In addition 
to jump detection functionality, the library provides several methods to aid in finding and measuring 
single-event mass adsorption events.


## Quickstart

### Installation

To install this package, first clone the repository:

`git clone https://github.com/NEMS-AI/jump-detection.git`

Change your current working directory to the Jump Detection directory:

`cd jump_detection`


Then, install the package in editable mode using pip:

`pip install --editable .`

### Import the Package
Next, you can import our package in your Python script:

```python
from jump_detection.processor import TimeSeriesProcessor
from jump_detection.plotting import *
```

### Load In Data
In the current implementation, code is expected to be formatted as _.csv_ whose columns
are in the format `time, mode 1, mode 2, ...`.

```python
processor = TimeSeriesProcessor(window_size=100, gap_size=10)
processor.load_data("data.csv")
```

### Step 1: Identify all frequency jumps

As outlined in `Nanomechanical mass measurements through feature-based time series clustering`, 
the first step in the outlined process is to find all frequency jumps in the time series.

```python
processor.process_data()
```

The results of which can be seen by running:

```python
print(r"A total of %s jumps were found: "%len(processor.segments))
```

### Step 2: Reduce to single-event jumps
In the following step, each jump event is summaries into a select number of features. These features are then 
clustered to identify the cluster of single-event jumps.

```python
jump_features = processor.get_all_features()
eps = get_eps(jump_features)
labels = get_class_labels(jump_features, eps)
plt.scatter(jump_features[:,1],jump_features[:,3], c = labels)
```

### Step 3: Output frequency shifts of single-event jumps
Once, the clusters/labels have been found, the output relative frequency shifts can be outputted and reported.

```python
diffs = processor.get_all_diffs()
plt.scatter(diffs[:,0],diffs[:,1], c = labels)
```

## Additinal Examples and Bootstrapping

For more examples and example usage, see the `notebooks` directory.

