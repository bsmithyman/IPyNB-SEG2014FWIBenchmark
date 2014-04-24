# IPython notebook for the Chevron SEG 2014 FWI Benchmark Dataset

This notebook is designed to make it straightforward to get started with the 2014 full-waveform inversion (FWI) blind-test dataset provided by the [SEG][] and [Chevron][]. Please feel free to remix and redistribute this notebook, under the terms of the license.

[View on IPython Notebook Viewer][VIEW]

## Requirements

- `IPython` with notebook support
- `matplotlib`
- `numpy`
- `scipy`
- [`pygeo`][pygeo] `segyread.SEGYFile` (included)

## Setup

Most Linux distributions should make the required packages available automatically. On [Ubuntu][], for example, run

    apt-get install ipython-notebook python-matplotlib python-scipy

To run the notebook, unpack the dataset to the `data` subdirectory (or you may symlink the data and/or modify the first code block of the notebook). Change to the project directory and run

    ipython notebook

To simplify installation on a Mac OS or Windows machine you may want to consider using [Enthought][]'s [Canopy][] distribution or [Continuum Analytics][]'s [Anaconda][] distribution, which contain all of the required packages.

## Static rendered versions

[Safari WebArchive file][ArSafari] | [Chrome HTML complete][ArChrome] | 2014-04-24

## Notes

This notebook should work on [Wakari][] with only minor modifications; you will need to remove the lines

    %pylab
    %matplotlib inline

These can be found in the third code box of the notebook, under **General Setup**. The full dataset fits within [Wakari][]'s free plan, and appears to run without difficulty.

## License information

This IPython notebook is licensed under [Creative Commons Attribution-ShareAlike 2.5 Canada][CCLic]. The SEG 2014 Chevron Full Waveform Inversion Synthetic dataset is licensed separately; see the license agreement included with the data. The SEG-Y library used is `segyread.SEGYFile` from [pygeo][], and is distributed under the [GNU Lesser GPL][LGPL].

[SEG]: http://www.seg.org/seg
[Chevron]: http://www.chevron.com/

[Ubuntu]: http://www.ubuntu.com/
[Enthought]: https://www.enthought.com/
[Canopy]: https://www.enthought.com/products/canopy/
[Continuum Analytics]: http://www.continuum.io/
[Anaconda]: https://store.continuum.io/cshop/anaconda/

[ArSafari]: https://www.dropbox.com/s/4mdpyx1n5cvjky0/ChevronNotebook.webarchive
[ArChrome]: https://www.dropbox.com/s/7jjh77s88htxc7r/ChevronNotebook.tar.gz

[Wakari]: https://www.wakari.io/

[CCLic]: http://creativecommons.org/licenses/by-sa/2.5/ca/
[pygeo]: https://github.com/bsmithyman/pygeo
[LGPL]: https://www.gnu.org/licenses/lgpl.html

[VIEW]: http://nbviewer.ipython.org/github/bsmithyman/IPyNB-SEG2014FWIBenchmark/blob/master/SEG%202014%20Chevron%20Initial.ipynb
