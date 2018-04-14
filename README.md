# Forked lightkurve

## Installation
To install,
```python
$ git clone git@github.com:jpdeleon/lightkurve2.git
$ cd lightkurve2
$ conda create -n test python=3.6 && source activate test
$ pip install .
```

## Quick test of scripts
To view tpf with stacked flux and superposed aperture mask
```shell
~ $ cd data
data $ make_mask Yu18candidates/ktwo211397844-unofficial-tpf.fits --r=4 --m='round' -s
```

To do aperture photometry with given mask and radii=[3,4,5] then save as fits file
```shell
data $ tpf2lc 'Yu18candidates/' 3 4 5 --mask='round' -s -v --o='reduced'
```
