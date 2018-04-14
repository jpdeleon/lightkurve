# Forked lightkurve

## Installation
To install,
```python
$ conda create -n test python=3.6 && source activate test
$ pip install .
```

## Quick test of scripts
To view tpf with stacked flux and superposed aperture mask
```shell
cd lightkurve2/scripts
./make_mask '../data/Yu18candidates/ktwo211397844-unofficial-tpf.fit' --r=4 --m='round' -s
```

To do aperture photometry with given mask and radii=[3,4,5] then save as fits file
```shell
$ ./tpf2lc '../data/Yu18candidates/' 3 4 5 --mask='round' -s -v --o='reduced'
```
