# Forked lightkurve

## Installation
To install,
```python
$ pip install -e git+git@github.com:jpdeleon/lightkurve2.git@master#egg=lightkurve
```

## Quick test of scripts
To view tpf with stacked flux and superposed aperture mask
```shell
fname = '../data/Yu18candidates/ktwo211397844-unofficial-tpf.fits'
./make_mask $fname --r=4 --m='round' -s
```

To do aperture photometry with given mask and save as fits file
```shell
$ ./tpf2lc $fname --mask='round' -s -v --o='reduced'
```
