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
alias fname='../data/Yu18candidates/ktwo211397844-unofficial-tpf.fits'
./make_mask $fname --r=4 --m='round' -s
```

To do aperture photometry with given mask and save as fits file
```shell
$ ./tpf2lc $fname --mask='round' -s -v --o='reduced'
```
