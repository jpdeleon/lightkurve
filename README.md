# Forked lightkurve

## Installation
To install,
```python
$ conda create -n test python=3.6 && source activate test
(test) $ git clone git@github.com:jpdeleon/lightkurve2.git
(test) $ cd lightkurve2
(test) $ pip install .
```

## Quick test of scripts
To view tpf with stacked flux and superposed aperture mask
```bash
(test) ~ $ cd data/Yu18candidates/
(test) Yu18candidates $ make_mask ktwo211397844-unofficial-tpf.fits --r=4 --m='round' -s
```

To do aperture photometry with given mask and radii=[2,3,4] then save as fits file
```bash
(test) Yu18candidates $ tpf2lc . 3 4 5 --mask='round' -s -v --o='reduced'
```

To check output showing two plots of (1) lightcurve and (2) mask:
```bash
$ cd reduced/
(test) reduced $ plot_lc 211397844_round.fits 3 -v -a -m -s=10
```
