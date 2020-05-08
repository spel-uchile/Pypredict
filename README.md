# Pypredict

[![GitHub tag](https://img.shields.io/github/tag/spel-uchile/Pypredict.svg)](https://github.com/spel-uchile/Pypredict/releases)
[![license](https://img.shields.io/github/license/spel-uchile/Pypredict)](https://www.gnu.org/licenses/gpl-3.0.en.html)

This software is a real-time satellite tracker and orbit propagator. It displays any satellite's position and orbital parameters in real time. It can also simulate satellite deployments of one satellite from another with different velocities, considering the mass of the satellites. In the near future it will be able to simulate range-based localization systems using TDOA.

Current version: 3.2.1

![](pypredict/img/Screenshot.png)

### Requirements:

* [cartopy>=0.17](https://github.com/SciTools/cartopy)
* [cython>=0.28](https://github.com/cython/cython)
* [Fiona>=1.8.13.post1](https://github.com/Toblerity/Fiona)
* [geos>=0.2.2](https://github.com/grst/geos)
* [matplotlib>=3.2.1](https://github.com/matplotlib/matplotlib)
* [numpy>=1.18.2](https://github.com/numpy/numpy)
* [Pillow>=7.1.1](https://github.com/python-pillow/Pillow)
* [pykdtree>=1.3.1](https://github.com/storpipfugl/pykdtree)
* [pymongo>=3.10.1](https://github.com/mongodb/mongo-python-driver)
* [pyorbital>=1.5](https://github.com/pytroll/pyorbital)
* [PyQt5>=5.14](https://pypi.org/project/PyQt5/)
* [pyshp>=1.1.4](https://github.com/GeospatialPython/pyshp)
* [sgp4>=2.7](https://github.com/brandon-rhodes/python-sgp4)
* [shapely>=1.5.6](https://github.com/simplegeo/shapely)

### Supported software:

* [Python>=3.6](https://www.python.org/downloads/)
* [Ubuntu>=19.10](https://ubuntu.com/download/desktop)

### Installation instructions (for GNU/Linux):

For Raspbian (Raspberry Pi) you may need to install this first:
```bash
sudo apt-get install python3-cairocffi
sudo apt-get install libatlas-base-dev
```

After this, the rest is the same.

1. Install dependencies:
```bash
sudo sh install.sh
pip3 install --upgrade pip
```
2. Clone this repository:
```bash
git clone https://github.com/spel-uchile/Pypredict.git
```
3. Move to the Pypredict folder:
```bash
cd Pypredict/
```
4. Install Pypredict:
```bash
pip3 install . --no-binary shapely
```

### Usage:

You can search "Pypredict" on your application window, or you can run the following command:
```bash
python3 -m pypredict
```

### Contact:

* Name: Matías Vidal Valladares
* e-mail: matias.vidal.v@gmail.com
