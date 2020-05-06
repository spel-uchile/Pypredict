# Pypredict

[![GitHub tag](https://img.shields.io/github/tag/spel-uchile/Pypredict.svg)](https://github.com/spel-uchile/Pypredict/releases)
[![license](https://img.shields.io/github/license/spel-uchile/Pypredict)](https://www.gnu.org/licenses/gpl-3.0.en.html)

This software is a real-time satellite tracker and orbit propagator. It displays any satellite's position and orbital parameters in real time. It can also simulate satellite deployments from one satellite to another with different velocities, considering the mass of the satellites. In the near future it will be able to simulate range-based localization systems using TDOA.

Current version: 3.2.1

![](pypredict/img/Screenshot.png)

### Supported software:

* Python 3.7.15
* Matplotlib 3.2.1
* Ubuntu 19.10

### Installation instructions:

1. Install dependencies:
```bash
sudo sh install.sh
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
pip install .
```

### Usage:

You can search "Pypredict" on your application window, or you can run the following command:
```bash
python3 -m pypredict
```

### Contact:

* Name: Matías Vidal Valladares
* e-mail: matias.vidal.v@gmail.com
