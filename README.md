Cinoas
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/RenkeHuang/cinoas/workflows/CI/badge.svg)](https://github.com/RenkeHuang/cinoas/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/RenkeHuang/cinoas/branch/main/graph/badge.svg)](https://codecov.io/gh/RenkeHuang/cinoas/branch/main)

The **Cinoas** package implements an active space selection scheme based on state-averaged configuration interaction singles natural orbital (CINO) occupations for multireference calculations of valence and core excited states.

You can learn more about theories and functionalities in the [documentation](under construction). 

Installation
------------

```bash
git clone https://github.com/RenkeHuang/cinoas.git
cd cinoas
pip install -e .
```

This procedure will register `cinoas` within `pip` and you should be able to see `cinoas` listed by calling
```bash
pip list
```

You can run tests to make sure the package is installed successfully
```bash
cd cinoas/tests
pytest -v
```

## Requirements
* numpy
* pytest
* [Psi4](https://github.com/psi4/psi4)

Getting Started
---------------

```python

```


### Copyright

Copyright (c) 2023, Renke Huang


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
