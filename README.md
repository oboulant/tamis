# tamis

Implementation in C, wrapped in Python of the exact solution by Block Coordinate Descent for fast detection of multiple change-points.

It follows the following paper :

* J.-P. Vert and K. Bleakley, "Fast detection of multiple change-points shared by many signals using group LARS", In J. Lafferty, C. K. I. Williams, J. Shawe-Taylor, R.S. Zemel and A. Culotta (Eds), Advances in Neural Information Processing Systems 23 (NIPS), p.2343-2351, 2010. [[paper]](https://members.cbio.mines-paristech.fr/~jvert/svn/ngs/Lasso/article/groupLARS/nips2010/nips2010.pdf) [[supplementary informations]](https://members.cbio.mines-paristech.fr/~jvert/svn/ngs/Lasso/article/groupLARS/nips2010/supplementary.pdf) [[poster]](https://members.cbio.mines-paristech.fr/~jvert/publi/nips2010poster/poster.pdf)

The Matlab implementation by the authors can be found [here](https://members.cbio.mines-paristech.fr/~jvert/svn/GFLseg/html/).

Details about the pseudo code by the same authors can be found [here](https://hal.archives-ouvertes.fr/hal-00602121).

This repository is the Python binding of the C implementation of [https://github.com/oboulant/block-coordinate-descent](https://github.com/oboulant/block-coordinate-descent).

## Local Build

If you are using `bash` :
```bash
> python -m pip install .[dev,display]
```

If you are using `zsh` :
```zsh
> python -m pip install .\[dev,display\]
```

## Contributing

If you wish to contribute, please install the pre-commit hooks. If you performed a previous local install with `python -m pip install .\[dev,display\]`, then the `pre-commit` package should already be installed. Then, you just have to tell `pre-commit` to install the hooks.

```
> pre-commit install
```

## Run example

```zsh
> python -m examples.example
```

## Run tests

From the top directory :

```
> make test
```
