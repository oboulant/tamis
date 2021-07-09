# tamis

Block Coordinate Descent, implemented in C, wrapped in Python, following Bleakley &amp; Vert 2011

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
