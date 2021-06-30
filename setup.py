from setuptools import Extension, setup

ext_modules = [
    Extension(
        "tamis._tamis.ebcd",
        sources=[
            "src/tamis/_tamis/ebcd.pyx",
            "src/tamis/_tamis/ebcd_computation.c",
        ],
    ),
]


if __name__ == "__main__":
    from Cython.Build import cythonize

    setup(
        ext_modules=cythonize(ext_modules, language_level="3"),
    )
