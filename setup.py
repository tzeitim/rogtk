import setuptools
from setuptools_rust import Binding, RustExtension

with open("README.md", "r") as fh:
    long_description = fh.read()

packages = setuptools.find_packages()
packages.append('rogtk')
setuptools.setup(
    name="rogtk",
    version="0.0.1",
    author="polivares",
    author_email="pedroy.final@gmail.com",
    description="general tools for genomics in Rust",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tzeitim/rogtk",
    #packages=setuptools.find_packages(),
    packages=packages,
    rust_extensions=[RustExtension("rogtk.rogtk", binding=Binding.PyO3 , debug=False)],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = ['setuptools_rust'],
    python_requires='>=3.6',
)

