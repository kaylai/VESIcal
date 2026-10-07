import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="VESIcal",
    version="1.2.12",
    author="Kayla Iacovino, Simon Matthews, Penny Wieser",
    author_email="kaylaiacovino@gmail.com",
    description=("A generalized python library for calculating and plotting various things "
                 "related to mixed volatile (H2O-CO2) solubility in silicate melts."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kaylai/VESIcal",
    packages=setuptools.find_packages(),
    install_requires=[
            'pandas>=2.1.4',
            'numpy>=1.26.3',
            'matplotlib>=3.10.7',
            'scipy>=1.11.4,!=1.15.*', # 1.15 breaks on macOS 27
            'sympy>=1.12.1',
            'openpyxl>=3.1.2'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
)
