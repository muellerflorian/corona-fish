import setuptools

setuptools.setup(
    name="covfish",
    version="0.0.1.1",
    author="Florian MUELLER",
    author_email="muellerf.research@gmail.com",
    description="smfish against cov-2",
    url="",
    packages=setuptools.find_packages(),
    install_requires=['biopython',
                      'pandas',
                      'seaborn',
                      'numpy',
                      'matplotlib',
                      'statsmodels'],
    include_package_data=True,
    classifiers=(
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ),
)
