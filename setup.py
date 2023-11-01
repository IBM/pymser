import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pymser",
    version="1.0.19",
    author="Felipe Lopes de Oliveira",
    author_email="felipe.lopes.oliveira@ibm.com",
    description="Library to apply the Marginal Standard Error Rule \
for transient regime detection and truncation on Grand Canonical \
Monte Carlo adsorption simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/IBM/pymser",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: BSD License",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.9",
    include_package_data=True,
    install_requires=['numpy',
                      'scipy',
                      'statsmodels'],
    license='BSD 3-Clause License'
)
