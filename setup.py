from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='SpatialViewPy',
      version='0.1.2',
      descriptiopn = 'Visualizing multi-sample spatial transcriptimic data using SpatialView',
      py_module='spatialviewpy/prepare_viz',
      package_dir = {'':'src'},
      classifiers=[
          "Programming Language :: Python :: 3.6",
          "Operating System :: OS Independent",
      ],
      long_description = long_description,
      long_description_content_type = "text/markdown",
      install_requires = [
           "requests >= 2.0.0",
            "numpy >= 1.23.0",
            "pandas >= 2.0.0",
            "scipy >= 1.10",
            "scanpy >= 1.9",
            "anndata >= 0.8",
            "pandas >= 1.5",
            "zipfile2 >= 0.0.12",
            "seaborn >= 0.12.2",
      ],
      url="https://github.com/kendziorski-lab/SpatialViewPy",
      author="Chitrasen Mohanty",
      author_email="mohantychitraen@gmail.com")