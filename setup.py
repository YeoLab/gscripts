from setuptools import setup, find_packages

with open("README") as file:
    long_description = file.read()

setup(
    name = "gscripts",
    long_description = long_description,
    version = "0.1.1",
    packages = find_packages(),
    
    
    
    install_requires = ['setuptools', 
                        'pysam >= 0.6',
                        'numpy >= 1.5.1 ',
                        'scipy >= 0.11.0',
                        'matplotlib >= 1.1.0',
                        'pybedtools >= 0.5',
                        'scikit-learn >= 0.13.0',
                        ],
      
    setup_requires = ["setuptools_git >= 0.3",],
    
    #metadata for upload to PyPI
    author = "Gabriel Pratt",
    author_email = "gpratt@ucsd.edu",
    description = "A set of scripts for analysis of high throughput data",
    license = "GPL2",
    keywords = "bioinformatics",
    url = "https://github.com/gpratt",
    
    #Other stuff I feel like including here
    include_package_data = True,
    zip_safe = False #True I think
)
