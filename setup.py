from setuptools import setup, find_packages


scripts = ['rnaseq/count_tags.py',
           'rnaseq/single_RPKM.py',
           'general/calculate_NRF.py',
           'general/normalize_bedGraph.py',
           'general/make_trackhubs.py',
           'general/count_aligned_from_sam.py',
           'general/negBedGraph.py',
           'general/parsers.py',
           'clipseq/perform_idr.py',
           'clipseq/run_piranha.py',
           #'riboseq/riboseq_coverage.py',
           'mapping/map_paired_with_STAR.py',
           'mapping/sam_to_bam_and_sort.py',
           'miso/submit_miso_pipeline.py',
           'miso/submit_index_gff.py',
           'rnaseq/oldsplice.py',
           'rnaseq/submit_oldsplice.py',
           'rnaseq/oldsplice_gff.py',
           'rnaseq/submit_oldsplice_gff.py',
           'rnaseq/parse_oldsplice.py',
           'rnaseq/submit_parse_oldsplice.py',
           'output_parsers/parseMiso.py']

scripts = map((lambda x: "gscripts/" + x), scripts)

with open("README.rst") as file:
    long_description = file.read()

setup(

    name="gscripts",
    long_description=long_description,
    version="0.1.5",
    packages=find_packages(),


    install_requires=['setuptools',
                      'pysam >= 0.6',
                      'numpy >= 1.5.1 ',
                      'scipy >= 0.11.0',
                      'matplotlib >= 1.1.0',
                      'pybedtools >= 0.5',
                      'scikit-learn >= 0.13.0',
    ],

    setup_requires=["setuptools_git >= 0.3", ],

    scripts=scripts,

    #metadata for upload to PyPI
    author="Gabriel Pratt",
    author_email="gpratt@ucsd.edu",
    description="A set of scripts for analysis of high throughput data",
    license="GPL2",
    keywords="bioinformatics",
    url="https://github.com/gpratt",

    #Other stuff I feel like including here
    include_package_data=True,
    zip_safe=False #True I think
)
