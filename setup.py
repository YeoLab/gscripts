from setuptools import setup, find_packages


scripts = [
    'annotations/gencode_annotate_utrs_5p_3p.pl',
    'mirna/miR_splitter.py',
    'clipseq/barcode_collapse.py',
    'clipseq/barcode_collapse_pe.py',
    'clipseq/fix_scores.py',
    'clipseq/run_kasey.py',
    'clipseq/kmer_extractor.py',
    'clipseq/perform_idr.py',
    'clipseq/remove_softclip.py',
    'clipseq/run_piranha.py',
    'clipseq/run_ripseeker.py',
    'clipseq/demultiplex_barcoded_fastq.py',
    'clipseq/demux_paired_end.py',
    'clipseq/sort_fastq.py',
    'editing/create_rna_editing_makefile.py',
    'general/biwwig_corr.py',
    'general/cat_files.py',
    'general/cat_biogem.py',
    'general/calculate_NRF.py',
    'general/count_aligned_from_sam.py',
    'general/make_bigwig_files.py',
    'general/make_trackhubs.py',
    'general/negBedGraph.py',
    'general/normalize_bedGraph.py',
    'general/parsers.py',
    'general/gsnap_index.py',
    'ipython_server/serve_ipython.py',
    'mapping/convert_sam.py',
    'mapping/downsample.py',
    'mapping/genome_generate.py',
    'mapping/sort_bam.py',
    'mapping/index_bam.py',
    'mapping/map_paired_or_single_with_STAR.py',
    'mapping/sam_to_bam_and_sort.py',
    'miso/submit_index_gff.py',
    'miso/submit_miso_pipeline.py',
    'miso/find_nan_miso_events.py',
    'output_parsers/parseMiso.py',
    'pwm/cisbp_to_meme.py',
    #'pwm/rbpdb_to_meme.py',
    'riboseq/riboseq_coverage.py',
    'riboseq/read_filter.py',
    'rnaseq/count_tags.py',
    'rnaseq/make_rnaseqc.py',
    'rnaseq/oldsplice.py', 
    'rnaseq/oldsplice_gff.py',
    'rnaseq/parse_oldsplice.py',
    'rnaseq/sailfish_index.py',
    'rnaseq/single_RPKM.py',
    'rnaseq/submit_oldsplice.py',
    'rnaseq/submit_oldsplice_gff.py',
    'rnaseq/submit_parse_oldsplice.py',
    'rnaseq/sailfish_index.py',
    'rnaseq/sailfish_quant.py'
]

scripts = map((lambda x: "gscripts/" + x), scripts)

with open("README.rst") as file:
    long_description = file.read()

setup(

    name="gscripts",
    long_description=long_description,
    version="0.1.6",
    packages=find_packages(),


    install_requires=['setuptools',
                      'pysam >= 0.6',
                      'numpy >= 1.5.1 ',
                      'scipy >= 0.11.0',
                      'matplotlib >= 1.1.0',
                      'pybedtools >= 0.5',
                      'scikit-learn >= 0.13.0',
                      'matplotlib_venn',
                      'clipper', 'HTSeq',
                      'screed >= 0.9',
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
