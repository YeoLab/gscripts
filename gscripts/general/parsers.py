"""

General parsers for QC output of pipeline file, generally pass a handle to the file you want to parse, returns a dict containing
all useful information

Currently this isn't standard

"""
import glob
import os

import pandas as pd
import pybedtools

def get_names(files, num_seps, sep):

    """ Given a list of files return that files base name and the path to that file
        files -- list of files
        num_seps -- int number of seperators in to call the real name
        sep -- str seperator to split on
    """

    return {sep.join(os.path.basename(file_name).split(sep)[0 : num_seps]) : file_name for file_name in files}

def rnaseq_metrics(analysis_dir, num_seps=1, sep="."):
    """

    Generates RNA-seq metrics

    analysis dir, directory to pull information from
    num_seps number of seperators to join back to get the name of the item
    sep: sperator to split / join on to get full name

    """
    
    nrf_files = glob.glob(os.path.join(analysis_dir, "*.NRF.metrics"))
    cutadapt_files = glob.glob(os.path.join(analysis_dir, "*.adapterTrim.metrics"))
    star_files = glob.glob(os.path.join(analysis_dir, "*.final.out"))
    

    
    nrf_names = get_names(nrf_files, num_seps, sep) 
    cutadapt_names = get_names(cutadapt_files, num_seps, sep)
    star_names = get_names(star_files, num_seps, sep) 
    
    nrf_df = pd.DataFrame({name : parse_nrf_file(nrf_file) for name, nrf_file in nrf_names.items()}).transpose()
    cutadapt_df = pd.DataFrame({name : parse_cutadapt_file(cutadapt_file) for name, cutadapt_file in cutadapt_names.items()}).transpose()
    star_df = pd.DataFrame({name : parse_star_file(star_file) for name, star_file in star_names.items()}).transpose()
    
    combined_df = pd.merge(cutadapt_df, star_df, left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, nrf_df, left_index=True, right_index=True, how="outer")

    return combined_df 
def clipseq_metrics(analysis_dir, iclip=False, sep="."):

    """
    
    Reports all clip-seq metrics in a given analysis directory (this is fragile for now, outputs must follow gabes naming clipseq pipeline /
    naming convetions

    """

    num_seps = 2 if iclip else 1
    
    rm_duped_files = glob.glob(os.path.join(analysis_dir, "*.rmDup.metrics"))
    peaks_files = glob.glob(os.path.join(analysis_dir, "*.rmDup.sorted.peaks.bed"))
    spot_files = glob.glob(os.path.join(analysis_dir, "*peaks.metrics"))

    rm_duped_names = get_names(rm_duped_files, num_seps, sep)
    peaks_names = get_names(peaks_files, num_seps, sep)
    spot_names = get_names(spot_files, num_seps, sep) 

    rm_duped_df = pd.DataFrame({name : parse_rm_duped_metrics_file(rm_duped_file) for name, rm_duped_file in rm_duped_names.items()}).transpose()
    spot_df = pd.DataFrame({name : parse_peak_metrics(spot_file) for name, spot_file in spot_names.items()}).transpose()
    peaks_df = pd.DataFrame({name : {"Num Peaks" : len(pybedtools.BedTool(peaks_file))} for name, peaks_file in peaks_names.items()}).transpose()
    combined_df = rnaseq_metrics(analysis_dir, num_seps, sep)
    combined_df = pd.merge(combined_df, rm_duped_df, left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, spot_df, left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, peaks_df, left_index=True, right_index=True, how="outer")
    return combined_df

def parse_star_file(star_file_name):
    with open(star_file_name) as star_file:
        star_dict = {}

        star_dict["Started job on"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Started mapping on"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Finished on"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Mapping speed, Million of reads per hour"] = star_file.next().strip().split("|")[1].strip()
        star_file.next()
        star_dict["Number of input reads"] = int(star_file.next().strip().split("|")[1].strip())
        star_dict["Average input read length"] = float(star_file.next().strip().split("|")[1].strip())
        star_file.next()
        star_dict["Uniquely mapped reads number"] = int(star_file.next().strip().split("|")[1].strip())
        star_dict["Uniquely mapped reads %"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Average mapped length"] = float(star_file.next().strip().split("|")[1].strip())
        star_dict["Number of splices: Total"] = int(star_file.next().strip().split("|")[1].strip())
        star_dict["Number of splices: Annotated (sjdb)"] = int(star_file.next().strip().split("|")[1].strip())
        star_dict["Number of splices: GT/AG"] = int(star_file.next().strip().split("|")[1].strip())
        star_dict["Number of splices: GC/AG"] = int(star_file.next().strip().split("|")[1].strip())
        star_dict["Number of splices: AT/AC"] = int(star_file.next().strip().split("|")[1].strip())
        star_dict["Number of splices: Non-canonical"] = int(star_file.next().strip().split("|")[1].strip())
        star_dict["Mismatch rate per base, percent"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Deletion rate per base"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Deletion average length"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Insertion rate per base"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Insertion average length"] = star_file.next().strip().split("|")[1].strip()
        star_file.next()
        star_dict["Number of reads mapped to multiple loci"] = star_file.next().strip().split("|")[1].strip()
        star_dict["% of reads mapped to multiple loci"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Number of reads mapped to too many loci"] = star_file.next().strip().split("|")[1].strip()
        star_dict["% of reads mapped to too many loci"] = star_file.next().strip().split("|")[1].strip()
        star_file.next()
        star_dict["% of reads unmapped: too many mismatches"] = star_file.next().strip().split("|")[1].strip()
        star_dict["% of reads unmapped: too short"] = star_file.next().strip().split("|")[1].strip()
        star_dict["% of reads unmapped: other"] = star_file.next().strip().split("|")[1].strip()
    return star_dict

def parse_peak_metrics(fn):
    with open(fn) as file_handle:
        file_handle.next()
        return {'spot' : float(file_handle.next())}
    
def parse_nrf_file(nrf_file):
    with open(nrf_file) as nrf_file:
        try:
            names = nrf_file.next().strip().split()
            values = nrf_file.next().strip().split()
            return {name : float(value) for name, value in zip(names, values)} 
        except:
            return {}
    
def parse_rm_duped_metrics_file(rmDup_file):
    try:
        total_count = 0
        removed_count = 0 
        df = pd.read_csv(rmDup_file, sep="\t")
        
        return {"total_count" : sum(df.total_count), 
                    "removed_count" : sum(df.removed_count), 
                    "Usable Reads" : sum(df.total_count) - sum(df.removed_count)}
    except Exception as e:
        print e
        return {"total_count" : None, 
                    "removed_count" : None, 
                    "Usable Reads" : None}

def parse_cutadapt_file(report):
    report_dir = {}
    try:
        with open(report) as report:
            report.next() #header
            report.next() #paramaters
            report.next() #max error rate
            report.next() #adapters (known already)
            processed_reads = [x.strip() for x in report.next().strip().split(":")]
            processed_bases = [x.strip() for x in report.next().strip().split(":")]
            trimmed_reads   = [x.strip() for x in report.next().strip().split(":")]
            quality_trimmed = [x.strip() for x in report.next().strip().split(":")]
            trimmed_bases   = [x.strip() for x in report.next().strip().split(":")]
            too_short_reads = [x.strip() for x in report.next().strip().split(":")]
            too_long_reads  = [x.strip() for x in report.next().strip().split(":")]
            total_time      = [x.strip() for x in report.next().strip().split(":")]
            time_pre_read   = [x.strip() for x in report.next().strip().split(":")]
            report_dir[processed_reads[0]] = int(processed_reads[1])
            report_dir[processed_bases[0]] = int(processed_bases[1].split()[0])
            report_dir[trimmed_reads[0]] = int(trimmed_reads[1].split()[0])
            report_dir[quality_trimmed[0]] = int(quality_trimmed[1].split()[0])
            report_dir[trimmed_bases[0]] = int(trimmed_bases[1].split()[0])
            report_dir[too_short_reads[0]] = int(too_short_reads[1].split()[0])
            report_dir[too_long_reads[0]] = int(too_long_reads[1].split()[0])
            report_dir[trimmed_bases[0]] = int(trimmed_bases[1].split()[0])
    except:
            print report
    return report_dir

def commas(x):
    try:
        return "{:,}".format(int(x))
    except:
        return "0"
