"""

General parsers for QC output of pipeline file, generally pass a handle to the file you want to parse, returns a dict
containing
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

    return {sep.join(os.path.basename(file_name).split(sep)[0: num_seps]): file_name for file_name in files}


def rnaseq_metrics(analysis_dir, num_seps=1, sep="."):
    """

    Generates RNA-seq metrics

    analysis dir, directory to pull information from
    num_seps number of seperators to join back to get the name of the item
    sep: sperator to split / join on to get full name

    """
    
    nrf_files = glob.glob(os.path.join(analysis_dir, "*.NRF.metrics"))
    cutadapt_files = glob.glob(os.path.join(analysis_dir, "*.adapterTrim.metrics"))
    rmrep_files = glob.glob(os.path.join(analysis_dir, "*rmRep.metrics"))
    
    star_files = glob.glob(os.path.join(analysis_dir, "*rmRep.bamLog.final.out"))
    if len(star_files) == 0:  #hack for old data
        star_files = glob.glob(os.path.join(analysis_dir, "*rmRep.samLog.final.out"))
    if len(star_files) == 0:  #Hack for new data
        star_files = glob.glob(os.path.join(analysis_dir, "*.bamLog.final.out"))
    nrf_names = get_names(nrf_files, num_seps, sep) 
    cutadapt_names = get_names(cutadapt_files, num_seps, sep)
    rmrep_names = get_names(rmrep_files, num_seps, sep)
    star_names = get_names(star_files, num_seps, sep) 
    
    nrf_df = pd.DataFrame({name: parse_nrf_file(nrf_file) for name, nrf_file in nrf_names.items()}).transpose()
    cutadapt_df = pd.DataFrame({name: parse_cutadapt_file(cutadapt_file) for name, cutadapt_file in cutadapt_names.items()}).transpose()
    rmrep_df = pd.DataFrame({name: parse_rmrep_file(rmrep_file) for name, rmrep_file in rmrep_names.items()}).transpose()
    star_df = pd.DataFrame({name: parse_star_file(star_file) for name, star_file in star_names.items()}).transpose()
    
    combined_df = pd.merge(cutadapt_df, star_df, left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, rmrep_df, left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, nrf_df, left_index=True, right_index=True, how="outer")

    #Rename columns to be useful
    combined_df = combined_df.rename(columns={"Processed bases": "Input Bases",
                                              "Processed reads": "Input Reads",
                                              "Number of input reads": "Reads Passing Quality Filter",
                                              "Uniquely mapped reads number": "Uniquely Mapped Reads",
           })

    #Make some useful stats
    try:
        combined_df["Percent Repetative"] = 1 - (combined_df['Reads Passing Quality Filter'] / combined_df['Input Reads'].astype(float))
    except ZeroDivisionError:
        pass
    except KeyError:
        print "cutadapt file maybe be broken, ignoring calculation"
        pass

    #combined_df["Repetative Reads"] = (combined_df['Input Reads'] - combined_df['Reads Passing Quality Filter']).astype(int)

    #Get Rid of worthless metrics
    combined_df = combined_df.drop(["Finished on",
                                    "Mapping speed, Million of reads per hour",
                                    "Started job on",
                                    "Started mapping on"], axis=1)

    return combined_df 


def clipseq_metrics(analysis_dir, iclip=False, single_end=False, num_seps=None, sep=".",
                    percent_usable=.01, number_usable=500000, frip=.05):

    """
    
    Reports all clip-seq metrics in a given analysis directory (this is fragile for now, outputs must follow gabes naming clipseq pipeline /
    naming convetions

    
    """
    if num_seps is None:
        num_seps = 2 if iclip else 1

    cutadapt_round2_files = glob.glob(os.path.join(analysis_dir, "*.adapterTrim.round2.metrics"))
    if single_end:
        rm_duped_files = glob.glob(os.path.join(analysis_dir, "*rmDup.bam.out"))
    else:
        rm_duped_files = glob.glob(os.path.join(analysis_dir, "*rmRep.rmDup.metrics"))
    
    rmRep_mapping_files = glob.glob(os.path.join(analysis_dir, "*.adapterTrim.round2.rep.bamLog.final.out"))

    if len(rm_duped_files) == 0:
        rm_duped_files = glob.glob(os.path.join(analysis_dir, "*r2.rmDup.metrics"))
        
    if len(rm_duped_files) == 0:
        rm_duped_files = glob.glob(os.path.join(analysis_dir, "*.rmDup.metrics"))
    
    peaks_files = glob.glob(os.path.join(analysis_dir, "*.peaks.bed"))
    spot_files = glob.glob(os.path.join(analysis_dir, "*peaks.metrics"))

    cutadapt_round2_names = get_names(cutadapt_round2_files, num_seps, sep)

    rm_duped_names = get_names(rm_duped_files, num_seps, sep)
    peaks_names = get_names(peaks_files, num_seps, sep)
    spot_names = get_names(spot_files, num_seps, sep)
    rmRep_mapping_names = get_names(rmRep_mapping_files, num_seps, sep)


    cutadapt_round2_df = pd.DataFrame({name: parse_cutadapt_file(cutadapt_file) for name, cutadapt_file in cutadapt_round2_names.items()}).transpose()
    cutadapt_round2_df.columns = ["{} Round 2".format(col) for col in cutadapt_round2_df.columns]

    rmRep_mapping_df = pd.DataFrame({name: parse_star_file(star_file) for name, star_file in rmRep_mapping_names.items()}).transpose()
    rmRep_mapping_df.columns = ["{} rmRep".format(col) for col in rmRep_mapping_df.columns]

    if single_end:
        rm_duped_df = pd.DataFrame({name: parse_se_umi(rm_duped_file) for name, rm_duped_file in rm_duped_names.items()}).transpose()
    else:
        rm_duped_df = pd.DataFrame({name: parse_rm_duped_metrics_file(rm_duped_file) for name, rm_duped_file in rm_duped_names.items()}).transpose()
    spot_df = pd.DataFrame({name: parse_peak_metrics(spot_file) for name, spot_file in spot_names.items()}).transpose()
    peaks_df = pd.DataFrame({name: {"Num Peaks": len(pybedtools.BedTool(peaks_file))} for name, peaks_file in peaks_names.items()}).transpose()
    combined_df = rnaseq_metrics(analysis_dir, num_seps, sep)

    combined_df = pd.merge(combined_df, cutadapt_round2_df, left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, rm_duped_df, left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, spot_df, left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, peaks_df, left_index=True, right_index=True, how="outer")
    combined_df = pd.merge(combined_df, rmRep_mapping_df, left_index=True, right_index=True, how="outer")

    combined_df['Uniquely Mapped Reads'] = combined_df['Uniquely Mapped Reads'].astype(float)
    try:
        combined_df["Percent Usable / Input"] = (combined_df['Usable Reads'] / combined_df['Uniquely Mapped Reads'])
        combined_df["Percent Usable / Mapped"] = (combined_df['Usable Reads'] / combined_df['Input Reads'])
        combined_df['Passed QC'] = (combined_df['Usable Reads'] > number_usable) & (combined_df['Percent Usable / Mapped'] > percent_usable)

    except ZeroDivisionError:
        pass

    return combined_df

def parse_rmrep_file(rmrep_file):
    try:
        df = pd.read_table(rmrep_file, header=None, sep=" ", index_col=0, names=["element", "repetitive_count"])
    except Exception as e:
        print rmrep_file
        raise e
    return df.sum().to_dict()


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
        return {'spot': float(file_handle.next())}
    
def parse_nrf_file(nrf_file):
    with open(nrf_file) as nrf_file:
        try:
            names = nrf_file.next().strip().split()
            values = nrf_file.next().strip().split()
            return {name: float(value) for name, value in zip(names, values)}
        except:
            return {}
    
def parse_rm_duped_metrics_file(rmDup_file):
    try:
        df = pd.read_csv(rmDup_file, sep="\t")
        
        return {"total_count": sum(df.total_count),
                    "removed_count": sum(df.removed_count),
                    "Usable Reads": sum(df.total_count) - sum(df.removed_count)}
    except Exception as e:
        print e
        return {"total_count": None,
                    "removed_count": None,
                    "Usable Reads": None}

def parse_se_umi(umi_file):
    input_reads = None
    usable_reads = None
    for line in open(umi_file):
        if "Input Reads" in line:
            input_reads = int(line.split(":")[-1].strip())
        if "Number of reads out" in line:
            usable_reads = int(line.split(":")[-1].strip())

    return {"total_count": input_reads,
            "Usable Reads": usable_reads}

def get_cutadapt_version(report):
    with open(report) as file_handle:
            version = file_handle.next()
    try:
        version = version.split()[-4]
    except:
        1
    return int(version.split(".")[1])


def parse_cutadapt_file(report):
    if os.path.getsize(report) == 0:
        return

    old_cutadapt = get_cutadapt_version(report) <= 8
    if old_cutadapt:
        return parse_old_cutadapt_file(report)
    else:
        return parse_new_cutadapt_file(report)


def parse_old_cutadapt_file(report):
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


def get_number_and_percent(line):
    line = [x.strip() for x in line.strip().split(":")]
    line = [line[0]] + line[1].split()
    line[2] = float(line[2][1:-2])
    line[1] = int(line[1].replace(",", ""))
    return line


def get_number(line):
    line = [x.strip() for x in line.strip().split(":")]
    line[1] = int(line[1].replace(",", ""))
    return line


def strip_bp(line):
    return line.replace("bp", "")


def remove_header(file_handle):
    """ for both SE and PE output removes header unifromly from cutadapt metrics file"""
    file_handle.next()
    file_handle.next()
    file_handle.next()
    file_handle.next()
    file_handle.next()
    file_handle.next()
    file_handle.next()
    #print foo.next()


def parse_new_cutadapt_file(report):
    report_dict = {}
    try:
        with open(report) as file_handle:
            remove_header(file_handle)
            processed_reads = get_number(file_handle.next())
            paired_file = processed_reads[0] == 'Total read pairs processed'
            if paired_file:
                r1_adapter = get_number_and_percent(file_handle.next())
                r2_adapter = get_number_and_percent(file_handle.next())
            else:
                adapter = get_number_and_percent(file_handle.next())

            too_short = get_number_and_percent(file_handle.next())
            written = get_number_and_percent(file_handle.next())
            file_handle.next()

            bp_processed = get_number(strip_bp(file_handle.next()))
            if paired_file:
                r1_bp_processed = get_number(strip_bp(file_handle.next()))
                r2_bp_processed = get_number(strip_bp(file_handle.next()))

            bp_quality_trimmed = get_number_and_percent(strip_bp(file_handle.next()))
            if paired_file:
                r1_bp_trimmed = get_number(strip_bp(file_handle.next()))
                r2_bp_trimmed = get_number(strip_bp(file_handle.next()))

            bp_written = get_number_and_percent(strip_bp(file_handle.next()))
            if paired_file:
                r1_bp_written = get_number(strip_bp(file_handle.next()))
                r2_bp_written = get_number(strip_bp(file_handle.next()))

    except Exception as e:
        print e
        print report
        return report_dict

    report_dict['Processed reads'] = processed_reads[1]
    if paired_file:
        report_dict["Read 1 with adapter"] = r1_adapter[1]
        report_dict["Read 1 with adapter percent"] = r1_adapter[2]
        report_dict["Read 2 with adapter"] = r2_adapter[1]
        report_dict["Read 2 with adapter percent"] = r2_adapter[2]
        report_dict['Read 1 basepairs processed'] = r1_bp_processed[1]
        report_dict['Read 2 basepairs processed'] = r2_bp_processed[1]
        report_dict['Read 1 Trimmed bases'] = r1_bp_trimmed[1]
        report_dict['Read 2 Trimmed bases'] = r2_bp_trimmed[1]
        report_dict['Read 1 {}'.format(bp_written[0])] = r1_bp_written[1]
        report_dict['Read 2 {}'.format(bp_written[0])] = r2_bp_written[1]
    else:
        report_dict['Reads with adapter'] = adapter[1]
        report_dict['Reads with adapter percent'] = adapter[2]


    report_dict['Too short reads'] = too_short[1]
    report_dict['Reads that were too short percent'] = too_short[2]
    report_dict['Reads Written'] = written[1]
    report_dict['Reads Written perccent'] = written[2]
    report_dict['Processed bases'] = bp_processed[1]
    report_dict['Trimmed bases'] = bp_quality_trimmed[1]
    report_dict['Trimmed bases percent'] = bp_quality_trimmed[2]
    report_dict[bp_written[0]] = bp_written[1]
    report_dict["{} percent".format(bp_written[0])] = bp_written[2]

    return report_dict


def commas(x):
    try:
        return "{:,}".format(int(x))
    except:
        return "0"
