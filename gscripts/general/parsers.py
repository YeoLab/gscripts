"""

General parsers for QC output of pipeline file, generally pass a handle to the file you want to parse, returns a dict containing
all useful information

Currently this isn't standard

"""

def parse_star_file(star_file_name):
    with open(star_file_name) as star_file:
        star_dict = {}

        star_dict["Started job on"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Started mapping on"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Finished on"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Mapping speed, Million of reads per hour"] = star_file.next().strip().split("|")[1].strip()
        star_file.next()
        star_dict["Number of input reads"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Average input read length"] = star_file.next().strip().split("|")[1].strip()
        star_file.next()
        star_dict["Uniquely mapped reads number"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Uniquely mapped reads %"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Average mapped length"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Number of splices: Total"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Number of splices: Annotated (sjdb)"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Number of splices: GT/AG"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Number of splices: GC/AG"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Number of splices: AT/AC"] = star_file.next().strip().split("|")[1].strip()
        star_dict["Number of splices: Non-canonical"] = star_file.next().strip().split("|")[1].strip()
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
        sum(df.total_count) - sum(df.removed_count)
            
        return {"total_count" : sum(df.total_count), 
                    "removed_count" : sum(df.removed_count), 
                    "Usable Reads" : sum(df.total_count) - sum(df.removed_count)}
    except:
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
