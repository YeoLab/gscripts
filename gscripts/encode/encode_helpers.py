import datetime
import functools 
import json
import os

import gspread
from oauth2client.client import SignedJwtAssertionCredentials
import pandas as pd
import pybedtools 
import pysam
import numpy as np

def _convert_color_to_list(rgb_string):
    return [float(color_value) for color_value in rgb_string[1:-1].split(",")]

def get_lab_manifest():
    json_key = json.load(open("/home/gpratt/ipython_notebook/public clip-588adbc137f3.json"))
    scope = ['https://spreadsheets.google.com/feeds']

    credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
    gc = gspread.authorize(credentials)
    
    sht1 = gc.open_by_url("https://docs.google.com/spreadsheets/d/1ZU2mQh54jentqvhR_oMnviLGWR8Nw_x338gULzKjNDI/edit#gid=0")
    ws = sht1.worksheet("Sheet1")
    list_of_lists = ws.get_all_values()
    manifest = pd.DataFrame(list_of_lists[1:], columns=list_of_lists[0])
    manifest.is_encode = manifest.is_encode == "TRUE"
    manifest.is_4000 = manifest.is_4000 == "TRUE"

    manifest.method_Paper_flag = manifest.Method_Paper_flag == "TRUE" 
    manifest['exp_id'] = manifest.ENCODE_ID.apply(lambda x: x.split("_")[0])
    manifest = manifest.drop(u'', axis=1) #Drops empty columns to try and get rid of a bug
    return manifest

def get_rbp_color_chooser():
    json_key = json.load(open("/home/gpratt/ipython_notebook/public clip-588adbc137f3.json"))
    scope = ['https://spreadsheets.google.com/feeds']

    credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
    gc = gspread.authorize(credentials)

    sht1 = gc.open_by_url("https://docs.google.com/spreadsheets/d/138x3BU5hRsMUGEooVmuLRy1HbZYYg8Z28SlTK-neVJI/edit#gid=0")
    ws = sht1.worksheet("Sheet1")
    list_of_lists = ws.get_all_values()
    manifest = pd.DataFrame(list_of_lists[1:], columns=list_of_lists[0])
    manifest = manifest.set_index("RBP")
    manifest['rgb'] = manifest.rgb.apply(_convert_color_to_list)
    

    return manifest

def get_mapped_reads(fn):
    try:
        return pysam.Samfile(fn).mapped
    except ValueError as e:
        print e, fn
        return np.nan
    except Exception as e:
        print e, fn
        return np.nan

def make_and_filter_clipper(fn, l2fc, pval, out_dir="/home/gpratt/projects/encode/analysis/peak_reanalysis_v14/", reset=False):
    bedtool = pybedtools.BedTool(fn)
    
    filter_data_inst = functools.partial(filter_data, l2fc=l2fc, pval=pval)
    out_file = os.path.join(out_dir, os.path.basename(bedtool.fn) + "l2fc_{}_pval_{}.clipper.bed".format(l2fc, pval))
    if not os.path.exists(out_file) or reset:
        bedtool = bedtool.filter(filter_data_inst).sort().saveas(out_file)

    return out_file

def make_clipper_ish(interval):
    interval.name = interval[7]
    interval[6] = interval.start
    interval[7] = interval.stop
    return interval

def filter_data(interval, l2fc, pval):
    #col4 is -log10 p-val
    #col5 is -log2 fold enrichment

    #This is the standard one 
    return (float(interval[3]) >= pval) and (float(interval[4]) >= l2fc)

def format_input_norm(row, input_norm_dir):
    uID = row.name[0]
    rep = row.name[-1]
    if rep == "rep1":
        rep_num = "01"
    elif rep == "rep2":
        rep_num = "02"
    else:
        print "error"
    #The one shitty part of this new approach is being careful how I deal with the input norm datasets
    return os.path.join(input_norm_dir, "{0}_{1}.basedon_{0}_{1}.peaks.l2inputnormnew.bed.compressed.bed.annotated".format(uID, rep_num))

def format_date(date):
    month, day, year = date.split("/")
    month = int(month)
    day = int(day)
    year = 2000 + int(year)
    return datetime.date(year, month, day)

def get_merged_data(input_norm_dir="/projects/ps-yeolab3/encode/analysis/Eric_Input_Norm/"):
    """
    Returns merged data from the google spreadsheet, currently only works for me because it assumes my gspead key
    """
    json_key = json.load(open("/home/gpratt/ipython_notebook/public clip-588adbc137f3.json"))
    scope = ['https://spreadsheets.google.com/feeds']
    
    credentials = SignedJwtAssertionCredentials(json_key['client_email'], json_key['private_key'], scope)
    gc = gspread.authorize(credentials)
    
    sht1 = gc.open_by_url("https://docs.google.com/spreadsheets/d/1NMjzbneXf8bGN13K9azoqK2YbMCrjQqfJBKFoCkkrtY/edit#gid=0")
    ws = sht1.worksheet("Sheet1")
    list_of_lists = ws.get_all_values()
    manifest = pd.DataFrame(list_of_lists[1:], columns=list_of_lists[0])
    merged_data = manifest.set_index("uID")
    
    non_clip_columns = list(merged_data.columns.difference(["CLIP_rep1", "CLIP_rep2"]))
    cols_to_put_back = list(non_clip_columns)
    cols_to_put_back.remove("RBP")
    cols_to_put_back.remove("Cell line")

    merged_data = merged_data.set_index(list(non_clip_columns), append=True)
    merged_data.columns = pd.MultiIndex.from_tuples([item.split("_") for item in merged_data.columns], names=['CLIP', "rep"])
    
    merged_data = merged_data.stack()

    for col in cols_to_put_back:
        merged_data[col] = merged_data.index.get_level_values(level=col)
        merged_data.index = merged_data.index.droplevel(col)

    merged_data.submitted = merged_data.submitted.apply(lambda x: False if x.strip() == "FALSE" else True)
    merged_data.peaks_submitable  = merged_data.peaks_submitable.apply(lambda x: False if x.strip() == "FALSE" else True)
    merged_data.family_mapping_submitable = merged_data.family_mapping_submitable.apply(lambda x: False if x.strip() == "FALSE" else True)
    merged_data['generally_submittable'] = merged_data.submitted | merged_data.family_mapping_submitable 
    
    format_input_norm_location = functools.partial(format_input_norm, input_norm_dir=input_norm_dir)
    merged_data['input_norm'] = merged_data.apply(format_input_norm_location, axis=1)
    merged_data = merged_data[merged_data.annotation != "not_qc"]
    
    merged_data['Submitted Date'] = merged_data['Submitted Date'].apply(format_date)
    merged_data.is_qced = merged_data.is_qced.apply(lambda x: False if x == "FALSE" else True)
    return merged_data

def confusion_numbers(threshold, threshold_col="FRiP", true_clasification_col='submitted', df=None):
    passed_frip = df[df[threshold_col]  >= threshold]
    failed_frip = df[df[threshold_col]  < threshold]
    
    tp = float(len(passed_frip[passed_frip[true_clasification_col]]))
    fp = float(len(passed_frip[~passed_frip[true_clasification_col]]))
    fn = float(len(failed_frip[failed_frip[true_clasification_col]]))
    tn = float(len(failed_frip[~failed_frip[true_clasification_col]]))
    
    return tp, fp, fn, tn

def get_best_f_score(threshold_col="FRiP", true_clasification_col='submitted', explore_range=np.arange(0,1, .001), df=None):
    true_positive_array = []
    false_positive_array = []
    threshold_array = []

    max_true_positive_rate = 0

    for threshold in explore_range:

        tp, fp, fn, tn = confusion_numbers(threshold, 
                                           threshold_col,
                                           true_clasification_col=true_clasification_col,
                                           df=df)
        try:
            true_positive_rate = tp / (tp + fn)
            false_positive_rate = fp / (fp + tn)

            precision = tp / (tp + fp)
            recall = tp / (tp + fn)

            f_score = 2 * ((precision * recall) / (precision + recall))
            #print threshold, f_score
        except ZeroDivisionError as e:
            print e
            break
        #lr = true_positive_rate / false_positive_rate
        if max_true_positive_rate < f_score:
            max_true_positive_rate = f_score
            best_threshold = threshold 

        true_positive_array.append(true_positive_rate)
        false_positive_array.append(false_positive_rate)
        threshold_array.append(threshold)
    
    
    print max_true_positive_rate
    print best_threshold
    return true_positive_array, false_positive_array, threshold_array, best_threshold
