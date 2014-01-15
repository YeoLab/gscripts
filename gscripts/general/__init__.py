from _read_sample_info_file import read_sample_info_file
import re


def check_text_in_filename(text, filename):
    return [re.search('failed', line) for line in open(filename)
            if re.search('failed', line)]