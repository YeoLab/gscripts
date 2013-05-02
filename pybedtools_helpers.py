def small_peaks(feature):
    """

    feature - pybedtools feature

    returns center of clipper called peak (the middle of the narrow start / stop)
    
    """
    feature.start = (int(feature[6]) + int(feature[7])) /  2
    feature.stop = (int(feature[6]) + int(feature[7])) /  2
    feature.name = feature.name.split("_")[0]
    return feature

def get_five_prime_end(feature):

    """

    gets 5' end interval in a strand intelgent way (pybedtools implementation of this doesn't work)

    """

    if feature.strand == "+":
        feature.stop = feature.start
    else:
        feature.start = feature.stop
    return feature

def get_three_prime_end(feature):

    """

    gets 3' end interval in a strand intelgent way (pybedtools implementation of this doesn't work)

    """

    if feature.strand == "+":
        feature.start = feature.stop
    else:
        feature.stop = feature.start
    return feature

