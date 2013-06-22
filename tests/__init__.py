import os

def get_test_dir():
    
    """
    Returns the data directory that contains example files for tests and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'data')


def get_file(fn):
    fn = os.path.join(get_test_dir(), fn)
    
    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % (fn))
    return fn