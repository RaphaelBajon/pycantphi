import importlib.resources
import pandas as pd
from typing import Union, TextIO
from io import StringIO


def get_package_data_file(filename: str) -> Union[str, TextIO]:
    """
    Get a data file from the package using importlib.resources
    
    Parameters
    ----------
    filename : str
        Name of the file to load from the data directory
        
    Returns
    -------
    Union[str, TextIO]
        File content or file object
    """
    try:
        # For Python >= 3.9
        with importlib.resources.files('pycantphi.data').joinpath(filename).open('r') as f:
            return f.read()
    except Exception:
        # Fallback for older Python versions
        with importlib.resources.path('pycantphi.data', filename) as path:
            return path

def read_sills_data() -> pd.DataFrame:
    """Read the sills data file from the package"""
    try:
        content = get_package_data_file('sills.csv')
        if isinstance(content, str):
            # If we got the content as a string
            return pd.read_csv(StringIO(content))
        else:
            # If we got a path
            return pd.read_csv(content)
    except Exception as e:
        raise RuntimeError(f"Error reading sills data: {str(e)}")