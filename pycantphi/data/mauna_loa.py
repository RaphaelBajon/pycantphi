import os 
import requests
import pandas as pd
import time
import datetime
import urllib
import datetime
import urllib.request
import pandas as pd
from pathlib import Path
import importlib.resources
from typing import Optional, Union, Dict
from functools import lru_cache

def load_maunaloa_pco2():
    """
    Load Mauna Loa annual pCO2 registry
    
    Returns:
    pandas DataFrame with year, pco2, and uncertainty columns
    """
    ml_fname = 'co2_annmean_mlo.txt'
    ml_url = 'ftp://ftp.cmdl.noaa.gov/ccg/co2/trends/' + ml_fname
    ml_alternative_url = 'ftp://aftp.cmdl.noaa.gov/ccg/co2/trends/' + ml_fname
    
    # Determine file directory (same directory as script)
    file_dir = os.path.dirname(os.path.abspath(os.getcwd()))
    file_path = os.path.join(file_dir, ml_fname)
    
    def download_ml_file():
        """Download Mauna Loa CO2 data file"""
        try:
            urllib.request.urlretrieve(ml_url, file_path)
        except Exception:
            try:
                urllib.request.urlretrieve(ml_alternative_url, file_path)
            except Exception as err:
                raise RuntimeError(f'Cannot download pCO2 Mauna Loa {ml_fname} file. Please download it manually.')
    
    # Check if file exists and is not empty
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        print(f'File with MaunaLoa record: {ml_fname} doesn\'t exist or is damaged,\n downloading it from: {ml_url}')
        download_ml_file()
    
    def process_ml_file():
        """Process the Mauna Loa CO2 data file"""
        try:
            # Read the file, skipping comment lines
            df = pd.read_csv(file_path, 
                             comment='#', 
                             delim_whitespace=True, 
                             names=['year', 'pco2', 'unc'], 
                             dtype={'year': int, 'pco2': float, 'unc': float})
            return df
        except Exception as e:
            raise RuntimeError(f'Error processing Mauna Loa file: {e}')
    
    # Process the file
    ret = process_ml_file()
    
    # Check if data is up to date (within last two years)
    current_year = datetime.datetime.now().year
    if ret['year'].max() < (current_year - 1):
        print(f'MaunaLoa record {ml_fname} could be obsolete. Trying to download it from: {ml_url}')
        download_ml_file()
        ret = process_ml_file()
    
    return ret

# TODO: need to improve this version!! multiple loads for the two versions to clear

class MaunaLoaData:
    """Handler for Mauna Loa CO2 data"""
    
    DEFAULT_FILENAME = 'co2_annmean_mlo.txt'
    BASE_URL = 'ftp://ftp.cmdl.noaa.gov/ccg/co2/trends/'
    ALT_URL = 'ftp://aftp.cmdl.noaa.gov/ccg/co2/trends/'
    
    def __init__(self):
        # Get the package's data directory
        self.data_dir = self._get_data_dir()
        
    def _get_data_dir(self) -> Path:
        """Get the package's data directory path using version-compatible approach"""
        try:
            # Try Python 3.9+ approach first
            with importlib.resources.files('pyphi.data') as data_path:
                return data_path
        except (AttributeError, ImportError):
            try:
                # Fallback for Python 3.7+
                return importlib.resources.path('pyphi.data', '').__enter__()
            except (AttributeError, ImportError):
                # Last resort: use package location
                import pyphi
                return os.path.dirname(pyphi.__file__) / 'data'
    
    @property
    def data_file(self) -> Path:
        """Path to data file in package directory"""
        return self.data_dir / self.DEFAULT_FILENAME
    
    def download_data(self) -> None:
        """Download latest data from Mauna Loa website"""
        urls = [
            self.BASE_URL + self.DEFAULT_FILENAME,
            self.ALT_URL + self.DEFAULT_FILENAME
        ]
        
        for url in urls:
            try:
                urllib.request.urlretrieve(url, self.data_file)
                return
            except Exception:
                continue
                
        raise RuntimeError(
            f'Cannot download Mauna Loa data. Please check your internet connection '
            f'or download {self.DEFAULT_FILENAME} manually.'
        )
    
    def process_data(self, file_path: Path) -> pd.DataFrame:
        """Process Mauna Loa data file"""
        try:
            df = pd.read_csv(
                file_path,
                comment='#',
                delim_whitespace=True,
                names=['year', 'pco2', 'unc'],
                dtype={'year': int, 'pco2': float, 'unc': float}
            )
            return df
        except Exception as e:
            raise RuntimeError(f'Error processing Mauna Loa file: {e}')
    
    def is_data_current(self, df: pd.DataFrame) -> bool:
        """Check if data is up to date"""
        current_year = datetime.datetime.now().year
        return df['year'].max() >= (current_year - 1)
    
    @lru_cache(maxsize=1)
    def load_data(self) -> pd.DataFrame:
        """Load Mauna Loa CO2 data"""
        # Check if file exists and process it
        if self.data_file.exists():
            df = self.process_data(self.data_file)
            
            # If data is not current, try to update it
            if not self.is_data_current(df):
                try:
                    print(f'MaunaLoa record could be obsolete. '
                          f'Trying to download it from: {self.BASE_URL}')
                    self.download_data()
                    df = self.process_data(self.data_file)
                except Exception as e:
                    print(f"Warning: Couldn't update data: {e}")
            
            return df
        else:
            # If file doesn't exist, try to download it
            try:
                print(f'File with MaunaLoa record does not exist, '
                      f'downloading it from: {self.BASE_URL}')
                self.download_data()
                return self.process_data(self.data_file)
            except Exception as e:
                raise RuntimeError(
                    f"Could not obtain Mauna Loa data: {str(e)}"
                )

def load_maunaloa_pco2_() -> pd.DataFrame:
    """Load Mauna Loa annual pCO2 registry"""
    return MaunaLoaData().load_data()