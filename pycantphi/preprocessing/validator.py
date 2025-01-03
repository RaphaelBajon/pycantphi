from typing import List, Set
import xarray as xr
from ..core.constants import MethodConstants


class DataValidator:
    """Validates input data requirements and consistency"""
    
    def __init__(self, ds: xr.Dataset):
        self.ds = ds
        self.validate_coordinates()
        self.validate_required_variables()
        
    def validate_required_variables(self, required: List[str] = None) -> None:
        """Check if all required variables are present"""
        if required is None:
            required = MethodConstants.REQUIRED_PARAMETERS
        missing = set(required) - set(self.ds.variables)
        if missing:
            raise ValueError(f"Missing required variables: {missing}")
            
    def validate_coordinates(self) -> None:
        """Validate coordinate variables"""
        required_coords = {'latitude', 'longitude', 'pressure'}
        missing_coords = required_coords - set(self.ds.coords)
        if missing_coords:
            raise ValueError(f"Missing required coordinates: {missing_coords}")
            