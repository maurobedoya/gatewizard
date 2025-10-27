"""
Optional dependency handling utilities for Gatewizard.

This module provides utilities for handling optional dependencies
in a consistent way across the application.
"""

import importlib
from typing import Optional, Dict, Any
from functools import wraps

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class OptionalDependencyError(Exception):
    """Exception raised when an optional dependency is missing."""
    pass

def require_optional_dependency(package_name: str, install_command: str = None):
    """
    Decorator to check for optional dependencies before function execution.
    
    Args:
        package_name: Name of the package to check
        install_command: Command to install the package (if different from pip install <package>)
    
    Raises:
        OptionalDependencyError: If the dependency is not available
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if not is_package_available(package_name):
                install_cmd = install_command or f"pip install {package_name}"
                raise OptionalDependencyError(
                    f"{package_name} is required for this operation. "
                    f"Install with: {install_cmd}"
                )
            return func(*args, **kwargs)
        return wrapper
    return decorator

def is_package_available(package_name: str) -> bool:
    """
    Check if a package is available for import.
    
    Args:
        package_name: Name of the package to check
        
    Returns:
        True if package is available, False otherwise
    """
    try:
        importlib.import_module(package_name)
        return True
    except ImportError:
        return False

def safe_import(package_name: str, alternative_name: str = None) -> Optional[Any]:
    """
    Safely import a package, returning None if not available.
    
    Args:
        package_name: Name of the package to import
        alternative_name: Alternative import name (e.g., 'parmed' for package, 'pmd' for import)
        
    Returns:
        The imported module or None if not available
    """
    try:
        module = importlib.import_module(package_name)
        logger.debug(f"Successfully imported {package_name}")
        return module
    except ImportError as e:
        logger.debug(f"Optional dependency {package_name} not available: {e}")
        return None

def get_optional_dependencies_status() -> Dict[str, bool]:
    """
    Get the status of all optional dependencies.
    
    Returns:
        Dictionary mapping package names to availability status
    """
    dependencies = {
        'parmed': 'ParmEd for NAMD conversion functionality',
        'biopython': 'BioPython for enhanced PDB file handling',
        'customtkinter': 'CustomTkinter for modern GUI',
        'matplotlib': 'Matplotlib for plotting and visualization',
    }
    
    status = {}
    for package, description in dependencies.items():
        status[package] = {
            'available': is_package_available(package),
            'description': description
        }
    
    return status

def check_and_warn_missing_dependencies():
    """Log warnings for missing optional dependencies."""
    status = get_optional_dependencies_status()
    
    missing = [pkg for pkg, info in status.items() if not info['available']]
    
    if missing:
        logger.warning(f"Missing optional dependencies: {', '.join(missing)}")
        logger.info("Some features may be limited. Install missing packages if needed.")
    else:
        logger.info("All optional dependencies are available")
