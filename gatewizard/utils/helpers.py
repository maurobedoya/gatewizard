# gatewizard/utils/helpers.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Helper utilities for Gatewizard.

This module contains various utility functions that are used
throughout the application.
"""

import re
import os
import sys
import time
from pathlib import Path
from typing import Union, Tuple, Optional, List
from datetime import datetime, timedelta

def format_time(seconds: float) -> str:
    """
    Format time duration in a human-readable way.
    
    Args:
        seconds: Time duration in seconds
        
    Returns:
        Formatted time string (e.g., "2h 15m 30s", "45.2s")
    """
    if seconds < 0:
        return "0s"
    
    if seconds < 60:
        return f"{seconds:.1f}s"
    
    minutes = seconds / 60
    if minutes < 60:
        return f"{minutes:.1f}m"
    
    hours = minutes / 60
    if hours < 24:
        if hours >= 1:
            mins = int((minutes % 60))
            return f"{int(hours)}h {mins}m" if mins > 0 else f"{int(hours)}h"
        return f"{minutes:.1f}m"
    
    days = hours / 24
    if days < 7:
        hrs = int((hours % 24))
        return f"{int(days)}d {hrs}h" if hrs > 0 else f"{int(days)}d"
    
    weeks = days / 7
    return f"{weeks:.1f}w"

def format_size(bytes_size: int) -> str:
    """
    Format file size in a human-readable way.
    
    Args:
        bytes_size: Size in bytes
        
    Returns:
        Formatted size string (e.g., "1.5 GB", "256 MB")
    """
    if bytes_size < 0:
        return "0 B"
    
    units = ["B", "KB", "MB", "GB", "TB", "PB"]
    size = float(bytes_size)
    unit_index = 0
    
    while size >= 1024 and unit_index < len(units) - 1:
        size /= 1024
        unit_index += 1
    
    if unit_index == 0:
        return f"{int(size)} {units[unit_index]}"
    else:
        return f"{size:.1f} {units[unit_index]}"

def parse_version(version_string: str) -> Tuple[int, int, int]:
    """
    Parse a version string into major, minor, patch components.
    
    Args:
        version_string: Version string (e.g., "1.2.3")
        
    Returns:
        Tuple of (major, minor, patch) integers
        
    Raises:
        ValueError: If version string is invalid
    """
    # Remove 'v' prefix if present
    if version_string.startswith('v'):
        version_string = version_string[1:]
    
    # Split on dots and take first 3 components
    parts = version_string.split('.')
    
    if len(parts) < 2:
        raise ValueError(f"Invalid version string: {version_string}")
    
    try:
        major = int(parts[0])
        minor = int(parts[1])
        patch = int(parts[2]) if len(parts) > 2 else 0
        
        return (major, minor, patch)
        
    except ValueError:
        raise ValueError(f"Invalid version string: {version_string}")

def validate_email(email: str) -> bool:
    """
    Validate an email address using a simple regex.
    
    Args:
        email: Email address to validate
        
    Returns:
        True if email appears valid, False otherwise
    """
    if not email or not isinstance(email, str):
        return False
    
    # Simple email validation regex
    pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
    return re.match(pattern, email.strip()) is not None

def sanitize_filename(filename: str, replacement: str = "_") -> str:
    """
    Sanitize a filename by removing or replacing invalid characters.
    
    Args:
        filename: Original filename
        replacement: Character to replace invalid characters with
        
    Returns:
        Sanitized filename
    """
    if not filename:
        return "untitled"
    
    # Remove or replace invalid characters for most filesystems
    invalid_chars = r'<>:"/\\|?*'
    sanitized = filename
    
    for char in invalid_chars:
        sanitized = sanitized.replace(char, replacement)
    
    # Remove control characters
    sanitized = re.sub(r'[\x00-\x1f\x7f]', replacement, sanitized)
    
    # Remove leading/trailing whitespace and dots
    sanitized = sanitized.strip(' .')
    
    # Ensure filename is not empty and not a reserved name
    if not sanitized or sanitized.upper() in [
        'CON', 'PRN', 'AUX', 'NUL', 'COM1', 'COM2', 'COM3', 'COM4', 
        'COM5', 'COM6', 'COM7', 'COM8', 'COM9', 'LPT1', 'LPT2', 
        'LPT3', 'LPT4', 'LPT5', 'LPT6', 'LPT7', 'LPT8', 'LPT9'
    ]:
        sanitized = "file"
    
    # Limit length (most filesystems support 255 characters)
    if len(sanitized) > 255:
        name, ext = os.path.splitext(sanitized)
        max_name_length = 255 - len(ext)
        sanitized = name[:max_name_length] + ext
    
    return sanitized

def get_platform_info() -> dict:
    """
    Get information about the current platform.
    
    Returns:
        Dictionary with platform information
    """
    import platform
    
    return {
        'system': platform.system(),
        'platform': platform.platform(),
        'processor': platform.processor(),
        'architecture': platform.architecture(),
        'python_version': platform.python_version(),
        'python_implementation': platform.python_implementation(),
    }

def is_running_in_gui() -> bool:
    """
    Check if the application is running in a GUI environment.
    
    Returns:
        True if GUI environment is detected, False otherwise
    """
    # Check for GUI environment variables
    if sys.platform.startswith('linux'):
        return bool(os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY'))
    elif sys.platform == 'darwin':  # macOS
        return True  # macOS always has GUI capabilities
    elif sys.platform.startswith('win'):  # Windows
        return True  # Windows always has GUI capabilities
    
    return False

def find_executable(name: str, paths: Optional[List[str]] = None) -> Optional[str]:
    """
    Find an executable in the system PATH or specified paths.
    
    Args:
        name: Name of the executable
        paths: Additional paths to search (optional)
        
    Returns:
        Full path to executable if found, None otherwise
    """
    import shutil
    
    # First try the standard PATH
    exe_path = shutil.which(name)
    if exe_path:
        return exe_path
    
    # Try additional paths if provided
    if paths:
        for path in paths:
            full_path = Path(path) / name
            if full_path.is_file() and os.access(full_path, os.X_OK):
                return str(full_path)
            
            # Also try with common executable extensions on Windows
            if sys.platform.startswith('win'):
                for ext in ['.exe', '.bat', '.cmd']:
                    full_path_with_ext = Path(path) / (name + ext)
                    if full_path_with_ext.is_file():
                        return str(full_path_with_ext)
    
    return None

def create_unique_filename(base_path: Union[str, Path], extension: str = "") -> Path:
    """
    Create a unique filename by appending a number if the file already exists.
    
    Args:
        base_path: Base file path (without extension)
        extension: File extension (with or without leading dot)
        
    Returns:
        Unique file path
    """
    base_path = Path(base_path)
    
    # Ensure extension starts with a dot
    if extension and not extension.startswith('.'):
        extension = '.' + extension
    
    # Try the original filename first
    full_path = base_path.with_suffix(extension)
    if not full_path.exists():
        return full_path
    
    # Find a unique filename by appending numbers
    counter = 1
    while True:
        name_with_counter = f"{base_path.stem}_{counter}"
        full_path = base_path.with_name(name_with_counter).with_suffix(extension)
        
        if not full_path.exists():
            return full_path
        
        counter += 1
        
        # Prevent infinite loop
        if counter > 9999:
            # Use timestamp as fallback
            timestamp = int(time.time())
            name_with_timestamp = f"{base_path.stem}_{timestamp}"
            return base_path.with_name(name_with_timestamp).with_suffix(extension)

def retry_on_failure(max_attempts: int = 3, delay: float = 1.0, backoff: float = 2.0):
    """
    Decorator to retry a function on failure.
    
    Args:
        max_attempts: Maximum number of attempts
        delay: Initial delay between attempts in seconds
        backoff: Multiplier for delay on each retry
        
    Returns:
        Decorated function
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            current_delay = delay
            last_exception = None
            
            for attempt in range(max_attempts):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    last_exception = e
                    if attempt < max_attempts - 1:  # Don't sleep on last attempt
                        time.sleep(current_delay)
                        current_delay *= backoff
            
            # If we get here, all attempts failed
            raise last_exception
        
        return wrapper
    return decorator

def safe_cast(value: any, target_type: type, default: any = None):
    """
    Safely cast a value to a target type.
    
    Args:
        value: Value to cast
        target_type: Target type to cast to
        default: Default value if casting fails
        
    Returns:
        Cast value or default if casting fails
    """
    try:
        return target_type(value)
    except (ValueError, TypeError):
        return default

def chunks(lst: List, chunk_size: int):
    """
    Yield successive chunks from a list.
    
    Args:
        lst: List to chunk
        chunk_size: Size of each chunk
        
    Yields:
        Chunks of the specified size
    """
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

def flatten_dict(d: dict, parent_key: str = '', sep: str = '.') -> dict:
    """
    Flatten a nested dictionary.
    
    Args:
        d: Dictionary to flatten
        parent_key: Parent key prefix
        sep: Separator for keys
        
    Returns:
        Flattened dictionary
    """
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def create_directory_robust(directory_path: Union[str, Path], max_retries: int = 3) -> None:
    """
    Robustly create a directory with multiple fallback strategies.
    
    This function handles edge cases that can occur on Windows filesystems,
    especially when directories are manually deleted and recreated quickly.
    
    Args:
        directory_path: Path for the directory to create
        max_retries: Maximum number of retry attempts
        
    Raises:
        OSError: If directory creation fails after all retries
        PermissionError: If insufficient permissions to create directory
    """
    import time
    import stat
    import shutil
    from gatewizard.utils.logger import get_logger
    
    directory_path = Path(directory_path)
    logger = get_logger(__name__)
    
    for attempt in range(max_retries + 1):
        try:
            # Strategy 1: Try normal mkdir with exist_ok=True
            if not directory_path.exists():
                directory_path.mkdir(parents=True, exist_ok=True)
            
            # Verify the directory was created and is writable
            if directory_path.exists() and directory_path.is_dir():
                # Test write access
                test_file = directory_path / ".test_write_access"
                try:
                    test_file.touch()
                    test_file.unlink()
                    logger.debug(f"Directory creation successful: {directory_path}")
                    return  # Success
                except (OSError, PermissionError):
                    logger.warning(f"Directory exists but not writable: {directory_path}")
                    raise PermissionError(f"Directory not writable: {directory_path}")
            
        except FileExistsError:
            # This shouldn't happen with exist_ok=True, but handle it anyway
            logger.warning(f"FileExistsError on attempt {attempt + 1}: {directory_path}")
            
            # Strategy 2: Check if it's actually a file instead of directory
            if directory_path.exists() and not directory_path.is_dir():
                logger.warning(f"Path exists as file, removing: {directory_path}")
                try:
                    # Make sure the file is not read-only
                    if directory_path.exists():
                        directory_path.chmod(0o777)
                    directory_path.unlink()
                    time.sleep(0.1)  # Give filesystem time to update
                except OSError as unlink_error:
                    logger.warning(f"Failed to unlink file: {unlink_error}")
                    # If we can't remove it, try to move it
                    try:
                        backup_path = directory_path.with_name(f"{directory_path.name}.backup.{int(time.time())}")
                        directory_path.rename(backup_path)
                        logger.info(f"Moved conflicting file to: {backup_path}")
                    except OSError as rename_error:
                        logger.error(f"Failed to move conflicting file: {rename_error}")
                        if attempt >= max_retries:
                            raise OSError(f"Cannot remove or move conflicting file: {directory_path}")
                continue  # Retry creation
            
            # Strategy 3: If it exists as directory, verify it's usable
            if directory_path.exists() and directory_path.is_dir():
                return  # It exists and is a directory, we're good
        
        except PermissionError as e:
            logger.warning(f"Permission error on attempt {attempt + 1}: {e}")
            if attempt < max_retries:
                time.sleep(0.1 * (attempt + 1))  # Progressive delay
                continue
            else:
                raise
        
        except OSError as e:
            logger.warning(f"OSError on attempt {attempt + 1}: {e}")
            
            # Strategy 4: Try to handle Windows-specific issues
            if "File exists" in str(e) or e.errno == 17:  # EEXIST
                logger.info(f"Handling 'File exists' error for: {directory_path}")
                
                # Wait a bit for filesystem consistency (Windows sometimes needs this)
                time.sleep(0.1)
                
                # Check current state
                if directory_path.exists():
                    if directory_path.is_dir():
                        logger.info(f"Directory already exists and is valid: {directory_path}")
                        return  # Success - directory exists
                    else:
                        # It's a file, remove it
                        logger.info(f"Removing file at directory path: {directory_path}")
                        try:
                            directory_path.unlink()
                        except OSError:
                            pass
                
                # Try again after cleanup
                if attempt < max_retries:
                    time.sleep(0.1)
                    continue
            
            # For other OSErrors, retry with delay
            if attempt < max_retries:
                time.sleep(0.1 * (attempt + 1))
                continue
            else:
                raise
    
    # If we get here, all retries failed
    raise OSError(f"Failed to create directory after {max_retries + 1} attempts: {directory_path}")