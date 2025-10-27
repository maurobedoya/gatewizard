# gatewizard/utils/logger.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Logging utilities for Gatewizard.

This module provides centralized logging configuration and utilities
for the entire application.
"""

import logging
import sys
from pathlib import Path
from typing import Optional, Dict
from logging.handlers import RotatingFileHandler

# Global logger configuration
_loggers = {}
_configured = False

def setup_logger(
    name: str = "gatewizard",
    level: int = logging.INFO,
    debug: bool = False,
    log_file: Optional[Path] = None,
    max_file_size: int = 10 * 1024 * 1024,  # 10MB
    backup_count: int = 5,
    config: Optional[Dict] = None
) -> logging.Logger:
    """
    Set up the main logger for the application.
    
    Args:
        name: Logger name
        level: Logging level
        debug: Enable debug mode
        log_file: Path to log file (optional, overrides config)
        max_file_size: Maximum log file size in bytes
        backup_count: Number of backup files to keep
        config: Logging configuration dict (from LoggingConfig)
        
    Returns:
        Configured logger instance
    """
    global _configured
    
    # Process configuration if provided
    if config:
        # Override parameters with config values
        if 'log_level' in config:
            level_name = config['log_level'].upper()
            level = getattr(logging, level_name, logging.INFO)
        if 'max_file_size_mb' in config:
            max_file_size = config['max_file_size_mb'] * 1024 * 1024
        if 'backup_count' in config:
            backup_count = config['backup_count']
    
    # Set level based on debug flag (overrides config)
    if debug:
        level = logging.DEBUG
    
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Clear existing handlers if reconfiguring
    if _configured:
        logger.handlers.clear()
    
    # Create formatters
    detailed_formatter = logging.Formatter(
        fmt='%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    simple_formatter = logging.Formatter(
        fmt='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.WARNING if not debug else logging.INFO)  # Only warnings/errors to console by default
    console_handler.setFormatter(simple_formatter)
    logger.addHandler(console_handler)
    
    # File handler - only if explicitly requested
    file_logging_enabled = False
    log_file_path = None
    
    # Check config for file logging settings
    if config and config.get('enable_file_logging', False):
        file_logging_enabled = True
        log_dir = config.get('log_file_path', '')
        log_name = config.get('log_file_name', 'gatewizard.log')
        
        if log_dir:
            log_file_path = Path(log_dir) / log_name
        else:
            log_file_path = Path.cwd() / log_name
    
    # Command line log_file parameter overrides config
    if log_file:
        file_logging_enabled = True
        log_file_path = log_file
    
    # Set up file logging if enabled
    if file_logging_enabled and log_file_path:
        try:
            # Ensure log directory exists
            log_file_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Get log level from config or use DEBUG for file logging
            file_log_level = logging.DEBUG  # File gets everything
            if debug:
                file_log_level = logging.DEBUG
            elif config and 'log_level' in config:
                level_str = config.get('log_level', 'DEBUG').upper()
                file_log_level = getattr(logging, level_str, logging.DEBUG)
            
            # Use larger file size from config if available
            if config and 'max_file_size_mb' in config:
                max_file_size = config.get('max_file_size_mb', 10) * 1024 * 1024
            if config and 'backup_count' in config:
                backup_count = config.get('backup_count', 5)
            
            file_handler = RotatingFileHandler(
                log_file_path,
                maxBytes=max_file_size,
                backupCount=backup_count,
                encoding='utf-8'
            )
            file_handler.setLevel(file_log_level)
            file_handler.setFormatter(detailed_formatter)
            logger.addHandler(file_handler)
            
            logger.info(f"File logging enabled: {log_file_path}")
            
        except Exception as e:
            logger.warning(f"Could not set up file logging: {e}")
    else:
        logger.info("File logging disabled")
    
    # Store in global registry
    _loggers[name] = logger
    _configured = True
    
    logger.info(f"Logger '{name}' configured (level: {logging.getLevelName(level)})")
    
    return logger

def get_logger(name: str = None) -> logging.Logger:
    """
    Get a logger instance.
    
    Args:
        name: Logger name (if None, uses the module name of the caller)
        
    Returns:
        Logger instance
    """
    if name is None:
        # Get the caller's module name
        import inspect
        frame = inspect.currentframe().f_back
        name = frame.f_globals.get('__name__', 'gatewizard')
    
    # Return existing logger or create a child logger
    if name in _loggers:
        return _loggers[name]
    
    # Create child logger from root gatewizard logger
    root_logger = _loggers.get('gatewizard')
    if root_logger:
        logger = root_logger.getChild(name.replace('gatewizard.', ''))
    else:
        # Fallback: create basic logger
        logger = logging.getLogger(name)
        logger.setLevel(logging.INFO)
        
        if not logger.handlers:
            handler = logging.StreamHandler(sys.stdout)
            handler.setFormatter(logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            ))
            logger.addHandler(handler)
    
    _loggers[name] = logger
    return logger

def set_log_level(level: int, logger_name: str = "gatewizard"):
    """
    Set the logging level for a specific logger.
    
    Args:
        level: New logging level
        logger_name: Name of the logger to modify
    """
    logger = _loggers.get(logger_name)
    if logger:
        logger.setLevel(level)
        for handler in logger.handlers:
            handler.setLevel(level)
        logger.info(f"Log level changed to {logging.getLevelName(level)}")

def add_file_handler(
    log_file: Path,
    logger_name: str = "gatewizard",
    level: int = logging.DEBUG
):
    """
    Add a file handler to an existing logger.
    
    Args:
        log_file: Path to log file
        logger_name: Name of the logger
        level: Logging level for the file handler
    """
    logger = _loggers.get(logger_name)
    if not logger:
        logger = get_logger(logger_name)
    
    try:
        # Ensure log directory exists
        log_file.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=10 * 1024 * 1024,  # 10MB
            backupCount=5,
            encoding='utf-8'
        )
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter(
            fmt='%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        ))
        
        logger.addHandler(file_handler)
        logger.info(f"File logging added: {log_file}")
        
    except Exception as e:
        logger.warning(f"Could not add file handler: {e}")

class LoggerContext:
    """Context manager for temporary logging configuration."""
    
    def __init__(self, logger_name: str = "gatewizard", level: int = None):
        self.logger_name = logger_name
        self.new_level = level
        self.old_level = None
        self.logger = None
    
    def __enter__(self):
        self.logger = get_logger(self.logger_name)
        if self.new_level is not None:
            self.old_level = self.logger.level
            self.logger.setLevel(self.new_level)
        return self.logger
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.old_level is not None:
            self.logger.setLevel(self.old_level)

# Convenience functions for common logging patterns
def log_function_call(func):
    """Decorator to log function calls."""
    def wrapper(*args, **kwargs):
        logger = get_logger(func.__module__)
        logger.debug(f"Calling {func.__name__} with args={args}, kwargs={kwargs}")
        try:
            result = func(*args, **kwargs)
            logger.debug(f"{func.__name__} completed successfully")
            return result
        except Exception as e:
            logger.error(f"{func.__name__} failed: {e}", exc_info=True)
            raise
    return wrapper

def log_performance(func):
    """Decorator to log function performance."""
    def wrapper(*args, **kwargs):
        import time
        logger = get_logger(func.__module__)
        
        start_time = time.time()
        try:
            result = func(*args, **kwargs)
            elapsed = time.time() - start_time
            logger.debug(f"{func.__name__} completed in {elapsed:.3f}s")
            return result
        except Exception as e:
            elapsed = time.time() - start_time
            logger.error(f"{func.__name__} failed after {elapsed:.3f}s: {e}")
            raise
    return wrapper