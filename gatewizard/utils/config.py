# gatewizard/utils/config.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Configuration management for Gatewizard.

This module handles application configuration, user preferences,
and settings persistence.
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional, Union
from dataclasses import dataclass, asdict, field

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

@dataclass
class GuiConfig:
    """GUI-specific configuration."""
    theme: str = "dark"
    color_scheme: str = "blue"
    window_width: int = 1250
    window_height: int = 800
    target_screen: int = 0
    remember_window_position: bool = True
    last_window_x: Optional[int] = None
    last_window_y: Optional[int] = None
    # Default font scale (Large by default)
    font_scale: float = 1.2

@dataclass 
class SystemBuilderConfig:
    """System builder configuration."""
    water_model: str = "opc"
    protein_ff: str = "ff19SB"
    lipid_ff: str = "lipid21"
    preoriented: bool = True
    parametrize: bool = True
    salt_concentration: float = 0.15
    default_cation: str = "K+"
    default_anion: str = "Cl-"
    default_ph: float = 7.0
    use_pdb2pqr: bool = False

@dataclass
class PropkaConfig:
    """Propka-specific configuration."""
    version: str = "3"
    default_ph: float = 7.0
    auto_apply_states: bool = False

@dataclass
class LoggingConfig:
    """Logging configuration."""
    enable_file_logging: bool = False
    log_file_path: str = ""
    log_file_name: str = "gatewizard.log"
    log_level: str = "INFO"
    max_file_size_mb: int = 10
    backup_count: int = 5
    auto_create_directories: bool = True

@dataclass
class PathsConfig:
    """Paths and directories configuration."""
    working_directory: str = ""
    last_pdb_directory: str = ""
    temp_directory: str = ""
    log_directory: str = ""

@dataclass
class Config:
    """Main configuration class."""
    gui: GuiConfig = field(default_factory=GuiConfig)
    system_builder: SystemBuilderConfig = field(default_factory=SystemBuilderConfig)
    propka: PropkaConfig = field(default_factory=PropkaConfig)
    paths: PathsConfig = field(default_factory=PathsConfig)
    logging: LoggingConfig = field(default_factory=LoggingConfig)
    version: str = "0.1.0"
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary with Path objects converted to strings."""
        def convert_paths(obj):
            """Recursively convert Path objects to strings."""
            if isinstance(obj, Path):
                return str(obj)
            elif isinstance(obj, dict):
                return {k: convert_paths(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_paths(item) for item in obj]
            elif hasattr(obj, '__dict__'):
                # Handle dataclass and other objects with __dict__
                return convert_paths(obj.__dict__)
            else:
                return obj
        
        # Convert to dict first, then handle Path objects
        result = asdict(self)
        
        # Exclude log_file_path from being saved (it should be auto-set on each launch)
        if 'logging' in result and 'log_file_path' in result['logging']:
            result['logging']['log_file_path'] = ""  # Reset to empty so it uses current directory
        
        return convert_paths(result)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Config':
        """Create configuration from dictionary."""
        try:
            # Handle nested structures
            gui_data = data.get('gui', {})
            system_builder_data = data.get('system_builder', {})
            propka_data = data.get('propka', {})
            paths_data = data.get('paths', {})
            logging_data = data.get('logging', {})
            
            return cls(
                gui=GuiConfig(**gui_data),
                system_builder=SystemBuilderConfig(**system_builder_data),
                propka=PropkaConfig(**propka_data),
                paths=PathsConfig(**paths_data),
                logging=LoggingConfig(**logging_data),
                version=data.get('version', '0.1.0')
            )
        except Exception as e:
            logger.warning(f"Error loading config from dict: {e}")
            return cls()  # Return default config

class ConfigManager:
    """Configuration manager for loading and saving settings."""
    
    def __init__(self, config_dir: Optional[Path] = None):
        """
        Initialize configuration manager.
        
        Args:
            config_dir: Directory for configuration files (uses default if None)
        """
        self.config_dir = config_dir or self._get_default_config_dir()
        self.config_file = self.config_dir / "config.json"
        self._ensure_config_dir()
    
    def _get_default_config_dir(self) -> Path:
        """Get default configuration directory."""
        if os.name == 'nt':  # Windows
            base_dir = Path(os.environ.get('APPDATA', Path.home()))
        else:  # Unix-like
            base_dir = Path.home() / '.config'
        
        return base_dir / 'gatewizard'
    
    def _ensure_config_dir(self):
        """Ensure configuration directory exists."""
        try:
            self.config_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logger.warning(f"Could not create config directory: {e}")
    
    def load_config(self) -> Config:
        """
        Load configuration from file.
        
        Returns:
            Configuration object (default if file doesn't exist or is invalid)
        """
        if not self.config_file.exists():
            logger.info("No config file found, using defaults")
            return Config()
        
        try:
            with open(self.config_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            config = Config.from_dict(data)
            logger.info(f"Configuration loaded from {self.config_file}")
            return config
            
        except (json.JSONDecodeError, IOError) as e:
            logger.warning(f"Error loading config file: {e}")
            return Config()
    
    def save_config(self, config: Config) -> bool:
        """
        Save configuration to file.
        
        Args:
            config: Configuration object to save
            
        Returns:
            True if successful, False otherwise
        """
        try:
            self._ensure_config_dir()
            
            with open(self.config_file, 'w', encoding='utf-8') as f:
                json.dump(config.to_dict(), f, indent=2, ensure_ascii=False)
            
            logger.info(f"Configuration saved to {self.config_file}")
            return True
            
        except Exception as e:
            logger.error(f"Error saving config file: {e}")
            return False
    
    def backup_config(self) -> Optional[Path]:
        """
        Create a backup of the current configuration.
        
        Returns:
            Path to backup file if successful, None otherwise
        """
        if not self.config_file.exists():
            return None
        
        try:
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_file = self.config_dir / f"config_backup_{timestamp}.json"
            
            import shutil
            shutil.copy2(self.config_file, backup_file)
            
            logger.info(f"Configuration backed up to {backup_file}")
            return backup_file
            
        except Exception as e:
            logger.error(f"Error creating config backup: {e}")
            return None
    
    def reset_config(self) -> bool:
        """
        Reset configuration to defaults.
        
        Returns:
            True if successful, False otherwise
        """
        try:
            # Backup current config if it exists
            if self.config_file.exists():
                self.backup_config()
            
            # Save default config
            default_config = Config()
            return self.save_config(default_config)
            
        except Exception as e:
            logger.error(f"Error resetting config: {e}")
            return False

# Global configuration manager instance
_config_manager: Optional[ConfigManager] = None
_current_config: Optional[Config] = None

def get_config_manager() -> ConfigManager:
    """Get the global configuration manager instance."""
    global _config_manager
    if _config_manager is None:
        _config_manager = ConfigManager()
    return _config_manager

def load_config() -> Config:
    """
    Load the global configuration.
    
    Returns:
        Current configuration
    """
    global _current_config
    if _current_config is None:
        manager = get_config_manager()
        _current_config = manager.load_config()
    return _current_config

def save_config(config: Optional[Config] = None) -> bool:
    """
    Save the global configuration.
    
    Args:
        config: Configuration to save (uses current if None)
        
    Returns:
        True if successful, False otherwise
    """
    global _current_config
    
    if config is None:
        config = _current_config
    
    if config is None:
        logger.warning("No configuration to save")
        return False
    
    manager = get_config_manager()
    success = manager.save_config(config)
    
    if success:
        _current_config = config
    
    return success

def update_config(**kwargs) -> bool:
    """
    Update specific configuration values.
    
    Args:
        **kwargs: Configuration values to update
        
    Returns:
        True if successful, False otherwise
    """
    config = load_config()
    
    # Update nested configuration sections
    for key, value in kwargs.items():
        if hasattr(config, key):
            if isinstance(value, dict):
                # Update nested object
                section = getattr(config, key)
                for nested_key, nested_value in value.items():
                    if hasattr(section, nested_key):
                        setattr(section, nested_key, nested_value)
            else:
                setattr(config, key, value)
    
    return save_config(config)

def get_config_value(key_path: str, default: Any = None) -> Any:
    """
    Get a configuration value using dot notation.
    
    Args:
        key_path: Dot-separated path to the value (e.g., 'gui.theme')
        default: Default value if key not found
        
    Returns:
        Configuration value or default
    """
    config = load_config()
    
    try:
        value = config
        for key in key_path.split('.'):
            value = getattr(value, key)
        return value
    except AttributeError:
        return default

def set_config_value(key_path: str, value: Any) -> bool:
    """
    Set a configuration value using dot notation.
    
    Args:
        key_path: Dot-separated path to the value (e.g., 'gui.theme')
        value: Value to set
        
    Returns:
        True if successful, False otherwise
    """
    config = load_config()
    
    try:
        keys = key_path.split('.')
        target = config
        
        # Navigate to the parent object
        for key in keys[:-1]:
            target = getattr(target, key)
        
        # Set the final value
        setattr(target, keys[-1], value)
        
        return save_config(config)
        
    except AttributeError as e:
        logger.error(f"Error setting config value {key_path}: {e}")
        return False

# Convenience functions for common configuration operations
def get_working_directory() -> str:
    """Get the current working directory from config."""
    return get_config_value('paths.working_directory', str(Path.cwd()))

def set_working_directory(directory: str) -> bool:
    """Set the working directory in config."""
    return set_config_value('paths.working_directory', directory)

def get_gui_theme() -> str:
    """Get the current GUI theme."""
    return get_config_value('gui.theme', 'dark')

def set_gui_theme(theme: str) -> bool:
    """Set the GUI theme."""
    return set_config_value('gui.theme', theme)

def get_default_force_fields() -> Dict[str, str]:
    """Get default force field settings."""
    config = load_config()
    return {
        'water': config.system_builder.water_model,
        'protein': config.system_builder.protein_ff,
        'lipid': config.system_builder.lipid_ff
    }