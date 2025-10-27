#!/usr/bin/env python3

"""
Main entry point for Gatewizard application.

This module provides the main function that serves as the entry point
for the Gatewizard GUI application.
"""

import sys
import argparse
from pathlib import Path

from gatewizard import GUI_AVAILABLE, __version__
from gatewizard.utils.logger import setup_logger

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Gatewizard - Membrane protein preparation tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  gatewizard                    # Launch GUI
  gatewizard --screen 1         # Launch GUI on screen 1
  gatewizard --version          # Show version
        """
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version=f"Gatewizard {__version__}"
    )
    
    parser.add_argument(
        "--screen",
        type=int,
        default=0,
        help="Target screen number for multi-monitor setups (default: 0)"
    )
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging"
    )
    
    parser.add_argument(
        "--log-file",
        type=Path,
        help="Path to log file (default: gatewizard.log in current directory)"
    )
    
    return parser.parse_args()

def check_dependencies():
    """Check if all required dependencies are available."""
    missing_deps = []
    
    # Check GUI dependencies
    if not GUI_AVAILABLE:
        try:
            import customtkinter
            import tkinter
        except ImportError as e:
            missing_deps.append(f"GUI dependencies: {e}")
    
    # Check NumPy
    try:
        import numpy
    except ImportError:
        missing_deps.append("numpy (required for numerical operations)")
    
    if missing_deps:
        print("Error: Missing required dependencies:")
        for dep in missing_deps:
            print(f"  - {dep}")
        print("\nPlease install missing dependencies and try again.")
        return False
    
    return True

def main():
    """Main entry point for the Gatewizard application."""
    args = parse_arguments()
    
    # Load configuration first
    from gatewizard.utils.config import ConfigManager
    config_manager = ConfigManager()
    config = config_manager.load_config()
    
    # Get the current working directory where the user opened the app
    # This should always be fresh and not saved in settings
    current_user_directory = str(Path.cwd())
    
    # Setup logging based on configuration
    # Command line arguments override configuration
    log_file = None
    if args.log_file:
        # Command line log file specified - always enable file logging
        log_file = Path(args.log_file)
    else:
        # Auto-configure log file in the current user directory
        # Use the actual current working directory, not saved settings
        config.logging.log_file_path = current_user_directory
        config.logging.enable_file_logging = True  # Enable file logging by default
        
        # Don't save this auto-configuration to avoid remembering previous directories
        # This ensures each launch creates/uses log in the directory where app was opened
    
    logger = setup_logger(
        debug=args.debug,
        log_file=log_file,
        config=config.logging.to_dict() if hasattr(config.logging, 'to_dict') else config.logging.__dict__
    )
    
    logger.info(f"Starting Gatewizard {__version__}")
    
    # Check dependencies
    if not check_dependencies():
        sys.exit(1)
    
    # Check if GUI is available
    if not GUI_AVAILABLE:
        logger.error("GUI dependencies not available")
        print("Error: GUI dependencies are not available.")
        print("Please install GUI dependencies with:")
        print("  pip install gatewizard")
        sys.exit(1)
    
    try:
        # Import and start GUI
        from gatewizard.gui.app import ProteinViewerApp
        
        logger.info(f"Launching GUI on screen {args.screen}")
        app = ProteinViewerApp(target_screen=args.screen)
        app.mainloop()
        
    except KeyboardInterrupt:
        logger.info("Application interrupted by user")
        print("\nApplication interrupted by user.")
        sys.exit(0)
        
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)
        print(f"Error: An unexpected error occurred: {e}")
        print("Check the log file for more details.")
        sys.exit(1)
    
    finally:
        logger.info("Gatewizard application ended")

if __name__ == "__main__":
    main()
