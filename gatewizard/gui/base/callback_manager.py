"""
Base callback manager for GUI components.

This module provides a base class that handles proper callback scheduling
and cleanup for all GUI components to prevent "invalid command name" errors
during application shutdown.
"""

import logging
from typing import List, Optional, Callable, Any, Protocol

logger = logging.getLogger(__name__)

class TkinterWidgetProtocol(Protocol):
    """Protocol defining the Tkinter widget interface needed by CallbackManager."""
    def after(self, ms: int, func: Callable, *args: Any) -> str: ...
    def after_cancel(self, id: str) -> None: ...

class CallbackManager:
    """
    Base class for managing Tkinter callbacks with proper cleanup.
    
    All GUI components should inherit from this class or use it as a mixin
    to ensure proper callback tracking and cleanup.
    """
    
    def __init__(self):
        """Initialize the callback manager."""
        self._callback_ids: List[str] = []
        self._is_destroyed = False
    
    def schedule_callback(self, delay: int, callback: Callable, *args: Any) -> Optional[str]:
        """
        Schedule a callback with proper tracking for cleanup.
        
        Args:
            delay: Delay in milliseconds
            callback: Function to call
            *args: Arguments to pass to the callback
            
        Returns:
            Callback ID or None if scheduling failed
        """
        if self._is_destroyed:
            return None
            
        try:
            # This method requires the class to be mixed with a Tkinter widget
            # that has the after() method
            if not hasattr(self, 'after'):
                raise AttributeError("CallbackManager must be used with a Tkinter widget class")
            
            callback_id = getattr(self, 'after')(delay, callback, *args)
            self._callback_ids.append(callback_id)
            return callback_id
        except Exception as e:
            logger.debug(f"Failed to schedule callback: {e}")
            return None
    
    def cancel_callback(self, callback_id: str) -> None:
        """
        Cancel a specific callback.
        
        Args:
            callback_id: ID of the callback to cancel
        """
        if callback_id and callback_id in self._callback_ids:
            try:
                if not hasattr(self, 'after_cancel'):
                    raise AttributeError("CallbackManager must be used with a Tkinter widget class")
                
                getattr(self, 'after_cancel')(callback_id)
                self._callback_ids.remove(callback_id)
            except Exception as e:
                logger.debug(f"Failed to cancel callback {callback_id}: {e}")
    
    def cleanup_callbacks(self) -> None:
        """Cancel all tracked callbacks."""
        self._is_destroyed = True
        
        for callback_id in self._callback_ids[:]:  # Copy list to avoid modification during iteration
            try:
                if hasattr(self, 'after_cancel'):
                    getattr(self, 'after_cancel')(callback_id)
            except Exception as e:
                logger.debug(f"Failed to cancel callback {callback_id}: {e}")
        
        self._callback_ids.clear()
        logger.debug(f"Cleaned up {len(self._callback_ids)} callbacks for {type(self).__name__}")
    
    def __del__(self):
        """Ensure cleanup on deletion."""
        if not self._is_destroyed:
            self.cleanup_callbacks()