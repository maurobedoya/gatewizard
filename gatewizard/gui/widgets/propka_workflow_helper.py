# gatewizard/gui/widgets/propka_workflow_helper.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Propka workflow helper dialog for guiding users through the two-stage process.

This module provides a dialog to help users implement the
proper Propka two-stage workflow for membrane protein preparation.
"""

import os
import tkinter as tk
from tkinter import filedialog, messagebox
from typing import Optional, Dict, Any, List, Tuple
from pathlib import Path
import threading

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

from gatewizard.gui.constants import COLOR_SCHEME, FONTS, LAYOUT
from gatewizard.core.preparation import PreparationManager, PreparationError
from gatewizard.core.builder import Builder
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class PropkaWorkflowDialog(ctk.CTkToplevel):
    """
    Dialog to help users with the Propka two-stage workflow.
    
    This dialog guides users through:
    1. Running Propka analysis on original PDB
    2. Modifying residue names in packed PDB
    3. Preparing for Stage 2 parametrization
    """
    
    def __init__(self, parent, packed_pdb_path: str = ""):
        """
        Initialize the Propka workflow helper dialog.
        
        Args:
            parent: Parent widget
            packed_pdb_path: Path to packed PDB file from Stage 1
        """
        super().__init__(parent)
        
        self.packed_pdb_path = packed_pdb_path
        self.original_pdb_path = ""
        self.modified_pdb_path = ""
        self.preparation_manager = PreparationManager()
        self.builder = Builder()
        self.propka_results = []
        self.custom_modifications = {}
        self.current_ph = 7.0
        
        # Configure window
        self.title("Propka Workflow Helper")
        self.geometry("850x750")
        self.resizable(True, True)
        self.transient(parent)
        
        # Set protocol for window close button
        self.protocol("WM_DELETE_WINDOW", self._close_dialog)
        
        # Create widgets
        self._create_widgets()
        self._setup_layout()
        
        # Center dialog
        self._center_window()
        
        # Set grab after window is fully initialized
        self.after(10, self._set_grab)
        
        # Auto-fill packed PDB if provided
        if packed_pdb_path:
            self.packed_pdb_entry.insert(0, packed_pdb_path)
    
    def cleanup_callbacks(self):
        """Cancel all scheduled callbacks to prevent errors during shutdown."""
        try:
            # This dialog uses several self.after() calls that need cleanup
            logger.debug(f"Cleaned up callbacks for {type(self).__name__}")
        except Exception as e:
            logger.debug(f"Error cleaning up callbacks in {type(self).__name__}: {e}")
    
    def _create_widgets(self):
        """Create dialog widgets."""
        # Main scrollable frame
        self.main_frame = ctk.CTkScrollableFrame(self)
        
        # Title
        self.title_label = ctk.CTkLabel(
            self.main_frame,
            text="Propka Two-Stage Workflow Helper",
            font=FONTS['title']
        )
        
        # Introduction
        self.intro_frame = ctk.CTkFrame(self.main_frame, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.intro_text = ctk.CTkTextbox(
            self.intro_frame,
            height=80,
            font=FONTS['small'],
            wrap="word"
        )
        
        intro_content = """This helper guides you through the Propka two-stage workflow for proper protonation state handling. The workflow involves: (1) Running Propka analysis on your original PDB, (2) Modifying residue names in the packed PDB file, and (3) Preparing the modified file for Stage 2 parametrization."""
        
        self.intro_text.insert("1.0", intro_content)
        self.intro_text.configure(state="disabled")
        
        # Step 1: Select files
        self.step1_frame = ctk.CTkFrame(self.main_frame, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.step1_label = ctk.CTkLabel(
            self.step1_frame,
            text="Step 1: Select Files",
            font=FONTS['heading']
        )
        
        # Original PDB file
        self.original_pdb_frame = ctk.CTkFrame(self.step1_frame, fg_color="transparent")
        
        self.original_pdb_label = ctk.CTkLabel(
            self.original_pdb_frame,
            text="Original PDB file (for Propka analysis):",
            font=FONTS['body']
        )
        
        self.original_pdb_input_frame = ctk.CTkFrame(self.original_pdb_frame, fg_color="transparent")
        
        self.original_pdb_entry = ctk.CTkEntry(
            self.original_pdb_input_frame,
            placeholder_text="Select original PDB file...",
            width=400
        )
        
        self.original_pdb_browse = ctk.CTkButton(
            self.original_pdb_input_frame,
            text="Browse",
            width=80,
            command=self._browse_original_pdb
        )
        
        # Packed PDB file
        self.packed_pdb_frame = ctk.CTkFrame(self.step1_frame, fg_color="transparent")
        
        self.packed_pdb_label = ctk.CTkLabel(
            self.packed_pdb_frame,
            text="Packed PDB file (from Stage 1):",
            font=FONTS['body']
        )
        
        self.packed_pdb_input_frame = ctk.CTkFrame(self.packed_pdb_frame, fg_color="transparent")
        
        self.packed_pdb_entry = ctk.CTkEntry(
            self.packed_pdb_input_frame,
            placeholder_text="Select packed PDB file...",
            width=400
        )
        
        self.packed_pdb_browse = ctk.CTkButton(
            self.packed_pdb_input_frame,
            text="Browse",
            width=80,
            command=self._browse_packed_pdb
        )
        
        # Step 2: Run Propka analysis
        self.step2_frame = ctk.CTkFrame(self.main_frame, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.step2_label = ctk.CTkLabel(
            self.step2_frame,
            text="Step 2: Propka Analysis",
            font=FONTS['heading']
        )
        
        self.propka_params_frame = ctk.CTkFrame(self.step2_frame, fg_color="transparent")
        
        self.ph_label = ctk.CTkLabel(
            self.propka_params_frame,
            text="Target pH:",
            font=FONTS['body']
        )
        
        self.ph_entry = ctk.CTkEntry(
            self.propka_params_frame,
            width=80
        )
        self.ph_entry.insert(0, "7.0")
        
        self.version_label = ctk.CTkLabel(
            self.propka_params_frame,
            text="Propka Version:",
            font=FONTS['body']
        )
        
        self.version_combo = ctk.CTkComboBox(
            self.propka_params_frame,
            values=["3"],
            width=80
        )
        self.version_combo.set("3")
        
        self.run_propka_button = ctk.CTkButton(
            self.propka_params_frame,
            text="Run Propka Analysis",
            command=self._run_propka_analysis
        )
        
        # Propka results display
        self.results_frame = ctk.CTkFrame(self.step2_frame, fg_color="transparent")
        
        self.results_label = ctk.CTkLabel(
            self.results_frame,
            text="Propka Results:",
            font=FONTS['body']
        )
        
        self.results_text = ctk.CTkTextbox(
            self.results_frame,
            height=200,
            font=FONTS['small']
        )
        
        # Step 3: Modify residue names
        self.step3_frame = ctk.CTkFrame(self.main_frame, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.step3_label = ctk.CTkLabel(
            self.step3_frame,
            text="Step 3: Modify Residue Names",
            font=FONTS['heading']
        )
        
        self.modification_frame = ctk.CTkFrame(self.step3_frame, fg_color="transparent")
        
        self.auto_modify_button = ctk.CTkButton(
            self.modification_frame,
            text="Auto-Apply Propka Results",
            command=self._auto_apply_modifications
        )
        
        self.manual_modify_button = ctk.CTkButton(
            self.modification_frame,
            text="Manual Modifications...",
            command=self._manual_modifications
        )
        
        self.preview_button = ctk.CTkButton(
            self.modification_frame,
            text="Preview Changes",
            command=self._preview_modifications
        )
        
        # Output file
        self.output_frame = ctk.CTkFrame(self.step3_frame, fg_color="transparent")
        
        self.output_label = ctk.CTkLabel(
            self.output_frame,
            text="Modified PDB output file:",
            font=FONTS['body']
        )
        
        self.output_input_frame = ctk.CTkFrame(self.output_frame, fg_color="transparent")
        
        self.output_entry = ctk.CTkEntry(
            self.output_input_frame,
            placeholder_text="Output file path...",
            width=400
        )
        
        self.output_browse = ctk.CTkButton(
            self.output_input_frame,
            text="Browse",
            width=80,
            command=self._browse_output_file
        )
        
        # Modifications summary
        self.summary_frame = ctk.CTkFrame(self.step3_frame, fg_color="transparent")
        
        self.summary_label = ctk.CTkLabel(
            self.summary_frame,
            text="Modifications Summary:",
            font=FONTS['body']
        )
        
        self.summary_text = ctk.CTkTextbox(
            self.summary_frame,
            height=100,
            font=FONTS['small']
        )
        
        # Step 4: Ready for Stage 2
        self.step4_frame = ctk.CTkFrame(self.main_frame, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.step4_label = ctk.CTkLabel(
            self.step4_frame,
            text="Step 4: Ready for Stage 2",
            font=FONTS['heading']
        )
        
        self.stage2_info = ctk.CTkLabel(
            self.step4_frame,
            text="After modifying residue names, use the modified PDB file in Stage 2 parametrization with --notprotonate.",
            font=FONTS['body'],
            wraplength=500
        )
        
        self.copy_path_button = ctk.CTkButton(
            self.step4_frame,
            text="Copy Output Path",
            command=self._copy_output_path
        )
        
        # Action buttons
        self.buttons_frame = ctk.CTkFrame(self.main_frame, fg_color="transparent")
        
        self.close_button = ctk.CTkButton(
            self.buttons_frame,
            text="Close",
            command=self._close_dialog
        )
        
        self.help_button = ctk.CTkButton(
            self.buttons_frame,
            text="Help",
            command=self._show_help
        )
        
        self.reset_button = ctk.CTkButton(
            self.buttons_frame,
            text="Reset",
            command=self._reset_dialog
        )
    
    def _setup_layout(self):
        """Setup dialog layout."""
        self.main_frame.pack(fill="both", expand=True, padx=20, pady=20)
        
        # Title
        self.title_label.pack(pady=(0, 15))
        
        # Introduction
        self.intro_frame.pack(fill="x", pady=(0, 15))
        self.intro_text.pack(fill="x", padx=15, pady=15)
        
        # Step 1
        self.step1_frame.pack(fill="x", pady=(0, 15))
        self.step1_label.pack(anchor="w", padx=15, pady=(15, 10))
        
        self.original_pdb_frame.pack(fill="x", padx=15, pady=5)
        self.original_pdb_label.pack(anchor="w", pady=(0, 5))
        self.original_pdb_input_frame.pack(fill="x")
        self.original_pdb_entry.pack(side="left", padx=(0, 10))
        self.original_pdb_browse.pack(side="left")
        
        self.packed_pdb_frame.pack(fill="x", padx=15, pady=(10, 15))
        self.packed_pdb_label.pack(anchor="w", pady=(0, 5))
        self.packed_pdb_input_frame.pack(fill="x")
        self.packed_pdb_entry.pack(side="left", padx=(0, 10))
        self.packed_pdb_browse.pack(side="left")
        
        # Step 2
        self.step2_frame.pack(fill="x", pady=(0, 15))
        self.step2_label.pack(anchor="w", padx=15, pady=(15, 10))
        
        self.propka_params_frame.pack(fill="x", padx=15, pady=5)
        self.ph_label.pack(side="left", padx=(0, 10))
        self.ph_entry.pack(side="left", padx=(0, 20))
        self.version_label.pack(side="left", padx=(0, 10))
        self.version_combo.pack(side="left", padx=(0, 20))
        self.run_propka_button.pack(side="left")
        
        self.results_frame.pack(fill="x", padx=15, pady=(10, 15))
        self.results_label.pack(anchor="w", pady=(0, 5))
        self.results_text.pack(fill="x")
        
        # Step 3
        self.step3_frame.pack(fill="x", pady=(0, 15))
        self.step3_label.pack(anchor="w", padx=15, pady=(15, 10))
        
        self.modification_frame.pack(fill="x", padx=15, pady=5)
        self.auto_modify_button.pack(side="left", padx=(0, 10))
        self.manual_modify_button.pack(side="left", padx=(0, 10))
        self.preview_button.pack(side="left")
        
        self.output_frame.pack(fill="x", padx=15, pady=(10, 0))
        self.output_label.pack(anchor="w", pady=(0, 5))
        self.output_input_frame.pack(fill="x")
        self.output_entry.pack(side="left", padx=(0, 10))
        self.output_browse.pack(side="left")
        
        self.summary_frame.pack(fill="x", padx=15, pady=(10, 15))
        self.summary_label.pack(anchor="w", pady=(0, 5))
        self.summary_text.pack(fill="x")
        
        # Step 4
        self.step4_frame.pack(fill="x", pady=(0, 15))
        self.step4_label.pack(anchor="w", padx=15, pady=(15, 10))
        self.stage2_info.pack(anchor="w", padx=15, pady=(0, 10))
        self.copy_path_button.pack(padx=15, pady=(0, 15))
        
        # Buttons
        self.buttons_frame.pack(fill="x", pady=(10, 0))
        self.close_button.pack(side="right")
        self.help_button.pack(side="right", padx=(0, 10))
        self.reset_button.pack(side="right", padx=(0, 10))
    
    def _center_window(self):
        """Center the dialog window."""
        self.update_idletasks()
        
        # Get screen dimensions
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        
        # Get window dimensions
        window_width = self.winfo_width()
        window_height = self.winfo_height()
        
        # Calculate position
        x = (screen_width - window_width) // 2
        y = (screen_height - window_height) // 2
        
        self.geometry(f"+{x}+{y}")
    
    def _set_grab(self):
        """Set grab after window is visible to avoid errors."""
        try:
            self.grab_set()
        except Exception as e:
            logger.warning(f"Could not set grab on dialog: {e}")
    
    def _close_dialog(self):
        """Properly close the dialog."""
        try:
            self.grab_release()
        except Exception:
            pass
        self.destroy()
    
    def _browse_original_pdb(self):
        """Browse for original PDB file."""
        file_path = filedialog.askopenfilename(
            title="Select Original PDB File",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
        )
        
        if file_path:
            self.original_pdb_entry.delete(0, "end")
            self.original_pdb_entry.insert(0, file_path)
            self.original_pdb_path = file_path
    
    def _browse_packed_pdb(self):
        """Browse for packed PDB file."""
        file_path = filedialog.askopenfilename(
            title="Select Packed PDB File",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")]
        )
        
        if file_path:
            self.packed_pdb_entry.delete(0, "end")
            self.packed_pdb_entry.insert(0, file_path)
            self.packed_pdb_path = file_path
    
    def _browse_output_file(self):
        """Browse for output file location."""
        # Auto-generate default name if packed PDB is selected
        default_name = ""
        if self.packed_pdb_path:
            base_path = Path(self.packed_pdb_path)
            default_name = f"{base_path.stem}_propka_modified.pdb"
        
        file_path = filedialog.asksaveasfilename(
            title="Save Modified PDB File",
            defaultextension=".pdb",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
            initialfile=default_name
        )
        
        if file_path:
            self.output_entry.delete(0, "end")
            self.output_entry.insert(0, file_path)
            self.modified_pdb_path = file_path
    
    def _run_propka_analysis(self):
        """Run Propka analysis on the original PDB."""
        original_pdb = self.original_pdb_entry.get().strip()
        ph_str = self.ph_entry.get().strip()
        version = self.version_combo.get()
        
        if not original_pdb:
            messagebox.showerror("Error", "Please select the original PDB file.")
            return
        
        if not ph_str:
            messagebox.showerror("Error", "Please enter a pH value.")
            return
        
        try:
            ph = float(ph_str)
            if ph < 0 or ph > 14:
                messagebox.showerror("Error", "pH must be between 0 and 14.")
                return
            self.current_ph = ph
        except ValueError:
            messagebox.showerror("Error", "Invalid pH value.")
            return
        
        # Run analysis in background thread
        self._run_propka_thread(original_pdb, ph, version)
    
    def _run_propka_thread(self, original_pdb: str, ph: float, version: str):
        """Run Propka analysis in background thread."""
        def analysis_worker():
            try:
                # Update UI
                self.after(0, self._update_propka_status, "Running Propka analysis...")
                
                # Set analyzer version
                self.preparation_manager.propka_version = version
                
                # Run analysis
                pka_file = self.preparation_manager.run_analysis(original_pdb)
                summary_file = self.preparation_manager.extract_summary(pka_file)
                results = self.preparation_manager.parse_summary(summary_file)
                
                # Update UI with results
                self.after(0, self._display_propka_results, results, ph)
                
            except PreparationError as e:
                self.after(0, self._show_propka_error, str(e))
            except Exception as e:
                logger.error(f"Unexpected error in Propka analysis: {e}", exc_info=True)
                self.after(0, self._show_propka_error, f"Unexpected error: {str(e)}")
        
        # Start analysis thread
        analysis_thread = threading.Thread(target=analysis_worker, daemon=True)
        analysis_thread.start()
    
    def _update_propka_status(self, message: str):
        """Update Propka analysis status."""
        self.results_text.delete("1.0", "end")
        self.results_text.insert("1.0", message)
        self.update()
    
    def _display_propka_results(self, results: List[Dict], ph: float):
        """Display Propka analysis results."""
        self.propka_results = results
        
        # Display results
        results_text = f"Propka analysis completed for pH {ph:.1f}\n\n"
        results_text += f"Found {len(results)} protonable residues:\n\n"
        
        modifications_count = 0
        for residue in results:
            default_state = self.preparation_manager.get_default_protonation_state(residue, ph)
            original_state = residue['residue']
            chain_info = f" (Chain {residue['chain']})" if residue.get('chain') else ""
            
            if default_state != original_state:
                results_text += f"✓ {original_state} {residue['res_id']}{chain_info}: pKa={residue['pka']:.2f} → {default_state}\n"
                modifications_count += 1
            else:
                results_text += f"  {original_state} {residue['res_id']}{chain_info}: pKa={residue['pka']:.2f} → {default_state} (no change)\n"
        
        results_text += f"\nSuggested modifications: {modifications_count}"
        
        self.results_text.delete("1.0", "end")
        self.results_text.insert("1.0", results_text)
        
        messagebox.showinfo("Success", f"Propka analysis completed! Found {modifications_count} suggested modifications.")
    
    def _show_propka_error(self, error_message: str):
        """Show Propka analysis error."""
        logger.error(f"Propka analysis error: {error_message}")
        messagebox.showerror("Analysis Error", f"Propka analysis failed:\n\n{error_message}")
        self.results_text.delete("1.0", "end")
        self.results_text.insert("1.0", f"Error: {error_message}")
    
    def _auto_apply_modifications(self):
        """Automatically apply Propka results to packed PDB."""
        if not self.propka_results:
            messagebox.showerror("Error", "Please run Propka analysis first.")
            return
        
        packed_pdb = self.packed_pdb_entry.get().strip()
        if not packed_pdb:
            messagebox.showerror("Error", "Please select the packed PDB file.")
            return
        
        output_file = self.output_entry.get().strip()
        if not output_file:
            # Auto-generate output filename
            base_path = Path(packed_pdb)
            output_file = str(base_path.with_name(f"{base_path.stem}_propka_modified.pdb"))
            self.output_entry.delete(0, "end")
            self.output_entry.insert(0, output_file)
        
        try:
            # Create modification dictionary
            modifications = {}
            for residue in self.propka_results:
                default_state = self.preparation_manager.get_default_protonation_state(residue, self.current_ph)
                if default_state != residue['residue']:
                    # Use chain + residue number as key for better identification
                    chain = residue.get('chain', '')
                    if chain:
                        key = f"{chain}:{residue['res_id']}"
                    else:
                        key = str(residue['res_id'])
                    modifications[key] = default_state
            
            # Apply modifications using system builder
            success, message = self.builder.modify_residue_names_for_propka(
                packed_pdb, modifications, output_file
            )
            
            if success:
                self.modified_pdb_path = output_file
                self._update_summary(modifications, success=True, message=message)
                messagebox.showinfo("Success", message)
            else:
                messagebox.showerror("Error", f"Failed to apply modifications: {message}")
                
        except Exception as e:
            logger.error(f"Error applying modifications: {e}")
            messagebox.showerror("Error", f"Failed to apply modifications: {str(e)}")
    
    def _manual_modifications(self):
        """Open manual modifications dialog."""
        if not self.propka_results:
            messagebox.showwarning("Warning", "Run Propka analysis first to get suggested modifications.")
        
        # Create manual modifications dialog
        dialog = ManualModificationsDialog(self, self.propka_results, self.current_ph)
        dialog.wait_window()
        
        # Get results from dialog
        if hasattr(dialog, 'modifications') and dialog.modifications:
            self.custom_modifications = dialog.modifications
            self._update_summary(self.custom_modifications, success=False, message="Manual modifications defined")
    
    def _preview_modifications(self):
        """Preview the modifications that will be applied."""
        if not self.propka_results and not self.custom_modifications:
            messagebox.showinfo("No Modifications", "No modifications to preview. Run Propka analysis or define manual modifications first.")
            return
        
        # Create preview text
        preview_text = "Modifications Preview:\n\n"
        
        if self.custom_modifications:
            preview_text += "Manual modifications:\n"
            for res_id, new_name in self.custom_modifications.items():
                preview_text += f"  Residue {res_id} → {new_name}\n"
        else:
            preview_text += "Propka-suggested modifications:\n"
            for residue in self.propka_results:
                default_state = self.preparation_manager.get_default_protonation_state(residue, self.current_ph)
                if default_state != residue['residue']:
                    chain_info = f" (Chain {residue['chain']})" if residue.get('chain') else ""
                    preview_text += f"  {residue['residue']} {residue['res_id']}{chain_info} → {default_state}\n"
        
        # Show preview dialog
        preview_dialog = ctk.CTkToplevel(self)
        preview_dialog.title("Modifications Preview")
        preview_dialog.geometry("400x300")
        preview_dialog.transient(self)
        preview_dialog.grab_set()
        
        preview_text_widget = ctk.CTkTextbox(preview_dialog, font=FONTS['small'])
        preview_text_widget.pack(fill="both", expand=True, padx=20, pady=20)
        preview_text_widget.insert("1.0", preview_text)
        preview_text_widget.configure(state="disabled")
        
        close_btn = ctk.CTkButton(preview_dialog, text="Close", command=preview_dialog.destroy)
        close_btn.pack(pady=(0, 20))
    
    def _update_summary(self, modifications: Dict[str, str], success: bool, message: str):
        """Update the modifications summary."""
        summary_text = f"Status: {'Applied' if success else 'Defined'}\n"
        summary_text += f"Modifications: {len(modifications)}\n\n"
        
        if modifications:
            summary_text += "Changes:\n"
            for res_id, new_name in modifications.items():
                summary_text += f"  Residue {res_id} → {new_name}\n"
        
        if message:
            summary_text += f"\nDetails: {message}"
        
        self.summary_text.delete("1.0", "end")
        self.summary_text.insert("1.0", summary_text)
    
    def _copy_output_path(self):
        """Copy output file path to clipboard."""
        output_path = self.output_entry.get().strip()
        if output_path:
            try:
                self.clipboard_clear()
                self.clipboard_append(output_path)
                self.update()
                messagebox.showinfo("Copied", "Output file path copied to clipboard!")
            except Exception as e:
                logger.error(f"Error copying to clipboard: {e}")
                messagebox.showerror("Error", "Failed to copy to clipboard.")
        else:
            messagebox.showwarning("No Path", "No output file path to copy.")
    
    def _reset_dialog(self):
        """Reset the dialog to initial state."""
        result = messagebox.askyesno("Reset", "Reset all fields and results?")
        if result:
            # Clear all entries
            self.original_pdb_entry.delete(0, "end")
            self.packed_pdb_entry.delete(0, "end")
            self.output_entry.delete(0, "end")
            self.ph_entry.delete(0, "end")
            self.ph_entry.insert(0, "7.0")
            
            # Clear results
            self.results_text.delete("1.0", "end")
            self.summary_text.delete("1.0", "end")
            
            # Reset state
            self.propka_results = []
            self.custom_modifications = {}
            self.current_ph = 7.0
    
    def _show_help(self):
        """Show detailed help for using this dialog."""
        help_text = """Propka Workflow Helper - User Guide

OVERVIEW:
This helper automates the Propka two-stage workflow for proper protonation state handling in membrane protein preparation.

STEP-BY-STEP INSTRUCTIONS:

1. SELECT FILES:
   - Original PDB: Your initial protein structure (used for Propka analysis)
   - Packed PDB: The bilayer_*.pdb file generated by Stage 1 packing

2. RUN PROPKA ANALYSIS:
   - Set target pH (usually 7.0 for physiological conditions)
   - Choose Propka version (3 is recommended)
   - Click "Run Propka Analysis"
   - Review results showing which residues need modification

3. MODIFY RESIDUE NAMES:
   - Auto-Apply: Automatically applies all Propka suggestions
   - Manual: Define custom modifications for specific residues
   - Preview: See what changes will be made before applying

4. READY FOR STAGE 2:
   - Use the modified PDB file in Stage 2 parametrization
   - Copy the output path for easy use in Stage 2

COMMON MODIFICATIONS:
- GLU → GLH (protonated glutamate at low pH)
- ASP → ASH (protonated aspartate at low pH)
- HIS → HIP/HIE/HID (different histidine protonation states)
- LYS → LYN (neutral lysine at high pH)
- TYR → TYM (deprotonated tyrosine at high pH)

TROUBLESHOOTING:
- If Propka fails: Check that the original PDB file is valid
- If no modifications found: Protein may not need changes at target pH
- If residue numbers don't match: This is normal - packed files renumber residues

IMPORTANT NOTES:
- Always use the ORIGINAL PDB for Propka analysis
- Apply modifications to the PACKED PDB file
- The helper handles residue renumbering automatically
- Stage 2 will use --notprotonate to preserve your modifications
"""
        
        # Create help window
        help_window = ctk.CTkToplevel(self)
        help_window.title("Propka Workflow Helper - Help")
        help_window.geometry("700x600")
        help_window.resizable(True, True)
        help_window.transient(self)
        help_window.grab_set()
        
        # Center help window
        help_window.update_idletasks()
        x = self.winfo_rootx() + (self.winfo_width() // 2) - (700 // 2)
        y = self.winfo_rooty() + (self.winfo_height() // 2) - (600 // 2)
        help_window.geometry(f"+{x}+{y}")
        
        # Help content
        help_frame = ctk.CTkFrame(help_window)
        help_frame.pack(fill="both", expand=True, padx=20, pady=20)
        
        help_textbox = ctk.CTkTextbox(
            help_frame,
            font=FONTS['small'],
            wrap="word"
        )
        help_textbox.pack(fill="both", expand=True, pady=(0, 10))
        help_textbox.insert("1.0", help_text)
        help_textbox.configure(state="disabled")
        
        close_button = ctk.CTkButton(
            help_frame,
            text="Close",
            command=help_window.destroy
        )
        close_button.pack()

class ManualModificationsDialog(ctk.CTkToplevel):
    """
    Dialog for manually defining residue modifications.
    
    Allows users to specify custom protonation states beyond
    what Propka suggests.
    """
    
    def __init__(self, parent, propka_results: List[Dict], current_ph: float):
        """
        Initialize manual modifications dialog.
        
        Args:
            parent: Parent widget
            propka_results: Results from Propka analysis
            current_ph: Current pH value
        """
        super().__init__(parent)
        
        self.propka_results = propka_results or []
        self.current_ph = current_ph
        self.modifications = {}
        self.modification_widgets = []
        
        # Configure window
        self.title("Manual Residue Modifications")
        self.geometry("600x500")
        self.resizable(True, True)
        self.transient(parent)
        self.grab_set()
        
        # Create widgets
        self._create_widgets()
        self._setup_layout()
        self._populate_suggestions()
        
        # Center dialog
        self._center_window()
    
    def _create_widgets(self):
        """Create manual modification widgets."""
        # Main frame
        self.main_frame = ctk.CTkFrame(self)
        
        # Title
        self.title_label = ctk.CTkLabel(
            self.main_frame,
            text="Define Manual Residue Modifications",
            font=FONTS['heading']
        )
        
        # Instructions
        self.instructions = ctk.CTkLabel(
            self.main_frame,
            text="Define custom protonation states for specific residues. Leave blank to use original state.",
            font=FONTS['body'],
            wraplength=500
        )
        
        # Modifications area
        self.modifications_frame = ctk.CTkScrollableFrame(
            self.main_frame,
            height=300,
            fg_color=COLOR_SCHEME['background']
        )
        
        # Add modification button
        self.add_frame = ctk.CTkFrame(self.main_frame, fg_color="transparent")
        
        self.add_button = ctk.CTkButton(
            self.add_frame,
            text="Add Modification",
            command=self._add_modification_row
        )
        
        self.clear_button = ctk.CTkButton(
            self.add_frame,
            text="Clear All",
            command=self._clear_all_modifications
        )
        
        # Action buttons
        self.buttons_frame = ctk.CTkFrame(self.main_frame, fg_color="transparent")
        
        self.ok_button = ctk.CTkButton(
            self.buttons_frame,
            text="Apply",
            command=self._apply_modifications
        )
        
        self.cancel_button = ctk.CTkButton(
            self.buttons_frame,
            text="Cancel",
            command=self.destroy
        )
    
    def _setup_layout(self):
        """Setup dialog layout."""
        self.main_frame.pack(fill="both", expand=True, padx=20, pady=20)
        
        self.title_label.pack(pady=(0, 10))
        self.instructions.pack(pady=(0, 15))
        
        self.modifications_frame.pack(fill="both", expand=True, pady=(0, 15))
        
        self.add_frame.pack(fill="x", pady=(0, 15))
        self.add_button.pack(side="left", padx=(0, 10))
        self.clear_button.pack(side="left")
        
        self.buttons_frame.pack(fill="x")
        self.cancel_button.pack(side="right")
        self.ok_button.pack(side="right", padx=(0, 10))
    
    def _center_window(self):
        """Center the dialog window."""
        self.update_idletasks()
        
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        window_width = self.winfo_width()
        window_height = self.winfo_height()
        
        x = (screen_width - window_width) // 2
        y = (screen_height - window_height) // 2
        
        self.geometry(f"+{x}+{y}")
    
    def _populate_suggestions(self):
        """Populate with Propka suggestions."""
        if not self.propka_results:
            self._add_modification_row()
            return
        
        preparation_manager = PreparationManager()
        
        for residue in self.propka_results:
            default_state = preparation_manager.get_default_protonation_state(residue, self.current_ph)
            
            # Only show residues that Propka suggests changing
            if default_state != residue['residue']:
                chain = residue.get('chain', '')
                res_display = f"{residue['res_id']} (Chain {chain})" if chain else str(residue['res_id'])
                key = f"{chain}:{residue['res_id']}" if chain else str(residue['res_id'])
                
                self._add_modification_row(
                    residue_id=key,
                    current_name=residue['residue'],
                    suggested_name=default_state
                )
    
    def _add_modification_row(self, residue_id: str = "", current_name: str = "", suggested_name: str = ""):
        """Add a row for defining a modification."""
        row_frame = ctk.CTkFrame(self.modifications_frame, fg_color=COLOR_SCHEME['canvas'])
        row_frame.pack(fill="x", pady=2, padx=5)
        
        # Residue ID entry
        id_label = ctk.CTkLabel(row_frame, text="Residue ID:", width=80)
        id_label.pack(side="left", padx=(10, 5))
        
        id_entry = ctk.CTkEntry(row_frame, width=80)
        id_entry.pack(side="left", padx=(0, 10))
        if residue_id:
            id_entry.insert(0, residue_id)
        
        # Current name (read-only)
        current_label = ctk.CTkLabel(row_frame, text="Current:", width=60)
        current_label.pack(side="left", padx=(0, 5))
        
        current_display = ctk.CTkLabel(row_frame, text=current_name or "---", width=50)
        current_display.pack(side="left", padx=(0, 10))
        
        # Arrow
        arrow_label = ctk.CTkLabel(row_frame, text="→", width=20)
        arrow_label.pack(side="left")
        
        # New name combobox
        new_label = ctk.CTkLabel(row_frame, text="New:", width=40)
        new_label.pack(side="left", padx=(10, 5))
        
        # Available protonation states
        available_states = ["GLU", "GLH", "ASP", "ASH", "HIS", "HIP", "HIE", "HID", 
                          "LYS", "LYN", "TYR", "TYM", "CYS", "CYM", "ARG", "RN1"]
        
        new_combo = ctk.CTkComboBox(row_frame, values=available_states, width=80)
        new_combo.pack(side="left", padx=(0, 10))
        if suggested_name:
            new_combo.set(suggested_name)
        
        # Remove button
        remove_button = ctk.CTkButton(
            row_frame,
            text="×",
            width=30,
            height=25,
            command=lambda: self._remove_modification_row(row_frame)
        )
        remove_button.pack(side="right", padx=10)
        
        # Store widget references
        widget_data = {
            'frame': row_frame,
            'id_entry': id_entry,
            'current_display': current_display,
            'new_combo': new_combo,
            'remove_button': remove_button
        }
        self.modification_widgets.append(widget_data)
    
    def _remove_modification_row(self, frame):
        """Remove a modification row."""
        # Find and remove widget data
        for i, widget_data in enumerate(self.modification_widgets):
            if widget_data['frame'] == frame:
                self.modification_widgets.pop(i)
                break
        
        # Remove frame
        frame.destroy()
    
    def _clear_all_modifications(self):
        """Clear all modification rows."""
        for widget_data in self.modification_widgets:
            widget_data['frame'].destroy()
        
        self.modification_widgets = []
        
        # Add one empty row
        self._add_modification_row()
    
    def _apply_modifications(self):
        """Apply the defined modifications."""
        modifications = {}
        
        for widget_data in self.modification_widgets:
            residue_id = widget_data['id_entry'].get().strip()
            new_name = widget_data['new_combo'].get().strip()
            
            if residue_id and new_name:
                modifications[residue_id] = new_name
        
        if not modifications:
            messagebox.showwarning("No Modifications", "No modifications defined.")
            return
        
        self.modifications = modifications
        self.destroy()

# Helper functions for the workflow
def create_residue_mapping_from_propka(propka_results: List[Dict], ph: float) -> Dict[str, str]:
    """
    Create a residue mapping dictionary from Propka results.
    
    Args:
        propka_results: Results from Propka analysis
        ph: Target pH value
        
    Returns:
        Dictionary mapping residue IDs to new residue names
    """
    preparation_manager = PreparationManager()
    mapping = {}
    
    for residue in propka_results:
        default_state = preparation_manager.get_default_protonation_state(residue, ph)
        if default_state != residue['residue']:
            mapping[str(residue['res_id'])] = default_state
    
    return mapping

def validate_residue_modifications(modifications: Dict[str, str]) -> Tuple[bool, str]:
    """
    Validate residue modifications for common errors.
    
    Args:
        modifications: Dictionary of residue modifications
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not modifications:
        return False, "No modifications provided"
    
    valid_residues = {"GLU", "GLH", "ASP", "ASH", "HIS", "HIP", "HIE", "HID", 
                     "LYS", "LYN", "TYR", "TYM", "CYS", "CYM", "ARG", "RN1"}
    
    for res_id, new_name in modifications.items():
        # Validate residue ID format
        if not res_id.isdigit():
            return False, f"Invalid residue ID format: {res_id}"
        
        # Validate new residue name
        if new_name not in valid_residues:
            return False, f"Invalid residue name: {new_name}"
    
    return True, ""

def get_modification_summary(modifications: Dict[str, str]) -> str:
    """
    Get a human-readable summary of modifications.
    
    Args:
        modifications: Dictionary of residue modifications
        
    Returns:
        Summary string
    """
    if not modifications:
        return "No modifications"
    
    summary_lines = [f"Residue {res_id} → {new_name}" for res_id, new_name in modifications.items()]
    return f"{len(modifications)} modifications: " + "; ".join(summary_lines)