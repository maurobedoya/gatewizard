# gatewizard/gui/widgets/ligand_param_widget.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Widget for ligand parametrization in the Builder frame.

Detects ligands from a PDB file, displays their information and 2D structures,
and provides controls for setting charge/method before parametrization.
Ligand cards are arranged horizontally (up to 4 per row, wrapping for more).
"""

import os
import re
import threading
from pathlib import Path
from typing import Optional, Callable, Dict, List, Any

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

from tkinter import messagebox

from gatewizard.gui.constants import (
    COLOR_SCHEME, FONTS, WIDGET_SIZES, LAYOUT
)
from gatewizard.tools.ligand_parametrization import (
    detect_ligands,
    extract_ligand_pdb,
    parametrize_ligand,
    get_ligand_2d_image,
    get_ligand_2d_image_from_pdb_lines,
    LigandInfo,
    LigandParametrizationError,
    CHARGE_METHODS,
    DEFAULT_CHARGE_METHOD,
    ATOM_TYPES,
    DEFAULT_ATOM_TYPE,
    NON_RECOMMENDED_COMBOS,
)
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

from PIL import Image as PILImage
from rdkit import Chem

# Maximum number of ligand columns per row
MAX_COLS = 4

# Default 2D image dimensions (pixels) – overridden dynamically
IMG_WIDTH = 200
IMG_HEIGHT = 160

# Adaptive size table: {num_ligands: (width, height)}
_ADAPTIVE_SIZES = {
    1: (520, 420),
    2: (380, 300),
    3: (300, 240),
}
_DEFAULT_SIZE = (IMG_WIDTH, IMG_HEIGHT)  # 4+ ligands


def _compute_image_size(num_ligands: int) -> tuple:
    """Return (width, height) adapted to how many ligands are shown."""
    return _ADAPTIVE_SIZES.get(num_ligands, _DEFAULT_SIZE)


class LigandParamWidget(ctk.CTkFrame):
    """
    Widget for detecting, viewing, and parametrizing ligands from a PDB file.

    Placed in the Builder frame between Output Folder Name and Lipid Composition.
    Ligand cards are arranged horizontally (up to 4 per row).
    """

    def __init__(
        self,
        parent,
        status_callback: Optional[Callable[[str], None]] = None,
        working_dir_callback: Optional[Callable[[], str]] = None,
    ):
        super().__init__(parent, fg_color=COLOR_SCHEME['content_inside_bg'])

        self.status_callback = status_callback
        self.working_dir_callback = working_dir_callback

        # State
        self.detected_ligands: List[LigandInfo] = []
        self.ligand_widgets: Dict[str, Dict[str, Any]] = {}
        self.parametrized_ligands: Dict[str, Dict[str, str]] = {}
        self.current_pdb_file: Optional[str] = None
        self._2d_images: Dict[str, Any] = {}  # Keep refs to prevent GC
        self._parametrizing = False

        self._create_widgets()
        self._setup_layout()

    # ------------------------------------------------------------------
    # Widget creation
    # ------------------------------------------------------------------

    def _create_widgets(self):
        """Create all sub-widgets."""
        # Title
        self.title_label = ctk.CTkLabel(
            self,
            text="Ligand Parametrization",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )

        # Description
        self.info_label = ctk.CTkLabel(
            self,
            text="Detect and parametrize non-standard residues (ligands) "
                 "for membrane system building. Uses AMBER/GAFF2.",
            font=FONTS['small'],
            text_color=COLOR_SCHEME['inactive'],
            wraplength=600,
            justify="left"
        )

        # Detect button
        self.detect_button = ctk.CTkButton(
            self,
            text="Detect Ligands from PDB",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._detect_ligands,
        )

        # Status label
        self.status_label = ctk.CTkLabel(
            self,
            text="No PDB file loaded",
            font=FONTS['small'],
            text_color=COLOR_SCHEME['inactive']
        )

        # Container for ligand cards (grid layout)
        self.ligands_container = ctk.CTkFrame(self, fg_color="transparent")

        # Global charge method selector
        self.global_frame = ctk.CTkFrame(self, fg_color="transparent")

        self.atom_type_label = ctk.CTkLabel(
            self.global_frame,
            text="Atom Type:",
            font=FONTS['body']
        )

        atom_type_values = [
            f"{k} - {v}" for k, v in ATOM_TYPES.items()
        ]
        self.atom_type_combo = ctk.CTkComboBox(
            self.global_frame,
            values=atom_type_values,
            width=140,
            height=WIDGET_SIZES['combobox_height'],
            command=self._on_combo_changed,
        )
        self.atom_type_combo.set(
            f"{DEFAULT_ATOM_TYPE} - {ATOM_TYPES[DEFAULT_ATOM_TYPE]}"
        )

        self.charge_method_label = ctk.CTkLabel(
            self.global_frame,
            text="Charge Method:",
            font=FONTS['body']
        )

        charge_method_values = [
            f"{k} - {v}" for k, v in CHARGE_METHODS.items()
        ]
        self.charge_method_combo = ctk.CTkComboBox(
            self.global_frame,
            values=charge_method_values,
            width=220,
            height=WIDGET_SIZES['combobox_height'],
            command=self._on_combo_changed,
        )
        self.charge_method_combo.set(
            f"{DEFAULT_CHARGE_METHOD} - "
            f"{CHARGE_METHODS[DEFAULT_CHARGE_METHOD]}"
        )

        # Warning label for non-recommended combos
        self.combo_warning_label = ctk.CTkLabel(
            self,
            text="",
            font=FONTS['small'],
            text_color="#e8a838",
            wraplength=600,
            justify="left",
        )

        # Parametrize all button
        self.parametrize_button = ctk.CTkButton(
            self,
            text="Parametrize All Ligands",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._parametrize_all,
            state="disabled",
        )

        # Progress label
        self.progress_label = ctk.CTkLabel(
            self,
            text="",
            font=FONTS['small'],
            text_color=COLOR_SCHEME['text']
        )

    def _setup_layout(self):
        """Layout all widgets."""
        pad_m = LAYOUT['padding_medium']
        pad_s = LAYOUT['padding_small']

        self.title_label.pack(anchor="w", padx=pad_m, pady=(pad_m, pad_s))
        self.info_label.pack(anchor="w", padx=pad_m, pady=(0, pad_s))
        self.detect_button.pack(anchor="w", padx=pad_m, pady=pad_s)
        self.status_label.pack(anchor="w", padx=pad_m, pady=(0, pad_s))

        self.ligands_container.pack(fill="x", padx=pad_m, pady=pad_s)

        self.global_frame.pack(fill="x", padx=pad_m, pady=pad_s)
        self.atom_type_label.pack(side="left", padx=(0, pad_s))
        self.atom_type_combo.pack(side="left", padx=pad_s)
        self.charge_method_label.pack(side="left", padx=(pad_m, pad_s))
        self.charge_method_combo.pack(side="left", padx=pad_s)

        self.combo_warning_label.pack(anchor="w", padx=pad_m, pady=0)

        self.parametrize_button.pack(anchor="w", padx=pad_m, pady=pad_s)
        self.progress_label.pack(anchor="w", padx=pad_m, pady=(0, pad_m))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_pdb_file(self, pdb_file: str):
        """Set the PDB file to analyze."""
        if pdb_file and Path(pdb_file).exists():
            self.current_pdb_file = pdb_file
            self.status_label.configure(
                text=f"PDB loaded: {Path(pdb_file).name}",
                text_color=COLOR_SCHEME['text']
            )
            # Auto-detect on file change
            self._detect_ligands()
        else:
            self.current_pdb_file = None
            self.status_label.configure(
                text="No PDB file loaded",
                text_color=COLOR_SCHEME['inactive']
            )

    def get_parametrized_ligands(self) -> Dict[str, Dict[str, str]]:
        """Return parametrized ligand file paths."""
        return dict(self.parametrized_ligands)

    def has_parametrized_ligands(self) -> bool:
        """Check if any ligands have been parametrized."""
        return len(self.parametrized_ligands) > 0

    def has_detected_ligands(self) -> bool:
        """Check if any ligands were detected."""
        return len(self.detected_ligands) > 0

    def update_fonts(self, scaled_fonts):
        """Update fonts when scaling changes."""
        try:
            if hasattr(self, 'title_label'):
                self.title_label.configure(
                    font=scaled_fonts.get('heading', FONTS['heading'])
                )
            if hasattr(self, 'info_label'):
                self.info_label.configure(
                    font=scaled_fonts.get('small', FONTS['small'])
                )
            if hasattr(self, 'status_label'):
                self.status_label.configure(
                    font=scaled_fonts.get('small', FONTS['small'])
                )
            if hasattr(self, 'progress_label'):
                self.progress_label.configure(
                    font=scaled_fonts.get('small', FONTS['small'])
                )
            if hasattr(self, 'atom_type_label'):
                self.atom_type_label.configure(
                    font=scaled_fonts.get('small', FONTS['small'])
                )
            if hasattr(self, 'combo_warning_label'):
                self.combo_warning_label.configure(
                    font=scaled_fonts.get('small', FONTS['small'])
                )
        except Exception as e:
            logger.warning(
                f"Error updating fonts in LigandParamWidget: {e}"
            )

    # ------------------------------------------------------------------
    # Detection
    # ------------------------------------------------------------------

    def _detect_ligands(self):
        """Detect ligands in the current PDB file."""
        if (not self.current_pdb_file
                or not Path(self.current_pdb_file).exists()):
            messagebox.showwarning(
                "No PDB File", "Please select a PDB file first."
            )
            return

        try:
            self.detected_ligands = detect_ligands(self.current_pdb_file)

            # Clear existing
            self._clear_ligand_widgets()
            self.parametrized_ligands.clear()

            if not self.detected_ligands:
                self.status_label.configure(
                    text="No ligands detected in PDB file",
                    text_color=COLOR_SCHEME['inactive']
                )
                self.parametrize_button.configure(state="disabled")
                if self.status_callback:
                    self.status_callback("No ligands detected")
                return

            self.status_label.configure(
                text=f"Found {len(self.detected_ligands)} ligand(s)",
                text_color="#10B981"
            )

            # Compute adaptive image size
            num = len(self.detected_ligands)
            self._current_img_size = _compute_image_size(num)

            # Configure grid columns
            cols = min(num, MAX_COLS)
            for c in range(cols):
                self.ligands_container.grid_columnconfigure(
                    c, weight=1, uniform="ligcol"
                )

            # Create cards in a horizontal grid
            for idx, ligand in enumerate(self.detected_ligands):
                row = idx // cols
                col = idx % cols
                self._create_ligand_card(ligand, row, col)

            self.parametrize_button.configure(state="normal")

            # Check for existing parametrization files
            self._check_existing_parametrization()

            # Re-bind scroll events so new widgets respond to mouse wheel
            self.after(50, self._rebind_scroll)

            if self.status_callback:
                names = ", ".join(l.name for l in self.detected_ligands)
                self.status_callback(f"Detected ligands: {names}")

        except Exception as e:
            logger.error(f"Error detecting ligands: {e}", exc_info=True)
            self.status_label.configure(
                text=f"Error: {str(e)[:80]}",
                text_color="#EF4444"
            )

    def _clear_ligand_widgets(self):
        """Remove all ligand card widgets."""
        for name, widgets in self.ligand_widgets.items():
            if 'frame' in widgets:
                widgets['frame'].destroy()
        self.ligand_widgets.clear()
        self._2d_images.clear()
        # Reset grid configuration
        for c in range(MAX_COLS):
            self.ligands_container.grid_columnconfigure(
                c, weight=0, uniform=""
            )

    # ------------------------------------------------------------------
    # Existing parametrization detection
    # ------------------------------------------------------------------

    def _check_existing_parametrization(self):
        """Check if previous parametrization files already exist and are valid.

        For each detected ligand, look for {working_dir}/ligand_params/{LIG}/
        containing {LIG}.frcmod, {LIG}.lib, and logs/tleap.log with
        ``Exiting LEaP: Errors = 0``.
        If ALL ligands have valid existing files, auto-populate
        ``parametrized_ligands`` and inform the user.
        """
        working_dir = None
        if self.working_dir_callback:
            working_dir = self.working_dir_callback()
        if not working_dir:
            working_dir = str(Path.cwd())

        found: Dict[str, Dict[str, str]] = {}
        missing: List[str] = []

        for ligand in self.detected_ligands:
            lig_dir = Path(working_dir) / "ligand_params" / ligand.name
            frcmod = lig_dir / f"{ligand.name}.frcmod"
            lib = lig_dir / f"{ligand.name}.lib"
            tleap_log = lig_dir / "logs" / "tleap.log"

            if frcmod.is_file() and lib.is_file() and tleap_log.is_file():
                if self._tleap_log_ok(tleap_log):
                    found[ligand.name] = {
                        'frcmod': str(frcmod),
                        'lib': str(lib),
                        'mol2': str(lig_dir / f"{ligand.name}.mol2"),
                        'prmtop': str(lig_dir / f"{ligand.name}.prmtop"),
                        'inpcrd': str(lig_dir / f"{ligand.name}.inpcrd"),
                    }
                else:
                    missing.append(ligand.name)
            else:
                missing.append(ligand.name)

        if not found:
            return  # nothing cached

        # Update card status labels for found ligands
        for name in found:
            w = self.ligand_widgets.get(name, {})
            sl = w.get('status_label')
            if sl:
                sl.configure(
                    text="Previous run found",
                    text_color="#60A5FA"  # blue
                )
            # Hide per-ligand button for already-parametrized ligands
            self._show_retry_button(name, show=False)

        if missing:
            # Some ligands lack valid parametrization
            names_ok = ", ".join(found.keys())
            names_miss = ", ".join(missing)
            self.progress_label.configure(
                text=(f"Found existing results for: {names_ok}. "
                      f"Missing: {names_miss}. "
                      f"Use per-ligand buttons or "
                      f"'Re-parametrize All' to re-run everything.")
            )
            # Pre-populate found ligands so they are reused
            self.parametrized_ligands.update(found)
            return

        # ALL ligands have valid cached results
        self.parametrized_ligands = found
        self.parametrize_button.configure(
            state="normal",
            text="Re-parametrize All Ligands",
            fg_color="#10B981"
        )

        names_str = ", ".join(found.keys())
        self.progress_label.configure(
            text=(f"All ligand(s) already parametrized ({names_str}). "
                  f"Ready to build. Click 'Re-parametrize All Ligands' "
                  f"to overwrite.")
        )

        if self.status_callback:
            self.status_callback(
                f"Loaded existing ligand parametrization: {names_str}"
            )

        logger.info(
            f"Loaded existing parametrization for: {names_str}"
        )

    @staticmethod
    def _tleap_log_ok(log_path: Path) -> bool:
        """Return True if tleap.log indicates zero errors.

        Looks for a line like:
            Exiting LEaP: Errors = 0; Warnings = 1; Notes = 0.
        Only the Errors count matters; warnings/notes are acceptable.
        """
        try:
            text = log_path.read_text(errors='replace')
            # Match "Exiting LEaP: Errors = <N>"
            m = re.search(
                r'Exiting LEaP:\s*Errors\s*=\s*(\d+)', text
            )
            if m:
                return int(m.group(1)) == 0
            return False
        except Exception:
            return False

    # ------------------------------------------------------------------
    # Card creation (horizontal grid)
    # ------------------------------------------------------------------

    def _create_ligand_card(
        self, ligand: LigandInfo, row: int, col: int
    ):
        """Create a compact card for a single ligand placed in a grid cell."""
        card = ctk.CTkFrame(
            self.ligands_container,
            fg_color=COLOR_SCHEME['canvas'],
            corner_radius=8,
        )
        card.grid(row=row, column=col, padx=4, pady=4, sticky="nsew")

        # ---- Header: name + status ----
        header = ctk.CTkFrame(card, fg_color="transparent")
        header.pack(fill="x", padx=8, pady=(6, 2))

        name_label = ctk.CTkLabel(
            header,
            text=ligand.name,
            font=FONTS['heading'],
            text_color="#60A5FA"
        )
        name_label.pack(side="left")

        status_label = ctk.CTkLabel(
            header,
            text="Not parametrized",
            font=FONTS['small'],
            text_color=COLOR_SCHEME['inactive']
        )
        status_label.pack(side="right", padx=4)

        # ---- Info line ----
        info_text = (
            f"{ligand.formula}  |  {ligand.num_atoms} atoms  |  "
            f"Chain {ligand.chain}  |  Res {ligand.res_id}"
        )
        info_label = ctk.CTkLabel(
            card,
            text=info_text,
            font=FONTS['small'],
            text_color=COLOR_SCHEME['text']
        )
        info_label.pack(anchor="w", padx=10, pady=(0, 4))

        # ---- 2D structure image (CTkLabel with CTkImage) ----
        iw, ih = getattr(self, '_current_img_size', _DEFAULT_SIZE)
        image_label = ctk.CTkLabel(
            card,
            text="Generating 2D structure...",
            font=FONTS['small'],
            text_color=COLOR_SCHEME['inactive'],
            width=iw,
            height=ih,
            fg_color="#1C1C1C",
            corner_radius=6,
        )
        image_label.pack(padx=8, pady=4)

        self._generate_2d_image_async(ligand, image_label)

        # ---- Controls row: charge + multiplicity side-by-side ----
        ctrl_frame = ctk.CTkFrame(card, fg_color="transparent")
        ctrl_frame.pack(fill="x", padx=8, pady=(4, 2))

        charge_lbl = ctk.CTkLabel(
            ctrl_frame, text="Charge:", font=FONTS['small']
        )
        charge_lbl.pack(side="left", padx=(0, 2))

        charge_entry = ctk.CTkEntry(
            ctrl_frame, width=50, height=WIDGET_SIZES['entry_height']
        )
        charge_entry.insert(0, "0")
        charge_entry.pack(side="left", padx=(0, 10))

        mult_lbl = ctk.CTkLabel(
            ctrl_frame, text="Mult:", font=FONTS['small']
        )
        mult_lbl.pack(side="left", padx=(0, 2))

        mult_entry = ctk.CTkEntry(
            ctrl_frame, width=50, height=WIDGET_SIZES['entry_height']
        )
        mult_entry.insert(0, "1")
        mult_entry.pack(side="left")

        # ---- Element summary ----
        elem_text = "  ".join(
            f"{e}: {c}" for e, c in sorted(ligand.elements.items())
        )
        elem_label = ctk.CTkLabel(
            card,
            text=f"Elements: {elem_text}",
            font=FONTS['small'],
            text_color=COLOR_SCHEME['inactive']
        )
        elem_label.pack(anchor="w", padx=10, pady=(2, 2))

        # ---- Per-ligand run/retry button (visible by default) ----
        retry_button = ctk.CTkButton(
            card,
            text=f"Parametrize {ligand.name}",
            font=FONTS['small'],
            height=24,
            fg_color="#3B82F6",
            hover_color="#2563EB",
            text_color="#FFFFFF",
            corner_radius=4,
            command=lambda n=ligand.name: self._parametrize_single(n),
        )
        retry_button.pack(padx=10, pady=(2, 6), anchor="w")

        # Store references
        self.ligand_widgets[ligand.name] = {
            'frame': card,
            'charge_entry': charge_entry,
            'multiplicity_entry': mult_entry,
            'status_label': status_label,
            'image_label': image_label,
            'retry_button': retry_button,
            'ligand': ligand,
        }

    # ------------------------------------------------------------------
    # 2D molecular image
    # ------------------------------------------------------------------

    def _generate_2d_image_async(self, ligand: LigandInfo, image_label):
        """Generate 2D image in a background thread and display it."""
        def _generate():
            try:
                import tempfile
                with tempfile.NamedTemporaryFile(
                    suffix='.png', delete=False
                ) as tmp:
                    tmp_path = tmp.name

                iw, ih = getattr(self, '_current_img_size', _DEFAULT_SIZE)
                result = get_ligand_2d_image_from_pdb_lines(
                    ligand.pdb_lines, tmp_path,
                    width=iw * 2, height=ih * 2,
                    remove_nonpolar_h=True,
                )

                if result and Path(tmp_path).exists():
                    image_label.after(
                        0,
                        lambda p=tmp_path: self._display_image(
                            image_label, p, ligand.name
                        )
                    )
                else:
                    image_label.after(
                        0,
                        lambda: image_label.configure(
                            text="2D structure\nnot available",
                            image=None
                        )
                    )
            except Exception as e:
                logger.warning(
                    f"Error generating 2D image for {ligand.name}: {e}"
                )
                try:
                    image_label.after(
                        0,
                        lambda: image_label.configure(
                            text="2D error", image=None
                        )
                    )
                except Exception:
                    pass

        thread = threading.Thread(target=_generate, daemon=True)
        thread.start()

    def _display_image(
        self, image_label, image_path: str, ligand_name: str
    ):
        """Display a 2D image on a CTkLabel (main-thread only)."""
        try:
            iw, ih = getattr(self, '_current_img_size', _DEFAULT_SIZE)
            img = PILImage.open(image_path)
            img = img.resize((iw, ih), PILImage.LANCZOS)
            ctk_img = ctk.CTkImage(
                dark_image=img, light_image=img,
                size=(iw, ih)
            )
            image_label.configure(image=ctk_img, text="")
            # Keep reference to prevent garbage collection
            self._2d_images[ligand_name] = ctk_img
        except Exception as e:
            logger.warning(f"Error displaying 2D image: {e}")
            image_label.configure(
                text="2D structure\nnot available", image=None
            )
        finally:
            try:
                import os
                os.unlink(image_path)
            except Exception:
                pass

    # ------------------------------------------------------------------
    # Mouse-wheel scroll propagation
    # ------------------------------------------------------------------

    def _rebind_scroll(self):
        """Bind mouse-wheel events on all children so scrolling works
        inside the ligand section.  Events are forwarded to the parent
        CTkScrollableFrame canvas."""
        scroll_canvas = self._find_scroll_canvas()
        if scroll_canvas is None:
            return

        def _on_mousewheel(event):
            try:
                scroll_canvas.yview_scroll(
                    int(-1 * (event.delta / 120)), "units"
                )
            except Exception:
                pass
            return "break"

        def _on_mousewheel_linux(event):
            try:
                if event.num == 4:
                    scroll_canvas.yview_scroll(-3, "units")
                elif event.num == 5:
                    scroll_canvas.yview_scroll(3, "units")
            except Exception:
                pass
            return "break"

        def _bind_recursive(widget):
            try:
                widget.bind("<MouseWheel>", _on_mousewheel)
                widget.bind("<Button-4>", _on_mousewheel_linux)
                widget.bind("<Button-5>", _on_mousewheel_linux)
                for child in widget.winfo_children():
                    _bind_recursive(child)
            except Exception:
                pass

        _bind_recursive(self)

    def _find_scroll_canvas(self):
        """Walk up the widget tree to find the parent CTkScrollableFrame
        and return its internal canvas."""
        parent = self.master
        while parent:
            if isinstance(parent, ctk.CTkScrollableFrame):
                if hasattr(parent, '_parent_canvas'):
                    return parent._parent_canvas
            parent = getattr(parent, 'master', None)
        return None

    # ------------------------------------------------------------------
    # Atom type / Charge method
    # ------------------------------------------------------------------

    def _get_selected_atom_type(self) -> str:
        """Get the selected atom type code."""
        value = self.atom_type_combo.get()
        if " - " in value:
            return value.split(" - ")[0].strip()
        return DEFAULT_ATOM_TYPE

    def _get_selected_charge_method(self) -> str:
        """Get the selected charge method code."""
        value = self.charge_method_combo.get()
        if " - " in value:
            return value.split(" - ")[0].strip()
        return DEFAULT_CHARGE_METHOD

    def _on_combo_changed(self, _event=None):
        """Check atom-type / charge-method pairing and show warning if needed."""
        at = self._get_selected_atom_type()
        cm = self._get_selected_charge_method()
        if (at, cm) in NON_RECOMMENDED_COMBOS:
            self.combo_warning_label.configure(
                text=f"Warning: {at}/{cm} is not recommended. "
                     f"Use gaff/bcc or gaff2/abcg2 instead (AMBER manual)."
            )
        else:
            self.combo_warning_label.configure(text="")

    # ------------------------------------------------------------------
    # Parametrization
    # ------------------------------------------------------------------

    def _show_retry_button(self, ligand_name: str, show: bool = True):
        """Show or hide the per-ligand button."""
        w = self.ligand_widgets.get(ligand_name, {})
        btn = w.get('retry_button')
        if not btn:
            return
        if show:
            # Determine text/color based on whether it already completed once
            if ligand_name in self.parametrized_ligands:
                btn.configure(
                    text=f"Re-parametrize {ligand_name}",
                    state="normal",
                    fg_color="#3B82F6",
                )
            else:
                btn.configure(
                    text=f"Parametrize {ligand_name}",
                    state="normal",
                    fg_color="#3B82F6",
                )
            btn.pack(padx=10, pady=(2, 6), anchor="w")
        else:
            btn.pack_forget()

    def _show_retry_button_failed(self, ligand_name: str):
        """Show per-ligand button in failed/retry style."""
        w = self.ligand_widgets.get(ligand_name, {})
        btn = w.get('retry_button')
        if not btn:
            return
        btn.configure(
            text=f"Retry {ligand_name}",
            state="normal",
            fg_color="#F59E0B",
            hover_color="#D97706",
        )
        btn.pack(padx=10, pady=(2, 6), anchor="w")

    def _parametrize_single(self, ligand_name: str):
        """Parametrize (or retry) a single ligand."""
        if self._parametrizing:
            messagebox.showinfo(
                "In Progress", "Parametrization is already running."
            )
            return

        # Find the ligand info
        ligand = None
        for lig in self.detected_ligands:
            if lig.name == ligand_name:
                ligand = lig
                break
        if ligand is None:
            return

        # Working directory
        working_dir = None
        if self.working_dir_callback:
            working_dir = self.working_dir_callback()
        if not working_dir:
            working_dir = str(Path.cwd())

        charge_method = self._get_selected_charge_method()
        atom_type = self._get_selected_atom_type()

        # Warn for non-recommended combos
        if (atom_type, charge_method) in NON_RECOMMENDED_COMBOS:
            proceed = messagebox.askyesno(
                "Non-Recommended Combination",
                f"The combination {atom_type}/{charge_method} is not "
                f"recommended by the AMBER manual.\n\n"
                f"Recommended pairings are gaff/bcc and gaff2/abcg2.\n\n"
                f"Do you want to proceed anyway?",
                icon="warning",
            )
            if not proceed:
                return

        w = self.ligand_widgets.get(ligand_name, {})
        try:
            charge = int(w['charge_entry'].get())
        except (KeyError, ValueError):
            charge = 0
        try:
            multiplicity = int(w['multiplicity_entry'].get())
        except (KeyError, ValueError):
            multiplicity = 1

        self._parametrizing = True
        # Hide retry button while running
        self._show_retry_button(ligand_name, show=False)
        sl = w.get('status_label')
        if sl:
            sl.configure(text="Running...", text_color="#F59E0B")
        self.progress_label.configure(
            text=f"Parametrizing {ligand_name}..."
        )

        def _run():
            try:
                lig_dir = str(
                    Path(working_dir) / "ligand_params" / ligand_name
                )
                lig_pdb = extract_ligand_pdb(
                    self.current_pdb_file, ligand_name, lig_dir
                )
                files = parametrize_ligand(
                    ligand_pdb=lig_pdb,
                    ligand_name=ligand_name,
                    output_dir=lig_dir,
                    charge=charge,
                    charge_method=charge_method,
                    atom_type=atom_type,
                    multiplicity=multiplicity,
                )
                self.parametrized_ligands[ligand_name] = files

                if sl:
                    sl.after(0, lambda s=sl: s.configure(
                        text="Completed", text_color="#10B981"
                    ))
                # Hide the per-ligand button after success
                self.after(
                    0,
                    lambda: self._show_retry_button(
                        ligand_name, show=False
                    )
                )

                # Check if all ligands are now parametrized
                all_done = all(
                    l.name in self.parametrized_ligands
                    for l in self.detected_ligands
                )
                if all_done:
                    total = len(self.detected_ligands)
                    self.progress_label.after(
                        0,
                        lambda: self.progress_label.configure(
                            text=f"All {total} ligand(s) parametrized "
                                 f"successfully"
                        )
                    )
                    self.parametrize_button.after(
                        0,
                        lambda: self.parametrize_button.configure(
                            state="normal",
                            text="Re-parametrize All Ligands",
                            fg_color="#10B981"
                        )
                    )
                else:
                    done = sum(
                        1 for l in self.detected_ligands
                        if l.name in self.parametrized_ligands
                    )
                    total = len(self.detected_ligands)
                    self.progress_label.after(
                        0,
                        lambda: self.progress_label.configure(
                            text=f"{ligand_name} completed "
                                 f"({done}/{total} done)"
                        )
                    )

                if self.status_callback:
                    self.status_callback(
                        f"Ligand {ligand_name} parametrized successfully"
                    )

            except (LigandParametrizationError, Exception) as e:
                err_msg = str(e)
                logger.error(
                    f"Parametrization error for {ligand_name}: {err_msg}",
                    exc_info=True
                )
                if sl:
                    sl.after(0, lambda s=sl: s.configure(
                        text="Failed", text_color="#EF4444"
                    ))
                self.after(
                    0,
                    lambda: self._show_retry_button_failed(ligand_name)
                )

                sync_hint = ""
                if self._is_cloud_sync_error(err_msg):
                    sync_hint = (
                        " -- Pause cloud sync and retry."
                    )
                display_msg = f"Error ({ligand_name}): " \
                              f"{err_msg[:80]}{sync_hint}"
                self.progress_label.after(
                    0,
                    lambda m=display_msg:
                        self.progress_label.configure(text=m)
                )
            finally:
                self._parametrizing = False

        thread = threading.Thread(target=_run, daemon=True)
        thread.start()

    def _parametrize_all(self):
        """Parametrize all detected ligands."""
        if self._parametrizing:
            messagebox.showinfo(
                "In Progress", "Parametrization is already running."
            )
            return

        if not self.detected_ligands:
            return

        # Working directory
        working_dir = None
        if self.working_dir_callback:
            working_dir = self.working_dir_callback()
        if not working_dir:
            working_dir = str(Path.cwd())

        charge_method = self._get_selected_charge_method()
        atom_type = self._get_selected_atom_type()

        # Warn for non-recommended atom-type / charge-method combos
        if (atom_type, charge_method) in NON_RECOMMENDED_COMBOS:
            proceed = messagebox.askyesno(
                "Non-Recommended Combination",
                f"The combination {atom_type}/{charge_method} is not "
                f"recommended by the AMBER manual.\n\n"
                f"Recommended pairings are gaff/bcc and gaff2/abcg2.\n\n"
                f"Do you want to proceed anyway?",
                icon="warning",
            )
            if not proceed:
                return

        # Collect per-ligand settings
        charges: Dict[str, int] = {}
        multiplicities: Dict[str, int] = {}
        for ligand in self.detected_ligands:
            w = self.ligand_widgets.get(ligand.name, {})
            try:
                charges[ligand.name] = int(w['charge_entry'].get())
            except (KeyError, ValueError):
                charges[ligand.name] = 0
            try:
                multiplicities[ligand.name] = int(
                    w['multiplicity_entry'].get()
                )
            except (KeyError, ValueError):
                multiplicities[ligand.name] = 1

        # Hide all retry buttons before starting
        for ligand in self.detected_ligands:
            self._show_retry_button(ligand.name, show=False)

        self._parametrizing = True
        self.parametrize_button.configure(
            state="disabled", text="Parametrizing..."
        )
        self.progress_label.configure(text="Starting parametrization...")

        # Background thread
        def _run():
            results: Dict[str, Dict[str, str]] = {}
            failed: List[str] = []
            total = len(self.detected_ligands)

            for i, ligand in enumerate(self.detected_ligands, 1):
                # Progress
                self.progress_label.after(
                    0,
                    lambda n=ligand.name, idx=i, t=total:
                        self.progress_label.configure(
                            text=f"Parametrizing {n} ({idx}/{t})..."
                        )
                )

                # Status -> Running
                w = self.ligand_widgets.get(ligand.name, {})
                sl = w.get('status_label')
                if sl:
                    sl.after(0, lambda s=sl: s.configure(
                        text="Running...", text_color="#F59E0B"
                    ))

                try:
                    lig_dir = str(
                        Path(working_dir) / "ligand_params" / ligand.name
                    )

                    lig_pdb = extract_ligand_pdb(
                        self.current_pdb_file, ligand.name, lig_dir
                    )

                    files = parametrize_ligand(
                        ligand_pdb=lig_pdb,
                        ligand_name=ligand.name,
                        output_dir=lig_dir,
                        charge=charges.get(ligand.name, 0),
                        charge_method=charge_method,
                        atom_type=atom_type,
                        multiplicity=multiplicities.get(ligand.name, 1),
                    )

                    results[ligand.name] = files

                    if sl:
                        sl.after(0, lambda s=sl: s.configure(
                            text="Completed", text_color="#10B981"
                        ))
                    # Hide per-ligand button on success
                    _ok_name = ligand.name
                    self.after(
                        0,
                        lambda n=_ok_name: self._show_retry_button(
                            n, show=False
                        )
                    )

                except (LigandParametrizationError, Exception) as e:
                    err_msg = str(e)
                    logger.error(
                        f"Parametrization error for {ligand.name}: "
                        f"{err_msg}", exc_info=True
                    )
                    failed.append(ligand.name)

                    if sl:
                        sl.after(0, lambda s=sl: s.configure(
                            text="Failed", text_color="#EF4444"
                        ))
                    # Show amber retry button for the failed ligand
                    _name = ligand.name
                    self.after(
                        0,
                        lambda n=_name: self._show_retry_button_failed(n)
                    )

            # Store partial (or full) results
            self.parametrized_ligands.update(results)

            if not failed:
                # All succeeded
                self.progress_label.after(
                    0,
                    lambda: self.progress_label.configure(
                        text=(f"All {total} ligand(s) parametrized "
                              f"successfully")
                    )
                )
                self.parametrize_button.after(
                    0,
                    lambda: self.parametrize_button.configure(
                        state="normal",
                        text="Re-parametrize All Ligands",
                        fg_color="#10B981"
                    )
                )
                if self.status_callback:
                    self.status_callback(
                        f"Ligand parametrization complete: "
                        f"{list(results.keys())}"
                    )
            else:
                # Some failed
                ok_count = len(results)
                fail_names = ", ".join(failed)
                summary = (
                    f"{ok_count}/{total} succeeded. "
                    f"Failed: {fail_names}. "
                    f"Use per-ligand Retry buttons or "
                    f"'Re-parametrize All' to re-run everything."
                )
                self.progress_label.after(
                    0,
                    lambda m=summary:
                        self.progress_label.configure(text=m)
                )
                btn_color = "#10B981" if ok_count > 0 else "#EF4444"
                self.parametrize_button.after(
                    0,
                    lambda c=btn_color: self.parametrize_button.configure(
                        state="normal",
                        text="Re-parametrize All Ligands",
                        fg_color=c
                    )
                )

                # Detect cloud-sync permission errors from last failure
                last_fail = failed[-1]
                for lig in self.detected_ligands:
                    if lig.name == last_fail:
                        break

                if self.status_callback:
                    self.status_callback(
                        f"Ligand parametrization: {ok_count}/{total} "
                        f"succeeded, failed: {fail_names}"
                    )

            self._parametrizing = False

        thread = threading.Thread(target=_run, daemon=True)
        thread.start()

    # ------------------------------------------------------------------
    # Cloud sync error detection
    # ------------------------------------------------------------------

    @staticmethod
    def _is_cloud_sync_error(error_msg: str) -> bool:
        """Detect if an error is likely caused by cloud sync file locking.

        Checks for PermissionError / Errno 13 patterns and common cloud
        service paths (Dropbox, OneDrive, Google Drive, iCloud).
        """
        msg = error_msg.lower()
        # Permission denied indicators
        perm_keywords = [
            'permission denied',
            'errno 13',
            '[errno 13]',
            'access is denied',
            'winerror 5',
            'winerror 32',      # sharing violation
            'being used by another process',
        ]
        if not any(kw in msg for kw in perm_keywords):
            return False
        # Cloud-sync path hints
        cloud_hints = [
            'dropbox', 'onedrive', 'google drive', 'googledrive',
            'icloud', 'box sync', 'mega', 'pcloud', 'syncthing',
            'nextcloud', 'seafile',
        ]
        return any(h in msg for h in cloud_hints)
