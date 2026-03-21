# gatewizard/gui/frames/visualize.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Visualization frame for protein structure viewing using VTK.

This module provides the GUI for loading and visualizing protein structures
using VTK-based 3D rendering with offscreen rendering to a tkinter Canvas.
"""

import json
import math
import os
import tkinter as tk
from tkinter import filedialog, messagebox, colorchooser
from typing import Optional, Callable, Dict, List
from pathlib import Path
from collections import defaultdict

import numpy as np

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

from gatewizard.gui.constants import COLOR_SCHEME, FONTS, FILE_FILTERS, LAYOUT
from gatewizard.core.file_manager import FileManager
from gatewizard.core.viewer import (
    parse_pdb, ProteinStructure, Selection, Residue,
    AA_NAMES, BACKBONE_NAMES, SS_COLORS, SS_LABELS, CHAIN_PALETTE,
    VDW_RADII, _assign_secondary_structure,
)
from gatewizard.gui.widgets.vtk_frame import VTKFrame
from gatewizard.gui.widgets.collapsible_section import CollapsibleSection


def _make_gear_image(size=16, color="white", teeth=8):
    """Draw a gear icon as a PIL image (no font/unicode needed)."""
    from PIL import Image, ImageDraw
    # Render at 2x for antialiased quality, then resize
    s2 = size * 2
    img = Image.new("RGBA", (s2, s2), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    cx, cy = s2 / 2, s2 / 2
    outer_r = s2 / 2 - 1
    inner_r = outer_r * 0.65
    tooth_half = math.pi / (teeth * 2) * 0.75
    points = []
    for i in range(teeth):
        ca = 2 * math.pi * i / teeth - math.pi / 2
        a1 = ca - tooth_half * 1.6
        points.append((cx + inner_r * math.cos(a1), cy + inner_r * math.sin(a1)))
        a2 = ca - tooth_half
        points.append((cx + outer_r * math.cos(a2), cy + outer_r * math.sin(a2)))
        a3 = ca + tooth_half
        points.append((cx + outer_r * math.cos(a3), cy + outer_r * math.sin(a3)))
        a4 = ca + tooth_half * 1.6
        points.append((cx + inner_r * math.cos(a4), cy + inner_r * math.sin(a4)))
    draw.polygon(points, fill=color)
    hole_r = inner_r * 0.45
    draw.ellipse([cx - hole_r, cy - hole_r, cx + hole_r, cy + hole_r],
                 fill=(0, 0, 0, 0))
    img = img.resize((size, size), Image.LANCZOS)
    return ctk.CTkImage(light_image=img, dark_image=img, size=(size, size))
from gatewizard.gui.rendering import (
    make_vdw_actor, make_ball_stick_actor, make_stick_actor,
    make_cartoon_actor, make_tube_ss_actor, make_backbone_actor,
    make_surface_actor, QUALITY_PRESETS, QUALITY_LABELS, MATERIAL_PRESETS,
)
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)


# -- UI constants ----------------------------------------------------------

REPRESENTATIONS_UI = ['VDW (Spacefill)', 'Ball & Stick', 'Sticks', 'Cartoon',
                      'Tube SS', 'Backbone', 'Surface']
COLOR_SCHEMES_UI = ['Element (CPK)', 'Chain', 'Secondary Structure',
                    'Residue Nature']
SELECTION_CRITERIA = ['All', 'Protein', 'Backbone', 'Sidechain', 'Water',
                      'Ligand', 'Chain...', 'Residue range...',
                      'Around selection...', 'MDAnalysis expression...']

REP_LABELS: Dict[str, str] = {
    'vdw': 'VDW', 'ball_stick': 'Ball&Stick', 'sticks': 'Sticks',
    'cartoon': 'Cartoon', 'tube_ss': 'Tube SS',
    'backbone': 'Backbone', 'surface': 'Surface',
}
REP_VALUES = list(REP_LABELS.values())
REP_KEYS = list(REP_LABELS.keys())


def _rep_key(label: str) -> str:
    m = {'VDW (Spacefill)': 'vdw', 'Ball & Stick': 'ball_stick',
         'Sticks': 'sticks', 'Cartoon': 'cartoon', 'Tube SS': 'tube_ss',
         'Backbone': 'backbone', 'Surface': 'surface'}
    return m.get(label, 'ball_stick')


def _cs_key(label: str) -> str:
    m = {'Element (CPK)': 'element', 'Chain': 'chain',
         'Secondary Structure': 'ss', 'Residue Nature': 'residue_nature',
         'Uniform': 'element',
         'element': 'element', 'chain': 'chain', 'ss': 'ss',
         'uniform': 'uniform', 'residue_nature': 'residue_nature'}
    return m.get(label, 'element')


def _sep(parent):
    ctk.CTkFrame(parent, height=2, fg_color="gray30").pack(fill="x", pady=8, padx=8)


# --------------------------------------------------------------------------
# Main frame
# --------------------------------------------------------------------------

class VisualizeFrame(ctk.CTkFrame):
    """VTK-based protein structure visualization frame."""

    def __init__(
        self,
        parent,
        pdb_changed_callback: Optional[Callable[[Optional[str]], None]] = None,
        status_callback: Optional[Callable[[str], None]] = None,
        initial_directory: Optional[str] = None,
    ):
        super().__init__(parent, fg_color=COLOR_SCHEME['content_bg'])

        self.pdb_changed_callback = pdb_changed_callback
        self.status_callback = status_callback
        self.initial_directory = initial_directory or str(Path.cwd())

        self.current_pdb_file = None
        self.file_manager = FileManager()
        self.structure: Optional[ProteinStructure] = None
        self.selections: List[Selection] = []
        self._pdb_filepath: Optional[str] = None

        self._ssao_enabled = False
        self._shadows_enabled = False
        self._shadow_light = None
        self._depth_cue_density = 0.0
        self._drag_data: dict = {'idx': None, 'start_y': 0, 'frames': []}
        self._ref_lines_visible = False
        self._axes_mode_current: Optional[str] = None

        self._build_ui()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        self.main_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.main_frame.pack(fill="both", expand=True, padx=6, pady=6)

        # --- Left panel: controls ---
        self.ctrl = ctk.CTkScrollableFrame(self.main_frame, width=220,
                                           label_text="Controls")
        self.ctrl.pack(side="left", fill="y", padx=(0, 4))

        # Load Structure (expanded by default)
        sec_load = CollapsibleSection(self.ctrl, "Load Structure", expanded=True)
        sec_load.pack(fill="x")
        ctk.CTkButton(sec_load.content, text="Open PDB File",
                      command=self._browse_pdb,
                      height=34).pack(pady=4, padx=8, fill="x")
        self.pdb_id_entry = ctk.CTkEntry(sec_load.content,
                                         placeholder_text="PDB ID (e.g. 1CRN)")
        self.pdb_id_entry.pack(pady=2, padx=8, fill="x")
        ctk.CTkButton(sec_load.content, text="Download PDB",
                      command=self._download_pdb,
                      height=34).pack(pady=4, padx=8, fill="x")

        # Edit Structure (collapsed by default)
        sec_edit = CollapsibleSection(self.ctrl, "Edit Structure", expanded=False)
        sec_edit.pack(fill="x")
        ctk.CTkButton(sec_edit.content, text="Rename Chain",
                      command=self._edit_chain_name,
                      height=30, fg_color="gray35").pack(pady=2, padx=8, fill="x")
        ctk.CTkButton(sec_edit.content, text="Rename Residues",
                      command=self._edit_residue_name,
                      height=30, fg_color="gray35").pack(pady=2, padx=8, fill="x")
        ctk.CTkButton(sec_edit.content, text="Renumber Residues",
                      command=self._edit_residue_numbering,
                      height=30, fg_color="gray35").pack(pady=2, padx=8, fill="x")
        ctk.CTkButton(sec_edit.content, text="Delete Selection Atoms",
                      command=self._delete_selection_atoms,
                      height=30, fg_color="#8a3a3a",
                      hover_color="#cc3333").pack(pady=2, padx=8, fill="x")

        # Save (collapsed by default)
        sec_save = CollapsibleSection(self.ctrl, "Save", expanded=False)
        sec_save.pack(fill="x")
        ctk.CTkButton(sec_save.content, text="Save as PDB",
                      command=self._save_pdb,
                      height=32, fg_color="#2a6e2a",
                      hover_color="#358535").pack(pady=3, padx=8, fill="x")
        ctk.CTkButton(sec_save.content, text="Save Image",
                      command=self._save_image,
                      height=32, fg_color="gray35").pack(pady=3, padx=8, fill="x")

        # Viewpoint (collapsed by default)
        sec_vp = CollapsibleSection(self.ctrl, "Viewpoint", expanded=False)
        sec_vp.pack(fill="x")
        ctk.CTkButton(sec_vp.content, text="Save Viewpoint",
                      command=self._save_viewpoint,
                      height=30, fg_color="gray35").pack(pady=2, padx=8, fill="x")
        ctk.CTkButton(sec_vp.content, text="Load Viewpoint",
                      command=self._load_viewpoint,
                      height=30, fg_color="gray35").pack(pady=2, padx=8, fill="x")

        _sep(self.ctrl)
        self.info_label = ctk.CTkLabel(self.ctrl,
                                       text="Load a structure to begin",
                                       wraplength=200)
        self.info_label.pack(pady=8, padx=8)

        # --- Centre: VTK 3D view ---
        self.vtk_frame = VTKFrame(self.main_frame, width=900, height=700)
        self.vtk_frame.pack(side="left", fill="both", expand=True, padx=4)
        self.vtk_frame._pick_callback = self._on_pick
        self.vtk_frame._right_click_callback = self._on_right_click_atom
        self.vtk_frame._resize_callback = self._on_vtk_resize
        self.vtk_frame._post_render_callback = self._post_render

        # Atom labels & measurements state
        self._atom_labels: list = []       # [(actor, coord, text), ...]
        self._label_font_size: int = 14
        self._label_color: tuple = (1.0, 1.0, 1.0)  # RGB 0-1
        self._measure_mode: Optional[str] = None   # 'distance'|'angle'|'dihedral'
        self._measure_picks: list = []     # accumulated picked atoms
        self._measure_actors: list = []    # [(actor, ...), ...]
        self._measure_highlight_actors: list = []  # yellow spheres during pick
        self._measure_border_active: bool = False   # whether border should show
        self._measure_font_size: int = 14
        self._measure_label_color: tuple = (1.0, 1.0, 0.0)  # yellow
        self._measure_line_color: tuple = (1.0, 1.0, 0.0)   # yellow
        self._measure_line_width: float = 2.0

        # --- Right panel ---
        right_frame = ctk.CTkFrame(self.main_frame, fg_color="transparent")
        right_frame.pack(side="right", fill="y", padx=(4, 0))

        # View section (expanded) — merges Quick Actions + View Helpers
        sec_view = CollapsibleSection(right_frame, "View", expanded=True)
        sec_view.pack(fill="x")
        btn_row1 = ctk.CTkFrame(sec_view.content, fg_color="transparent")
        btn_row1.pack(fill="x", pady=(2, 2), padx=2)
        ctk.CTkButton(btn_row1, text="Auto-Detect", width=85, height=26,
                      font=("", 11), fg_color="gray35",
                      command=self._auto_detect_molecules).pack(
                          side="left", padx=2)
        ctk.CTkButton(btn_row1, text="Reset View", width=75, height=26,
                      font=("", 11), fg_color="gray35",
                      command=self._reset_view).pack(side="left", padx=2)
        btn_row2 = ctk.CTkFrame(sec_view.content, fg_color="transparent")
        btn_row2.pack(fill="x", pady=(0, 2), padx=2)
        ctk.CTkButton(btn_row2, text="BG Color", width=70, height=26,
                      font=("", 11), fg_color="gray35",
                      command=self._change_bg_color).pack(side="left", padx=2)
        self._axes_btn = ctk.CTkButton(btn_row2, text="XYZ Axes", width=80,
                      height=26, font=("", 11), fg_color="gray35",
                      command=self._axes_dialog)
        self._axes_btn.pack(side="left", padx=2)
        self._ref_lines_btn = ctk.CTkButton(btn_row2, text="Ref. Lines",
                      width=80, height=26, font=("", 11), fg_color="gray35",
                      command=self._toggle_ref_lines)
        self._ref_lines_btn.pack(side="left", padx=2)

        # Measurement section (collapsed)
        sec_meas = CollapsibleSection(right_frame, "Measurement", expanded=False)
        sec_meas.pack(fill="x")
        meas_row1 = ctk.CTkFrame(sec_meas.content, fg_color="transparent")
        meas_row1.pack(fill="x", pady=(2, 1), padx=2)
        self._measure_var = ctk.StringVar(value="Measure")
        self._measure_menu = ctk.CTkOptionMenu(
            meas_row1, variable=self._measure_var,
            values=["Distance", "Angle", "Dihedral"],
            width=110, height=26, font=("", 11),
            fg_color="gray35", button_color="gray45",
            command=self._on_measure_selected)
        self._measure_menu.pack(side="left", padx=2)
        ctk.CTkButton(meas_row1, text="x", width=28, height=26,
                      font=("", 11), fg_color="#8a3a3a",
                      hover_color="#cc3333",
                      command=self._remove_all_measurements
                      ).pack(side="left", padx=2)
        meas_row2 = ctk.CTkFrame(sec_meas.content, fg_color="transparent")
        meas_row2.pack(fill="x", pady=(1, 2), padx=2)
        ctk.CTkLabel(meas_row2, text="Size:", font=("", 10),
                     anchor="w").pack(side="left", padx=(4, 0))
        self._meas_size_var = ctk.StringVar(value="14")
        ctk.CTkOptionMenu(
            meas_row2, variable=self._meas_size_var,
            values=[str(s) for s in [8, 10, 12, 14, 16, 18, 20, 24, 28]],
            width=55, height=24, font=("", 10),
            fg_color="gray30", button_color="gray40",
            command=self._on_meas_size_change
        ).pack(side="left", padx=2)
        self._meas_color_btn = ctk.CTkButton(
            meas_row2, text="", width=24, height=24,
            corner_radius=4, fg_color="#ffff00",
            hover_color="#dddd00", border_width=1,
            border_color="gray50",
            command=self._on_meas_color_pick)
        self._meas_color_btn.pack(side="left", padx=2)
        meas_row3 = ctk.CTkFrame(sec_meas.content, fg_color="transparent")
        meas_row3.pack(fill="x", pady=(1, 2), padx=2)
        ctk.CTkLabel(meas_row3, text="Line:", font=("", 10),
                     anchor="w").pack(side="left", padx=(4, 0))
        self._meas_width_var = ctk.StringVar(value="2")
        ctk.CTkOptionMenu(
            meas_row3, variable=self._meas_width_var,
            values=["1", "2", "3", "4", "5", "6"],
            width=50, height=24, font=("", 10),
            fg_color="gray30", button_color="gray40",
            command=self._on_meas_width_change
        ).pack(side="left", padx=2)
        self._meas_line_color_btn = ctk.CTkButton(
            meas_row3, text="", width=24, height=24,
            corner_radius=4, fg_color="#ffff00",
            hover_color="#dddd00", border_width=1,
            border_color="gray50",
            command=self._on_meas_line_color_pick)
        self._meas_line_color_btn.pack(side="left", padx=2)

        # Labels section (collapsed)
        sec_labels = CollapsibleSection(right_frame, "Labels", expanded=False)
        sec_labels.pack(fill="x")
        lbl_frame = ctk.CTkFrame(sec_labels.content, fg_color="transparent")
        lbl_frame.pack(fill="x", pady=(2, 2), padx=2)
        self._label_size_var = ctk.StringVar(value="14")
        ctk.CTkLabel(lbl_frame, text="Size:", font=("", 10),
                     anchor="w").pack(side="left", padx=(4, 0))
        size_spin = ctk.CTkOptionMenu(
            lbl_frame, variable=self._label_size_var,
            values=[str(s) for s in [8, 10, 12, 14, 16, 18, 20, 24, 28]],
            width=55, height=24, font=("", 10),
            fg_color="gray30", button_color="gray40",
            command=self._on_label_size_change)
        size_spin.pack(side="left", padx=2)
        self._label_color_btn = ctk.CTkButton(
            lbl_frame, text="", width=24, height=24,
            corner_radius=4, fg_color="#ffffff",
            hover_color="#dddddd", border_width=1,
            border_color="gray50",
            command=self._on_label_color_pick)
        self._label_color_btn.pack(side="left", padx=2)

        # Rendering section (collapsed)
        sec_render = CollapsibleSection(right_frame, "Rendering", expanded=False)
        sec_render.pack(fill="x")
        render_content = sec_render.content

        row_ao = ctk.CTkFrame(render_content, fg_color="transparent")
        row_ao.pack(fill="x", padx=6, pady=1)
        ctk.CTkLabel(row_ao, text="Ambient Occlusion", font=("", 11),
                     width=130, anchor="w").pack(side="left")
        self._ssao_var = ctk.StringVar(value="off")
        self._ssao_switch = ctk.CTkSwitch(
            row_ao, text="", variable=self._ssao_var,
            onvalue="on", offvalue="off", width=40,
            command=self._update_render_passes)
        self._ssao_switch.pack(side="right", padx=4)

        row_sh = ctk.CTkFrame(render_content, fg_color="transparent")
        row_sh.pack(fill="x", padx=6, pady=1)
        ctk.CTkLabel(row_sh, text="Shadows", font=("", 11),
                     width=130, anchor="w").pack(side="left")
        self._shadows_var = ctk.StringVar(value="off")
        self._shadows_switch = ctk.CTkSwitch(
            row_sh, text="", variable=self._shadows_var,
            onvalue="on", offvalue="off", width=40,
            command=self._update_render_passes)
        self._shadows_switch.pack(side="right", padx=4)

        row_dc = ctk.CTkFrame(render_content, fg_color="transparent")
        row_dc.pack(fill="x", padx=6, pady=1)
        ctk.CTkLabel(row_dc, text="Depth Cueing", font=("", 11),
                     width=130, anchor="w").pack(side="left")
        self._depth_cue_var = ctk.DoubleVar(value=0.0)
        self._depth_cue_slider = ctk.CTkSlider(
            row_dc, from_=0.0, to=1.0, variable=self._depth_cue_var,
            width=100, command=self._on_depth_cue_change)
        self._depth_cue_slider.pack(side="right", padx=4)

        row_persp = ctk.CTkFrame(render_content, fg_color="transparent")
        row_persp.pack(fill="x", padx=6, pady=(1, 4))
        ctk.CTkLabel(row_persp, text="Perspective", font=("", 11),
                     width=130, anchor="w").pack(side="left")
        self._persp_var = ctk.DoubleVar(value=0.0)
        self._persp_slider = ctk.CTkSlider(
            row_persp, from_=0.0, to=1.0, variable=self._persp_var,
            width=100, command=self._on_perspective_change)
        self._persp_slider.pack(side="right", padx=4)

        # Selections section (expanded by default, collapsible)
        sec_sel = CollapsibleSection(right_frame, "Selections", expanded=True,
                                      fill_vertical=True)
        sec_sel.pack(fill="both", expand=True)
        self.sel_panel = ctk.CTkScrollableFrame(sec_sel.content, width=280,
                                                fg_color="transparent")
        self.sel_panel.pack(fill="both", expand=True)
        ctk.CTkButton(self.sel_panel, text="+ New Selection",
                      command=self._new_selection,
                      height=30, fg_color="#2a6e2a",
                      hover_color="#358535").pack(pady=(4, 4), padx=6, fill="x")
        self.sel_list_frame = ctk.CTkFrame(self.sel_panel, fg_color="transparent")
        self.sel_list_frame.pack(fill="both", expand=True)

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def _browse_pdb(self):
        path = filedialog.askopenfilename(
            title="Select PDB File",
            initialdir=self.initial_directory,
            filetypes=[("PDB files", "*.pdb"), ("All", "*.*")])
        if path:
            self._load(path)

    def _download_pdb(self):
        pid = self.pdb_id_entry.get().strip().upper()
        if not pid:
            messagebox.showwarning("Warning", "Enter a PDB ID")
            return
        try:
            import requests
            url = f"https://files.rcsb.org/download/{pid}.pdb"
            resp = requests.get(url, timeout=15)
            resp.raise_for_status()
            out = Path(self.initial_directory) / f"{pid}.pdb"
            out.write_text(resp.text)
            self._load(str(out))
        except Exception as e:
            messagebox.showerror("Error", f"Download failed: {e}")

    def _load(self, path: str):
        try:
            self.structure = parse_pdb(path)
            self._pdb_filepath = path
            s = self.structure
            title = s.title[:60] if s.title else os.path.basename(path)
            self.info_label.configure(
                text=f"{title}\n{len(s.atoms)} atoms  {len(s.residues)} residues\n"
                     f"{len(s.chains)} chain(s)  {len(s.bonds)} bonds")
            self._remove_all_labels()
            self._remove_all_measurements()
            self.selections.clear()
            self.selections.append(
                Selection("All", list(range(len(s.atoms))),
                          representation='vdw', color_scheme='element'))
            self._refresh_sel_ui()
            self._rebuild_full()
            # Notify other frames
            self.current_pdb_file = path
            if self.pdb_changed_callback:
                self.pdb_changed_callback(path)
            if self.status_callback:
                self.status_callback(f"Loaded {os.path.basename(path)}")
        except Exception as e:
            messagebox.showerror("Error", f"Parse failed: {e}")
            logger.error(f"PDB parse error: {e}")

    # ------------------------------------------------------------------
    # Scene rebuild
    # ------------------------------------------------------------------

    def _rebuild(self, *_):
        if self.structure is None:
            return
        # Build SS lookup once for all selections
        ss_map = {}
        for res in self.structure.residues:
            ss_map[(res.chain_id, res.seq_id)] = res.ss
        self.vtk_frame.clear_actors()
        for sel in self.selections:
            sel.actors.clear()
            if not sel.visible:
                continue
            atoms = [self.structure.atoms[i] for i in sel.atom_indices]
            bonds = self._local_bonds(sel.atom_indices)
            cs = sel.color_scheme
            uc = sel.uniform_color
            cc = sel.carbon_color
            rep = sel.representation
            q = sel.quality
            if rep == 'vdw':
                a = make_vdw_actor(atoms, cs, uc, scale=sel.atom_scale,
                                   carbon_color=cc, quality=q, ss_map=ss_map)
                if a:
                    sel.actors.append(a)
            elif rep == 'ball_stick':
                sel.actors.extend(make_ball_stick_actor(
                    atoms, bonds, cs, uc, ball_scale=sel.ball_scale,
                    stick_radius=sel.bond_radius, carbon_color=cc, quality=q,
                    ss_map=ss_map))
            elif rep == 'sticks':
                sel.actors.extend(make_stick_actor(
                    atoms, bonds, cs, uc, radius=sel.stick_radius,
                    carbon_color=cc, quality=q, ss_map=ss_map))
            elif rep == 'cartoon':
                sub = self._sub_struct(sel.atom_indices)
                sel.actors.extend(make_cartoon_actor(
                    sub, color_scheme=cs, uniform_color=uc, quality=q,
                    helix_w=sel.helix_width, sheet_w=sel.sheet_width,
                    coil_w=sel.coil_width, ss_colors=sel.ss_colors))
            elif rep == 'tube_ss':
                sub = self._sub_struct(sel.atom_indices)
                sel.actors.extend(make_tube_ss_actor(
                    sub, color_scheme=cs, uniform_color=uc, quality=q,
                    helix_w=sel.helix_width, sheet_w=sel.sheet_width,
                    coil_w=sel.coil_width, ss_colors=sel.ss_colors))
            elif rep == 'backbone':
                sub = self._sub_struct(sel.atom_indices)
                sel.actors.extend(make_backbone_actor(
                    sub, color_scheme=cs, uniform_color=uc, quality=q,
                    radius=sel.backbone_radius))
            elif rep == 'surface':
                c = uc if uc else (200, 200, 255)
                a = make_surface_actor(atoms, color=c, opacity=sel.opacity,
                                       resolution=sel.surface_resolution,
                                       radius=sel.surface_radius,
                                       color_scheme=cs, uniform_color=uc,
                                       carbon_color=cc, ss_map=ss_map)
                if a:
                    sel.actors.append(a)
            for actor in sel.actors:
                prop = actor.GetProperty()
                prop.SetAmbient(sel.ambient)
                prop.SetDiffuse(sel.diffuse)
                prop.SetSpecular(sel.specular)
                prop.SetSpecularPower(sel.specular_power)
                self.vtk_frame.add_actor(actor)
        # Measurement line actors and text labels both live in _label_renderer
        # (always on top) and are not cleared by clear_actors().
        self.vtk_frame.renderer.ResetCameraClippingRange()
        self._apply_depth_cueing()
        self.vtk_frame.render()

    def _rebuild_full(self):
        self._rebuild()
        self.vtk_frame.reset_camera()
        # Apply perspective setting from slider
        cam = self.vtk_frame.renderer.GetActiveCamera()
        t = self._persp_var.get()
        if t < 0.01:
            cam.SetParallelProjection(1)
        else:
            import math as _m
            angle = 5.0 + (t ** 1.5) * 55.0
            cam.SetParallelProjection(0)
            cam.SetViewAngle(angle)
        self._apply_depth_cueing()
        self._update_axes_position()
        self.vtk_frame.render()

    def _safe_refresh_rebuild(self, refresh=True):
        """Refresh UI + rebuild scene while preserving exact camera state."""
        cam = self.vtk_frame.renderer.GetActiveCamera()
        state = (cam.GetPosition(), cam.GetFocalPoint(),
                 cam.GetViewUp(), cam.GetClippingRange(),
                 cam.GetParallelScale())
        self.vtk_frame.lock_camera()
        if refresh:
            self._refresh_sel_ui()
        self._rebuild()
        cam.SetPosition(state[0])
        cam.SetFocalPoint(state[1])
        cam.SetViewUp(state[2])
        cam.SetClippingRange(state[3])
        cam.SetParallelScale(state[4])
        self._update_axes_position()
        self.vtk_frame.render()
        self.after(200, self.vtk_frame.unlock_camera)

    def _update_axes_position(self):
        """Recalculate and reapply axes position when mode is 'center'."""
        if self._axes_mode_current != "center" or not self.structure:
            return
        import numpy as _np
        vis_indices = set()
        for sel in self.selections:
            if sel.visible:
                vis_indices.update(sel.atom_indices)
        if vis_indices:
            coords = _np.array(
                [self.structure.atoms[i].coord for i in vis_indices])
        else:
            coords = _np.array([a.coord for a in self.structure.atoms])
        center = tuple(coords.mean(axis=0).tolist())
        self.vtk_frame.set_axes(mode="center", center=center)

    def _local_bonds(self, indices):
        idx_set = set(indices)
        idx_map = {g: l for l, g in enumerate(indices)}
        return [(idx_map[a], idx_map[b])
                for a, b in self.structure.bonds
                if a in idx_set and b in idx_set]

    def _sub_struct(self, indices):
        atom_ids = {id(self.structure.atoms[i]) for i in indices}
        sub = ProteinStructure()
        sub.atoms = [self.structure.atoms[i] for i in indices]
        seen = set()
        for res in self.structure.residues:
            for a in res.atoms:
                if id(a) in atom_ids:
                    key = (res.chain_id, res.seq_id)
                    if key not in seen:
                        seen.add(key)
                        sub.residues.append(res)
                        sub.chains.setdefault(res.chain_id, []).append(res)
                    break
        return sub

    # ------------------------------------------------------------------
    # Picking
    # ------------------------------------------------------------------

    def _on_pick(self, cx, cy):
        if self.structure is None:
            return
        best_atom = self._pick_atom_at(cx, cy)
        # If a measurement mode is active, delegate to measurement handler
        if self._measure_mode and best_atom is not None:
            self._handle_measure_pick(best_atom)
            return
        if best_atom is not None:
            self.info_label.configure(
                text=f"Picked: {best_atom.name} ({best_atom.res_name} "
                     f"{best_atom.res_id} {best_atom.chain_id})")
        else:
            if self._measure_mode:
                return  # don't reset view while measuring

    def _center_on_selection(self, sel):
        if not sel.atom_indices:
            return
        coords = np.array([self.structure.atoms[i].coord for i in sel.atom_indices])
        center = coords.mean(axis=0)
        cam = self.vtk_frame.renderer.GetActiveCamera()
        pos = cam.GetPosition()
        fp = cam.GetFocalPoint()
        cur_dist = math.sqrt(sum((pos[i] - fp[i]) ** 2 for i in range(3)))
        extent = np.linalg.norm(coords.max(axis=0) - coords.min(axis=0))
        dist = max(cur_dist * 0.8, extent * 1.5, 20.0)
        self.vtk_frame.focus_on_point(center, distance=dist)
        self.vtk_frame.renderer.ResetCameraClippingRange()

    def _reset_view(self):
        if self.structure:
            self.vtk_frame.reset_camera()
            self.vtk_frame.render()

    # ------------------------------------------------------------------
    # Atom picking helper (shared by pick & right-click)
    # ------------------------------------------------------------------

    def _pick_atom_at(self, cx, cy, threshold=400):
        """Return the nearest *visible* atom to canvas coords *(cx, cy)*.

        Only atoms belonging to a visible selection whose representation
        actually renders individual atoms (vdw, ball_stick, sticks) are
        considered.  Returns *None* when no suitable atom is nearby.
        """
        if self.structure is None:
            return None
        w = self.vtk_frame.render_window.GetSize()[0]
        h = self.vtk_frame.render_window.GetSize()[1]
        if w < 1 or h < 1:
            return None
        # Build set of atom indices that are actually shown
        _ATOM_REPS = {'vdw', 'ball_stick', 'sticks'}
        visible_indices = set()
        for sel in self.selections:
            if sel.visible and sel.representation in _ATOM_REPS:
                visible_indices.update(sel.atom_indices)
        if not visible_indices:
            return None
        renderer = self.vtk_frame.renderer
        best_dist = float('inf')
        best_atom = None
        atoms = self.structure.atoms
        for idx in visible_indices:
            atom = atoms[idx]
            renderer.SetWorldPoint(*atom.coord, 1.0)
            renderer.WorldToDisplay()
            dp = renderer.GetDisplayPoint()
            dx = dp[0] - cx
            dy = dp[1] - (h - cy)
            d = dx * dx + dy * dy
            if d < best_dist:
                best_dist = d
                best_atom = atom
        if best_atom is not None and best_dist < threshold:
            return best_atom
        return None

    # ------------------------------------------------------------------
    # Right-click atom label context menu
    # ------------------------------------------------------------------

    def _on_right_click_atom(self, cx, cy, event):
        """Show a context menu when right-clicking near an atom."""
        atom = self._pick_atom_at(cx, cy)

        menu = tk.Menu(self, tearoff=0)

        # Label placement submenu (only if an atom was picked)
        if atom is not None:
            label_menu = tk.Menu(menu, tearoff=0)
            formats = [
                f"{atom.chain_id}:{atom.res_id}:{atom.name}",
                f"{atom.chain_id}:{atom.res_name}{atom.res_id}",
                f"{atom.res_name}{atom.res_id}:{atom.name}",
                f"{atom.chain_id}:{atom.res_id}",
                f"{atom.res_name}{atom.res_id}",
            ]
            for text in formats:
                label_menu.add_command(
                    label=text,
                    command=lambda t=text, a=atom: self._place_atom_label(a, t))
            label_menu.add_separator()
            label_menu.add_command(
                label="Custom...",
                command=lambda a=atom: self._custom_label_dialog(a))
            menu.add_cascade(label="Add label", menu=label_menu)
            menu.add_separator()

        # Per-label removal submenu
        if self._atom_labels:
            rm_menu = tk.Menu(menu, tearoff=0)
            for idx, (_actor, _coord, text) in enumerate(self._atom_labels):
                rm_menu.add_command(
                    label=text,
                    command=lambda i=idx: self._remove_label(i))
            rm_menu.add_separator()
            rm_menu.add_command(label="Remove ALL labels",
                               command=self._remove_all_labels)
            menu.add_cascade(label="Remove label", menu=rm_menu)

        # Per-measurement removal submenu
        if self._measure_actors:
            rm_meas = tk.Menu(menu, tearoff=0)
            for idx, item in enumerate(self._measure_actors):
                rm_meas.add_command(
                    label=item[-1],  # stored description
                    command=lambda i=idx: self._remove_measurement(i))
            rm_meas.add_separator()
            rm_meas.add_command(label="Remove ALL measurements",
                                command=self._remove_all_measurements)
            menu.add_cascade(label="Remove measurement", menu=rm_meas)

        # Don't show an empty menu
        if menu.index("end") is None:
            return

        rx = self.vtk_frame.canvas.winfo_rootx() + event.x
        ry = self.vtk_frame.canvas.winfo_rooty() + event.y
        menu.tk_popup(rx, ry)

    def _place_atom_label(self, atom, text):
        """Add a 3-D text label at *atom*'s position using a vtkBillboardTextActor3D."""
        from vtkmodules.vtkRenderingCore import vtkBillboardTextActor3D
        actor = vtkBillboardTextActor3D()
        actor.SetInput(text)
        actor.SetPosition(*atom.coord)
        tp = actor.GetTextProperty()
        tp.SetFontSize(self._label_font_size)
        tp.SetColor(*self._label_color)
        tp.SetJustificationToCentered()
        tp.SetBold(True)
        tp.SetFontFamilyToCourier()
        actor.SetDisplayOffset(0, 12)
        # Add to overlay label renderer so labels are never occluded
        self.vtk_frame._label_renderer.AddActor(actor)
        self._atom_labels.append((actor, atom.coord, text))
        self.vtk_frame.render()

    def _custom_label_dialog(self, atom):
        """Prompt user for custom label text, then place it on *atom*."""
        d = ctk.CTkToplevel(self)
        d.title("Custom Label")
        d.geometry("280x120")
        d.transient(self)
        d.attributes('-topmost', True)
        d.after(100, d.grab_set)
        ctk.CTkLabel(d, text="Label text:", font=("", 12)).pack(
            pady=(10, 2), padx=12, anchor="w")
        entry = ctk.CTkEntry(d, height=28)
        entry.pack(padx=12, fill="x")
        entry.focus_set()

        def _ok(_event=None):
            txt = entry.get().strip()
            if txt:
                self._place_atom_label(atom, txt)
            d.destroy()

        entry.bind("<Return>", _ok)
        ctk.CTkButton(d, text="OK", height=28,
                      command=_ok).pack(pady=(8, 4), padx=12, fill="x")

    def _on_label_size_change(self, val):
        """Update font size on all existing labels."""
        self._label_font_size = int(val)
        for actor, _coord, _text in self._atom_labels:
            actor.GetTextProperty().SetFontSize(self._label_font_size)
        self.vtk_frame.render()

    def _on_label_color_pick(self):
        """Open color chooser for label color."""
        cur = tuple(int(c * 255) for c in self._label_color)
        result = colorchooser.askcolor(
            initialcolor=f"#{cur[0]:02x}{cur[1]:02x}{cur[2]:02x}",
            title="Label Color")
        if result[0] is None:
            return
        r, g, b = result[0]
        self._label_color = (r / 255.0, g / 255.0, b / 255.0)
        hex_color = result[1]
        self._label_color_btn.configure(fg_color=hex_color,
                                         hover_color=hex_color)
        for actor, _coord, _text in self._atom_labels:
            actor.GetTextProperty().SetColor(*self._label_color)
        self.vtk_frame.render()

    def _on_meas_size_change(self, val):
        """Update font size on all existing measurement labels."""
        self._measure_font_size = int(val)
        for _lines, txt, _desc in self._measure_actors:
            txt.GetTextProperty().SetFontSize(self._measure_font_size)
        self.vtk_frame.render()

    def _on_meas_color_pick(self):
        """Open color chooser for measurement label color."""
        cur = tuple(int(c * 255) for c in self._measure_label_color)
        result = colorchooser.askcolor(
            initialcolor=f"#{cur[0]:02x}{cur[1]:02x}{cur[2]:02x}",
            title="Measurement Label Color")
        if result[0] is None:
            return
        r, g, b = result[0]
        self._measure_label_color = (r / 255.0, g / 255.0, b / 255.0)
        hex_color = result[1]
        self._meas_color_btn.configure(fg_color=hex_color,
                                        hover_color=hex_color)
        for _line_actors, txt, _desc in self._measure_actors:
            txt.GetTextProperty().SetColor(*self._measure_label_color)
        self.vtk_frame.render()

    def _on_meas_line_color_pick(self):
        """Open color chooser for measurement line color."""
        cur = tuple(int(c * 255) for c in self._measure_line_color)
        result = colorchooser.askcolor(
            initialcolor=f"#{cur[0]:02x}{cur[1]:02x}{cur[2]:02x}",
            title="Measurement Line Color")
        if result[0] is None:
            return
        r, g, b = result[0]
        self._measure_line_color = (r / 255.0, g / 255.0, b / 255.0)
        hex_color = result[1]
        self._meas_line_color_btn.configure(fg_color=hex_color,
                                             hover_color=hex_color)
        for line_actors, _txt, _desc in self._measure_actors:
            for la in line_actors:
                la.GetProperty().SetColor(*self._measure_line_color)
        self.vtk_frame.render()

    def _on_meas_width_change(self, val):
        """Update line width on all existing measurement lines."""
        self._measure_line_width = float(val)
        for line_actors, _txt, _desc in self._measure_actors:
            for la in line_actors:
                la.GetProperty().SetLineWidth(self._measure_line_width)
        self.vtk_frame.render()

    def _remove_label(self, idx):
        """Remove a single atom label by index."""
        if 0 <= idx < len(self._atom_labels):
            actor, _coord, _text = self._atom_labels.pop(idx)
            self.vtk_frame._label_renderer.RemoveActor(actor)
            self.vtk_frame.render()

    def _remove_all_labels(self):
        """Remove every atom label actor."""
        for actor, _coord, _text in self._atom_labels:
            self.vtk_frame._label_renderer.RemoveActor(actor)
        self._atom_labels.clear()
        self.vtk_frame.render()

    # ------------------------------------------------------------------
    # Measurement tools (distance / angle / dihedral)
    # ------------------------------------------------------------------

    def _on_measure_selected(self, choice: str):
        """Called when user selects from the Measure dropdown."""
        self._start_measure(choice.lower())

    def _start_measure(self, mode: str):
        """Enter measurement picking mode: 'distance', 'angle', or 'dihedral'."""
        self._measure_mode = mode
        self._measure_picks.clear()
        self._clear_measure_highlights()
        need = {'distance': 2, 'angle': 3, 'dihedral': 4}[mode]
        self.info_label.configure(
            text=f"Measurement: click {need} atoms ({mode})")
        self._show_measure_border(True)

    def _cancel_measure(self):
        self._clear_measure_highlights()
        self._show_measure_border(False)
        self._measure_mode = None
        self._measure_picks.clear()
        self._measure_var.set("Measure")
        if self.structure:
            s = self.structure
            title = s.title[:60] if s.title else ''
            self.info_label.configure(
                text=f"{title}\n{len(s.atoms)} atoms  {len(s.residues)} residues\n"
                     f"{len(s.chains)} chain(s)  {len(s.bonds)} bonds")

    def _show_measure_border(self, show: bool):
        """Enable or disable the yellow border overlay."""
        self._measure_border_active = show
        # Redraw immediately via the post-render path
        self._draw_measure_border()

    def _draw_measure_border(self):
        """Actually draw/remove the border rectangle on the canvas."""
        canvas = self.vtk_frame.canvas
        canvas.delete("measure_border")
        if self._measure_border_active:
            w = canvas.winfo_width()
            h = canvas.winfo_height()
            pad = 3
            canvas.create_rectangle(
                pad, pad, w - pad, h - pad,
                outline="#FFD700", width=3, tags="measure_border")

    def _post_render(self):
        """Called after every VTK render blit to re-draw canvas overlays."""
        self._draw_measure_border()

    def _on_vtk_resize(self, w, h):
        """Redraw measure border when canvas is resized."""
        self._draw_measure_border()

    def _get_atom_display_radius(self, atom):
        """Return the visual radius of *atom* based on its current representation."""
        vdw_r = VDW_RADII.get(atom.element, VDW_RADII['DEFAULT'])
        cpk_fallback = vdw_r * 0.3  # ball_stick default ball_scale
        if not self.structure:
            return cpk_fallback
        # Find atom index
        atom_idx = None
        for i, a in enumerate(self.structure.atoms):
            if a is atom:
                atom_idx = i
                break
        if atom_idx is None:
            return cpk_fallback
        # Check visible selections containing this atom
        for sel in self.selections:
            if not sel.visible or atom_idx not in sel.atom_indices:
                continue
            rep = sel.representation
            if rep == 'vdw':
                return vdw_r * sel.atom_scale
            elif rep == 'ball_stick':
                return vdw_r * sel.ball_scale
            elif rep == 'sticks':
                return sel.stick_radius
            else:
                # cartoon, tube_ss, backbone, surface – atom not shown as sphere
                return cpk_fallback
        return cpk_fallback

    def _add_measure_highlight(self, atom):
        """Add a yellow highlight sphere around *atom* sized to its representation."""
        from vtkmodules.vtkFiltersSources import vtkSphereSource
        from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper
        radius = self._get_atom_display_radius(atom) * 1.10
        src = vtkSphereSource()
        src.SetCenter(*atom.coord)
        src.SetRadius(radius)
        src.SetPhiResolution(16)
        src.SetThetaResolution(16)
        src.Update()
        mapper = vtkPolyDataMapper()
        mapper.SetInputData(src.GetOutput())
        actor = vtkActor()
        actor.SetMapper(mapper)
        prop = actor.GetProperty()
        prop.SetColor(1.0, 0.85, 0.0)
        prop.SetOpacity(0.30)
        prop.SetAmbient(0.6)
        prop.SetDiffuse(0.4)
        prop.SetSpecular(0.0)
        prop.LightingOn()
        self.vtk_frame.renderer.AddActor(actor)
        self._measure_highlight_actors.append(actor)
        self.vtk_frame.render()

    def _clear_measure_highlights(self):
        """Remove all temporary measurement highlight spheres."""
        for actor in self._measure_highlight_actors:
            self.vtk_frame.renderer.RemoveActor(actor)
        self._measure_highlight_actors.clear()

    def _handle_measure_pick(self, atom):
        """Accumulate picks and compute once enough atoms are selected."""
        self._measure_picks.append(atom)
        self._add_measure_highlight(atom)
        need = {'distance': 2, 'angle': 3, 'dihedral': 4}[self._measure_mode]
        remaining = need - len(self._measure_picks)
        if remaining > 0:
            names = " -> ".join(
                f"{a.chain_id}:{a.res_name}{a.res_id}:{a.name}"
                for a in self._measure_picks)
            self.info_label.configure(
                text=f"{names}\nPick {remaining} more atom(s)")
            return
        coords = [np.array(a.coord) for a in self._measure_picks]
        mode = self._measure_mode
        if mode == 'distance':
            value = np.linalg.norm(coords[1] - coords[0])
            label = f"{value:.2f} Å"
        elif mode == 'angle':
            v1 = coords[0] - coords[1]
            v2 = coords[2] - coords[1]
            cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
            value = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))
            label = f"{value:.1f}°"
        else:  # dihedral
            b1 = coords[1] - coords[0]
            b2 = coords[2] - coords[1]
            b3 = coords[3] - coords[2]
            n1 = np.cross(b1, b2)
            n2 = np.cross(b2, b3)
            n1_norm = np.linalg.norm(n1) + 1e-12
            n2_norm = np.linalg.norm(n2) + 1e-12
            n1 = n1 / n1_norm
            n2 = n2 / n2_norm
            m1 = np.cross(n1, b2 / (np.linalg.norm(b2) + 1e-12))
            x = np.dot(n1, n2)
            y = np.dot(m1, n2)
            value = np.degrees(np.arctan2(y, x))
            label = f"{value:.1f}°"

        names = " -> ".join(
            f"{a.chain_id}:{a.res_name}{a.res_id}:{a.name}"
            for a in self._measure_picks)
        desc = f"{mode}: {names} = {label}"
        self._draw_measurement(coords, label, mode, desc)
        self.info_label.configure(text=f"{names}\n{label}")
        # Deactivate measurement mode
        self._clear_measure_highlights()
        self._show_measure_border(False)
        self._measure_mode = None
        self._measure_picks.clear()
        self._measure_var.set("Measure")

    @staticmethod
    def _make_dashed_polydata(p1, p2, dash_len=0.3, gap_len=0.20):
        """Return a vtkPolyData with dashed-line segments between *p1* and *p2*."""
        from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkCellArray
        from vtkmodules.vtkCommonCore import vtkPoints
        direction = p2 - p1
        total_len = np.linalg.norm(direction)
        if total_len < 1e-8:
            return None
        unit = direction / total_len
        pts = vtkPoints()
        cells = vtkCellArray()
        t = 0.0
        while t < total_len:
            t_end = min(t + dash_len, total_len)
            id0 = pts.InsertNextPoint(*(p1 + unit * t))
            id1 = pts.InsertNextPoint(*(p1 + unit * t_end))
            cells.InsertNextCell(2)
            cells.InsertCellPoint(id0)
            cells.InsertCellPoint(id1)
            t = t_end + gap_len
        pd = vtkPolyData()
        pd.SetPoints(pts)
        pd.SetLines(cells)
        return pd

    def _draw_measurement(self, coords, label, mode, desc=""):
        """Draw dashed lines between picked atoms and a text label.

        Lines are added to the label overlay renderer so they are
        always visible on top of the geometry.
        Each measurement is stored as a single entry:
        ``(line_actors_list, text_actor, description_string)``.
        """
        from vtkmodules.vtkRenderingCore import (
            vtkActor, vtkPolyDataMapper, vtkBillboardTextActor3D,
        )

        line_actors = []
        for i in range(len(coords) - 1):
            pd = self._make_dashed_polydata(coords[i], coords[i + 1])
            if pd is None:
                continue
            mapper = vtkPolyDataMapper()
            mapper.SetInputData(pd)
            line_actor = vtkActor()
            line_actor.SetMapper(mapper)
            line_actor.GetProperty().SetColor(*self._measure_line_color)
            line_actor.GetProperty().SetLineWidth(self._measure_line_width)
            self.vtk_frame._label_renderer.AddActor(line_actor)
            line_actors.append(line_actor)

        mid = sum(coords) / len(coords)
        txt = vtkBillboardTextActor3D()
        txt.SetInput(label)
        txt.SetPosition(*mid)
        tp = txt.GetTextProperty()
        tp.SetFontSize(self._measure_font_size)
        tp.SetColor(*self._measure_label_color)
        tp.SetBold(True)
        tp.SetFontFamilyToCourier()
        txt.SetDisplayOffset(0, 14)
        # Add to label overlay renderer (always on top)
        self.vtk_frame._label_renderer.AddActor(txt)

        # Store as one grouped entry: (lines, text_actor, description)
        self._measure_actors.append((line_actors, txt, desc))
        self.vtk_frame.render()

    def _remove_measurement(self, idx):
        """Remove a single measurement by index."""
        if 0 <= idx < len(self._measure_actors):
            line_actors, txt, _desc = self._measure_actors.pop(idx)
            for la in line_actors:
                self.vtk_frame._label_renderer.RemoveActor(la)
            self.vtk_frame._label_renderer.RemoveActor(txt)
            self.vtk_frame.render()

    def _remove_all_measurements(self):
        """Remove all measurement actors."""
        for line_actors, txt, _desc in self._measure_actors:
            for la in line_actors:
                self.vtk_frame._label_renderer.RemoveActor(la)
            self.vtk_frame._label_renderer.RemoveActor(txt)
        self._measure_actors.clear()
        self._cancel_measure()
        self.vtk_frame.render()

    # ------------------------------------------------------------------
    # Auto-detect molecules
    # ------------------------------------------------------------------

    def _auto_detect_molecules(self):
        if not self.structure:
            messagebox.showwarning("Warning", "Load a structure first")
            return
        res_groups: Dict[str, list] = defaultdict(list)
        for i, a in enumerate(self.structure.atoms):
            if a.res_name in AA_NAMES:
                res_groups['Protein'].append(i)
            elif a.res_name in ('HOH', 'WAT', 'TIP'):
                res_groups['Water'].append(i)
            else:
                res_groups[a.res_name].append(i)
        self.selections.clear()
        color_idx = 0
        for name, indices in res_groups.items():
            if name == 'Protein':
                sel = Selection(name, indices, representation='tube_ss',
                                color_scheme='ss', criteria='Protein')
            elif name == 'Water':
                sel = Selection(name, indices, representation='vdw',
                                color_scheme='element', visible=False,
                                criteria='Water')
            else:
                c = CHAIN_PALETTE[color_idx % len(CHAIN_PALETTE)]
                color_idx += 1
                sel = Selection(name, indices, representation='vdw',
                                color_scheme='element', carbon_color=c,
                                criteria='MDAnalysis expression...',
                                criteria_extra=f'resname {name}')
            self.selections.append(sel)
        self._refresh_sel_ui()
        self._rebuild_full()

    # ------------------------------------------------------------------
    # Selection management
    # ------------------------------------------------------------------

    def _new_selection(self):
        if not self.structure:
            messagebox.showwarning("Warning", "Load a structure first")
            return
        dlg = SelectionDialog(self, self.structure, self.selections)
        self.wait_window(dlg)
        if dlg.result:
            self.selections.append(dlg.result)
            self._safe_refresh_rebuild()

    def _refresh_sel_ui(self):
        self._drag_data = {'idx': None, 'start_y': 0, 'frames': []}
        for w in self.sel_list_frame.winfo_children():
            w.destroy()
        frames = []
        for idx, sel in enumerate(self.selections):
            f = ctk.CTkFrame(self.sel_list_frame, fg_color="gray20",
                             corner_radius=8, border_width=1,
                             border_color="gray45")
            f.pack(fill="x", pady=(0, 5), padx=4)
            f._sel_idx = idx
            frames.append(f)

            # Row 1: checkbox | name | edit | [X]
            row1 = ctk.CTkFrame(f, fg_color="transparent")
            row1.pack(fill="x", padx=4, pady=(4, 0))
            vv = ctk.BooleanVar(value=sel.visible)
            ctk.CTkCheckBox(row1, text="", variable=vv, width=24,
                            command=lambda s=sel, v=vv: self._toggle(s, v)
                            ).pack(side="left")
            name_lbl = ctk.CTkLabel(row1, text=sel.name, font=("", 13, "bold"),
                                    cursor="hand2")
            name_lbl.pack(side="left", padx=4)
            name_lbl.bind("<ButtonPress-1>",
                          lambda e, i=idx: self._drag_start(e, i))
            name_lbl.bind("<B1-Motion>", self._drag_motion)
            name_lbl.bind("<ButtonRelease-1>", self._drag_end)
            name_lbl.bind("<Double-Button-1>",
                          lambda e, s=sel: self._dbl_click_center(s))
            if idx > 0:
                ctk.CTkButton(row1, text="x", width=24, height=22,
                              fg_color="#8a3a3a", hover_color="#cc3333",
                              font=("", 11),
                              command=lambda s=sel: self._del_sel(s)
                              ).pack(side="right")
            ctk.CTkButton(row1, text="Edit", width=36, height=22,
                          fg_color="gray35", font=("", 11),
                          command=lambda s=sel: self._edit_selection(s)
                          ).pack(side="right", padx=(0, 2))

            # Row 2: rep dropdown | color swatch | gear | atom count
            row2 = ctk.CTkFrame(f, fg_color="transparent")
            row2.pack(fill="x", padx=4, pady=(2, 4))

            rep_label = REP_LABELS.get(sel.representation, 'VDW')

            def _on_rep(val, s=sel):
                k = REP_KEYS[REP_VALUES.index(val)]
                s.representation = k
                self._safe_refresh_rebuild()

            om = ctk.CTkOptionMenu(row2, values=REP_VALUES, width=100,
                                   height=24, font=("", 11), command=_on_rep,
                                   dynamic_resizing=False)
            om.set(rep_label)
            om.pack(side="left", padx=(0, 3))

            swatch = self._make_swatch(row2, sel)
            swatch.pack(side="left", padx=(0, 3))

            if not hasattr(self, '_gear_img'):
                self._gear_img = _make_gear_image(size=16, color="white")
            ctk.CTkButton(row2, text="", image=self._gear_img,
                          width=24, height=24,
                          fg_color="gray35",
                          command=lambda s=sel: self._sel_settings(s)
                          ).pack(side="left", padx=(0, 3))

            ctk.CTkLabel(row2, text=f"{len(sel.atom_indices)} atoms",
                         font=("", 10), text_color="gray55"
                         ).pack(side="right", padx=2)

        self._drag_data['frames'] = frames

    def _sel_swatch_color(self, sel):
        if sel.uniform_color:
            return "#%02X%02X%02X" % sel.uniform_color
        if sel.carbon_color:
            return "#%02X%02X%02X" % sel.carbon_color
        cs_colors = {'element': '#909090', 'chain': '#E6194B',
                     'ss': '#FF0000', 'uniform': '#909090',
                     'residue_nature': '#3CB44B'}
        return cs_colors.get(sel.color_scheme, '#909090')

    # Swatch stripe palettes for multi-color schemes
    _SWATCH_STRIPES = {
        'element': [(144, 144, 144), (255, 13, 13), (48, 80, 248), (255, 255, 255)],
        'chain':   [(230, 25, 75), (60, 180, 75), (255, 225, 25), (0, 130, 200)],
        'ss':      [(180, 141, 218), (33, 150, 166), (232, 232, 232), (181, 213, 200)],
        'residue_nature': [(220, 60, 60), (70, 100, 220), (60, 180, 75), (230, 200, 50)],
    }

    def _make_swatch(self, parent, sel):
        """Create a 24x24 canvas swatch with vertical stripes or solid color."""
        cs = sel.color_scheme
        if cs == 'uniform' and sel.uniform_color:
            colors = [sel.uniform_color] * 4
        elif sel.carbon_color and cs == 'element':
            colors = [sel.carbon_color, (255, 13, 13),
                      (48, 80, 248), (255, 255, 255)]
        else:
            colors = self._SWATCH_STRIPES.get(
                cs, self._SWATCH_STRIPES['element'])
        size = 24
        sw = tk.Canvas(parent, width=size, height=size,
                       highlightthickness=1, highlightbackground="#808080",
                       bd=0, cursor="hand2")
        stripe_w = size // len(colors)
        for i, c in enumerate(colors):
            hx = "#%02X%02X%02X" % tuple(c)
            x0 = i * stripe_w
            x1 = (i + 1) * stripe_w if i < len(colors) - 1 else size
            sw.create_rectangle(x0, 0, x1, size, fill=hx, outline='')
        sw.bind("<Button-1>", lambda e, s=sel: self._chg_cs(s))
        return sw

    def _dbl_click_center(self, sel):
        self._drag_data['idx'] = None
        self._drag_data['_dragged'] = False
        self._center_on_selection(sel)

    def _drag_start(self, event, idx):
        self._drag_data['idx'] = idx
        self._drag_data['start_y'] = event.y_root
        self._drag_data['_dragged'] = False
        f = self._drag_data['frames'][idx]
        f.configure(border_color="#5599ff")

    def _drag_motion(self, event):
        drag = self._drag_data
        if drag['idx'] is None:
            return
        idx = drag['idx']
        frames = drag['frames']
        f = frames[idx]
        f_h = f.winfo_height() + 6
        dy = event.y_root - drag['start_y']
        if dy > f_h * 0.5 and idx < len(self.selections) - 1:
            drag['_dragged'] = True
            self.selections[idx], self.selections[idx + 1] = \
                self.selections[idx + 1], self.selections[idx]
            drag['idx'] = idx + 1
            drag['start_y'] = event.y_root
            self._drag_repack(drag['idx'])
        elif dy < -f_h * 0.5 and idx > 0:
            drag['_dragged'] = True
            self.selections[idx], self.selections[idx - 1] = \
                self.selections[idx - 1], self.selections[idx]
            drag['idx'] = idx - 1
            drag['start_y'] = event.y_root
            self._drag_repack(drag['idx'])

    def _drag_repack(self, highlight_idx):
        self._drag_data['idx'] = None
        self._refresh_sel_ui()
        if highlight_idx < len(self._drag_data['frames']):
            self._drag_data['frames'][highlight_idx].configure(
                border_color="#5599ff")
            self._drag_data['idx'] = highlight_idx

    def _drag_end(self, event):
        drag = self._drag_data
        if drag['idx'] is not None and drag['idx'] < len(drag['frames']):
            drag['frames'][drag['idx']].configure(border_color="gray45")
        dragged = drag.get('_dragged', False)
        drag['idx'] = None
        if dragged:
            self._rebuild()

    def _toggle(self, sel, var):
        sel.visible = var.get()
        self._rebuild()

    def _del_sel(self, sel):
        self.selections.remove(sel)
        self._safe_refresh_rebuild()

    # ------------------------------------------------------------------
    # Edit selection
    # ------------------------------------------------------------------

    def _edit_selection(self, sel):
        if not self.structure:
            return
        d = ctk.CTkToplevel(self)
        d.title(f"Edit: {sel.name}")
        d.geometry("360x420")
        d.transient(self)
        d.attributes('-topmost', True)

        ctk.CTkLabel(d, text="Name:", font=("", 12)).pack(
            pady=(10, 1), padx=12, anchor="w")
        name_entry = ctk.CTkEntry(d, height=28)
        name_entry.pack(padx=12, fill="x")
        name_entry.insert(0, sel.name)

        # Show how the selection was created
        if sel.criteria:
            info = sel.criteria
            if sel.criteria_extra:
                info += f"  ({sel.criteria_extra})"
            ctk.CTkLabel(d, text=f"Current: {info}",
                         font=("", 11), text_color="gray60").pack(
                pady=(6, 0), padx=12, anchor="w")

        ctk.CTkLabel(d, text="Re-select atoms:", font=("", 12)).pack(
            pady=(6, 1), padx=12, anchor="w")
        crit_menu = ctk.CTkOptionMenu(d,
                                      values=['(keep current)'] + SELECTION_CRITERIA,
                                      height=28)
        crit_menu.set('(keep current)')
        crit_menu.pack(padx=12, fill="x")

        extra_frame = ctk.CTkFrame(d, fg_color="transparent")
        _extra = [ctk.CTkEntry(extra_frame)]

        def _show_extra(val, prefill=''):
            for w in extra_frame.winfo_children():
                w.destroy()
            extra_frame.pack_forget()
            if val == 'Chain...':
                chains = ', '.join(sorted(self.structure.chains.keys()))
                ctk.CTkLabel(extra_frame,
                             text=f"Available: {chains}").pack(anchor="w", pady=2)
                _extra[0] = ctk.CTkEntry(extra_frame,
                                         placeholder_text="Chain IDs (e.g. A B)")
                _extra[0].pack(fill="x", pady=2)
                if prefill:
                    _extra[0].insert(0, prefill)
                extra_frame.pack(fill="x", padx=12, after=crit_menu)
            elif val == 'Residue range...':
                _extra[0] = ctk.CTkEntry(
                    extra_frame,
                    placeholder_text="e.g. A:50 or A:10-50,B:1-20")
                _extra[0].pack(fill="x", pady=2)
                if prefill:
                    _extra[0].insert(0, prefill)
                extra_frame.pack(fill="x", padx=12, after=crit_menu)
            elif val == 'MDAnalysis expression...':
                hdr = ctk.CTkFrame(extra_frame, fg_color="transparent")
                hdr.pack(fill="x")
                ctk.CTkLabel(hdr, text="MDAnalysis selection:").pack(
                    side="left", anchor="w")
                ctk.CTkButton(
                    hdr, text="?", width=28, height=24,
                    fg_color="#5a5a8a", hover_color="#7a7aaa",
                    command=lambda: _show_selection_help(d)
                ).pack(side="right")
                _extra[0] = ctk.CTkEntry(
                    extra_frame,
                    placeholder_text="e.g. protein and name CA")
                _extra[0].pack(fill="x", pady=2)
                if prefill:
                    _extra[0].insert(0, prefill)
                extra_frame.pack(fill="x", padx=12, after=crit_menu)

        def _on_crit(val):
            _show_extra(val)

        crit_menu.configure(command=_on_crit)

        # Pre-fill with the stored criteria if available
        if sel.criteria and sel.criteria in SELECTION_CRITERIA:
            crit_menu.set(sel.criteria)
            _show_extra(sel.criteria, sel.criteria_extra)

        def _apply():
            new_name = name_entry.get().strip()
            if new_name:
                sel.name = new_name
            crit = crit_menu.get()
            if crit != '(keep current)':
                extra_text = _extra[0].get().strip() \
                    if crit in ('Chain...', 'Residue range...',
                                'MDAnalysis expression...') else ''
                indices = self._resolve_sel_criteria(crit, extra_text)
                if indices:
                    sel.atom_indices = indices
                    sel.criteria = crit
                    sel.criteria_extra = extra_text
                else:
                    messagebox.showwarning("Warning", "No atoms matched",
                                           parent=d)
                    return
            d.destroy()
            self._safe_refresh_rebuild()

        ctk.CTkButton(d, text="Apply", height=30, fg_color="#2a6e2a",
                      hover_color="#358535", command=_apply
                      ).pack(pady=(12, 4), padx=12, fill="x")
        ctk.CTkButton(d, text="Cancel", height=28, fg_color="gray40",
                      command=d.destroy).pack(padx=12, fill="x", pady=(0, 8))

    def _resolve_sel_criteria(self, crit, extra=''):
        atoms = self.structure.atoms
        if crit == 'All':
            return list(range(len(atoms)))
        if crit == 'Protein':
            return [i for i, a in enumerate(atoms) if a.res_name in AA_NAMES]
        if crit == 'Backbone':
            return [i for i, a in enumerate(atoms)
                    if a.res_name in AA_NAMES and a.name in BACKBONE_NAMES]
        if crit == 'Sidechain':
            return [i for i, a in enumerate(atoms)
                    if a.res_name in AA_NAMES and a.name not in BACKBONE_NAMES]
        if crit == 'Water':
            return [i for i, a in enumerate(atoms) if a.res_name in ('HOH', 'WAT', 'TIP')]
        if crit == 'Ligand':
            return [i for i, a in enumerate(atoms)
                    if a.res_name not in AA_NAMES
                    and a.res_name not in ('HOH', 'WAT', 'TIP')]
        if crit == 'Chain...':
            chains = {c.strip().upper() for c in extra.replace(',', ' ').split() if c.strip()}
            return [i for i, a in enumerate(atoms) if a.chain_id in chains]
        if crit == 'Residue range...':
            return self._parse_range_text(extra)
        if crit == 'Around selection...':
            return []
        if crit == 'MDAnalysis expression...':
            return _resolve_mda_expression(self.structure, extra, parent=self)
        return []

    def _parse_range_text(self, text):
        indices = []
        for part in text.split(','):
            part = part.strip()
            if ':' not in part:
                continue
            ch, rng = part.split(':', 1)
            ch = ch.strip().upper()
            try:
                if '-' in rng:
                    lo, hi = int(rng.split('-')[0]), int(rng.split('-')[1])
                else:
                    lo = hi = int(rng.strip())
            except ValueError:
                continue
            for i, a in enumerate(self.structure.atoms):
                if a.chain_id == ch and lo <= a.res_id <= hi:
                    indices.append(i)
        return indices

    # ------------------------------------------------------------------
    # Color scheme dialog
    # ------------------------------------------------------------------

    def _chg_cs(self, sel):
        d = ctk.CTkToplevel(self)
        d.title("Color")
        d.geometry("260x360")
        d.transient(self)
        d.attributes('-topmost', True)
        d.after(100, d.grab_set)
        ctk.CTkLabel(d, text="Color Scheme",
                     font=("", 13, "bold")).pack(pady=(8, 4), padx=12)
        cs = {'Element (CPK)': 'element', 'Chain': 'chain',
              'Secondary Structure': 'ss', 'Residue Nature': 'residue_nature'}
        for label, key in cs.items():
            ctk.CTkButton(d, text=label, height=30,
                          command=lambda k=key, dd=d: (
                              setattr(sel, 'color_scheme', k),
                              setattr(sel, 'uniform_color', None),
                              dd.destroy(),
                              self._safe_refresh_rebuild())
                          ).pack(pady=2, padx=12, fill="x")
        _sep(d)

        def _pick_uniform(dd=d):
            dd.destroy()
            init = ("#%02X%02X%02X" % sel.uniform_color
                    if sel.uniform_color else None)
            c = colorchooser.askcolor(initialcolor=init,
                                     title="Pick uniform color")
            if c[0]:
                sel.uniform_color = (int(c[0][0]), int(c[0][1]), int(c[0][2]))
                sel.color_scheme = 'uniform'
                self._safe_refresh_rebuild()

        ctk.CTkButton(d, text="Pick Uniform Color...", height=30,
                      fg_color="#5a3a6e", hover_color="#7a4a9e",
                      command=_pick_uniform).pack(pady=2, padx=12, fill="x")

        def _pick_carbon(dd=d):
            dd.destroy()
            init = ("#%02X%02X%02X" % sel.carbon_color
                    if sel.carbon_color else None)
            c = colorchooser.askcolor(initialcolor=init,
                                     title="Pick carbon color")
            if c[0]:
                sel.carbon_color = (int(c[0][0]), int(c[0][1]), int(c[0][2]))
                sel.uniform_color = None
                if sel.color_scheme == 'uniform':
                    sel.color_scheme = 'element'
                self._safe_refresh_rebuild()

        ctk.CTkButton(d, text="Pick Carbon Color...", height=30,
                      fg_color="#3a5a6e", hover_color="#4a7a9e",
                      command=_pick_carbon).pack(pady=2, padx=12, fill="x")
        if sel.carbon_color:
            ctk.CTkButton(d, text="Reset Carbon Color", height=28,
                          fg_color="gray40",
                          command=lambda dd=d: (
                              setattr(sel, 'carbon_color', None),
                              dd.destroy(),
                              self._safe_refresh_rebuild())
                          ).pack(pady=2, padx=12, fill="x")

    # ------------------------------------------------------------------
    # Per-selection settings (quality, sizes, material, SS colors)
    # ------------------------------------------------------------------

    def _sel_settings(self, sel):
        rep = sel.representation
        _attrs = ('quality', 'atom_scale', 'bond_radius', 'ball_scale',
                  'stick_radius', 'backbone_radius', 'helix_width',
                  'sheet_width', 'coil_width', 'opacity',
                  'surface_resolution', 'surface_radius',
                  'ambient', 'diffuse', 'specular', 'specular_power')
        _orig = {a: getattr(sel, a) for a in _attrs}
        _orig['ss_colors'] = dict(sel.ss_colors)
        _preview_pending = [False]

        def _schedule_preview(*_):
            if not _preview_pending[0]:
                _preview_pending[0] = True
                d.after(80, _do_preview)

        def _do_preview():
            _preview_pending[0] = False
            _write_sel()
            self._rebuild()

        def _write_sel():
            sel.quality = int(q_sl.get())
            for attr, sl in size_sliders.items():
                val = sl.get()
                if attr in ('surface_resolution',):
                    val = int(val)
                else:
                    val = round(val, 3)
                setattr(sel, attr, val)
            sel.ambient = round(amb_sl.get(), 2)
            sel.diffuse = round(dif_sl.get(), 2)
            sel.specular = round(spc_sl.get(), 2)
            sel.specular_power = round(spp_sl.get(), 1)

        size_h = {'vdw': 55, 'ball_stick': 110, 'sticks': 55,
                  'cartoon': 165, 'tube_ss': 165, 'backbone': 55,
                  'surface': 165}.get(rep, 0)
        ss_h = 210 if rep in ('cartoon', 'tube_ss', 'backbone') else 0
        h = 240 + size_h + ss_h + 260

        d = ctk.CTkToplevel(self)
        d.title(f"Settings: {sel.name}")
        d.geometry(f"340x{h}")
        d.transient(self)
        d.attributes('-topmost', True)
        d.after(100, d.grab_set)

        def _mk(txt, lo, hi, steps, val, fmt=".2f"):
            lbl = ctk.CTkLabel(d, text=f"{txt}: {val:{fmt}}", font=("", 11))
            lbl.pack(pady=(4, 0), padx=12, anchor="w")
            sl = ctk.CTkSlider(d, from_=lo, to=hi, number_of_steps=steps)
            sl.set(val)

            def _on_change(v, l=lbl, t=txt, f=fmt):
                l.configure(text=f"{t}: {v:{f}}")
                _schedule_preview()

            sl.configure(command=_on_change)
            sl.pack(padx=12, fill="x")
            return sl, lbl

        # Quality
        q_lbl = ctk.CTkLabel(d, text=f"Quality: {QUALITY_LABELS[sel.quality - 1]}",
                             font=("", 11))
        q_lbl.pack(pady=(4, 0), padx=12, anchor="w")
        q_sl = ctk.CTkSlider(d, from_=1, to=5, number_of_steps=4)
        q_sl.set(sel.quality)

        def _on_q(v):
            q_lbl.configure(text=f"Quality: {QUALITY_LABELS[int(v) - 1]}")
            _schedule_preview()

        q_sl.configure(command=_on_q)
        q_sl.pack(padx=12, fill="x")

        _sep(d)
        ctk.CTkLabel(d, text="Size Parameters",
                     font=("", 12, "bold")).pack(pady=(2, 2), padx=12, anchor="w")
        size_sliders: dict = {}
        if rep == 'vdw':
            size_sliders['atom_scale'] = _mk("Atom Scale", 0.1, 3.0, 29, sel.atom_scale)[0]
        elif rep == 'ball_stick':
            size_sliders['ball_scale'] = _mk("Ball Scale", 0.05, 1.5, 29, sel.ball_scale)[0]
            size_sliders['bond_radius'] = _mk("Bond Radius", 0.02, 0.5, 24, sel.bond_radius)[0]
        elif rep == 'sticks':
            size_sliders['stick_radius'] = _mk("Stick Radius", 0.05, 0.6, 22, sel.stick_radius)[0]
        elif rep in ('cartoon', 'tube_ss'):
            size_sliders['helix_width'] = _mk("Helix Width", 0.5, 5.0, 18, sel.helix_width)[0]
            size_sliders['sheet_width'] = _mk("Sheet Width", 0.5, 5.0, 18, sel.sheet_width)[0]
            size_sliders['coil_width'] = _mk("Coil Width", 0.1, 3.0, 14, sel.coil_width)[0]
        elif rep == 'backbone':
            size_sliders['backbone_radius'] = _mk("Backbone Radius", 0.05, 1.5, 14, sel.backbone_radius)[0]
        elif rep == 'surface':
            size_sliders['opacity'] = _mk("Opacity", 0.1, 1.0, 18, sel.opacity)[0]
            size_sliders['surface_resolution'] = _mk(
                "Resolution", 16, 128, 14,
                float(sel.surface_resolution), ".0f")[0]
            size_sliders['surface_radius'] = _mk(
                "Surface Radius", 0.02, 0.40, 19,
                sel.surface_radius, ".3f")[0]

        # SS Colors
        if rep in ('cartoon', 'tube_ss', 'backbone'):
            _sep(d)
            ctk.CTkLabel(d, text="Secondary Structure Colors",
                         font=("", 12, "bold")).pack(
                             pady=(2, 2), padx=12, anchor="w")
            ss_frame = ctk.CTkFrame(d, fg_color="transparent")
            ss_frame.pack(fill="x", padx=12)
            for ss_key, ss_label in SS_LABELS.items():
                row = ctk.CTkFrame(ss_frame, fg_color="transparent")
                row.pack(fill="x", pady=1)
                ctk.CTkLabel(row, text=ss_label, font=("", 11),
                             width=90).pack(side="left")
                c = sel.ss_colors.get(ss_key,
                                      SS_COLORS.get(ss_key, (200, 200, 200)))
                hex_c = "#%02X%02X%02X" % tuple(c)
                sw = ctk.CTkButton(row, text="", width=28, height=22,
                                   fg_color=hex_c, hover_color=hex_c,
                                   corner_radius=4, border_width=1,
                                   border_color="gray50")

                def _pick_ss(k=ss_key, btn=sw):
                    cur = sel.ss_colors.get(k,
                                            SS_COLORS.get(k, (200, 200, 200)))
                    cur_hex = "#%02X%02X%02X" % tuple(cur)
                    cc = colorchooser.askcolor(
                        initialcolor=cur_hex,
                        title=f"Pick {SS_LABELS[k]} color", parent=d)
                    if cc[0]:
                        rgb = (int(cc[0][0]), int(cc[0][1]), int(cc[0][2]))
                        sel.ss_colors[k] = rgb
                        hx = "#%02X%02X%02X" % rgb
                        btn.configure(fg_color=hx, hover_color=hx)
                        _schedule_preview()

                sw.configure(command=_pick_ss)
                sw.pack(side="left", padx=4)

                def _reset_ss(k=ss_key, btn=sw):
                    sel.ss_colors[k] = SS_COLORS.get(k, (200, 200, 200))
                    hx = "#%02X%02X%02X" % tuple(sel.ss_colors[k])
                    btn.configure(fg_color=hx, hover_color=hx)
                    _schedule_preview()

                ctk.CTkButton(row, text="Reset", width=46, height=20,
                              font=("", 9), fg_color="gray35",
                              command=_reset_ss).pack(side="left", padx=2)

        # Material
        _sep(d)
        ctk.CTkLabel(d, text="Material",
                     font=("", 12, "bold")).pack(pady=(2, 2), padx=12, anchor="w")
        amb_sl, amb_lbl = _mk("Ambient", 0.0, 1.0, 20, sel.ambient)
        dif_sl, dif_lbl = _mk("Diffuse", 0.0, 1.0, 20, sel.diffuse)
        spc_sl, spc_lbl = _mk("Specular", 0.0, 1.0, 20, sel.specular)
        spp_sl, spp_lbl = _mk("Spec Power", 1, 100, 99,
                               sel.specular_power, ".0f")

        pf = ctk.CTkFrame(d, fg_color="transparent")
        pf.pack(fill="x", padx=12, pady=(4, 0))
        for pname, (a, df, sp, spp) in MATERIAL_PRESETS.items():
            def _setp(a=a, df=df, sp=sp, spp=spp):
                amb_sl.set(a)
                dif_sl.set(df)
                spc_sl.set(sp)
                spp_sl.set(spp)
                amb_lbl.configure(text=f"Ambient: {a:.2f}")
                dif_lbl.configure(text=f"Diffuse: {df:.2f}")
                spc_lbl.configure(text=f"Specular: {sp:.2f}")
                spp_lbl.configure(text=f"Spec Power: {spp:.0f}")
                _schedule_preview()

            ctk.CTkButton(pf, text=pname, width=60, height=22,
                          font=("", 10), fg_color="gray35",
                          command=_setp
                          ).pack(side="left", padx=1, expand=True)

        def _save_cam():
            cam = self.vtk_frame.renderer.GetActiveCamera()
            return (cam.GetPosition(), cam.GetFocalPoint(),
                    cam.GetViewUp(), cam.GetClippingRange(),
                    cam.GetParallelScale())

        def _restore_cam(state):
            cam = self.vtk_frame.renderer.GetActiveCamera()
            cam.SetPosition(state[0])
            cam.SetFocalPoint(state[1])
            cam.SetViewUp(state[2])
            cam.SetClippingRange(state[3])
            cam.SetParallelScale(state[4])
            self.vtk_frame.render()

        def _ok():
            _write_sel()
            state = _save_cam()
            self.vtk_frame.lock_camera()
            d.destroy()
            self._refresh_sel_ui()
            self._rebuild()
            _restore_cam(state)
            self.after(200, self.vtk_frame.unlock_camera)

        def _cancel():
            for attr in _attrs:
                setattr(sel, attr, _orig[attr])
            sel.ss_colors = _orig['ss_colors']
            state = _save_cam()
            self.vtk_frame.lock_camera()
            d.destroy()
            self._rebuild()
            _restore_cam(state)
            self.after(200, self.vtk_frame.unlock_camera)

        ctk.CTkButton(d, text="OK", height=32, fg_color="#2a6e2a",
                      hover_color="#358535", command=_ok
                      ).pack(pady=(10, 4), padx=12, fill="x")
        ctk.CTkButton(d, text="Cancel", height=30, fg_color="gray40",
                      command=_cancel).pack(padx=12, fill="x")

    # ------------------------------------------------------------------
    # Edit operations
    # ------------------------------------------------------------------

    def _edit_chain_name(self):
        if not self.structure:
            messagebox.showwarning("Warning", "Load a structure first")
            return
        d = ctk.CTkToplevel(self)
        d.title("Rename Chain")
        d.geometry("320x240")
        d.transient(self)
        d.attributes('-topmost', True)
        d.after(100, d.grab_set)

        chains = sorted(self.structure.chains.keys())
        ctk.CTkLabel(d, text=f"Available chains: {', '.join(chains)}").pack(
            pady=(10, 4), padx=12)
        ctk.CTkLabel(d, text="Old chain ID:").pack(pady=(4, 1), padx=12, anchor="w")
        old_entry = ctk.CTkEntry(d, placeholder_text="e.g. A")
        old_entry.pack(padx=12, fill="x")
        ctk.CTkLabel(d, text="New chain ID:").pack(pady=(4, 1), padx=12, anchor="w")
        new_entry = ctk.CTkEntry(d, placeholder_text="e.g. B")
        new_entry.pack(padx=12, fill="x")

        def _apply():
            old_ch = old_entry.get().strip()
            new_ch = new_entry.get().strip()
            if not old_ch or not new_ch:
                messagebox.showwarning("Warning", "Fill both fields", parent=d)
                return
            if len(new_ch) != 1:
                messagebox.showwarning("Warning", "Chain ID must be 1 character",
                                       parent=d)
                return
            count = 0
            for atom in self.structure.atoms:
                if atom.chain_id == old_ch:
                    atom.chain_id = new_ch
                    count += 1
            for res in self.structure.residues:
                if res.chain_id == old_ch:
                    res.chain_id = new_ch
            self.structure.chains.clear()
            for res in self.structure.residues:
                self.structure.chains.setdefault(res.chain_id, []).append(res)
            d.destroy()
            messagebox.showinfo("Done",
                                f"Renamed {count} atoms from chain {old_ch} to {new_ch}")
            self._safe_refresh_rebuild()

        ctk.CTkButton(d, text="Apply", height=30, fg_color="#2a6e2a",
                      command=_apply).pack(pady=(10, 4), padx=12, fill="x")

    def _edit_residue_name(self):
        if not self.structure:
            messagebox.showwarning("Warning", "Load a structure first")
            return
        d = ctk.CTkToplevel(self)
        d.title("Rename Residues")
        d.geometry("340x280")
        d.transient(self)
        d.attributes('-topmost', True)
        d.after(100, d.grab_set)

        ctk.CTkLabel(d, text="Chain:").pack(pady=(10, 1), padx=12, anchor="w")
        ch_entry = ctk.CTkEntry(d, placeholder_text="e.g. A")
        ch_entry.pack(padx=12, fill="x")
        ctk.CTkLabel(d, text="Residue range (e.g. 10-20, or 10 for single):").pack(
            pady=(4, 1), padx=12, anchor="w")
        rng_entry = ctk.CTkEntry(d, placeholder_text="e.g. 10-20")
        rng_entry.pack(padx=12, fill="x")
        ctk.CTkLabel(d, text="New residue name:").pack(
            pady=(4, 1), padx=12, anchor="w")
        name_entry = ctk.CTkEntry(d, placeholder_text="e.g. LIG")
        name_entry.pack(padx=12, fill="x")

        def _apply():
            ch = ch_entry.get().strip().upper()
            rng = rng_entry.get().strip()
            new_name = name_entry.get().strip().upper()
            if not ch or not rng or not new_name:
                messagebox.showwarning("Warning", "Fill all fields", parent=d)
                return
            try:
                if '-' in rng:
                    lo, hi = int(rng.split('-')[0]), int(rng.split('-')[1])
                else:
                    lo = hi = int(rng)
            except ValueError:
                messagebox.showwarning("Warning", "Invalid range", parent=d)
                return
            count = 0
            for atom in self.structure.atoms:
                if atom.chain_id == ch and lo <= atom.res_id <= hi:
                    atom.res_name = new_name
                    count += 1
            for res in self.structure.residues:
                if res.chain_id == ch and lo <= res.seq_id <= hi:
                    res.name = new_name
            d.destroy()
            messagebox.showinfo("Done", f"Renamed {count} atoms to {new_name}")
            self._safe_refresh_rebuild()

        ctk.CTkButton(d, text="Apply", height=30, fg_color="#2a6e2a",
                      command=_apply).pack(pady=(10, 4), padx=12, fill="x")

    def _edit_residue_numbering(self):
        if not self.structure:
            messagebox.showwarning("Warning", "Load a structure first")
            return
        d = ctk.CTkToplevel(self)
        d.title("Renumber Residues")
        d.geometry("340x280")
        d.transient(self)
        d.attributes('-topmost', True)
        d.after(100, d.grab_set)

        ctk.CTkLabel(d, text="Chain:").pack(pady=(10, 1), padx=12, anchor="w")
        ch_entry = ctk.CTkEntry(d, placeholder_text="e.g. A")
        ch_entry.pack(padx=12, fill="x")
        ctk.CTkLabel(d, text="Current range (e.g. 100-200):").pack(
            pady=(4, 1), padx=12, anchor="w")
        rng_entry = ctk.CTkEntry(d, placeholder_text="e.g. 100-200")
        rng_entry.pack(padx=12, fill="x")
        ctk.CTkLabel(d, text="New start number:").pack(
            pady=(4, 1), padx=12, anchor="w")
        start_entry = ctk.CTkEntry(d, placeholder_text="e.g. 1")
        start_entry.pack(padx=12, fill="x")

        def _apply():
            ch = ch_entry.get().strip().upper()
            rng = rng_entry.get().strip()
            try:
                new_start = int(start_entry.get().strip())
                if '-' in rng:
                    lo, hi = int(rng.split('-')[0]), int(rng.split('-')[1])
                else:
                    lo = hi = int(rng)
            except ValueError:
                messagebox.showwarning("Warning", "Invalid numbers", parent=d)
                return
            old_ids = sorted(set(a.res_id for a in self.structure.atoms
                                 if a.chain_id == ch and lo <= a.res_id <= hi))
            remap = {old: new_start + i for i, old in enumerate(old_ids)}
            count = 0
            for atom in self.structure.atoms:
                if atom.chain_id == ch and atom.res_id in remap:
                    atom.res_id = remap[atom.res_id]
                    count += 1
            for res in self.structure.residues:
                if res.chain_id == ch and res.seq_id in remap:
                    res.seq_id = remap[res.seq_id]
            d.destroy()
            messagebox.showinfo(
                "Done",
                f"Renumbered {count} atoms ({len(remap)} residues)")
            self._safe_refresh_rebuild()

        ctk.CTkButton(d, text="Apply", height=30, fg_color="#2a6e2a",
                      command=_apply).pack(pady=(10, 4), padx=12, fill="x")

    def _delete_selection_atoms(self):
        if not self.structure:
            messagebox.showwarning("Warning", "No structure loaded")
            return
        d = ctk.CTkToplevel(self)
        d.title("Delete Atoms")
        d.geometry("400x380")
        d.transient(self)
        d.attributes('-topmost', True)

        # Track preview highlight actors so we can remove them
        _highlight_actors = []

        def _clear_highlight():
            for a in _highlight_actors:
                self.vtk_frame.remove_actor(a)
            _highlight_actors.clear()
            self.vtk_frame.render()

        ctk.CTkLabel(d, text="Delete Atoms",
                     font=("", 14, "bold")).pack(pady=(10, 4), padx=12)

        # --- Mode selector ---
        mode_var = ctk.StringVar(value="criteria")

        # Container that holds the mode-specific widgets (one at a time)
        body = ctk.CTkFrame(d, fg_color="transparent")
        body.pack(fill="both", padx=12, pady=(0, 4), expand=True)

        # -- Panel A: By criteria --
        panel_crit = ctk.CTkFrame(body, fg_color="transparent")
        crit_menu = ctk.CTkOptionMenu(panel_crit, values=SELECTION_CRITERIA,
                                      height=28)
        crit_menu.pack(fill="x", pady=(4, 2))
        crit_extra_frame = ctk.CTkFrame(panel_crit, fg_color="transparent")
        _crit_extra = [ctk.CTkEntry(crit_extra_frame)]

        def _on_del_crit(val):
            for w in crit_extra_frame.winfo_children():
                w.destroy()
            crit_extra_frame.pack_forget()
            if val == 'Chain...':
                chains = ', '.join(sorted(self.structure.chains.keys()))
                ctk.CTkLabel(crit_extra_frame,
                             text=f"Available: {chains}").pack(
                    anchor="w", pady=2)
                _crit_extra[0] = ctk.CTkEntry(
                    crit_extra_frame, placeholder_text="Chain IDs (e.g. A B)")
                _crit_extra[0].pack(fill="x", pady=2)
                crit_extra_frame.pack(fill="x")
            elif val == 'Residue range...':
                _crit_extra[0] = ctk.CTkEntry(
                    crit_extra_frame,
                    placeholder_text="e.g. A:50 or A:10-50,B:1-20")
                _crit_extra[0].pack(fill="x", pady=2)
                crit_extra_frame.pack(fill="x")
            elif val == 'MDAnalysis expression...':
                hdr = ctk.CTkFrame(crit_extra_frame, fg_color="transparent")
                hdr.pack(fill="x")
                ctk.CTkLabel(hdr, text="Expression:").pack(
                    side="left", anchor="w")
                ctk.CTkButton(
                    hdr, text="?", width=28, height=24,
                    fg_color="#5a5a8a", hover_color="#7a7aaa",
                    command=lambda: _show_selection_help(d)
                ).pack(side="right")
                _crit_extra[0] = ctk.CTkEntry(
                    crit_extra_frame,
                    placeholder_text="e.g. resname LIG and name CA")
                _crit_extra[0].pack(fill="x", pady=2)
                crit_extra_frame.pack(fill="x")

        crit_menu.configure(command=_on_del_crit)

        # -- Panel B: Existing selection --
        panel_exist = ctk.CTkFrame(body, fg_color="transparent")
        names = ([s.name for s in self.selections[1:]]
                 if len(self.selections) > 1 else ["(no selections)"])
        exist_menu = ctk.CTkOptionMenu(panel_exist, values=names, height=28)
        exist_menu.pack(fill="x", pady=(4, 2))

        # -- Panel C: MDAnalysis expression --
        panel_mda = ctk.CTkFrame(body, fg_color="transparent")
        mda_hdr = ctk.CTkFrame(panel_mda, fg_color="transparent")
        mda_hdr.pack(fill="x", pady=(4, 0))
        ctk.CTkLabel(mda_hdr, text="MDAnalysis selection:").pack(
            side="left", anchor="w")
        ctk.CTkButton(
            mda_hdr, text="?", width=28, height=24,
            fg_color="#5a5a8a", hover_color="#7a7aaa",
            command=lambda: _show_selection_help(d)
        ).pack(side="right")
        mda_entry = ctk.CTkEntry(panel_mda,
                                 placeholder_text="e.g. protein and not backbone")
        mda_entry.pack(fill="x", pady=2)

        panels = {"criteria": panel_crit,
                  "existing": panel_exist,
                  "mda": panel_mda}

        def _switch_mode(*_):
            for p in panels.values():
                p.pack_forget()
            panels[mode_var.get()].pack(fill="x")

        # Mode radio buttons
        modes_frame = ctk.CTkFrame(body, fg_color="transparent")
        modes_frame.pack(fill="x", pady=(0, 2))
        for label, val in [("By criteria", "criteria"),
                           ("Existing selection", "existing"),
                           ("MDAnalysis expression", "mda")]:
            ctk.CTkRadioButton(
                modes_frame, text=label, variable=mode_var, value=val,
                command=_switch_mode
            ).pack(anchor="w", pady=1)

        _switch_mode()  # show initial panel

        # --- preview label ---
        preview_lbl = ctk.CTkLabel(d, text="", font=("", 11))
        preview_lbl.pack(pady=(4, 0), padx=12, anchor="w")

        def _gather_indices():
            mode = mode_var.get()
            if mode == "criteria":
                crit = crit_menu.get()
                if crit == 'MDAnalysis expression...':
                    expr = _crit_extra[0].get().strip()
                    return _resolve_mda_expression(
                        self.structure, expr, parent=d)
                extra = (_crit_extra[0].get().strip()
                         if crit in ('Chain...', 'Residue range...') else '')
                return self._resolve_sel_criteria(crit, extra)
            elif mode == "existing":
                chosen = exist_menu.get()
                if chosen == "(no selections)":
                    return []
                sel = next((s for s in self.selections
                            if s.name == chosen), None)
                return list(sel.atom_indices) if sel else []
            else:  # mda
                expr = mda_entry.get().strip()
                return _resolve_mda_expression(
                    self.structure, expr, parent=d)

        def _preview():
            _clear_highlight()
            indices = _gather_indices()
            if not indices:
                preview_lbl.configure(text="No atoms matched")
                return
            preview_lbl.configure(text=f"Matched: {len(indices)} atom(s)")
            # Build highlight actors (keep current camera view)
            actors = self._make_highlight_actors(indices)
            for a in actors:
                self.vtk_frame.add_actor(a)
                _highlight_actors.append(a)
            self.vtk_frame.render()

        btn_frame = ctk.CTkFrame(d, fg_color="transparent")
        btn_frame.pack(fill="x", padx=12, pady=(6, 4))
        ctk.CTkButton(btn_frame, text="Preview", height=28,
                      fg_color="#5a5a6e", hover_color="#7a7a9e",
                      command=_preview
                      ).pack(fill="x", pady=(0, 4))

        def _do_delete():
            _clear_highlight()
            indices = _gather_indices()
            if not indices:
                messagebox.showwarning("Warning", "No atoms matched",
                                       parent=d)
                return
            to_remove = set(indices)
            if not messagebox.askyesno(
                    "Confirm",
                    f"Delete {len(to_remove)} atom(s)?\n"
                    "This cannot be undone.", parent=d):
                return
            # Build old→new index mapping
            idx_map = {}
            new_idx = 0
            for old_idx in range(len(self.structure.atoms)):
                if old_idx not in to_remove:
                    idx_map[old_idx] = new_idx
                    new_idx += 1
            new_atoms = [a for i, a in enumerate(self.structure.atoms)
                         if i not in to_remove]
            self.structure.atoms = new_atoms
            self.structure._rebuild_residues_and_chains()
            self.structure.build_bonds()
            # Remap selections: update indices, drop empty ones
            kept = []
            for sel in self.selections:
                new_indices = [idx_map[i] for i in sel.atom_indices
                               if i in idx_map]
                if sel.name == "All":
                    sel.atom_indices = list(range(len(self.structure.atoms)))
                    kept.append(sel)
                elif new_indices:
                    sel.atom_indices = new_indices
                    kept.append(sel)
                # else: selection is now empty → drop it
            self.selections.clear()
            self.selections.extend(kept)
            # Re-assign secondary structure (psique / heuristic)
            import tempfile, os
            fd, tmp = tempfile.mkstemp(suffix='.pdb')
            os.close(fd)
            try:
                self.structure.write_pdb(tmp)
                _assign_secondary_structure(self.structure, filepath=tmp)
            finally:
                try:
                    os.unlink(tmp)
                except OSError:
                    pass
            d.destroy()
            s = self.structure
            self.info_label.configure(
                text=f"{len(s.atoms)} atoms  {len(s.residues)} residues\n"
                     f"{len(s.chains)} chain(s)")
            self._safe_refresh_rebuild()

        def _on_cancel():
            _clear_highlight()
            d.destroy()

        ctk.CTkButton(btn_frame, text="Delete", height=30,
                      fg_color="#8a3a3a", hover_color="#cc3333",
                      command=_do_delete).pack(fill="x", pady=(0, 2))
        ctk.CTkButton(btn_frame, text="Cancel", height=30,
                      fg_color="gray40",
                      command=_on_cancel).pack(fill="x")

        d.protocol("WM_DELETE_WINDOW", _on_cancel)

    # ------------------------------------------------------------------
    # Highlight helpers for preview
    # ------------------------------------------------------------------

    def _make_highlight_actors(self, indices):
        """Create semi-transparent glow spheres around the given atom indices.

        Each atom gets a translucent sphere slightly larger than its VDW
        radius, rendered in yellow/orange with smooth shading.  This gives
        a soft "aura" that clearly marks the selection without obscuring
        the underlying representation.
        """
        from vtkmodules.vtkFiltersSources import vtkSphereSource
        from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper
        from vtkmodules.vtkFiltersCore import vtkAppendPolyData

        atoms = self.structure.atoms
        appender = vtkAppendPolyData()
        for i in indices:
            a = atoms[i]
            vdw_r = VDW_RADII.get(a.element, VDW_RADII['DEFAULT'])
            src = vtkSphereSource()
            src.SetCenter(*a.coord)
            src.SetRadius(vdw_r * 1.25)  # ~25% larger than atom
            src.SetPhiResolution(16)
            src.SetThetaResolution(16)
            src.Update()
            appender.AddInputData(src.GetOutput())

        appender.Update()
        mapper = vtkPolyDataMapper()
        mapper.SetInputData(appender.GetOutput())

        actor = vtkActor()
        actor.SetMapper(mapper)
        prop = actor.GetProperty()
        prop.SetColor(1.0, 0.85, 0.0)     # warm yellow-orange
        prop.SetOpacity(0.20)
        prop.SetAmbient(0.6)
        prop.SetDiffuse(0.4)
        prop.SetSpecular(0.0)
        prop.LightingOn()
        return [actor]

    def _focus_camera_on(self, indices):
        """Move the camera to focus on the given atom indices."""
        atoms = self.structure.atoms
        coords = np.array([atoms[i].coord for i in indices])
        center = coords.mean(axis=0)
        extent = coords.max(axis=0) - coords.min(axis=0)
        radius = max(np.linalg.norm(extent) * 0.7, 5.0)
        cam = self.vtk_frame.renderer.GetActiveCamera()
        cam.SetFocalPoint(*center)
        cam.SetPosition(center[0], center[1], center[2] + radius * 2.5)
        cam.SetViewUp(0, 1, 0)
        self.vtk_frame.renderer.ResetCameraClippingRange()

    # ------------------------------------------------------------------
    # Background, rendering passes, depth cueing
    # ------------------------------------------------------------------

    def _change_bg_color(self):
        d = ctk.CTkToplevel(self)
        d.title("Background Color")
        d.geometry("220x120")
        d.transient(self)
        d.attributes('-topmost', True)
        d.after(100, d.grab_set)

        def _pick():
            d.destroy()
            c = colorchooser.askcolor(title="Background color")
            if c[0]:
                r, g, b = c[0][0] / 255, c[0][1] / 255, c[0][2] / 255
                self.vtk_frame.renderer.SetBackground(r, g, b)
                self.vtk_frame.render()

        def _reset():
            d.destroy()
            self.vtk_frame.renderer.SetBackground(0.18, 0.18, 0.18)
            self.vtk_frame.render()

        ctk.CTkButton(d, text="Pick Color...", height=30,
                      command=_pick).pack(pady=(12, 4), padx=12, fill="x")
        ctk.CTkButton(d, text="Reset to Default", height=30,
                      fg_color="gray40",
                      command=_reset).pack(pady=4, padx=12, fill="x")

    def _update_render_passes(self, *_):
        ren = self.vtk_frame.renderer
        ssao_on = self._ssao_var.get() == "on"
        shadows_on = self._shadows_var.get() == "on"
        self._ssao_enabled = ssao_on
        self._shadows_enabled = shadows_on
        if not ssao_on and not shadows_on:
            ren.SetPass(None)
            if self._shadow_light:
                ren.RemoveLight(self._shadow_light)
                self._shadow_light = None
            self.vtk_frame.render()
            return
        try:
            from vtkmodules.vtkRenderingOpenGL2 import (
                vtkSSAOPass, vtkRenderStepsPass, vtkShadowMapPass,
            )
            from vtkmodules.vtkRenderingCore import vtkLight
        except ImportError:
            self.vtk_frame.render()
            return

        basic = vtkRenderStepsPass()
        top_pass = basic

        if shadows_on:
            if not self._shadow_light:
                light = vtkLight()
                light.SetPositional(True)
                light.SetPosition(50, 80, 100)
                light.SetFocalPoint(0, 0, 0)
                light.SetConeAngle(60)
                light.SetIntensity(0.8)
                ren.AddLight(light)
                self._shadow_light = light
            shadow = vtkShadowMapPass()
            shadow.SetOpaqueSequence(basic)
            baker = shadow.GetShadowMapBakerPass()
            baker.SetResolution(2048)
            top_pass = shadow
        else:
            if self._shadow_light:
                ren.RemoveLight(self._shadow_light)
                self._shadow_light = None

        if ssao_on:
            ssao = vtkSSAOPass()
            ssao.SetDelegatePass(top_pass)
            bounds = ren.ComputeVisiblePropBounds()
            dx = bounds[1] - bounds[0]
            dy = bounds[3] - bounds[2]
            dz = bounds[5] - bounds[4]
            scene_size = (dx * dx + dy * dy + dz * dz) ** 0.5
            radius = max(1.0, scene_size * 0.05)
            ssao.SetRadius(radius)
            ssao.SetKernelSize(128)
            ssao.SetBias(0.001)
            ssao.BlurOn()
            top_pass = ssao

        ren.SetPass(top_pass)
        self.vtk_frame.render()

    def _on_depth_cue_change(self, val):
        self._depth_cue_density = float(val)
        self.vtk_frame.fog_density = self._depth_cue_density
        self.vtk_frame.render()

    def _apply_depth_cueing(self):
        self.vtk_frame.fog_density = self._depth_cue_density

    # ------------------------------------------------------------------
    # Perspective control
    # ------------------------------------------------------------------

    def _on_perspective_change(self, val):
        """Slider 0→1: 0 = orthographic, 1 = full perspective (view-angle 60°).

        Zoom now uses dolly (camera distance) in both modes, so the view
        angle is exclusively controlled by this slider — no size jumps.
        """
        import math as _m
        cam = self.vtk_frame.renderer.GetActiveCamera()
        t = float(val)
        if t < 0.01:
            # Switch to orthographic
            if not cam.GetParallelProjection():
                # Save distance-based scale so ortho matches current view
                cam.SetParallelScale(
                    cam.GetDistance()
                    * _m.tan(_m.radians(cam.GetViewAngle() / 2.0)))
            cam.SetParallelProjection(1)
        else:
            # Perspective: map slider with power curve for smooth onset
            angle = 5.0 + (t ** 1.5) * 55.0  # 5° … 60°
            if cam.GetParallelProjection():
                # Transition from ortho → perspective: match apparent size
                ps = cam.GetParallelScale()
                new_dist = ps / _m.tan(_m.radians(angle / 2.0))
                fp = cam.GetFocalPoint()
                pos = cam.GetPosition()
                d = [pos[i] - fp[i] for i in range(3)]
                norm = _m.sqrt(sum(x * x for x in d)) or 1.0
                d = [x / norm for x in d]
                cam.SetPosition(*[fp[i] + d[i] * new_dist for i in range(3)])
            cam.SetParallelProjection(0)
            cam.SetViewAngle(angle)
        self.vtk_frame.renderer.ResetCameraClippingRange()
        self.vtk_frame.render()

    # ------------------------------------------------------------------
    # Axes widget
    # ------------------------------------------------------------------

    def _axes_dialog(self):
        """Open small dialog to choose axis display mode."""
        d = ctk.CTkToplevel(self)
        d.title("XYZ Axes")
        d.geometry("230x210")
        d.transient(self)
        d.attributes('-topmost', True)
        d.after(100, d.grab_set)

        mode_var = ctk.StringVar(
            value=self._axes_mode_current or "off")

        for label, val in [("Off", "off"),
                           ("Corner (bottom-right)", "corner"),
                           ("Visible atoms center", "center"),
                           ("Origin (0, 0, 0)", "origin")]:
            ctk.CTkRadioButton(
                d, text=label, variable=mode_var, value=val,
            ).pack(anchor="w", padx=16, pady=3)

        def _apply():
            m = mode_var.get()
            self._axes_mode_current = m if m != "off" else None
            center = None
            if m == "center" and self.structure:
                import numpy as _np
                # Compute center from visible selections only
                vis_indices = set()
                for sel in self.selections:
                    if sel.visible:
                        vis_indices.update(sel.atom_indices)
                if vis_indices:
                    coords = _np.array(
                        [self.structure.atoms[i].coord for i in vis_indices])
                else:
                    coords = _np.array(
                        [a.coord for a in self.structure.atoms])
                center = tuple(coords.mean(axis=0).tolist())
            self.vtk_frame.set_axes(
                mode=self._axes_mode_current, center=center)
            self.vtk_frame.render()
            c = "#4a7a4a" if self._axes_mode_current else "gray35"
            self._axes_btn.configure(fg_color=c)
            d.destroy()

        ctk.CTkButton(d, text="Apply", height=30,
                      command=_apply).pack(pady=(10, 4), padx=12, fill="x")
        ctk.CTkButton(d, text="Cancel", height=30,
                      fg_color="gray40",
                      command=d.destroy).pack(pady=4, padx=12, fill="x")

    # ------------------------------------------------------------------
    # Reference lines
    # ------------------------------------------------------------------

    def _toggle_ref_lines(self):
        """Toggle X/Y/Z reference lines through the origin."""
        self._ref_lines_visible = not self._ref_lines_visible
        self.vtk_frame.set_reference_lines(self._ref_lines_visible)
        self.vtk_frame.render()
        c = "#4a7a4a" if self._ref_lines_visible else "gray35"
        self._ref_lines_btn.configure(fg_color=c)

    # ------------------------------------------------------------------
    # Viewpoint save / load
    # ------------------------------------------------------------------

    def _get_camera_state(self):
        cam = self.vtk_frame.renderer.GetActiveCamera()
        return {
            'position': list(cam.GetPosition()),
            'focal_point': list(cam.GetFocalPoint()),
            'view_up': list(cam.GetViewUp()),
            'clipping_range': list(cam.GetClippingRange()),
            'parallel_scale': cam.GetParallelScale(),
            'view_angle': cam.GetViewAngle(),
            'parallel_projection': bool(cam.GetParallelProjection()),
        }

    def _set_camera_state(self, state):
        cam = self.vtk_frame.renderer.GetActiveCamera()
        cam.SetPosition(*state['position'])
        cam.SetFocalPoint(*state['focal_point'])
        cam.SetViewUp(*state['view_up'])
        cam.SetClippingRange(*state['clipping_range'])
        cam.SetParallelScale(state.get('parallel_scale', 1.0))
        if 'view_angle' in state:
            cam.SetViewAngle(state['view_angle'])
        if 'parallel_projection' in state:
            cam.SetParallelProjection(state['parallel_projection'])

    def _sel_to_dict(self, sel):
        return {
            'name': sel.name,
            'atom_indices': sel.atom_indices,
            'representation': sel.representation,
            'color_scheme': sel.color_scheme,
            'uniform_color': list(sel.uniform_color) if sel.uniform_color else None,
            'visible': sel.visible,
            'surface_resolution': sel.surface_resolution,
            'surface_radius': sel.surface_radius,
            'opacity': sel.opacity,
            'carbon_color': list(sel.carbon_color) if sel.carbon_color else None,
            'quality': sel.quality,
            'atom_scale': sel.atom_scale,
            'bond_radius': sel.bond_radius,
            'ball_scale': sel.ball_scale,
            'stick_radius': sel.stick_radius,
            'backbone_radius': sel.backbone_radius,
            'helix_width': sel.helix_width,
            'sheet_width': sel.sheet_width,
            'coil_width': sel.coil_width,
            'ambient': sel.ambient,
            'diffuse': sel.diffuse,
            'specular': sel.specular,
            'specular_power': sel.specular_power,
            'ss_colors': {k: list(v) for k, v in sel.ss_colors.items()},
        }

    def _dict_to_sel(self, d):
        uc = tuple(d['uniform_color']) if d.get('uniform_color') else None
        cc = tuple(d['carbon_color']) if d.get('carbon_color') else None
        ss_c = {k: tuple(v) for k, v in d.get('ss_colors', {}).items()} or None
        return Selection(
            name=d['name'], atom_indices=d['atom_indices'],
            representation=d.get('representation', 'ball_stick'),
            color_scheme=d.get('color_scheme', 'element'),
            uniform_color=uc, visible=d.get('visible', True),
            surface_resolution=d.get('surface_resolution', 64),
            surface_radius=d.get('surface_radius', 0.12),
            opacity=d.get('opacity', 0.5), carbon_color=cc,
            quality=d.get('quality', 3),
            atom_scale=d.get('atom_scale', 1.0),
            bond_radius=d.get('bond_radius', 0.15),
            ball_scale=d.get('ball_scale', 0.3),
            stick_radius=d.get('stick_radius', 0.2),
            backbone_radius=d.get('backbone_radius', 0.3),
            helix_width=d.get('helix_width', 3.25),
            sheet_width=d.get('sheet_width', 2.5),
            coil_width=d.get('coil_width', 0.5),
            ambient=d.get('ambient', 0.2),
            diffuse=d.get('diffuse', 0.8),
            specular=d.get('specular', 0.05),
            specular_power=d.get('specular_power', 1.0),
            ss_colors=ss_c,
        )

    def _save_viewpoint(self):
        if not self.structure:
            messagebox.showwarning("Warning", "No structure loaded")
            return
        bg = list(self.vtk_frame.renderer.GetBackground())
        state = {
            'pdb_file': os.path.abspath(self._pdb_filepath) if self._pdb_filepath else '',
            'camera': self._get_camera_state(),
            'background': bg,
            'selections': [self._sel_to_dict(s) for s in self.selections],
            'ssao': self._ssao_enabled,
            'shadows': self._shadows_enabled,
            'depth_cue': self._depth_cue_density,
        }
        path = filedialog.asksaveasfilename(
            title="Save Viewpoint",
            initialdir=self.initial_directory,
            defaultextension=".json",
            filetypes=[("Viewpoint JSON", "*.json"), ("All", "*.*")])
        if path:
            try:
                with open(path, 'w') as f:
                    json.dump(state, f, indent=2)
                messagebox.showinfo("Saved",
                                    f"Viewpoint saved to {os.path.basename(path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Save failed: {e}")

    def _load_viewpoint(self):
        path = filedialog.askopenfilename(
            title="Load Viewpoint",
            initialdir=self.initial_directory,
            filetypes=[("Viewpoint JSON", "*.json"), ("All", "*.*")])
        if not path:
            return
        try:
            with open(path, 'r') as f:
                state = json.load(f)
            if not self.structure:
                pdb_path = state.get('pdb_file', '')
                if pdb_path and os.path.isfile(pdb_path):
                    self._load(pdb_path)
                else:
                    messagebox.showwarning(
                        "Warning",
                        f"PDB file not found:\n{pdb_path}\n\n"
                        "Please load the structure manually first.")
                    return
            n_atoms = len(self.structure.atoms)
            self._remove_all_labels()
            self._remove_all_measurements()
            self.selections.clear()
            for sd in state.get('selections', []):
                sd['atom_indices'] = [i for i in sd['atom_indices'] if i < n_atoms]
                if sd['atom_indices']:
                    self.selections.append(self._dict_to_sel(sd))
            if not self.selections:
                self.selections.append(
                    Selection("All", list(range(n_atoms)),
                              representation='vdw', color_scheme='element'))
            bg = state.get('background', [0.18, 0.18, 0.18])
            self.vtk_frame.renderer.SetBackground(*bg)
            if state.get('ssao'):
                self._ssao_var.set('on')
                self._ssao_switch.select()
            else:
                self._ssao_var.set('off')
                self._ssao_switch.deselect()
            if state.get('shadows'):
                self._shadows_var.set('on')
                self._shadows_switch.select()
            else:
                self._shadows_var.set('off')
                self._shadows_switch.deselect()
            dc = state.get('depth_cue', 0.0)
            self._depth_cue_density = dc
            self._depth_cue_var.set(dc)
            self._depth_cue_slider.set(dc)
            self.vtk_frame.fog_density = dc
            self._safe_refresh_rebuild()
            self._update_render_passes()
            if 'camera' in state:
                self._set_camera_state(state['camera'])
            self.vtk_frame.render()
            messagebox.showinfo("Loaded",
                                f"Viewpoint loaded from {os.path.basename(path)}")
        except Exception as e:
            messagebox.showerror("Error", f"Load failed: {e}")
            logger.error(f"Viewpoint load error: {e}")

    # ------------------------------------------------------------------
    # Save image / PDB
    # ------------------------------------------------------------------

    def _save_image(self):
        if not self.structure:
            messagebox.showwarning("Warning", "No structure loaded")
            return
        d = ctk.CTkToplevel(self)
        d.title("Save Image")
        d.geometry("320x300")
        d.transient(self)
        d.attributes('-topmost', True)
        d.after(100, d.grab_set)

        ctk.CTkLabel(d, text="Resolution Scale:", font=("", 12)).pack(
            pady=(10, 2), padx=12, anchor="w")
        scale_menu = ctk.CTkOptionMenu(d,
                                       values=["1x", "2x", "3x", "4x"],
                                       height=28)
        scale_menu.set("1x")
        scale_menu.pack(padx=12, fill="x")

        ctk.CTkLabel(d, text="Background:", font=("", 12)).pack(
            pady=(8, 2), padx=12, anchor="w")
        bg_menu = ctk.CTkOptionMenu(d,
                                    values=["As displayed", "Transparent"],
                                    height=28)
        bg_menu.set("As displayed")
        bg_menu.pack(padx=12, fill="x")

        ctk.CTkLabel(d, text="Format:", font=("", 12)).pack(
            pady=(8, 2), padx=12, anchor="w")
        fmt_menu = ctk.CTkOptionMenu(d,
                                     values=["PNG", "JPEG", "TIFF", "BMP"],
                                     height=28)
        fmt_menu.set("PNG")
        fmt_menu.pack(padx=12, fill="x")

        def _on_bg(val):
            if val == "Transparent":
                fmt_menu.configure(values=["PNG", "TIFF"])
                if fmt_menu.get() not in ("PNG", "TIFF"):
                    fmt_menu.set("PNG")
            else:
                fmt_menu.configure(values=["PNG", "JPEG", "TIFF", "BMP"])

        bg_menu.configure(command=_on_bg)

        def _do_save():
            scale = int(scale_menu.get()[0])
            transparent = bg_menu.get() == "Transparent"
            fmt = fmt_menu.get()
            ext_map = {"PNG": ".png", "JPEG": ".jpg",
                       "TIFF": ".tiff", "BMP": ".bmp"}
            ext = ext_map[fmt]
            d.destroy()
            save_path = filedialog.asksaveasfilename(
                title="Save Image",
                initialdir=self.initial_directory,
                defaultextension=ext,
                filetypes=[(f"{fmt} files", f"*{ext}"), ("All", "*.*")])
            if not save_path:
                return
            try:
                img = self.vtk_frame.render_to_image(scale=scale,
                                                     transparent=transparent)
                if fmt == "JPEG" and img.mode == "RGBA":
                    img = img.convert("RGB")
                img.save(save_path)
                messagebox.showinfo("Saved",
                                    f"Image saved to {os.path.basename(save_path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Save failed: {e}")

        ctk.CTkButton(d, text="Save...", height=32, fg_color="#2a6e2a",
                      hover_color="#358535",
                      command=_do_save).pack(pady=(12, 4), padx=12, fill="x")
        ctk.CTkButton(d, text="Cancel", height=28, fg_color="gray40",
                      command=d.destroy).pack(padx=12, fill="x", pady=(0, 8))

    def _save_pdb(self):
        if not self.structure:
            messagebox.showwarning("Warning", "No structure loaded")
            return
        path = filedialog.asksaveasfilename(
            title="Save PDB",
            initialdir=self.initial_directory,
            defaultextension=".pdb",
            filetypes=[("PDB files", "*.pdb"), ("All", "*.*")])
        if path:
            try:
                self.structure.write_pdb(path)
                messagebox.showinfo("Saved",
                                    f"Saved to {os.path.basename(path)}")
            except Exception as e:
                messagebox.showerror("Error", f"Save failed: {e}")

    # ------------------------------------------------------------------
    # Frame interface methods (required by app.py)
    # ------------------------------------------------------------------

    def on_stage_shown(self):
        """Called when this stage becomes active."""
        pass

    def on_pdb_changed(self, pdb_file: Optional[str]):
        """Called when PDB file changes in another frame."""
        if pdb_file != self.current_pdb_file:
            self.current_pdb_file = pdb_file
            if pdb_file:
                self._load(pdb_file)
            else:
                self.structure = None
                self.selections.clear()
                self.vtk_frame.clear_actors()
                self.vtk_frame.render()
                self.info_label.configure(text="No structure loaded")

    def cleanup(self):
        """Cleanup resources when frame is destroyed."""
        try:
            self.vtk_frame.clear_actors()
        except Exception as e:
            logger.error(f"Error during VTK cleanup: {e}")

    def update_fonts(self, scaled_fonts):
        """Update fonts in the visualization frame."""
        try:
            if hasattr(self, 'info_label'):
                self.info_label.configure(font=scaled_fonts.get('body', FONTS['body']))
        except Exception:
            pass


# --------------------------------------------------------------------------
# MDAnalysis selection helpers
# --------------------------------------------------------------------------

def _resolve_mda_expression(structure, expression, parent=None):
    """Resolve an MDAnalysis selection expression to atom indices.

    Writes the structure to a temporary PDB, creates an MDAnalysis Universe,
    runs the selection, and returns matching atom indices.
    """
    if not expression:
        return []
    try:
        import MDAnalysis as mda
    except ImportError:
        messagebox.showerror(
            "MDAnalysis Required",
            "MDAnalysis is not installed.\n"
            "Install it with: pip install MDAnalysis",
            parent=parent)
        return []
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
    try:
        tmp_path = tmp.name
        tmp.close()
        structure.write_pdb(tmp_path)
        u = mda.Universe(tmp_path)
        ag = u.select_atoms(expression)
        return list(ag.indices)
    except Exception as exc:
        messagebox.showerror(
            "Selection Error",
            f"Invalid MDAnalysis expression:\n{exc}",
            parent=parent)
        return []
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass


def _show_selection_help(parent):
    """Show help window for MDAnalysis atom selection syntax."""
    hw = ctk.CTkToplevel(parent)
    hw.title("MDAnalysis Selection Help")
    hw.geometry("620x520")
    hw.transient(parent)
    hw.attributes('-topmost', True)
    hw.after(100, hw.grab_set)

    ctk.CTkLabel(
        hw, text="MDAnalysis Atom Selection Syntax",
        font=("", 15, "bold")
    ).pack(pady=(10, 4), padx=10)

    tb = ctk.CTkTextbox(hw, wrap="word")
    tb.pack(fill="both", expand=True, padx=10, pady=(4, 8))

    content = """\
MDAnalysis uses a powerful and flexible atom selection language.

BASIC SELECTIONS:
  protein          -- All protein atoms
  backbone         -- Protein backbone atoms (CA, C, N, O)
  name CA          -- All alpha-carbon atoms
  all              -- All atoms
  water            -- All water molecules
  resname ALA      -- All alanine residues

COMBINING SELECTIONS:
  protein and backbone        -- Protein backbone only
  protein and not name H*     -- Protein without hydrogens
  name CA or name CB          -- Alpha and beta carbons
  resid 1:50                  -- Residues 1 through 50
  protein and resid 10:100    -- Protein residues 10-100

RESIDUE SELECTIONS:
  resname ALA GLY VAL         -- Specific amino acids
  resid 1 5 10                -- Specific residue numbers
  resid 1:50 and name CA      -- CA atoms in residues 1-50

CHAIN / SEGMENT SELECTIONS:
  segid A            -- Segment (chain) A
  segid A B          -- Segments A and B

ATOM PROPERTIES:
  type CA            -- Atoms of type CA
  mass > 12          -- Atoms with mass > 12

SPATIAL SELECTIONS:
  around 5.0 protein -- Atoms within 5 A of protein
  sphzone 10 name CA -- Atoms in sphere of 10 A around CA

EXAMPLES FOR VISUALIZATION:
  "protein"
      Select all protein atoms

  "name CA"
      Select alpha-carbon atoms only

  "resname LIG"
      Select ligand residues named LIG

  "protein and resid 50:150"
      Select specific protein region

  "not (protein or water)"
      Select everything except protein and water

  "segid A and resid 1:100 and backbone"
      Chain A backbone for residues 1-100

For full documentation see:
https://docs.mdanalysis.org/stable/documentation_pages/selections.html
"""
    tb.insert("0.0", content)
    tb.configure(state="disabled")

    ctk.CTkButton(hw, text="Close", width=100,
                  command=hw.destroy).pack(pady=(0, 10))


# --------------------------------------------------------------------------
# Selection creation dialog
# --------------------------------------------------------------------------

class SelectionDialog(ctk.CTkToplevel):
    def __init__(self, parent, structure, existing):
        super().__init__(parent)
        self.title("New Selection")
        self.geometry("360x560")
        self.transient(parent)
        self.attributes('-topmost', True)
        self.after(100, self.grab_set)
        self.structure = structure
        self.existing = existing
        self.result = None
        self._picked_color = None

        ctk.CTkLabel(self, text="Selection Name:").pack(
            pady=(8, 1), padx=10, anchor="w")
        self.name_entry = ctk.CTkEntry(self,
                                       placeholder_text="e.g. Active site",
                                       height=28)
        self.name_entry.pack(pady=(0, 4), padx=10, fill="x")

        ctk.CTkLabel(self, text="Select:").pack(
            pady=(4, 1), padx=10, anchor="w")
        self.crit_menu = ctk.CTkOptionMenu(self, values=SELECTION_CRITERIA,
                                           command=self._on_crit, height=28)
        self.crit_menu.pack(pady=(0, 2), padx=10, fill="x")

        self.extra_frame = ctk.CTkFrame(self, fg_color="transparent")
        self.extra_entry = ctk.CTkEntry(self.extra_frame)
        self.around_frame = ctk.CTkFrame(self, fg_color="transparent")

        self.rep_label = ctk.CTkLabel(self, text="Representation:")
        self.rep_label.pack(pady=(4, 1), padx=10, anchor="w")
        self.rep_menu = ctk.CTkOptionMenu(self, values=REPRESENTATIONS_UI,
                                          height=28)
        self.rep_menu.pack(pady=(0, 2), padx=10, fill="x")
        self.rep_menu.set("Ball & Stick")

        ctk.CTkLabel(self, text="Color Scheme:").pack(
            pady=(4, 1), padx=10, anchor="w")
        self.cs_menu = ctk.CTkOptionMenu(self, values=COLOR_SCHEMES_UI,
                                         height=28)
        self.cs_menu.pack(pady=(0, 2), padx=10, fill="x")
        self.cs_menu.set("Element (CPK)")

        color_frame = ctk.CTkFrame(self, fg_color="transparent")
        color_frame.pack(fill="x", padx=10, pady=(4, 0))
        self._color_btn = ctk.CTkButton(
            color_frame, text="Pick Specific Color",
            height=28, fg_color="#5a3a6e", hover_color="#7a4a9e",
            command=self._pick_specific_color)
        self._color_btn.pack(side="left", fill="x", expand=True)
        self._color_preview = ctk.CTkLabel(color_frame, text="  ", width=28,
                                           fg_color="gray40", corner_radius=4)
        self._color_preview.pack(side="right", padx=(4, 0))

        self._q_label = ctk.CTkLabel(self, text="Quality: Medium")
        self._q_label.pack(pady=(6, 1), padx=10, anchor="w")
        self.q_slider = ctk.CTkSlider(self, from_=1, to=5, number_of_steps=4)
        self.q_slider.set(3)
        self.q_slider.configure(command=self._on_quality)
        self.q_slider.pack(pady=(0, 2), padx=10, fill="x")

        bf = ctk.CTkFrame(self, fg_color="transparent")
        bf.pack(fill="x", padx=10, pady=(8, 8))
        ctk.CTkButton(bf, text="Create", command=self._create, height=30,
                      fg_color="#2a6e2a",
                      hover_color="#358535").pack(side="left", expand=True, padx=4)
        ctk.CTkButton(bf, text="Cancel", command=self.destroy, height=30,
                      fg_color="gray40").pack(side="left", expand=True, padx=4)

    def _pick_specific_color(self):
        c = colorchooser.askcolor(title="Pick color for this selection",
                                  parent=self)
        if c[0]:
            self._picked_color = (int(c[0][0]), int(c[0][1]), int(c[0][2]))
            hex_c = "#%02X%02X%02X" % self._picked_color
            self._color_preview.configure(fg_color=hex_c)
            self._color_btn.configure(text=f"Color: {hex_c}")
        else:
            self._picked_color = None
            self._color_preview.configure(fg_color="gray40")
            self._color_btn.configure(text="Pick Specific Color")

    def _on_quality(self, val):
        self._q_label.configure(text=f"Quality: {QUALITY_LABELS[int(val) - 1]}")

    def _on_crit(self, val):
        for w in self.extra_frame.winfo_children():
            w.destroy()
        for w in self.around_frame.winfo_children():
            w.destroy()
        self.extra_frame.pack_forget()
        self.around_frame.pack_forget()
        if hasattr(self, '_mda_frame'):
            self._mda_frame.pack_forget()
        if val == 'Chain...':
            chains = ', '.join(sorted(self.structure.chains.keys()))
            ctk.CTkLabel(self.extra_frame,
                         text=f"Available: {chains}").pack(anchor="w", pady=2)
            self.extra_entry = ctk.CTkEntry(
                self.extra_frame,
                placeholder_text="Chain IDs (e.g. A B or A,B)")
            self.extra_entry.pack(fill="x", pady=2)
            self.extra_frame.pack(fill="x", padx=10, before=self.rep_label)
        elif val == 'Residue range...':
            self.extra_entry = ctk.CTkEntry(
                self.extra_frame,
                placeholder_text="e.g. A:50 or A:10-50,B:1-20")
            self.extra_entry.pack(fill="x", pady=2)
            self.extra_frame.pack(fill="x", padx=10, before=self.rep_label)
        elif val == 'Around selection...':
            if self.existing:
                names = [s.name for s in self.existing]
                ctk.CTkLabel(self.around_frame,
                             text="Base selection:").pack(anchor="w", pady=2)
                self.around_sel = ctk.CTkOptionMenu(self.around_frame,
                                                    values=names)
                self.around_sel.pack(fill="x", pady=2)
                ctk.CTkLabel(self.around_frame,
                             text="Distance (A):").pack(anchor="w", pady=2)
                self.around_dist = ctk.CTkEntry(
                    self.around_frame, placeholder_text="5.0")
                self.around_dist.pack(fill="x", pady=2)
                self.around_frame.pack(fill="x", padx=10,
                                       before=self.rep_label)
        elif val == 'MDAnalysis expression...':
            if not hasattr(self, '_mda_frame'):
                self._mda_frame = ctk.CTkFrame(self, fg_color="transparent")
            else:
                for w in self._mda_frame.winfo_children():
                    w.destroy()
            hdr = ctk.CTkFrame(self._mda_frame, fg_color="transparent")
            hdr.pack(fill="x", pady=(2, 0))
            ctk.CTkLabel(hdr, text="MDAnalysis selection:").pack(
                side="left", anchor="w")
            ctk.CTkButton(
                hdr, text="?", width=28, height=24,
                fg_color="#5a5a8a", hover_color="#7a7aaa",
                command=lambda: _show_selection_help(self)
            ).pack(side="right")
            self.mda_entry = ctk.CTkEntry(
                self._mda_frame,
                placeholder_text="e.g. protein and name CA")
            self.mda_entry.pack(fill="x", pady=2)
            self._mda_frame.pack(fill="x", padx=10, before=self.rep_label)

    def _create(self):
        name = self.name_entry.get().strip()
        if not name:
            messagebox.showwarning("Warning", "Enter a name", parent=self)
            return
        crit = self.crit_menu.get()
        indices = self._resolve(crit)
        if not indices:
            messagebox.showwarning("Warning", "No atoms matched", parent=self)
            return
        uc = self._picked_color
        cs = _cs_key(self.cs_menu.get())
        if uc:
            cs = 'uniform'
        # Determine extra text for criteria that need it
        crit_extra = ''
        if crit in ('Chain...', 'Residue range...'):
            crit_extra = self.extra_entry.get().strip()
        elif crit == 'MDAnalysis expression...':
            mda = getattr(self, 'mda_entry', None)
            crit_extra = mda.get().strip() if mda else ''
        self.result = Selection(
            name, indices,
            representation=_rep_key(self.rep_menu.get()),
            color_scheme=cs,
            uniform_color=uc,
            quality=int(self.q_slider.get()),
            criteria=crit,
            criteria_extra=crit_extra)
        self.destroy()

    def _resolve(self, crit):
        atoms = self.structure.atoms
        if crit == 'All':
            return list(range(len(atoms)))
        if crit == 'Protein':
            return [i for i, a in enumerate(atoms) if a.res_name in AA_NAMES]
        if crit == 'Backbone':
            return [i for i, a in enumerate(atoms)
                    if a.res_name in AA_NAMES and a.name in BACKBONE_NAMES]
        if crit == 'Sidechain':
            return [i for i, a in enumerate(atoms)
                    if a.res_name in AA_NAMES and a.name not in BACKBONE_NAMES]
        if crit == 'Water':
            return [i for i, a in enumerate(atoms) if a.res_name in ('HOH', 'WAT', 'TIP')]
        if crit == 'Ligand':
            return [i for i, a in enumerate(atoms)
                    if a.res_name not in AA_NAMES
                    and a.res_name not in ('HOH', 'WAT', 'TIP')]
        if crit == 'Chain...':
            raw = self.extra_entry.get().strip().upper()
            chains = {c.strip() for c in raw.replace(',', ' ').split() if c.strip()}
            return [i for i, a in enumerate(atoms) if a.chain_id in chains]
        if crit == 'Residue range...':
            return self._parse_range(self.extra_entry.get().strip())
        if crit == 'Around selection...':
            return self._around()
        if crit == 'MDAnalysis expression...':
            expr = getattr(self, 'mda_entry', None)
            if expr is None:
                return []
            return _resolve_mda_expression(self.structure, expr.get().strip(),
                                           parent=self)
        return []

    def _parse_range(self, text):
        indices = []
        for part in text.split(','):
            part = part.strip()
            if ':' not in part:
                continue
            ch, rng = part.split(':', 1)
            ch = ch.strip().upper()
            try:
                if '-' in rng:
                    lo, hi = int(rng.split('-')[0]), int(rng.split('-')[1])
                else:
                    lo = hi = int(rng.strip())
            except ValueError:
                continue
            for i, a in enumerate(self.structure.atoms):
                if a.chain_id == ch and lo <= a.res_id <= hi:
                    indices.append(i)
        return indices

    def _around(self):
        try:
            base_name = self.around_sel.get()
            dist = float(self.around_dist.get().strip())
        except (AttributeError, ValueError):
            return []
        base = next((s for s in self.existing if s.name == base_name), None)
        if not base:
            return []
        bc = np.array([self.structure.atoms[i].coord for i in base.atom_indices])
        bs = set(base.atom_indices)
        result = []
        for i, a in enumerate(self.structure.atoms):
            if i in bs:
                continue
            if np.min(np.linalg.norm(bc - a.coord, axis=1)) <= dist:
                result.append(i)
        return result
