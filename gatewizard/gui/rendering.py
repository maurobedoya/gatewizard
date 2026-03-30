# gatewizard/gui/rendering.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
VTK actor builders for molecular visualization.

Builds VTK actors for various molecular representations (VDW, Ball & Stick,
Sticks, Cartoon, Tube SS, Backbone, Surface) from the data model objects
in ``gatewizard.core.viewer``.
"""

import math
from typing import Dict, List, Optional, Tuple

import numpy as np

from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkProperty
from vtkmodules.vtkFiltersSources import vtkSphereSource
from vtkmodules.vtkFiltersCore import (
    vtkAppendPolyData,
    vtkGlyph3D,
    vtkTubeFilter,
    vtkPolyDataNormals,
)
from vtkmodules.vtkCommonCore import (
    vtkPoints,
    vtkFloatArray,
    vtkUnsignedCharArray,
    vtkIdList,
)
from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkCellArray, vtkLine

from gatewizard.core.viewer import (
    VDW_RADII,
    COVALENT_RADII,
    ELEMENT_COLORS,
    SS_COLORS,
    CHAIN_PALETTE,
    HELIX_SS_TYPES,
    ProteinStructure,
    RESIDUE_NATURE,
    RESIDUE_NATURE_COLORS,
)

# Quality presets: (sphere_theta, sphere_phi, tube_sides, smooth_factor)
QUALITY_PRESETS = {
    1: (6, 6, 4, 2),
    2: (10, 10, 6, 3),
    3: (16, 16, 8, 4),
    4: (24, 24, 12, 6),
    5: (48, 48, 32, 16),
}
QUALITY_LABELS = ["Low", "Medium-Low", "Medium", "High", "Ultra"]

MATERIAL_PRESETS = {
    "Default": (0.1, 0.7, 0.3, 20.0),
    "Shiny": (0.05, 0.6, 0.8, 40.0),
    "Matte": (0.2, 0.8, 0.05, 1.0),
    "Metallic": (0.15, 0.45, 0.7, 60.0),
    "Plastic": (0.1, 0.6, 0.5, 30.0),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _rgb_f(rgb):
    return (rgb[0] / 255.0, rgb[1] / 255.0, rgb[2] / 255.0)


def _atom_color(atom, scheme="element", carbon_color=None, ss_map=None):
    if carbon_color and atom.element == "C":
        return carbon_color
    if scheme == "element":
        return ELEMENT_COLORS.get(atom.element, ELEMENT_COLORS["DEFAULT"])
    if scheme == "chain":
        idx = (
            list("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789").index(atom.chain_id)
            if atom.chain_id in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
            else 0
        )
        return CHAIN_PALETTE[idx % len(CHAIN_PALETTE)]
    if scheme == "residue_nature":
        cat = RESIDUE_NATURE.get(atom.res_name, "other")
        return RESIDUE_NATURE_COLORS[cat]
    if scheme == "ss":
        if ss_map:
            ss = ss_map.get((atom.chain_id, atom.res_id), "C")
            return SS_COLORS.get(ss, SS_COLORS["DEFAULT"])
        return ELEMENT_COLORS.get(atom.element, ELEMENT_COLORS["DEFAULT"])
    return ELEMENT_COLORS.get(atom.element, ELEMENT_COLORS["DEFAULT"])


def _res_color(res, scheme="ss", ss_colors=None):
    colors = ss_colors if ss_colors else SS_COLORS
    if scheme == "ss":
        return colors.get(res.ss, colors.get("DEFAULT", SS_COLORS["DEFAULT"]))
    if scheme == "chain":
        idx = (
            list("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789").index(res.chain_id)
            if res.chain_id in "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
            else 0
        )
        return CHAIN_PALETTE[idx % len(CHAIN_PALETTE)]
    if scheme == "residue_nature":
        cat = RESIDUE_NATURE.get(res.name, "other")
        return RESIDUE_NATURE_COLORS[cat]
    return colors.get(res.ss, colors.get("DEFAULT", SS_COLORS["DEFAULT"]))


# ---------------------------------------------------------------------------
# VDW / Spacefill
# ---------------------------------------------------------------------------


def make_vdw_actor(
    atoms,
    color_scheme="element",
    uniform_color=None,
    scale=1.0,
    carbon_color=None,
    quality=3,
    ss_map=None,
):
    if not atoms:
        return None
    qp = QUALITY_PRESETS.get(quality, QUALITY_PRESETS[3])
    points = vtkPoints()
    colors = vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    colors.SetName("Colors")
    radii = vtkFloatArray()
    radii.SetName("Radii")
    for a in atoms:
        points.InsertNextPoint(*a.coord)
        r = VDW_RADII.get(a.element, VDW_RADII["DEFAULT"]) * scale
        radii.InsertNextValue(r)
        c = (
            uniform_color
            if uniform_color
            else _atom_color(a, color_scheme, carbon_color, ss_map=ss_map)
        )
        colors.InsertNextTuple3(*c)
    pd = vtkPolyData()
    pd.SetPoints(points)
    pd.GetPointData().SetScalars(colors)
    pd.GetPointData().AddArray(radii)
    sphere = vtkSphereSource()
    sphere.SetThetaResolution(qp[0])
    sphere.SetPhiResolution(qp[1])
    sphere.SetRadius(1.0)
    glyph = vtkGlyph3D()
    glyph.SetInputData(pd)
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.SetScaleModeToScaleByScalar()
    glyph.SetColorModeToColorByScalar()
    glyph.SetScaleFactor(1.0)
    glyph.SetInputArrayToProcess(0, 0, 0, 0, "Radii")
    glyph.Update()
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())
    mapper.ScalarVisibilityOn()
    actor = vtkActor()
    actor.SetMapper(mapper)
    return actor


# ---------------------------------------------------------------------------
# Sticks (bonds as tubes)
# ---------------------------------------------------------------------------


def _make_sticks(
    atoms,
    bonds,
    color_scheme="element",
    uniform_color=None,
    radius=0.15,
    carbon_color=None,
    quality=3,
    ss_map=None,
):
    if not bonds:
        return None
    points = vtkPoints()
    lines = vtkCellArray()
    colors = vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    colors.SetName("Colors")
    pt_idx = 0
    atom_set = set(range(len(atoms)))
    for ai, aj in bonds:
        if ai not in atom_set or aj not in atom_set:
            continue
        a1, a2 = atoms[ai], atoms[aj]
        mid = (a1.coord + a2.coord) / 2.0
        c1 = (
            uniform_color
            if uniform_color
            else _atom_color(a1, color_scheme, carbon_color, ss_map=ss_map)
        )
        c2 = (
            uniform_color
            if uniform_color
            else _atom_color(a2, color_scheme, carbon_color, ss_map=ss_map)
        )
        p1 = pt_idx
        points.InsertNextPoint(*a1.coord)
        colors.InsertNextTuple3(*c1)
        pt_idx += 1
        p2 = pt_idx
        points.InsertNextPoint(*mid)
        colors.InsertNextTuple3(*c1)
        pt_idx += 1
        ln = vtkLine()
        ln.GetPointIds().SetId(0, p1)
        ln.GetPointIds().SetId(1, p2)
        lines.InsertNextCell(ln)
        p3 = pt_idx
        points.InsertNextPoint(*mid)
        colors.InsertNextTuple3(*c2)
        pt_idx += 1
        p4 = pt_idx
        points.InsertNextPoint(*a2.coord)
        colors.InsertNextTuple3(*c2)
        pt_idx += 1
        ln2 = vtkLine()
        ln2.GetPointIds().SetId(0, p3)
        ln2.GetPointIds().SetId(1, p4)
        lines.InsertNextCell(ln2)
    pd = vtkPolyData()
    pd.SetPoints(points)
    pd.SetLines(lines)
    pd.GetPointData().SetScalars(colors)
    tube = vtkTubeFilter()
    qp = QUALITY_PRESETS.get(quality, QUALITY_PRESETS[3])
    tube.SetInputData(pd)
    tube.SetRadius(radius)
    tube.SetNumberOfSides(qp[2])
    tube.CappingOn()
    tube.Update()
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())
    mapper.ScalarVisibilityOn()
    actor = vtkActor()
    actor.SetMapper(mapper)
    return actor


# ---------------------------------------------------------------------------
# Ball & Stick
# ---------------------------------------------------------------------------


def make_ball_stick_actor(
    atoms,
    bonds,
    color_scheme="element",
    uniform_color=None,
    ball_scale=0.3,
    stick_radius=0.15,
    carbon_color=None,
    quality=3,
    ss_map=None,
):
    actors = []
    ba = make_vdw_actor(
        atoms,
        color_scheme,
        uniform_color,
        scale=ball_scale,
        carbon_color=carbon_color,
        quality=quality,
        ss_map=ss_map,
    )
    if ba:
        actors.append(ba)
    sa = _make_sticks(
        atoms,
        bonds,
        color_scheme,
        uniform_color,
        stick_radius,
        carbon_color=carbon_color,
        quality=quality,
        ss_map=ss_map,
    )
    if sa:
        actors.append(sa)
    return actors


# ---------------------------------------------------------------------------
# Sticks only
# ---------------------------------------------------------------------------


def make_stick_actor(
    atoms,
    bonds,
    color_scheme="element",
    uniform_color=None,
    radius=0.2,
    carbon_color=None,
    quality=3,
    ss_map=None,
):
    actors = []
    sa = _make_sticks(
        atoms,
        bonds,
        color_scheme,
        uniform_color,
        radius,
        carbon_color=carbon_color,
        quality=quality,
        ss_map=ss_map,
    )
    if sa:
        actors.append(sa)
    return actors


# ---------------------------------------------------------------------------
# Cartoon
# ---------------------------------------------------------------------------


def make_cartoon_actor(
    structure,
    chains=None,
    color_scheme="ss",
    uniform_color=None,
    helix_w=2.0,
    sheet_w=2.5,
    coil_w=0.5,
    quality=3,
    ss_colors=None,
):
    """Cartoon: continuous flat ribbon per chain, width/color by SS."""
    actors = []
    chain_ids = chains or list(structure.chains.keys())
    qp = QUALITY_PRESETS.get(quality, QUALITY_PRESETS[3])
    smooth_factor = qp[3] * 2
    cs_sides = max(8, qp[2])
    half_thick = 0.15
    for cid in chain_ids:
        residues = structure.chains.get(cid, [])
        valid_res, ca_all, o_all = [], [], []
        for r in residues:
            if r.ca_coord is not None:
                valid_res.append(r)
                ca_all.append(r.ca_coord)
                o_all.append(r.o_coord)
        nr = len(ca_all)
        if nr < 2:
            continue
        chain_normals = _compute_ribbon_normals(ca_all, o_all)
        sm_pts, sm_nms = _smooth_coords_and_normals(
            ca_all, chain_normals, smooth_factor
        )
        n_sm = len(sm_pts)
        ss_per_pt, res_per_pt = [], []
        for ri in range(nr):
            count = smooth_factor if ri < nr - 1 else 1
            for _ in range(count):
                ss_per_pt.append(valid_res[ri].ss)
                res_per_pt.append(ri)
        segments = []
        seg_start = 0
        for j in range(1, n_sm):
            if ss_per_pt[j] != ss_per_pt[j - 1]:
                segments.append((seg_start, j - 1, ss_per_pt[seg_start]))
                seg_start = j
        segments.append((seg_start, n_sm - 1, ss_per_pt[seg_start]))
        hw = helix_w * 0.35
        sw = sheet_w * 0.35
        aw = sheet_w * 0.55
        cw = coil_w * 0.25
        widths = np.empty(n_sm)
        thicknesses = np.empty(n_sm)
        for i in range(n_sm):
            ss = ss_per_pt[i]
            if ss in HELIX_SS_TYPES:
                widths[i] = hw
                thicknesses[i] = half_thick
            elif ss == "E":
                widths[i] = sw
                thicknesses[i] = half_thick
            else:
                widths[i] = cw
                thicknesses[i] = cw
        for s_s, s_e, ss in segments:
            if ss != "E":
                continue
            seg_len = s_e - s_s + 1
            arrow_start = s_s + max(1, int(seg_len * 0.70))
            for i in range(arrow_start, s_e + 1):
                frac = (i - arrow_start) / max(1, s_e - arrow_start)
                widths[i] = aw * (1.0 - frac)
        transition = max(3, smooth_factor)
        for seg_idx in range(len(segments) - 1):
            bnd = segments[seg_idx][1] + 1
            for d in range(transition):
                t = 1.0 - d / transition
                alpha = 0.5 * t
                ib = bnd - 1 - d
                ia = bnd + d
                if 0 <= ib < n_sm and 0 <= ia < n_sm:
                    avg_w = 0.5 * (widths[ib] + widths[ia])
                    widths[ib] += alpha * (avg_w - widths[ib])
                    widths[ia] += alpha * (avg_w - widths[ia])
                    avg_t = 0.5 * (thicknesses[ib] + thicknesses[ia])
                    thicknesses[ib] += alpha * (avg_t - thicknesses[ib])
                    thicknesses[ia] += alpha * (avg_t - thicknesses[ia])
        colors_list = []
        for i in range(n_sm):
            ri = res_per_pt[i]
            if uniform_color:
                colors_list.append(uniform_color)
            else:
                colors_list.append(
                    _res_color(valid_res[ri], color_scheme, ss_colors=ss_colors)
                )
        pd = _build_colored_ribbon_pd(
            sm_pts,
            sm_nms,
            widths.tolist(),
            colors_list,
            half_thickness=thicknesses.tolist(),
            cs_sides=cs_sides,
        )
        if pd:
            actors.append(_colored_ribbon_actor(pd))
    return actors


# ---------------------------------------------------------------------------
# Tube SS
# ---------------------------------------------------------------------------


def make_tube_ss_actor(
    structure,
    chains=None,
    color_scheme="ss",
    uniform_color=None,
    helix_w=2.0,
    sheet_w=2.5,
    coil_w=0.5,
    quality=3,
    ss_colors=None,
):
    """One continuous tube per chain varying radius/color by SS type."""
    actors = []
    chain_ids = chains or list(structure.chains.keys())
    qp = QUALITY_PRESETS.get(quality, QUALITY_PRESETS[3])
    smooth_factor = qp[3] * 2
    tube_sides = qp[2]
    for cid in chain_ids:
        residues = structure.chains.get(cid, [])
        valid_res, ca_all = [], []
        for r in residues:
            if r.ca_coord is not None:
                valid_res.append(r)
                ca_all.append(r.ca_coord)
        nr = len(ca_all)
        if nr < 2:
            continue
        sm_pts = _smooth_coords(ca_all, smooth_factor)
        n_sm = len(sm_pts)
        ss_per_pt, res_per_pt = [], []
        for ri in range(nr):
            count = smooth_factor if ri < nr - 1 else 1
            for _ in range(count):
                ss_per_pt.append(valid_res[ri].ss)
                res_per_pt.append(ri)
        segments = []
        seg_start = 0
        for j in range(1, n_sm):
            if ss_per_pt[j] != ss_per_pt[j - 1]:
                segments.append((seg_start, j - 1, ss_per_pt[seg_start]))
                seg_start = j
        segments.append((seg_start, n_sm - 1, ss_per_pt[seg_start]))
        radii = np.empty(n_sm)
        for i in range(n_sm):
            ss = ss_per_pt[i]
            if ss in HELIX_SS_TYPES:
                radii[i] = helix_w * 0.35
            elif ss == "E":
                radii[i] = sheet_w * 0.30
            else:
                radii[i] = coil_w * 0.4
        arrow_r = sheet_w * 0.50
        for s_s, s_e, ss in segments:
            if ss != "E":
                continue
            seg_len = s_e - s_s + 1
            arrow_start = s_s + max(1, int(seg_len * 0.70))
            for i in range(arrow_start, s_e + 1):
                frac = (i - arrow_start) / max(1, s_e - arrow_start)
                radii[i] = arrow_r * (1.0 - frac)
        transition = max(3, smooth_factor)
        boundaries = [j for j in range(1, n_sm) if ss_per_pt[j] != ss_per_pt[j - 1]]
        for bnd in boundaries:
            for d in range(transition):
                t = 1.0 - d / transition
                alpha = 0.5 * t
                ib = bnd - 1 - d
                ia = bnd + d
                if 0 <= ib < n_sm and 0 <= ia < n_sm:
                    avg = 0.5 * (radii[ib] + radii[ia])
                    radii[ib] += alpha * (avg - radii[ib])
                    radii[ia] += alpha * (avg - radii[ia])
        colors_list = []
        for i in range(n_sm):
            ri = res_per_pt[i]
            if uniform_color:
                colors_list.append(uniform_color)
            else:
                colors_list.append(
                    _res_color(valid_res[ri], color_scheme, ss_colors=ss_colors)
                )
        a = _make_colored_tube_actor(
            sm_pts, radii.tolist(), colors_list, tube_sides=tube_sides
        )
        if a:
            actors.append(a)
    return actors


# ---------------------------------------------------------------------------
# Backbone
# ---------------------------------------------------------------------------


def make_backbone_actor(
    structure,
    chains=None,
    color_scheme="chain",
    uniform_color=None,
    radius=0.3,
    quality=3,
):
    actors = []
    chain_ids = chains or list(structure.chains.keys())
    qp = QUALITY_PRESETS.get(quality, QUALITY_PRESETS[3])
    for idx, cid in enumerate(chain_ids):
        residues = structure.chains.get(cid, [])
        ca = [r.ca_coord for r in residues if r.ca_coord is not None]
        valid_res = [r for r in residues if r.ca_coord is not None]
        if len(ca) < 2:
            continue
        if uniform_color:
            cf = _rgb_f(uniform_color)
            a = _make_tube_actor(ca, radius, cf, tube_sides=qp[2])
            if a:
                actors.append(a)
        elif color_scheme in ("ss", "residue_nature"):
            # Per-residue colored segments
            colors_list = []
            for r in valid_res:
                if color_scheme == "ss":
                    colors_list.append(SS_COLORS.get(r.ss, SS_COLORS["DEFAULT"]))
                else:
                    cat = RESIDUE_NATURE.get(r.name, "other")
                    colors_list.append(RESIDUE_NATURE_COLORS[cat])
            a = _make_colored_tube_actor(
                ca, [radius] * len(ca), colors_list, tube_sides=qp[2]
            )
            if a:
                actors.append(a)
        elif color_scheme == "chain":
            cf = _rgb_f(CHAIN_PALETTE[idx % len(CHAIN_PALETTE)])
            a = _make_tube_actor(ca, radius, cf, tube_sides=qp[2])
            if a:
                actors.append(a)
        elif color_scheme == "element":
            # Use CA element color (mostly gray)
            cf = _rgb_f(ELEMENT_COLORS.get("C", ELEMENT_COLORS["DEFAULT"]))
            a = _make_tube_actor(ca, radius, cf, tube_sides=qp[2])
            if a:
                actors.append(a)
        else:
            cf = (1, 1, 1)
            a = _make_tube_actor(ca, radius, cf, tube_sides=qp[2])
            if a:
                actors.append(a)
    return actors


# ---------------------------------------------------------------------------
# Surface
# ---------------------------------------------------------------------------


def make_surface_actor(
    atoms,
    color=(200, 200, 255),
    opacity=0.5,
    resolution=64,
    radius=0.12,
    color_scheme=None,
    uniform_color=None,
    carbon_color=None,
    ss_map=None,
):
    try:
        from vtkmodules.vtkFiltersCore import vtkContourFilter
        from vtkmodules.vtkImagingHybrid import vtkGaussianSplatter
    except ImportError:
        return None
    if not atoms:
        return None
    res = max(16, min(128, int(resolution)))
    rad = max(0.02, min(0.5, float(radius)))
    pts = vtkPoints()
    for a in atoms:
        pts.InsertNextPoint(*a.coord)
    pd = vtkPolyData()
    pd.SetPoints(pts)
    splat = vtkGaussianSplatter()
    splat.SetInputData(pd)
    splat.SetSampleDimensions(res, res, res)
    splat.SetRadius(rad)
    splat.SetExponentFactor(-5)
    splat.ScalarWarpingOff()
    splat.Update()
    contour = vtkContourFilter()
    contour.SetInputConnection(splat.GetOutputPort())
    contour.SetValue(0, 0.01)
    contour.Update()
    surf_pd = contour.GetOutput()
    use_scheme = color_scheme and color_scheme != "uniform"
    if use_scheme and surf_pd.GetNumberOfPoints() > 0:
        # Color surface vertices by nearest atom
        coords = np.array([a.coord for a in atoms])
        n_verts = surf_pd.GetNumberOfPoints()
        vert_colors = vtkUnsignedCharArray()
        vert_colors.SetNumberOfComponents(3)
        vert_colors.SetName("Colors")
        for vi in range(n_verts):
            vp = np.array(surf_pd.GetPoint(vi))
            dists = np.sum((coords - vp) ** 2, axis=1)
            nearest = atoms[int(np.argmin(dists))]
            c = _atom_color(nearest, color_scheme, carbon_color, ss_map=ss_map)
            vert_colors.InsertNextTuple3(*c)
        surf_pd.GetPointData().SetScalars(vert_colors)
    mapper = vtkPolyDataMapper()
    mapper.SetInputData(surf_pd)
    if use_scheme and surf_pd.GetNumberOfPoints() > 0:
        mapper.ScalarVisibilityOn()
    else:
        mapper.ScalarVisibilityOff()
    actor = vtkActor()
    actor.SetMapper(mapper)
    p = actor.GetProperty()
    if not use_scheme:
        p.SetColor(*_rgb_f(color))
    p.SetOpacity(opacity)
    p.SetSpecular(0.3)
    p.SetSpecularPower(20)
    return actor


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------


def _smooth_coords(coords, factor=3):
    pts = np.array(coords)
    if len(pts) < 3:
        return pts
    result = []
    for i in range(len(pts) - 1):
        p0 = pts[max(i - 1, 0)]
        p1 = pts[i]
        p2 = pts[min(i + 1, len(pts) - 1)]
        p3 = pts[min(i + 2, len(pts) - 1)]
        for t_idx in range(factor):
            t = t_idx / factor
            tt = t * t
            ttt = tt * t
            q = 0.5 * (
                (2 * p1)
                + (-p0 + p2) * t
                + (2 * p0 - 5 * p1 + 4 * p2 - p3) * tt
                + (-p0 + 3 * p1 - 3 * p2 + p3) * ttt
            )
            result.append(q)
    result.append(pts[-1])
    return result


def _make_tube_actor(coords, radius, color_float, tube_sides=12):
    pts = vtkPoints()
    n = len(coords)
    for c in coords:
        pts.InsertNextPoint(*c)
    ids = vtkIdList()
    for i in range(n):
        ids.InsertNextId(i)
    cells = vtkCellArray()
    cells.InsertNextCell(ids)
    pd = vtkPolyData()
    pd.SetPoints(pts)
    pd.SetLines(cells)
    tube = vtkTubeFilter()
    tube.SetInputData(pd)
    tube.SetRadius(radius)
    tube.SetNumberOfSides(tube_sides)
    tube.CappingOn()
    tube.Update()
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())
    mapper.ScalarVisibilityOff()
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color_float)
    actor.GetProperty().SetSpecular(0.3)
    actor.GetProperty().SetSpecularPower(20)
    return actor


def _compute_ribbon_normals(ca_coords, o_coords):
    n = len(ca_coords)
    normals = []
    for i in range(n):
        if i == 0:
            t = np.asarray(ca_coords[1]) - np.asarray(ca_coords[0])
        elif i == n - 1:
            t = np.asarray(ca_coords[-1]) - np.asarray(ca_coords[-2])
        else:
            t = np.asarray(ca_coords[i + 1]) - np.asarray(ca_coords[i - 1])
        tl = np.linalg.norm(t)
        t = t / tl if tl > 1e-12 else np.array([0.0, 0.0, 1.0])
        if o_coords[i] is not None:
            co = np.asarray(o_coords[i]) - np.asarray(ca_coords[i])
        else:
            co = (
                np.array([1.0, 0.0, 0.0])
                if abs(t[0]) < 0.9
                else np.array([0.0, 1.0, 0.0])
            )
        co = co - np.dot(co, t) * t
        cl = np.linalg.norm(co)
        if cl < 1e-6:
            co = (
                np.array([1.0, 0.0, 0.0])
                if abs(t[0]) < 0.9
                else np.array([0.0, 1.0, 0.0])
            )
            co = co - np.dot(co, t) * t
            cl = np.linalg.norm(co)
        normals.append(co / (cl + 1e-12))
    for i in range(1, n):
        if np.dot(normals[i], normals[i - 1]) < 0:
            normals[i] = -normals[i]
    return normals


def _smooth_coords_and_normals(coords, normals, factor=3):
    pts = np.array(coords)
    nms = np.array(normals)
    if len(pts) < 3:
        return list(pts), list(nms)
    r_pts, r_nms = [], []
    for i in range(len(pts) - 1):
        p0, p1 = pts[max(i - 1, 0)], pts[i]
        p2 = pts[min(i + 1, len(pts) - 1)]
        p3 = pts[min(i + 2, len(pts) - 1)]
        n0, n1 = nms[max(i - 1, 0)], nms[i]
        n2 = nms[min(i + 1, len(nms) - 1)]
        n3 = nms[min(i + 2, len(nms) - 1)]
        for t_idx in range(factor):
            t = t_idx / factor
            tt = t * t
            ttt = tt * t
            q = 0.5 * (
                (2 * p1)
                + (-p0 + p2) * t
                + (2 * p0 - 5 * p1 + 4 * p2 - p3) * tt
                + (-p0 + 3 * p1 - 3 * p2 + p3) * ttt
            )
            r_pts.append(q)
            nm = 0.5 * (
                (2 * n1)
                + (-n0 + n2) * t
                + (2 * n0 - 5 * n1 + 4 * n2 - n3) * tt
                + (-n0 + 3 * n1 - 3 * n2 + n3) * ttt
            )
            nm_len = np.linalg.norm(nm)
            r_nms.append(nm / (nm_len + 1e-12))
    r_pts.append(pts[-1])
    r_nms.append(nms[-1])
    return r_pts, r_nms


def _build_colored_ribbon_pd(
    points, normals, widths, colors, half_thickness=0.15, cs_sides=10
):
    n = len(points)
    if n < 2:
        return None
    if isinstance(half_thickness, (int, float)):
        thicknesses = [half_thickness] * n
    else:
        thicknesses = half_thickness
    vtk_pts = vtkPoints()
    vtk_colors = vtkUnsignedCharArray()
    vtk_colors.SetNumberOfComponents(3)
    vtk_colors.SetName("Colors")
    polys = vtkCellArray()
    angles = [2.0 * math.pi * k / cs_sides for k in range(cs_sides)]
    cos_a = [math.cos(a) for a in angles]
    sin_a = [math.sin(a) for a in angles]
    for i in range(n):
        p = np.asarray(points[i], dtype=np.float64)
        N = np.asarray(normals[i], dtype=np.float64)
        if i == 0:
            T = np.asarray(points[1]) - p
        elif i == n - 1:
            T = p - np.asarray(points[i - 1])
        else:
            T = np.asarray(points[i + 1]) - np.asarray(points[i - 1])
        tl = np.linalg.norm(T)
        T = T / tl if tl > 1e-12 else np.array([0.0, 0.0, 1.0])
        B = np.cross(T, N)
        bl = np.linalg.norm(B)
        B = B / bl if bl > 1e-12 else np.array([0.0, 1.0, 0.0])
        w = widths[i]
        th = thicknesses[i]
        r, g, b = int(colors[i][0]), int(colors[i][1]), int(colors[i][2])
        if w < 0.001 and th < 0.001:
            for _ in range(cs_sides):
                vtk_pts.InsertNextPoint(*p)
                vtk_colors.InsertNextTuple3(r, g, b)
        else:
            for k in range(cs_sides):
                pt = p + w * cos_a[k] * N + th * sin_a[k] * B
                vtk_pts.InsertNextPoint(*pt)
                vtk_colors.InsertNextTuple3(r, g, b)
    S = cs_sides
    for i in range(n - 1):
        base = i * S
        nb = (i + 1) * S
        for k in range(S):
            k1 = (k + 1) % S
            polys.InsertNextCell(3)
            polys.InsertCellPoint(base + k)
            polys.InsertCellPoint(nb + k)
            polys.InsertCellPoint(nb + k1)
            polys.InsertNextCell(3)
            polys.InsertCellPoint(base + k)
            polys.InsertCellPoint(nb + k1)
            polys.InsertCellPoint(base + k1)
    pd = vtkPolyData()
    pd.SetPoints(vtk_pts)
    pd.SetPolys(polys)
    pd.GetPointData().SetScalars(vtk_colors)
    return pd


def _colored_ribbon_actor(pd):
    norm_f = vtkPolyDataNormals()
    norm_f.SetInputData(pd)
    norm_f.ComputePointNormalsOn()
    norm_f.ConsistencyOn()
    norm_f.AutoOrientNormalsOn()
    norm_f.SplittingOff()
    norm_f.Update()
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(norm_f.GetOutputPort())
    mapper.ScalarVisibilityOn()
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetSpecular(0.3)
    actor.GetProperty().SetSpecularPower(20)
    return actor


def _make_colored_tube_actor(coords, radii, colors, tube_sides=12):
    n = len(coords)
    if n < 2:
        return None
    pts = vtkPoints()
    vtk_radii = vtkFloatArray()
    vtk_radii.SetName("Radii")
    vtk_radii.SetNumberOfComponents(1)
    vtk_colors = vtkUnsignedCharArray()
    vtk_colors.SetNumberOfComponents(3)
    vtk_colors.SetName("Colors")
    for i in range(n):
        pts.InsertNextPoint(*coords[i])
        vtk_radii.InsertNextValue(radii[i])
        r, g, b = int(colors[i][0]), int(colors[i][1]), int(colors[i][2])
        vtk_colors.InsertNextTuple3(r, g, b)
    ids = vtkIdList()
    for i in range(n):
        ids.InsertNextId(i)
    cells = vtkCellArray()
    cells.InsertNextCell(ids)
    pd = vtkPolyData()
    pd.SetPoints(pts)
    pd.SetLines(cells)
    pd.GetPointData().SetScalars(vtk_radii)
    pd.GetPointData().AddArray(vtk_colors)
    tube = vtkTubeFilter()
    tube.SetInputData(pd)
    tube.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tube.SetNumberOfSides(tube_sides)
    tube.CappingOn()
    tube.Update()
    output = tube.GetOutput()
    output.GetPointData().SetActiveScalars("Colors")
    mapper = vtkPolyDataMapper()
    mapper.SetInputData(output)
    mapper.ScalarVisibilityOn()
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetSpecular(0.3)
    actor.GetProperty().SetSpecularPower(20)
    return actor
