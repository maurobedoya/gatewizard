# gatewizard/gui/widgets/vtk_frame.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
VTK offscreen rendering widget for tkinter / customtkinter.

Renders a VTK scene offscreen and blits it to a tkinter Canvas via PIL.
Supports mouse-driven rotation, panning, zooming, depth cueing (fog),
and high-resolution image export.
"""

import math
import tkinter as tk
from typing import Callable, Optional

import numpy as np

try:
    from vtkmodules.vtkRenderingCore import (
        vtkRenderer,
        vtkRenderWindow,
        vtkActor,
        vtkWindowToImageFilter,
    )
    from vtkmodules.vtkRenderingOpenGL2 import (
        vtkOpenGLRenderer,
    )  # noqa: force GL backend
    from vtkmodules.util.numpy_support import vtk_to_numpy

    VTK_AVAILABLE = True
except ImportError:
    VTK_AVAILABLE = False

try:
    from PIL import Image, ImageTk

    PIL_AVAILABLE = True
except ImportError:
    PIL_AVAILABLE = False


class VTKFrame(tk.Frame):
    """Renders a VTK scene offscreen and blits it to a tkinter Canvas."""

    def __init__(self, master, width: int = 900, height: int = 700, **kw):
        super().__init__(master, bg="#212121", **kw)
        self._vw = width
        self._vh = height

        self.renderer = vtkRenderer()
        self.renderer.SetBackground(0.18, 0.18, 0.18)
        self.renderer.GradientBackgroundOff()

        self.render_window = vtkRenderWindow()
        self.render_window.SetOffScreenRendering(1)
        self.render_window.SetSize(width, height)
        self.render_window.AddRenderer(self.renderer)
        self.render_window.SetMultiSamples(4)

        # Overlay renderer for corner axes widget
        self._axes_renderer = vtkRenderer()
        self._axes_renderer.SetViewport(0.8, 0.0, 1.0, 0.2)
        self._axes_renderer.SetBackground(0, 0, 0)
        self._axes_renderer.SetBackgroundAlpha(0.0)
        self._axes_renderer.SetLayer(1)
        self._axes_renderer.InteractiveOff()

        # Overlay renderer for labels (always on top, full viewport)
        self._label_renderer = vtkRenderer()
        self._label_renderer.SetViewport(0.0, 0.0, 1.0, 1.0)
        self._label_renderer.SetBackgroundAlpha(0.0)
        self._label_renderer.SetLayer(2)
        self._label_renderer.InteractiveOff()
        # Share the main camera so labels track the same view
        self._label_renderer.SetActiveCamera(self.renderer.GetActiveCamera())

        self.render_window.SetNumberOfLayers(3)
        self.render_window.AddRenderer(self._axes_renderer)
        self.render_window.AddRenderer(self._label_renderer)
        self._axes_actor = None  # vtkAxesActor when visible
        self._axes_mode = None  # None, 'corner', 'center', 'origin'

        self.w2i = vtkWindowToImageFilter()
        self.w2i.SetInput(self.render_window)
        self.w2i.SetInputBufferTypeToRGB()

        self.w2i_depth = vtkWindowToImageFilter()
        self.w2i_depth.SetInput(self.render_window)
        self.w2i_depth.SetInputBufferTypeToZBuffer()

        self.fog_density: float = 0.0

        self.canvas = tk.Canvas(
            self, width=width, height=height, bg="#212121", highlightthickness=0
        )
        self.canvas.pack(fill="both", expand=True)

        self._last_x = 0
        self._last_y = 0
        self._render_pending = False
        self._drag_started = False
        self._pick_callback: Optional[Callable] = None
        self._right_click_callback: Optional[Callable] = None
        self._resize_callback: Optional[Callable] = None
        self._post_render_callback: Optional[Callable] = None
        self._camera_locked = False

        self.canvas.bind("<ButtonPress-1>", self._on_press)
        self.canvas.bind("<B1-Motion>", self._on_rotate)
        self.canvas.bind("<ButtonRelease-1>", self._on_click)
        self.canvas.bind("<ButtonPress-3>", self._on_press)
        self.canvas.bind("<B3-Motion>", self._on_pan)
        self.canvas.bind("<ButtonRelease-3>", self._on_right_click)
        self.canvas.bind("<ButtonPress-2>", self._on_press)
        self.canvas.bind("<B2-Motion>", self._on_pan)
        self.canvas.bind("<MouseWheel>", self._on_zoom)
        self.canvas.bind("<Button-4>", self._on_zoom)
        self.canvas.bind("<Button-5>", self._on_zoom)
        self.canvas.bind("<Configure>", self._on_resize)

    # -- public helpers ------------------------------------------------

    def add_actor(self, actor):
        self.renderer.AddActor(actor)

    def remove_actor(self, actor):
        self.renderer.RemoveActor(actor)

    def clear_actors(self):
        self.renderer.RemoveAllViewProps()
        # Re-add persistent overlays (axes in center/origin, ref-lines)
        if self._axes_mode in ("center", "origin") and self._axes_actor:
            self.renderer.AddActor(self._axes_actor)
        for a in getattr(self, "_ref_line_actors", []):
            self.renderer.AddActor(a)

    def lock_camera(self):
        """Block mouse-driven camera changes (rotate/pan/zoom)."""
        self._camera_locked = True

    def unlock_camera(self):
        """Re-enable mouse-driven camera changes."""
        self._camera_locked = False

    def reset_camera(self):
        self.renderer.ResetCamera()
        self.renderer.GetActiveCamera().Zoom(1.0)

    # -- Axes widget ---------------------------------------------------

    def set_axes(self, mode=None, center=None):
        """Show/hide orientation axes.

        *mode*: ``'corner'``, ``'center'``, ``'origin'`` or ``None`` (hide).
        *center*: (x, y, z) centre of the protein (used when mode='center').
        """
        from vtkmodules.vtkRenderingAnnotation import vtkAxesActor

        # Clean up old actor
        if self._axes_actor:
            if self._axes_mode == "corner":
                self._axes_renderer.RemoveActor(self._axes_actor)
            else:
                self.renderer.RemoveActor(self._axes_actor)
            self._axes_actor = None
        self._axes_mode = mode
        if mode is None:
            return

        ax = vtkAxesActor()
        ax.SetShaftTypeToCylinder()
        ax.SetTipTypeToCone()
        ax.SetConeResolution(50)  # smooth cone tips
        ax.SetCylinderResolution(20)  # smooth cylinder shafts

        def _configure_labels(ax, font_size):
            for getter in (
                ax.GetXAxisCaptionActor2D,
                ax.GetYAxisCaptionActor2D,
                ax.GetZAxisCaptionActor2D,
            ):
                cap = getter()
                tp = cap.GetCaptionTextProperty()
                tp.SetFontSize(font_size)
                tp.BoldOn()
                tp.ItalicOff()
                tp.ShadowOn()
                # Fixed screen-size text, not scaled by distance
                cap.SetVisibility(1)

        if mode == "corner":
            ax.SetTotalLength(1.5, 1.5, 1.5)
            ax.SetCylinderRadius(0.08)
            ax.SetConeRadius(0.40)
            ax.SetNormalizedLabelPosition(1.3, 1.3, 1.3)
            _configure_labels(ax, 14)
            self._axes_renderer.RemoveAllViewProps()
            self._axes_renderer.AddActor(ax)
            self._sync_axes_camera()
        elif mode in ("center", "origin"):
            length = 8.0
            ax.SetTotalLength(length, length, length)
            ax.SetCylinderRadius(0.05)
            ax.SetConeRadius(0.40)
            ax.SetNormalizedLabelPosition(1.25, 1.25, 1.25)
            _configure_labels(ax, 10)
            for getter in (
                ax.GetXAxisCaptionActor2D,
                ax.GetYAxisCaptionActor2D,
                ax.GetZAxisCaptionActor2D,
            ):
                cap = getter()
                cap.LeaderOff()
                cap.SetAttachmentPoint(0, 0, 0)
                cap.GetTextActor().SetTextScaleModeToNone()
            # Position via transform so the whole actor (shafts + labels) moves
            from vtkmodules.vtkCommonTransforms import vtkTransform

            t = vtkTransform()
            if mode == "center" and center is not None:
                t.Translate(*center)
            else:
                t.Translate(0, 0, 0)
            ax.SetUserTransform(t)
            self.renderer.AddActor(ax)

        self._axes_actor = ax

    def _sync_axes_camera(self):
        """Copy main camera orientation to corner-axes camera."""
        if self._axes_mode != "corner" or not self._axes_actor:
            return
        main_cam = self.renderer.GetActiveCamera()
        ax_cam = self._axes_renderer.GetActiveCamera()
        ax_cam.SetViewUp(main_cam.GetViewUp())
        # Compute direction from main camera
        fp = main_cam.GetFocalPoint()
        pos = main_cam.GetPosition()
        d = [pos[i] - fp[i] for i in range(3)]
        norm = math.sqrt(sum(x * x for x in d)) or 1.0
        d = [x / norm for x in d]
        # Place axes camera looking at origin from same direction
        ax_cam.SetFocalPoint(0, 0, 0)
        ax_cam.SetPosition(d[0] * 5, d[1] * 5, d[2] * 5)
        ax_cam.SetParallelProjection(1)
        ax_cam.SetParallelScale(2.0)

    # -- Reference grid lines ------------------------------------------

    def set_reference_lines(self, show, length=200.0):
        """Show/hide X/Y/Z reference lines through the origin."""
        from vtkmodules.vtkFiltersSources import vtkLineSource
        from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper

        # Remove existing
        if hasattr(self, "_ref_line_actors"):
            for a in self._ref_line_actors:
                self.renderer.RemoveActor(a)
        self._ref_line_actors = []
        if not show:
            return

        half = length / 2.0
        colors = [(0.9, 0.2, 0.2), (0.2, 0.9, 0.2), (0.2, 0.2, 0.9)]
        endpoints = [
            ((-half, 0, 0), (half, 0, 0)),
            ((0, -half, 0), (0, half, 0)),
            ((0, 0, -half), (0, 0, half)),
        ]
        for (p1, p2), color in zip(endpoints, colors):
            src = vtkLineSource()
            src.SetPoint1(*p1)
            src.SetPoint2(*p2)
            mapper = vtkPolyDataMapper()
            mapper.SetInputConnection(src.GetOutputPort())
            actor = vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(*color)
            actor.GetProperty().SetLineWidth(1.5)
            actor.GetProperty().SetOpacity(0.6)
            self.renderer.AddActor(actor)
            self._ref_line_actors.append(actor)

    def render(self):
        self._sync_axes_camera()
        need_fog = self.fog_density > 0.01
        if need_fog:
            # Capture Z-buffer from main scene only (overlay layers off)
            # so that overlay actors don't overwrite the depth buffer.
            self._axes_renderer.DrawOff()
            self._label_renderer.DrawOff()
            self.render_window.Render()
            self.w2i_depth.Modified()
            self.w2i_depth.Update()
            self._axes_renderer.DrawOn()
            self._label_renderer.DrawOn()
        # Full render with all layers for the RGB image
        self.render_window.Render()
        self.w2i.Modified()
        self.w2i.Update()
        vtk_image = self.w2i.GetOutput()
        dims = vtk_image.GetDimensions()
        if dims[0] < 1 or dims[1] < 1:
            return
        scalars = vtk_image.GetPointData().GetScalars()
        if scalars is None:
            return
        arr = vtk_to_numpy(scalars).reshape(dims[1], dims[0], 3)
        arr = np.flipud(arr)

        if need_fog:
            arr = self._apply_fog(arr, dims)

        img = Image.fromarray(arr, "RGB")
        self._photo = ImageTk.PhotoImage(img)
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, anchor="nw", image=self._photo)
        if self._post_render_callback:
            self._post_render_callback()

    def focus_on_point(self, point, distance: float = 20.0):
        cam = self.renderer.GetActiveCamera()
        cam.SetFocalPoint(*point)
        pos = list(cam.GetPosition())
        fp = list(cam.GetFocalPoint())
        direction = [pos[i] - fp[i] for i in range(3)]
        norm = math.sqrt(sum(d * d for d in direction)) or 1
        direction = [d / norm for d in direction]
        cam.SetPosition(*[fp[i] + direction[i] * distance for i in range(3)])
        cam.OrthogonalizeViewUp()
        self.renderer.ResetCameraClippingRange()
        self._schedule_render()

    def render_to_image(self, scale: int = 1, transparent: bool = False):
        """Render the current scene to a PIL Image."""
        w, h = int(self._vw * scale), int(self._vh * scale)
        old_size = self.render_window.GetSize()
        self.render_window.SetSize(w, h)
        if transparent:
            old_bg = self.renderer.GetBackground()
            self.render_window.SetAlphaBitPlanes(1)
            self.renderer.SetBackgroundAlpha(0.0)
        self.render_window.Render()
        w2i = vtkWindowToImageFilter()
        w2i.SetInput(self.render_window)
        if transparent:
            w2i.SetInputBufferTypeToRGBA()
        else:
            w2i.SetInputBufferTypeToRGB()
        w2i.Update()
        vtk_img = w2i.GetOutput()
        dims = vtk_img.GetDimensions()
        scalars = vtk_img.GetPointData().GetScalars()
        if transparent:
            arr = vtk_to_numpy(scalars).reshape(dims[1], dims[0], 4)
            img = Image.fromarray(np.flipud(arr), "RGBA")
            self.render_window.SetAlphaBitPlanes(0)
            self.renderer.SetBackgroundAlpha(1.0)
            self.renderer.SetBackground(*old_bg)
        else:
            arr = vtk_to_numpy(scalars).reshape(dims[1], dims[0], 3)
            img = Image.fromarray(np.flipud(arr), "RGB")
        self.render_window.SetSize(*old_size)
        self.render()
        return img

    # -- internal ------------------------------------------------------

    def _apply_fog(self, arr, dims):
        # Z-buffer was already captured in render() from main scene only
        z_img = self.w2i_depth.GetOutput()
        z_scalars = z_img.GetPointData().GetScalars()
        if z_scalars is None:
            return arr
        z_arr = vtk_to_numpy(z_scalars).reshape(dims[1], dims[0])
        z_arr = np.flipud(z_arr)
        bg = self.renderer.GetBackground()
        bg_rgb = np.array([bg[0] * 255, bg[1] * 255, bg[2] * 255], dtype=np.float32)
        obj_mask = z_arr < 0.999
        if not np.any(obj_mask):
            return arr
        z_obj = z_arr[obj_mask]
        z_min, z_max = z_obj.min(), z_obj.max()
        if z_max - z_min > 1e-6:
            depth_norm = (z_arr - z_min) / (z_max - z_min)
        else:
            depth_norm = np.zeros_like(z_arr)
        depth_norm = np.clip(depth_norm, 0, 1)
        depth_norm[~obj_mask] = 0.0
        fog = (depth_norm * self.fog_density)[..., np.newaxis]
        arr = arr.astype(np.float32)
        arr = arr * (1.0 - fog) + bg_rgb * fog
        return np.clip(arr, 0, 255).astype(np.uint8)

    def _schedule_render(self):
        if not self._render_pending:
            self._render_pending = True
            self.after(16, self._do_render)

    def _do_render(self):
        self._render_pending = False
        self.render()

    def _on_press(self, event):
        self._last_x = event.x
        self._last_y = event.y
        self._drag_started = False

    def _on_click(self, event):
        if self._drag_started:
            return
        if self._pick_callback:
            self._pick_callback(event.x, event.y)

    def _on_right_click(self, event):
        """Right-click without drag → context menu callback."""
        if self._drag_started:
            return
        if hasattr(self, "_right_click_callback") and self._right_click_callback:
            self._right_click_callback(event.x, event.y, event)

    def _on_rotate(self, event):
        if self._camera_locked:
            return
        self._drag_started = True
        dx = event.x - self._last_x
        dy = event.y - self._last_y
        cam = self.renderer.GetActiveCamera()
        cam.Azimuth(-dx * 0.4)
        cam.Elevation(dy * 0.4)
        cam.OrthogonalizeViewUp()
        self._last_x = event.x
        self._last_y = event.y
        self._schedule_render()

    def _on_pan(self, event):
        if self._camera_locked:
            return
        dx = event.x - self._last_x
        dy = event.y - self._last_y
        cam = self.renderer.GetActiveCamera()
        fp = list(cam.GetFocalPoint())
        pos = list(cam.GetPosition())
        dist = cam.GetDistance()
        scale = dist * 0.0008
        up = list(cam.GetViewUp())
        view_dir = [fp[i] - pos[i] for i in range(3)]
        right = [
            view_dir[1] * up[2] - view_dir[2] * up[1],
            view_dir[2] * up[0] - view_dir[0] * up[2],
            view_dir[0] * up[1] - view_dir[1] * up[0],
        ]
        norm = math.sqrt(sum(r * r for r in right)) or 1
        right = [r / norm for r in right]
        move = [(-dx * right[i] + dy * up[i]) * scale for i in range(3)]
        cam.SetFocalPoint(*[fp[i] + move[i] for i in range(3)])
        cam.SetPosition(*[pos[i] + move[i] for i in range(3)])
        self._last_x = event.x
        self._last_y = event.y
        self._schedule_render()

    def _on_zoom(self, event):
        if self._camera_locked:
            return
        if event.delta:
            factor = 1.15 if event.delta > 0 else 1 / 1.15
        else:
            factor = 1.15 if event.num == 4 else 1 / 1.15
        cam = self.renderer.GetActiveCamera()
        if cam.GetParallelProjection():
            cam.SetParallelScale(cam.GetParallelScale() / factor)
        else:
            # Dolly: move camera closer/farther — keeps view angle untouched
            # so the perspective slider stays consistent.
            fp = cam.GetFocalPoint()
            pos = cam.GetPosition()
            d = [pos[i] - fp[i] for i in range(3)]
            cam.SetPosition(*[fp[i] + d[i] / factor for i in range(3)])
        self.renderer.ResetCameraClippingRange()
        self._schedule_render()

    def _on_resize(self, event):
        w, h = event.width, event.height
        if w < 10 or h < 10:
            return
        self._vw = w
        self._vh = h
        self.render_window.SetSize(w, h)
        if self._resize_callback:
            self._resize_callback(w, h)
        self._schedule_render()
