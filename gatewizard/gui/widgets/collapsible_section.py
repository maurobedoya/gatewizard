"""Collapsible section widget for CustomTkinter."""

from PIL import Image, ImageDraw

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")


def _make_arrow_image(direction="right", size=12, color="white"):
    """Draw a tiny triangle arrow as a PIL image (no font needed)."""
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    m = 2  # margin
    if direction == "down":
        draw.polygon([(m, m), (size - m, m), (size // 2, size - m)], fill=color)
    else:  # right
        draw.polygon([(m, m), (size - m, size // 2), (m, size - m)], fill=color)
    return ctk.CTkImage(light_image=img, dark_image=img, size=(size, size))


class CollapsibleSection(ctk.CTkFrame):
    """A section with a clickable header that expands/collapses its content."""

    def __init__(
        self,
        parent,
        title: str,
        expanded: bool = True,
        fill_vertical: bool = False,
        **kwargs,
    ):
        super().__init__(parent, fg_color="transparent", **kwargs)

        self._expanded = expanded
        self._title = title
        self._fill_vertical = fill_vertical
        self._content_pack_kw = (
            dict(fill="both", expand=True, padx=2)
            if fill_vertical
            else dict(fill="x", padx=2)
        )

        # Pre-build both arrow images so toggling is instant
        self._img_down = _make_arrow_image("down")
        self._img_right = _make_arrow_image("right")

        # Header row (clickable)
        self._header = ctk.CTkFrame(self, fg_color="gray25", corner_radius=4, height=28)
        self._header.pack(fill="x", pady=(4, 0), padx=2)
        self._header.pack_propagate(False)

        self._arrow_label = ctk.CTkLabel(
            self._header,
            text="",
            image=self._img_down if expanded else self._img_right,
            width=16,
            anchor="w",
        )
        self._arrow_label.pack(side="left", padx=(6, 0))

        self._title_label = ctk.CTkLabel(
            self._header, text=title, font=("", 12, "bold"), anchor="w"
        )
        self._title_label.pack(side="left", padx=(2, 6))

        # Content frame
        self._content = ctk.CTkFrame(self, fg_color="transparent")
        if expanded:
            self._content.pack(**self._content_pack_kw)

        # Bind clicks on header elements
        for w in (self._header, self._arrow_label, self._title_label):
            w.bind("<Button-1>", self._toggle)

    @property
    def content(self) -> ctk.CTkFrame:
        """The frame where child widgets should be placed."""
        return self._content

    @property
    def expanded(self) -> bool:
        return self._expanded

    def _toggle(self, event=None):
        if self._expanded:
            self.collapse()
        else:
            self.expand()

    def expand(self):
        if not self._expanded:
            self._expanded = True
            self._arrow_label.configure(image=self._img_down)
            self._content.pack(**self._content_pack_kw)

    def collapse(self):
        if self._expanded:
            self._expanded = False
            self._arrow_label.configure(image=self._img_right)
            self._content.pack_forget()
