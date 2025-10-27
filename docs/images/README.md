# Documentation Images

This directory contains all images used in the GateWizard documentation.

## Directory Structure

```
docs/images/
├── api/              # Images for API reference documentation
│   ├── propka_titration_curves.png
│   ├── propka_pka_distribution.png
│   └── ...
├── user-guide/       # Images for user guide
│   └── ...
├── installation/     # Images for installation guide
│   └── ...
└── troubleshooting/  # Images for troubleshooting
    └── ...
```

## Guidelines

### File Naming Convention
- Use lowercase with underscores: `module_feature_description.png`
- Be descriptive: `propka_titration_curves.png` instead of `figure1.png`
- Include module/section prefix for organization

### Image Formats
- **PNG**: Preferred for plots, screenshots, diagrams (lossless compression)
- **SVG**: For vector graphics, logos, diagrams (scalable)
- **JPG**: For photos only (lossy compression)

### Image Dimensions
- **Plots/Figures**: 1200-2400px wide (for high-DPI displays)
- **Screenshots**: Actual size or scaled to fit documentation width
- **Icons**: 64px, 128px, or 256px square

### DPI Requirements
- Save matplotlib figures at **300 DPI** for publication quality
- Use `plt.savefig('filename.png', dpi=300, bbox_inches='tight')`

## Referencing Images in Markdown

### Relative Path Format
From `docs/api/preparation.md`:
```markdown
![Alt Text](../images/api/propka_titration_curves.png)
```

From `docs/user-guide.md`:
```markdown
![Alt Text](images/user-guide/example.png)
```

### With Caption
```markdown
![Titration Curves](../images/api/propka_titration_curves.png)

*Figure: Description of the figure with relevant details.*
```

### Inline with Width Control
```markdown
<img src="../images/api/propka_titration_curves.png" alt="Titration Curves" width="800">
```

## Image Optimization

Before committing images:
- Compress PNG files: `optipng -o7 *.png` or `pngcrush`
- Keep file sizes < 500KB when possible
- For large plots, consider linking to full-resolution version

## Adding New Images

1. Generate your image/figure
2. Save to appropriate subdirectory in `docs/images/`
3. Follow naming convention
4. Reference in documentation using relative path
5. Add caption/description below image
6. Commit both image and documentation update
