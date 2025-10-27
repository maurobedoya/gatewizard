# Image Organization Guide for GateWizard Documentation

## ✅ Directory Structure Created

```
docs/
├── images/
│   ├── README.md                    # Complete guidelines
│   ├── api/                         # API reference figures
│   │   └── PLACEHOLDER_propka_titration_curves.txt
│   ├── user-guide/                  # User guide screenshots
│   ├── installation/                # Installation process images
│   ├── troubleshooting/             # Troubleshooting diagrams
│   └── examples/                    # Example outputs
```

## 📝 How to Add Figures to Documentation

### 1. Generate Your Figure

When creating figures with matplotlib, save them directly to the docs/images directory:

```python
import matplotlib.pyplot as plt

# ... your plotting code ...

# Save to appropriate subdirectory
plt.savefig('docs/images/api/propka_titration_curves.png', 
            dpi=300, bbox_inches='tight')
```

### 2. Reference in Markdown

In your markdown file (e.g., `docs/api/preparation.md`), add:

```markdown
![Titration Curves](../images/api/propka_titration_curves.png)

*Figure: Description of the figure.*
```

**Path Notes:**
- From `docs/api/preparation.md` → use `../images/api/filename.png`
- From `docs/user-guide.md` → use `images/user-guide/filename.png`
- From `docs/index.md` → use `images/filename.png`

### 3. File Naming Convention

✅ **Good:**
- `propka_titration_curves.png`
- `system_builder_workflow.png`
- `equilibration_namd_setup.png`

❌ **Avoid:**
- `figure1.png`
- `Screenshot 2025-01-15.png`
- `Untitled.png`

## 🎨 Image Best Practices

### For Matplotlib Figures
```python
fig, ax = plt.subplots(figsize=(10, 6))
# ... plotting code ...
plt.tight_layout()
plt.savefig('docs/images/api/figure_name.png', 
            dpi=300,                      # High resolution
            bbox_inches='tight',          # No whitespace
            facecolor='white',            # White background
            edgecolor='none')             # No border
```

### Recommended Dimensions
- **Standard plot:** 10×6 inches at 300 DPI = 3000×1800 pixels
- **Wide plot:** 12×6 inches at 300 DPI = 3600×1800 pixels
- **Square plot:** 8×8 inches at 300 DPI = 2400×2400 pixels

### Image Optimization
After generating, optionally compress:
```bash
# PNG optimization (lossless)
optipng -o7 docs/images/api/*.png

# Or using pngcrush
pngcrush -brute input.png output.png
```

## 📂 Organization by Section

### API Reference (`docs/images/api/`)
- Method outputs
- Class diagrams
- Algorithm visualizations
- Example results

### User Guide (`docs/images/user-guide/`)
- GUI screenshots
- Workflow diagrams
- Tutorial images

### Installation (`docs/images/installation/`)
- Installation steps
- Environment setup
- Dependency diagrams

### Troubleshooting (`docs/images/troubleshooting/`)
- Error screenshots
- Solution diagrams
- Before/after comparisons

### Examples (`docs/images/examples/`)
- Complete workflow outputs
- Sample results
- Use case demonstrations

## 🔗 Advanced Usage

### With Captions and Alt Text
```markdown
![Titration curves for ASP12, GLU13, and CYS38](../images/api/propka_titration_curves.png)

*Figure 1: Protonation state transitions across pH 2-12 for selected residues. 
The vertical dashed line indicates physiological pH (7.4). ASP12 and GLU13 
transition at low pH, while CYS38 remains deprotonated across the range.*
```

### Controlling Size with HTML
```markdown
<p align="center">
  <img src="../images/api/propka_titration_curves.png" 
       alt="Titration Curves" 
       width="800">
</p>
<p align="center"><em>Figure 1: Protein titration curves</em></p>
```

### Side-by-Side Images
```markdown
| Before | After |
|--------|-------|
| ![Before](../images/examples/before.png) | ![After](../images/examples/after.png) |
```

## 📋 Quick Checklist

When adding a new figure:

- [ ] Image saved in appropriate `docs/images/` subdirectory
- [ ] Filename uses lowercase with underscores
- [ ] Saved at 300 DPI for quality
- [ ] Referenced in markdown with correct relative path
- [ ] Alt text provided for accessibility
- [ ] Caption/description added below image
- [ ] File size reasonable (< 1MB preferred)
- [ ] Image added to git and committed

## 🔍 Current Example

Updated `docs/api/preparation.md` with:
- Image reference: `![Titration Curves](../images/api/propka_titration_curves.png)`
- Location: Line ~743
- Expected image: `docs/images/api/propka_titration_curves.png`

To create the actual image, run the plotting code from the example and save with:
```python
plt.savefig('docs/images/api/propka_titration_curves.png', dpi=300, bbox_inches='tight')
```
