# Quick Reference: Adding Figures to Documentation

## File Structure
```
docs/
├── api/preparation.md                      (your markdown file)
└── images/
    └── api/
        └── propka_titration_curves.png     (your image)
```

## In Your Python Code
```python
# Generate figure
plt.savefig('docs/images/api/propka_titration_curves.png', 
            dpi=300, bbox_inches='tight')
```

## In Your Markdown (docs/api/preparation.md)
```markdown
![Titration Curves](../images/api/propka_titration_curves.png)

*Figure: Description here.*
```

## Path Rules
- From `docs/api/file.md` → use `../images/api/image.png` (go up one level)
- From `docs/file.md` → use `images/api/image.png` (same level)
- From `docs/section/file.md` → use `../images/api/image.png` (go up one level)

## That's it! 🎉
