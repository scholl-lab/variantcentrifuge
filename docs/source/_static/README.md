# Static Files Directory

Place your logo file (e.g., `logo.png` or `logo.svg`) in this directory.

To configure the logo in Sphinx documentation, add this line to `conf.py`:

```python
html_logo = "_static/logo.png"
```

For the cohort report template, the logo needs to be embedded as Base64. Convert your logo to Base64 and add it to the HTML template:

```html
<img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgA..." alt="Logo" style="height: 40px;">
```