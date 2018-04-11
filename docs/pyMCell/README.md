# pyMCell Sphinx Autodoc

The documentation for pyMCell is automatically generated based off of the
docstrings using Sphinx.

To build the HTML, you'll want to activate a python virtual environment:

    python3 -m venv sphinx_env
    source sphinx_env/bin/activate

Then install the requirements:

    pip3 install -r requirements.txt

Run make:

    make html

View docs:

    firefox _build/html/index.html

Obviously, you can use a different browser if you want.
