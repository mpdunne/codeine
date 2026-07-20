import os

from importlib.metadata import PackageNotFoundError, version as get_version

project = 'Codeine'
author = 'Michael Dunne'
html_title = 'Codeine documentation'

try:
    release = version = get_version('codeine')
except PackageNotFoundError:
    release = version = 'x.x.x'

rtd_version = os.environ.get('READTHEDOCS_VERSION')

if rtd_version == 'latest':
    docs_version = 'development'
else:
    docs_version = release

rst_epilog = f'''
.. |docs_version| replace:: {docs_version}
'''

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

html_theme = 'furo'

html_theme_options = {
    'sidebar_hide_name': False,
    'navigation_with_keys': True,
}

html_sidebars = {
    '**': [
        'sidebar/brand.html',
        'sidebar/search.html',
        'sidebar/scroll-start.html',
        'sidebar/navigation.html',
        'sidebar/ethical-ads.html',
        'sidebar/scroll-end.html',
    ]
}

templates_path = ['_templates']

autosummary_generate = True
autosummary_generate_overwrite = True

autoclass_content = 'both'

autodoc_default_options = {
    'exclude-members': '__init__',
}