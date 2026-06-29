project = 'Codeine'
author = 'Michael Dunne'
html_title = 'Codeine documentation'

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