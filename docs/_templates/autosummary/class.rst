{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :exclude-members: __init__
   :show-inheritance:

{% if attributes %}


Properties
~~~~~~~~~~

.. autosummary::

{% for item in attributes %}
   ~{{ name }}.{{ item }}
{% endfor %}
{% endif %}

{% if methods %}


Methods
~~~~~~~~~~

.. autosummary::

{% for item in methods %}
{% if item != "__init__" %}
   ~{{ name }}.{{ item }}
{% endif %}
{% endfor %}


{% for item in methods %}
{% if item != "__init__" %}
.. automethod:: {{ name }}.{{ item }}

{% endif %}
{% endfor %}
{% endif %}