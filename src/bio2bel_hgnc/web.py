# -*- coding: utf-8 -*-

"""This module builds a :mod:`Flask` application for interacting with the underlying database. When installing,
use the web extra like:

.. source-code:: sh

    pip install bio2bel_hgnc[web]
"""

import logging

from bio2bel_hgnc.manager import Manager

log = logging.getLogger(__name__)

if __name__ == '__main__':
    manager = Manager()
    app_ = manager.get_flask_admin_app()
    app_.run(debug=True, host='0.0.0.0', port=5000)
