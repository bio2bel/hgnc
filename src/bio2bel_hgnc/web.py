# -*- coding: utf-8 -*-

""" This module contains the flask application to visualize the db

when pip installing

.. source-code:: sh

    pip install bio2bel_expasy[web]

"""

import flask_admin
from flask import Flask
from flask_admin.contrib.sqla import ModelView

from .manager import Manager
from .models import GeneFamily, HGNC


def add_admin(app, session, **kwargs):
    admin = flask_admin.Admin(app, **kwargs)
    admin.add_view(ModelView(HGNC, session))
    admin.add_view(ModelView(GeneFamily, session))
    return admin


def create_app(connection=None, url=None):
    """Creates a Flask application

    :type connection: Optional[str]
    :type url: Optional[str]
    :rtype: flask.Flask
    """
    app = Flask(__name__)
    manager = Manager(connection=connection)
    add_admin(app, manager.session, url=url)
    return app


if __name__ == '__main__':
    app = create_app()
    app.run(debug=True, host='0.0.0.0', port=5000)
