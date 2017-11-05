# -*- coding: utf-8 -*-

from six import string_types

from pyhgnc import QueryManager


def ensure_pyhgnc_manager(connection=None):
    """Ensures there's a query manager

    :param connection: Either a :class:`QueryManager` or arguments to build one.
    :type connection: None or str or QueryManager
    :rtype: QueryManager
    """
    if isinstance(connection, QueryManager):
        return connection

    if connection is None or isinstance(connection, string_types):
        return QueryManager(connection=connection)

    raise TypeError
