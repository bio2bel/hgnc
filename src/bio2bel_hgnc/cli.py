# -*- coding: utf-8 -*-

"""Command line interface for Bio2BEL HGNC."""

from .manager import Manager

main = Manager.get_cli()

if __name__ == '__main__':
    main()
