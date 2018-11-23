try:
    import aiida
except ImportError:
    pass
else:
    c = get_config()
    c.InteractiveShellApp.extensions = [
          'aiida.common.ipython.ipython_magics'
    ]
