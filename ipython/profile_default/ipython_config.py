print('~/.ipython/profile_default/ipython_config.py start ...')
c = get_config()
print('loading autoreload (reloading modules when executing script)')
c.InteractiveShellApp.extensions = ['autoreload']
# lines of code to run at IPython startup.
c.TerminalIPythonApp.exec_lines = ['%autoreload 2']

#%load_ext autoreload
#%autoreload 2
try:
    import aiida
except ImportError:
    pass
else:
    c = get_config()
    c.InteractiveShellApp.extensions = [
          'aiida.common.ipython.ipython_magics'
    ]
