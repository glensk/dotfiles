
==> Pouring octave-4.0.3_4.el_capitan.bottle.tar.gz
Warning: homebrew/science/octave dependency gcc was built with a different C++ standard
library (libstdc++ from clang). This may cause problems at runtime.
==> Caveats

Octave's graphical user interface is disabled; compile Octave with the option
--with-gui to enable it.


Several graphics toolkit are available. You can select them by using the command
'graphics_toolkit' in Octave.  Individual Gnuplot terminals can be chosen by setting
the environment variable GNUTERM and building gnuplot with the following options:

  setenv('GNUTERM','qt')    # Requires QT; install gnuplot --with-qt5
  setenv('GNUTERM','x11')   # Requires XQuartz; install gnuplot --with-x11
  setenv('GNUTERM','wxt')   # Requires wxmac; install gnuplot --with-wxmac
  setenv('GNUTERM','aqua')  # Requires AquaTerm; install gnuplot --with-aquaterm

You may also set this variable from within Octave. For printing the cairo backend
is recommended, i.e., install gnuplot with --with-cairo, and use

  print -dpdfcairo figure.pdf


When using the native Qt or fltk toolkits then invisible figures do not work because
osmesa is incompatible with Mac's OpenGL. The usage of gnuplot is recommended.


==> Summary
🍺  /usr/local/Cellar/octave/4.0.3_4: 2,184 files, 47.3M
==> Upgrading open-mpi
