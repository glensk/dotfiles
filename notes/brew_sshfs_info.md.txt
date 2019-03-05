run in shell:
$brew info fuse4x-kext


fuse4x-kext: stable 0.9.2
http://fuse4x.github.com
/usr/local/Cellar/fuse4x-kext/0.9.2 (6 files, 276K) *
  Built from source
From: https://github.com/mxcl/homebrew/commits/master/Library/Formula/fuse4x-kext.rb
==> Caveats
In order for FUSE-based filesystems to work, the fuse4x kernel extension
must be installed by the root user:

  sudo /bin/cp -rfX /usr/local/Cellar/fuse4x-kext/0.9.2/Library/Extensions/fuse4x.kext /Library/Extensions
  sudo chmod +s /Library/Extensions/fuse4x.kext/Support/load_fuse4x

If upgrading from a previous version of Fuse4x, the old kernel extension
will need to be unloaded before performing the steps listed above. First,
check that no FUSE-based filesystems are running:

  mount -t fuse4x

Unmount all FUSE filesystems and then unload the kernel extension:

  sudo kextunload -b org.fuse4x.kext.fuse4x
