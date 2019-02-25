brew link osxfuse
brew link --overwrite --dry-run osxfuse

%brew link --overwrite --dry-run osxfuse
Would remove:
/usr/local/lib/pkgconfig/fuse.pc -> /usr/local/lib/pkgconfig/osxfuse.pc
/usr/local/lib/libosxfuse_i64.dylib -> /usr/local/lib/libosxfuse_i64.2.dylib
/usr/local/lib/libosxfuse_i64.2.dylib
/usr/local/lib/libosxfuse_i32.dylib -> /usr/local/lib/libosxfuse_i32.2.dylib
/usr/local/lib/libosxfuse_i32.2.dylib
/usr/local/lib/libosxfuse.dylib -> /usr/local/lib/libosxfuse_i64.dylib
/usr/local/lib/libosxfuse.2.dylib -> /usr/local/lib/libosxfuse_i64.2.dylib

