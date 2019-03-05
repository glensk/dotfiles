yum install inotify-tools

python pyinotify.py -a -e IN_ATTRIB,IN_CLOSE_WRITE,IN_MODIFY,IN_CREATE,IN_DELETE,IN_DELETE_SELF,IN_MOVE_SELF,IN_MOVED_FROM,IN_MOVED_TO -r /home/glensk/talks

~/scripts/pyinotify/my_new_env/lib/python2.7/site-packages/pyinotify.py
