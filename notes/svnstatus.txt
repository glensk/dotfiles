' '
No modifications.

'A'
Item is scheduled for Addition.

'D'
Item is scheduled for Deletion.

'M'
Item has been modified.

'C'
Item is in conflict with updates received from the repository.

'X'
Item is related to an externals definition.

'I'
Item is being ignored (e.g. with the svn:ignore property).

'?'
Item is not under version control.

'!'
Item is missing (e.g. you moved or deleted it without using svn). This also indicates that a directory is incomplete (a checkout or update was interrupted).

'~'
Item is versioned as a directory, but has been replaced by a file, or vice versa.

The second column tells the status of a file's or directory's properties.

' '
No modifications.

'M'
Properties for this item have been modified.

'C'
Properties for this item are in conflict with property updates received from the repository.

The third column is populated only if the working copy directory is locked.

' '
Item is not locked.

'L'
Item is locked.

The fourth column is populated only if the item is scheduled for addition-with-history.

' '
No history scheduled with commit.

'+'
History scheduled with commit.

The fifth column is populated only if the item is switched relative to its parent (see the section called “Switching a Working Copy”).

' '
Item is a child of its parent directory.

'S'
Item is switched.

The out-of-date information appears in the eighth column (only if you pass the --show-updates switch).

' '
The item in your working copy is up-to-date.

'*'
A newer revision of the item exists on the server.

The remaining fields are variable width and delimited by spaces. The working revision is the next field if the --show-updates or --verbose switches are passed.

If the --verbose switch is passed, the last committed revision and last committed author are displayed next.

The working copy path is always the final field, so it can include spaces.


