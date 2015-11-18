
" ================ Keystrokes to remember / great shortcuts/ muscle memory{{{
"  BASICS:
"       i 	insert before character under cursor
"       a 	insert after cursor
"       I 	insert at beginning of current line
"       A 	insert at end of the line
"       o 	starts insert mode in a new line below current one
"       O 	insert in a new line above current one 
"  MOVING IN VIM BETWEEN SPLITS:
"           capslock + w
"           splits, move between splits (control||capslock + hjkl) ###############
"  MOVING IN VIM:
"       → 	← 	description
"       / 	? 	search for a pattern of text, jump to it by hitting Enter (<CR>)
"       * 	# 	search for the word under cursor
"       n 	N 	jump to the next match for the previous search
"       $ 	^ 	position cursor at end of current line
"       f 	F 	position cursor on the character in the same line that matches the next keystroke
"       t 	T 	position cursor before the next character that matches the keystroke
"       ; 	, 	repeat the last f, F, t, or T
"       w 	b 	move to start of next word
"       W 	B 	move to start of next "WORD" (sequence of non-blank characters)
"       } 	{ 	move down one paragraph (block of text separated by blank lines)
"       gg 	jump to first line of document
"       G 	jump to end of document 
"       ''  jump to last plase before jumping
"       '.  jump to last plase before editing
"  NAVIGATE:
"       ll   go to next l character;  fk   go to next k character
"       Fl   go to last l character;  Fk   go to last k character
"       )   Jump forward one sentence. 
"       (   Jump backward one sentence. 
"       }   Jump forward one paragraph. 
"       {  Jump backward one paragraph. 
"       cw   to delete next word and go into insert mode
"       <ctrl + w> deletes last word also in insert mode
"       <ctrl + u> deltes to the beginning of line deletes line in insert hode
"       <ctrl + h> backspace: deletes last char also in insert mode
"       <ctr. + n> or <ctrl +p> go up/down tab vorschlag
"       command + hjkl to move left down up right also in insert mode ( needs to have have Report Terminal Type: linux under iterm2 to be set
"   NAVIGATE IN COMMADN MODE:
"       0   beginning of line
"       $   end of line
"       shift i     first non whitespace char
"       shift a     append to end of line
"       o           create a new line below
"       shift o     create a new line above
"       C (shift c)     change entire line
"       x           deletes
"       shift x     deletes back
"       .           repeates the last command
"       shift h     high
"       shift m     middle
"       shift l     low
"       shift j     fuege nexte zeile an die aktuelle
"  IMPORTANT MOVEMENTSS:
"       control-o   jump to the last visited places
"       TAB         brings you back to rute with control-o
"  YANKING:
"       Y 	yank current line; prepend with number to yank that many lines
"       y} 	yank until end of paragraph
"       dd 	delete current line and yank it too (think "cut")
"       d3d 	delete 3 lines starting from current one
"       p 	paste yanked text at cursor; prepend number to paste that many times
"       P 	paste before cursor 
"  UPPERCASE LOWERCASE: toggle:
"  ( http://vim.wikia.com/wiki/Switching_case_of_characters ) 
"  mark correspondign word:
"    Toggle case "HellO" to "hELLo" with g~ then a movement. 
"    Uppercase "HellO" to "HELLO" with gU then a movement. 
"    Lowercase "HellO" to "hello" with gu then a movement. 
"  HIGHLIGHT INNER BRACKETS OR QUOTES:
"       vi( || vi) || vi" || vi' || vi{ all select inner quotes (vib) to mark whole
"       parentheses block (hallo wie )
"       yi( to to to beginning of bracket (yi" yi' yi[)
"  CHANGE:
"       ciw 	("change inner word") change word under cursor
"       ci" 	change double-quoted string (but keep the quotes)
"       ci( 	change text between matching parentheses, also works with brackets
"       cc 	    change whole line 
"  NEARDCOMMENTER:
"  <leader>ci to toggle comment of line
"  <leader>c<space> to toggle comment of line
"
"  NEARDTREE:
"  shift + T : opens file in new tab
"
"  SEARCH:
"  ggn : jump to fist search result
"  GN  : jump to last search result
"  <space> : unhighlight search results

"  VIMLATEX:
"  <leader>ll    -> comile source
"  <leader>ls    -> go to skim an show line
"  :retab     -> Change all existing tab characters to match current tab settings

"  REPLACE: ( http://www.guckes.net/vi/substitute.html )
"  :s/pattern/replacement/g             to replace stuff in current line
"  :%s/pattern/replacement/g            to replace in whole document
"  :%s/pattern/replacement/gc           to be asked at every instance
"  :%s/foo/bar/g    Find each occurrence of 'foo' (in all lines), and replace it with 'bar'. 
"  :s/foo/bar/g     Find each occurrence of 'foo' (in the current line only), and replace it with 'bar'. 
"  :%s/foo/bar/gc   Change each 'foo' to 'bar', but ask for confirmation first. 
"  :%s/\<foo\>/bar/gc   Change only whole words exactly matching 'foo' to 'bar'; ask for confirmation. 
"  :%s/foo/bar/gci      Change each 'foo' (case insensitive) to 'bar'; ask for confirmation. 
"                       This may be wanted after using :set noignorecase to make searches case sensitive (the default). 
"  :%s/foo/bar/gcI      Change each 'foo' (case sensitive) to 'bar'; ask for confirmation. 
"    This may be wanted after using :set ignorecase to make searches case insensitive. 
"
"  ,ev        -> open ~/.vimrc
"
"  REPLACE VISUAL SELECTION:
"  Select the first block: ctrl-v (move and select your block) "ay
"  Select the second block: ctrl-v (move and select block to change) c ctrl-o "aP <Esc>
"
"   MISCELLANIOUS:
"       :echo has('clipboard') if 1 we will directly paste to system clipboard
"       otherwise not  (run vim -v to get access to system clipboard)
"       highlight text :!fmt    -> formats text
" }}}
