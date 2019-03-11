tell application "Firefox"
	open location "omnifocus:///task/jNBEsQBVHjp"
end tell
(*
tell application "System Events" to tell process "firefox"
	set frontmost to true
	repeat with w from 1 to count of windows
		perform action "AXRaise" of window w
		set startTab to name of window 1
		repeat
			if name of window 1 contains "omnifocus:///task/jNBEsQBVHjp" then
				keystroke "w" using command down
			else
				keystroke "}" using command down
			end if
			delay 0.2
			if name of window 1 is startTab then exit repeat
		end repeat
	end repeat
end tell
*)