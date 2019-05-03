#------------------------------------------------------------------------------
# PromptManager configuration
#------------------------------------------------------------------------------
	# attributes of IPython.utils.coloransi.InputTermColors
	#{color.Black}           {color.Green}
	#{color.BlinkBlack}      {color.LightBlue}
	#{color.BlinkBlue}       {color.LightCyan}
	#{color.BlinkCyan}       {color.LightGray}
	#{color.BlinkGreen}      {color.LightGreen}
	#{color.BlinkLightGray}  {color.LightPurple}
	#{color.BlinkPurple}     {color.LightRed}
	#{color.BlinkRed}        {color.Purple}
	#{color.BlinkYellow}     {color.Red}
	#{color.Blue}            {color.White}
	#{color.Brown}           {color.Yellow}
	#{color.Cyan}            {color.Normal} (resets to terminal default)
	#{color.DarkGray}                       (try {color.normal} if above doesn't work)

	# \#
	#     the prompt/history count number. This escape is automatically
	#     wrapped in the coloring codes for the currently active color scheme.
	# \N  the 'naked' prompt/history count number: this is just the number
	#     itself, without any coloring applied to it. This lets you produce
	#     numbered prompts with your own colors.
	# \D  the prompt/history count, with the actual digits replaced by dots.
	#     Used mainly in continuation prompts (prompt_in2)
	# \w  the current working directory
	# \T  the current time
	# \W  the basename of current working directory
	# \Xn
	#     where $n=0\ldots5.$ The current working directory, with $HOME
	#     replaced by ~, and filtered out to contain only $n$ path elements
	# \Yn
	#     Similar to \Xn, but with the $n+1$ element included if it is ~ (this
	#     is similar to the behavior of the %cn escapes in tcsh)
	# \u
	#     the username of the current user
	# \$
	#     if the effective UID is 0, a #, otherwise a $
	# \h
	#     the hostname up to the first '.'
	# \H
	#     the hostname
	# \n
	#     a newline
	# \r
	#     a carriage return
	# \v
	#     IPython version string
	
	# Output prompt. '\#' will be transformed to the prompt number
	# c.PromptManager.out_template = 'Out[\\#]: '
	
	# Continuation prompt.
	# c.PromptManager.in2_template = '   .\\D.: '
	
	# If True (default), each prompt will be right-aligned with the preceding one.
	# c.PromptManager.justify = True
	
	# Input prompt.  '\#' will be transformed to the prompt number
#c.PromptManager.in_template = r'\T {color.LightGreen}\u@\h{color.LightBlue} {color.LightCyan}\w{color.LightBlue}{color.Green}\n[\#]'
c.PromptManager.in_template = r'\T {color.LightGreen}\u@\h \w\n{color.LightPurple}[\#{color.LightPurple}]'
#c.PromptManager.in_template = "[{_foobar}]\n[\#]>>> "
	
	# 
	# c.PromptManager.color_scheme = 'Linux'
