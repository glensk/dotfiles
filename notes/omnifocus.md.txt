⌃ + ⌘ + Arrow keys to move item up or down
⌘ + [  nest item entry 
⌘ + ]  unnest nest item entry
⌘ + }  Add Child
⌘ + {  Add Aunt
⌘ + shift + L  Add Flagg status
⌘ + '  To open add attachment of task
Show in Projects/Contexts (Option-Command-R)

################################################################
# general
################################################################
mpie p and mpie c are for planning
MPIE P and MPIE C are just the flagged and due actions


nice icons: git clone https://github.com/deaghean/omnifocus-perspective-icons   (see scripts/mac_tools/icons

my perspectives:
    select mpie and home in projects view
    Today: (Icon Ical)
        Group Actions by: Project
        Sort actons by: Due
        Fiter status by: Any status
        Filter by availibility: Available
        Fiter by duration: Any duration
        Filter Context: Remaining
        Focus: 
            mpie (folder)


how to make link to folder:
    --> drag the folder / file from finder to the omnifocus address/bar / attachment place -> click on Alias

how to make link from finder/terminal to omnifocus:
    --> open  ~/scripts/commands/omnifocus.applescript
    --> in omnifocus control click on the project and Copy as Link to get the link
    --> in the opened applescript change the link(s) (one to open and one to close)
    --> File/Export; File Format: Application --> Save  (maybe filename __link_omnifocus.app)
    --> to open omnifocus: open ~/proj/results_/__2013.10_dilute_limit_al_si/__link_omnifocus.app/


how to add a link to a perspective
    --> omnifocus:///perspective/MissingTimeEstimate
    --> Right Klick on the link and choose from the menue "Edit Link"

Do define if you see Projects (which usually dont have contexts) under contexts surf to following site
omnifocus:///change-setting?ContextModeShowsParents=false


go in browser to:
omnifocus:///change-preference?ContentLayout=compact
omnifocus:///change-preference?ContentLayout to get back to old (usual) view

-- from http://braintags.com/blog/2014/05/omnifocus-2-configuration/
Hide empty contexts in the main outline
omnifocus:///change-preference?MainOutlineIncludesEmptyContexts=false


Hide On Hold projects in forecast view
omnifocus:///change-preference?ForecastIncludesProjectsOnHold=false
omnifocus:///change-preference?ForecastIncludesProjectsOnHold

to sync every 10 minutes (instead of every hour)
omnifocus:///change-preference?MaximumTimeBetweenSync=600



##################################################################
# scripting
##################################################################

To place scripts on your OmniFocus 2 toolbar, make sure those scripts are located in Library/Application Scripts/com.omnigroup.OmniFocus2. (Sandboxing rules prevent OmniFocus 2 from running scripts which are located elsewhere.)

P.S. — In r206599, you can get to this folder by choosing Open Scripts Folder from the Help menu.
