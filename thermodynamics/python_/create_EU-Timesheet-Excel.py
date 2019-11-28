from decimal import Decimal
import calendar
from xlwt import *
import glob
import sys

zoom = 75

person       = sys.argv[1]
workHours    = sys.argv[2]
organisation = 'Max-Planck-Institut fur Eisenforschung GmbH'

book = Workbook()
book.set_colour_RGB(0x36,255,255,204) # replaces blue_gray     -> yellow
book.set_colour_RGB(0x0B,0  ,102,255) # replaces bright_green  -> blue
book.set_colour_RGB(0x1D,141,180,226) # replaces coral         -> lightBlue
book.set_colour_RGB(0x0F,0  ,176,80 ) # replaces cyan_ega      -> green
book.set_colour_RGB(0x12,217,217,217) # replaces dark_blue_ega -> gray

col1 = [ 'yellow', 'blue', 'lightBlue', 'green', 'gray' ]
col2 = [  0x36,     0x0B,   0x1D,        0x0F,    0x12  ]

files = glob.glob('_tmp3_[0-9]*_[0-9]*_[a-zA-Z]*')
files.sort()

for f in files:
  month        = f.split('_')[4]
  year         = int(f.split('_')[2])

  # we need to know the days for the given year and month
  m = 0
  for i in range(1,13):
    if calendar.month_name[i] == month:
      m = i
      break
  if m == 0:
    raise Exception("given month: "+month+" invalid")

  # get number of days for month
  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0: n = i[0]

  sheet = book.add_sheet(str(year)+" "+month,cell_overwrite_ok=True)
  sheet.portrait = 0       # set to Landscape
  sheet.print_scaling = 61 # printer scaling to 61%
  sheet.normal_magn = zoom
  for c in range(n+1):
    sheet.col(c).width = 1460
  sheet.col(0).width  = 6408
  sheet.col(n+1).width = 2218
  sheet.col(n+2).width = 2218

  # header
  sheet.write(1,12, 'Timesheet'       ,easyxf('font: height 520, bold True;'))
  sheet.write(3,0,  'Organisation: '  ,easyxf('font: height 280;' 'alignment: horiz right;'))
  sheet.write(3,1,  str(organisation) ,easyxf('font: height 280, bold True;'))
  sheet.write(4,0,  'Person: '        ,easyxf('font: height 280;' 'alignment: horiz right;'))
  sheet.write(4,1,  str(person)       ,easyxf('font: height 280, bold True;'))
  sheet.write(4,25, 'Number of hours envisaged i.e. according to the employment contract: '
                                      ,easyxf('font: height 280;' 'alignment: horiz right;'))
  sheet.write(4,26, str(workHours)    ,easyxf('font: height 280, bold True;'))
  sheet.write(6,0,  year              ,easyxf('font: height 280, bold True;' 'alignment: horiz right;'))
  sheet.write(6,2,  str(month)        ,easyxf('font: height 280, bold True;'))

  # first column of table
  p1 = ['SMARTMET','yellow']; p2 = ['Project y','yellow']; p3 = ['Project z','yellow'];

  ll = [ ['Date'], ['Day'], ['EU-Projects','bold','blue'], ['RTD Activities','lightBlue'], p1, p2, p3,
         ['Total RTD','right'], ['Demonstration','lightBlue'], p1, p2, p3, ['Total Demonstration','right'],
         ['Management','lightBlue'], p1, p2, p3, ['Total Management','right'], ['Other Activities','lightBlue'],
         p1, p2, p3, ['Total Other','right'], ['Internal and National Projects','bold','green'],
         ['Teaching','yellow'], ['B','yellow'], ['C','yellow'], ['Total','right'],
         ['Absences and activities not to be part of productive hours','bold','green'], ['Annual Leave'],
         ['Special Leave'], ['Illness'], ['Training / internal meetings'], ['Total','right'], [''],
         ['Total productive hours','right'], [''], ['Total hours','right'] ]

  bor = Borders()
  bor.left   = Borders.THIN
  bor.right  = Borders.THIN
  bor.top    = Borders.THIN
  bor.bottom = Borders.THIN

  aliCenter = Alignment()
  aliRight  = Alignment()
  aliCenter.horz = 2
  aliRight.horz  = 3

  fntBold = Font()
  fntBold.bold = True

  for i in range(len(ll)):
    style = XFStyle()
    style.borders = bor
    for j in ll[i][1:]:
      if j == 'bold':  style.font = fntBold
      if j == 'right': style.alignment = aliRight
      if j in col1:
        pat = Pattern()
        pat.pattern = Pattern.SOLID_PATTERN
        pat.pattern_fore_colour = col2[col1.index(j)]
        style.pattern = pat
    sheet.write(i+9,0,ll[i][0],style)

  # first two rows of table
  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0:
      style = XFStyle()
      style.borders = bor
      style.alignment = aliCenter
      pat = Pattern()
      pat.pattern = Pattern.SOLID_PATTERN
      if i[1] >= 5: # only for Sa and So here
        pat.pattern_fore_colour = col2[col1.index('gray')]
      else:
        pat.pattern_fore_colour = 0x09
      style.pattern = pat
      sheet.write( 9,i[0],i[0],style)
      style2 = XFStyle()
      pat2 = Pattern()
      pat2.pattern = Pattern.SOLID_PATTERN
      pat2.pattern_fore_colour = col2[col1.index('gray')]
      style2.alignment = aliCenter
      style2.borders = bor
      style2.pattern = pat2
      sheet.write(10,i[0],calendar.day_abbr[i[1]],style2)

  style = XFStyle()
  style.borders = bor
  style.alignment = aliRight
  sheet.write(9,n+1,'Total',style)
  sheet.write(9,n+2,'Notes',style)

  # last but one column containing the sum formulas
  style = XFStyle()
  styleYellow = XFStyle()
  style.num_format_str = '0.00'
  styleYellow.num_format_str = '0.00'
  style.borders = bor
  styleYellow.borders = bor
  pat = Pattern()
  pat.pattern = Pattern.SOLID_PATTERN
  pat.pattern_fore_colour = 0x36
  styleYellow.pattern = pat
  for i in range(37):
    sheet.write(10+i,n+1,None,style)
    sheet.write(10+i,n+2,None,styleYellow)

  if (n==28): column = "AC"
  if (n==29): column = "AD"
  if (n==30): column = "AE"
  if (n==31): column = "AF"

  for i in range(6):
    for j in range(4):
      sheet.write(13+5*i+j,n+1,Formula('SUM(B%d:%s%d)' % (13+5*i+j+1,column,13+5*i+j+1)), style)

  styleBold = XFStyle()
  styleBold.borders = bor
  styleBold.font = fntBold
  styleBold.num_format_str = '0.00'

  sheet.write(42,n+1,Formula('SUM(AG39:AG42)'),style)
  sheet.write(43,n+1,None,style)
  sheet.write(44,n+1,Formula('AG17+AG27+AG32+AG37'),styleBold)
  sheet.write(45,n+1,None,style)
  sheet.write(46,n+1,Formula('SUM(AG17+AG27+AG32+AG37+AG43)'),style)

  # interior of table
  styleBlue = XFStyle()
  styleBlue.pattern.pattern = 1
  styleBlue.pattern.pattern_fore_colour = col2[col1.index('blue')]

  styleLightBlue = XFStyle()
  styleLightBlue.pattern.pattern = 1
  styleLightBlue.pattern.pattern_fore_colour = col2[col1.index('lightBlue')]
  styleLightBlue.borders = bor

  styleYellow = XFStyle()
  styleYellow.pattern.pattern = 1
  styleYellow.pattern.pattern_fore_colour = col2[col1.index('yellow')]
  styleYellow.borders = bor

  styleGray = XFStyle()
  styleGray.pattern.pattern = 1
  styleGray.pattern.pattern_fore_colour = col2[col1.index('gray')]
  styleGray.borders = bor

  styleGreen = XFStyle()
  styleGreen.pattern.pattern = 1
  styleGreen.pattern.pattern_fore_colour = col2[col1.index('green')]

  styleBlank = XFStyle()
  styleBlank.borders = bor

  styleBlank2 = XFStyle()
  styleBlank2.borders.top = 1
  styleBlank2.borders.bottom = 1

  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0:
      sheet.write(11,i[0],None,styleBlue)

  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0:
      for k in range(6):
        sheet.write(12+5*k,i[0],None,styleLightBlue)
        if i[1]>=5:
          for j in range(4):
            sheet.write(13+5*k+j,i[0],None,styleGray)
        else:
          for j in range(3):
            sheet.write(13+5*k+j,i[0],None,styleYellow)
          sheet.write(13+5*k+j+1,i[0],None,styleBlank)

  # overwrite the two green rows
  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0:
      sheet.write(32,i[0],None,styleGreen)
      sheet.write(37,i[0],None,styleGreen)

  # last five rows of table have an unperiodic formatting
  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0:
      if i[1]>=5:
        for j in range(4):
          sheet.write(41+j,i[0],None,styleGray)
      else:
        sheet.write(41,i[0],None,styleYellow)
        for j in range(3):
          sheet.write(42+j,i[0],None,styleBlank)
      sheet.write(45,i[0],None,styleBlank2)
      sheet.write(46,i[0],None,styleBlank2)

  # footer
  sheet.write(49,0 ,'Signed:')
  sheet.write(49,10,'Approved:')
  sheet.write(52,0 ,'_____________________________')
  sheet.write(52,10,'_____________________________')
  sheet.write(53,0 ,'Date and Signature')
  sheet.write(53,10,'Date and Signature')

  nn = n - 7
  sheet.write(48,nn,'Productive hours per project:',easyxf('borders: left thin, top thin;'))
  sheet.write(49,nn,None,easyxf('borders: left thin;'))
  sheet.write(50,nn,None,easyxf('borders: left thin, bottom thin;'))
  nn = n - 6
  for i in range(6):
    sheet.write(48,nn+i,None,easyxf('borders: top thin;'))
    sheet.write(50,nn+i,None,easyxf('borders: bottom thin;'))
  sheet.write(48,nn+i+1,Formula('A14'),easyxf('alignment: horz right;' 'borders: top thin;'))
  sheet.write(49,nn+i+1,Formula('A15'),easyxf('alignment: horz right;'))
  sheet.write(50,nn+i+1,Formula('A16'),easyxf('alignment: horz right;' 'borders: bottom thin;'))
  sheet.write(48,nn+i+2,Formula('AG14+AG19+AG24+AG29'),easyxf('borders: right thin, top thin;'))
  sheet.write(49,nn+i+2,Formula('AG15+AG20+AG25+AG30'),easyxf('borders: right thin;'))
  sheet.write(50,nn+i+2,Formula('AG16+AG21+AG26+AG31'),easyxf('borders: right thin, bottom thin;'))
  nn = n - 8
  sheet.write(52,nn,None,easyxf('borders: left thin, top thin;'))
  sheet.write(53,nn,None,easyxf('borders: left thin;'))
  sheet.write(54,nn,None,easyxf('borders: left thin;'))
  nn = n - 7
  for i in range(7):
    sheet.write(52,nn+i,None,easyxf('borders: top thin;'))
  sheet.write(52,nn+i+1,'Time of illness and training / internal meetings:',easyxf('alignment: horz right;' 'borders: top thin;'))
  sheet.write(53,nn+i+1,'Hours:',easyxf('alignment: horz right;'))
  sheet.write(54,nn+i+1,'Days (here 8 hours are one day):',easyxf('alignment: horz right;'))
  sheet.write(52,nn+i+2,None,easyxf('borders: right thin, top thin;'))
  sheet.write(53,nn+i+2,Formula('AG41+AG42'),easyxf('borders: right thin;'))
  sheet.write(54,nn+i+2,Formula('AG54/8'),easyxf('borders: right thin;'))
  sheet.write_merge(55,57,n-8,n+1,'If the days are higher than 15 for a duration of one year the time '
                           'should be duly justified (Financial Guide, Version 02/04/2009, page 45).'
                           ,easyxf('borders: left thin, bottom thin, right thin;'
                                   'alignment: wrap True, indent 1, vertical center'))

  # finally images
  sheet.insert_bitmap('seventh-framework.bmp', 0, 0, scale_x = 1.18, scale_y = 1.44)
  sheet.insert_bitmap('EUB.bmp', 0, 21, scale_x = 1.18, scale_y = 1.44)
  sheet.insert_bitmap('NKS.bmp', 1, 27, scale_x = 1.18, scale_y = 1.44)

book.save('EU-Timesheet.xls')

