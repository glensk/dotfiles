# _*_ coding: utf-8

from decimal import Decimal
import calendar
from xlsxwriter import *
import glob
import sys

zoom = 100

person       = sys.argv[1].decode('UTF-8')
workHours    = sys.argv[2]
trips        = sys.argv[3]
nrProjects   = int(sys.argv[4])
projName     = sys.argv[5]
organisation = 'Max-Planck-Institut fur Eisenforschung GmbH'

book = Workbook('EU-Timesheet.xlsx')

col1 = [ 'yellow',   'blue',    'lightBlue', 'green',   'gray'    , 'red']
col2 = [ '#ffffcc',  '#0066ff', '#8db4e2',   '#00b050', '#d9d9d9' , '#ff0000']

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

  sheet = book.add_worksheet(str(year)+" "+month)
  sheet.set_print_scale(61)  # printer scaling to 61%
  sheet.set_landscape()      # set to Landscape
  sheet.fit_to_pages(1,1)
  sheet.set_zoom(zoom)
  sheet.set_column(0  ,n+1, 4.8)
  sheet.set_column(0  ,0  ,24.3)
  sheet.set_column(n+1,n+1, 8.0)
  sheet.set_column(n+2,n+2, 8.0)

  h1 = 12; h2 = 20; h3 = 32
  sheet.set_row(0,h1)
  sheet.set_row(1,h3)
  sheet.set_row(2,h1)
  sheet.set_row(3,h2)
  sheet.set_row(4,h2)
  sheet.set_row(5,h1)
  sheet.set_row(6,h2)
  for i in range(7,70): sheet.set_row(i,h1)

  # header
  big         = book.add_format({'font_name': 'Arial', 'font_size': 26, 'bold': True})
  bold14      = book.add_format({'font_name': 'Arial', 'font_size': 14, 'bold': True})
  right14     = book.add_format({'font_name': 'Arial', 'font_size': 14, 'align': 'right'})
  boldRight14 = book.add_format({'font_name': 'Arial', 'font_size': 14, 'align': 'right', 'bold': True})

  numberOfHours = 'Number of hours envisaged i.e. according to the employment contract: '

  sheet.write(1,12, 'Timesheet'      ,big)
  sheet.write(3,0,  'Organisation: ' ,right14)
  sheet.write(3,1,  str(organisation),bold14)
  sheet.write(4,0,  'Person: '       ,right14)
  #sheet.write(4,1,  str(person)      ,bold14)
  sheet.write(4,1,  unicode(person)      ,bold14)
  sheet.write(4,25, numberOfHours    ,right14)
  sheet.write(4,26, str(workHours)   ,bold14)
  sheet.write(6,0,  year             ,boldRight14)
  sheet.write(6,2,  str(month)       ,bold14)

  # first column of table
  p1 = [projName,'yellow'];
  if nrProjects == 1:
    p2 = ['Other projects','yellow']; p3 = ['Project z','yellow'];
  else:
    p2 = ['TIME-BRIDGE','yellow']; p3 = ['Other projects','yellow'];

  ll = [ ['Date'], ['Day'], ['EU-Projects','bold','blue'], ['RTD Activities','lightBlue'], p1, p2, p3,
         ['Total RTD','right'], ['Demonstration','lightBlue'], p1, p2, p3, ['Total Demonstration','right'],
         ['Management','lightBlue'], p1, p2, p3, ['Total Management','right'], ['Other Activities','lightBlue'],
         p1, p2, p3, ['Total Other','right'], ['Internal and National Projects','bold','green'],
         ['Teaching','yellow'], ['B','yellow'], ['C','yellow'], ['Total','right'],
         ['Absences and activities not to be part of productive hours','bold','green'], ['Annual Leave'],
         ['Special Leave'], ['Illness'], ['Training / internal meetings'], ['Total','right'], [''],
         ['Total productive hours','right'], [''], ['Total hours','right'] ]

  for i in range(len(ll)):
    form = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True})
    for j in ll[i][1:]:
      if j == 'bold':  form.set_bold()
      if j == 'right': form.set_align(j)
      if j in col1: form.set_bg_color(col2[col1.index(j)])
    sheet.write(i+9,0,ll[i][0],form)

  # first two rows of table
  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0:
      form  = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'align': 'center'})
      form2 = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'align': 'center', 'bg_color': col2[col1.index('gray')]})
      if i[1] >= 5: form.set_bg_color(col2[col1.index('gray')]) # only for Sa and So
      sheet.write( 9,i[0],i[0],form)
      sheet.write(10,i[0],calendar.day_abbr[i[1]],form2)

  formRight = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'align': 'right'})
  sheet.write(9,n+1,'Total',formRight)
  sheet.write(9,n+2,'Notes',formRight)

  # last but one column containing the sum formulas
  formNum       = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'num_format': '0.00'})
  formNumBold   = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'num_format': '0.00', 'bold': True})
  formNumYellow = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'num_format': '0.00', 'bg_color': col2[col1.index('yellow')]})

  for i in range(37):
    sheet.write(10+i,n+1,None,formNum)
    sheet.write(10+i,n+2,None,formNumYellow)

  if (n==28): c1 = "AC"; c2 = "AD"
  if (n==29): c1 = "AD"; c2 = "AE"
  if (n==30): c1 = "AE"; c2 = "AF"
  if (n==31): c1 = "AF"; c2 = "AG"

  for i in range(6):
    for j in range(4):
      if j<3 or i==5:
        sheet.write(13+5*i+j,n+1,'=SUM(B%d:%s%d)'  % (13+5*i+j+1,c1,13+5*i+j+1), formNum)
      else:
        sheet.write(13+5*i+j,n+1,'=SUM(%s%d:%s%d)' % (c2,13+5*i+1,c2,13+5*i+j), formNum)

  sheet.write(42,n+1,'=SUM(%s39:%s42)'                % (c2,c2)          ,formNum)
  sheet.write(43,n+1, None                                               ,formNum)
  sheet.write(44,n+1,'=%s17+%s27+%s32+%s37'           % (c2,c2,c2,c2)    ,formNumBold)
  sheet.write(45,n+1, None                                               ,formNum)
  sheet.write(46,n+1,'=SUM(%s17+%s27+%s32+%s37+%s43)' % (c2,c2,c2,c2,c2) ,formNum)

  # interior of table
  formLightBlue = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'bg_color': col2[col1.index('lightBlue')]})
  formYellow    = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'bg_color': col2[col1.index('yellow')]})
  formYellowNum = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'bg_color': col2[col1.index('yellow')], 'num_format': '0.00'})
  formError     = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'bg_color': col2[col1.index('red')], 'bold': True})
  formGreen     = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'bg_color': col2[col1.index('green')]})
  formGray      = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True, 'bg_color': col2[col1.index('gray')]})
  formBlue      = book.add_format({'font_name': 'Arial', 'font_size': 10,                 'bg_color': col2[col1.index('blue')]})
  formBlank     = book.add_format({'font_name': 'Arial', 'font_size': 10, 'border': True})
  formBlank2    = book.add_format({'font_name': 'Arial', 'font_size': 10, 'top': True, 'bottom': True})

  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0: sheet.write(11,i[0],None,formBlue)

  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0:
      for k in range(6):
        sheet.write(12+5*k,i[0],None,formLightBlue)
        if i[1]>=5:
          for j in range(4): sheet.write(13+5*k+j,i[0],None,formGray)
        else:
          for j in range(3): sheet.write(13+5*k+j,i[0],None,formYellow)
          sheet.write(13+5*k+j+1,i[0],None,formBlank)

  # overwrite the two green rows
  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0:
      sheet.write(32,i[0],None,formGreen)
      sheet.write(37,i[0],None,formGreen)

  # last five rows of table have an unperiodic formatting
  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0:
      if i[1]>=5:
        for j in range(4): sheet.write(41+j,i[0],None,formGray)
      else:
        sheet.write(41,i[0],None,formYellow)
        for j in range(3): sheet.write(42+j,i[0],None,formBlank)
      sheet.write(45,i[0],None,formBlank2)
      sheet.write(46,i[0],None,formBlank2)

  ff = open(f,'r')
  for i in calendar.Calendar().itermonthdays2(year,m):
    if i[0] != 0:
      line = ff.readline()
      if line != '':
        day = line.split()[0]
        if day != calendar.day_abbr[i[1]]:
          raise Exception("days in nova time and python calendar do not match;\n"
                          " month: "+month+" nova day: ",day," python day: ",calendar.day_abbr[i[1]])
        if len(line.split()) == 1: continue
        last = line.split()[-1]
        if last == "F": # Feiertag = Holiday
          for k in range(6):
            for j in range(4):
              sheet.write(13+5*k+j,i[0],None,formGray)
          for j in range(3):
            sheet.write(42+j,i[0],None,formGray)
          continue
        if last == "BA": # Dienstreise = Business trip
          soll = line.split()[1]
          sheet.write(13,i[0],last,formYellowNum)
          if trips == "relevant":
            sheet.write(28,i[0],soll,formYellowNum)
          else:
            if nrProjects == 1:
              sheet.write(29,i[0],soll,formYellowNum)
            else:
              sheet.write(30,i[0],soll,formYellowNum)
          continue
        if last == "U": # Urlaub = Annual leave
          soll = line.split()[1]
          sheet.write(13,i[0],last,formYellow)
          sheet.write(38,i[0],soll,formYellowNum)
          continue
        if last == "SU": # Sonderurlaub = Special leave
          soll = line.split()[1]
          sheet.write(13,i[0],last,formYellow)
          sheet.write(39,i[0],soll,formYellowNum)
          continue
        if last == "KR": # Krankheit = Illness
          soll = line.split()[1]
          sheet.write(13,i[0],last,formYellow)
          sheet.write(39,i[0],soll,formYellowNum)
          continue
        if last == "GL": # Gleitzeit
          sheet.write(13,i[0],last,formYellow)
          continue
        if last == "!!!!!!!!": # Unentschuldigt or some other problem
          sheet.write(13,i[0],last,formError)
          continue
        ist = line.split()[1]
        if float(ist)!=0:
          if nrProjects == 1:
            sheet.write(13,i[0],ist,formYellowNum)
          else:
            sheet.write(13,i[0],float(ist)/2.,formYellowNum)
            sheet.write(14,i[0],float(ist)/2.,formYellowNum)
  ff.close()

  # footer
  form = book.add_format({'font_name': 'Arial', 'font_size': 10})
  sheet.write(49,0 ,'Signed:'                      ,form)
  sheet.write(49,10,'Approved:'                    ,form)
  sheet.write(52,0 ,'_____________________________',form)
  sheet.write(52,10,'_____________________________',form)
  sheet.write(53,0 ,'Date and Signature'           ,form)
  sheet.write(53,10,'Date and Signature'           ,form)

  formRightBottom = book.add_format({'font_name': 'Arial', 'font_size': 10, 'right': True, 'bottom': True})
  formRightTop    = book.add_format({'font_name': 'Arial', 'font_size': 10, 'right': True, 'top': True})
  formLeftBottom  = book.add_format({'font_name': 'Arial', 'font_size': 10, 'left': True, 'bottom': True})
  formLeftTop     = book.add_format({'font_name': 'Arial', 'font_size': 10, 'left': True, 'top': True})
  formBottom      = book.add_format({'font_name': 'Arial', 'font_size': 10, 'bottom': True})
  formRight       = book.add_format({'font_name': 'Arial', 'font_size': 10, 'right': True})
  formLeft        = book.add_format({'font_name': 'Arial', 'font_size': 10, 'left': True})
  formTop         = book.add_format({'font_name': 'Arial', 'font_size': 10, 'top': True})

  formBottom2 = book.add_format({'font_name': 'Arial', 'font_size': 10, 'align': 'right', 'bottom': True})
  formTop2    = book.add_format({'font_name': 'Arial', 'font_size': 10, 'align': 'right', 'top': True})
  form2       = book.add_format({'font_name': 'Arial', 'font_size': 10, 'align': 'right'})
  formSpecial = book.add_format({'font_name': 'Arial', 'font_size': 10, 'valign': 'center', 'text_wrap': True,
                                 'indent': 1, 'left': True, 'bottom': True, 'right': True})

  nn = n - 7
  sheet.write(48,nn,'Productive hours per project:',formLeftTop)
  sheet.write(49,nn,None,formLeft)
  sheet.write(50,nn,None,formLeftBottom)
  nn = n - 6
  for i in range(6):
    sheet.write(48,nn+i,None,formTop)
    sheet.write(50,nn+i,None,formBottom)
  sheet.write(48,nn+i+1,'=A14',formTop2)
  sheet.write(49,nn+i+1,'=A15',form2)
  sheet.write(50,nn+i+1,'=A16',formBottom2)
  sheet.write(48,nn+i+2,'=%s14+%s19+%s24+%s29' % (c2,c2,c2,c2),formRightTop)
  sheet.write(49,nn+i+2,'=%s15+%s20+%s25+%s30' % (c2,c2,c2,c2),formRight)
  sheet.write(50,nn+i+2,'=%s16+%s21+%s26+%s31' % (c2,c2,c2,c2),formRightBottom)
  nn = n - 8
  sheet.write(52,nn,None,formLeftTop)
  sheet.write(53,nn,None,formLeft)
  sheet.write(54,nn,None,formLeft)
  nn = n - 7
  for i in range(7): sheet.write(52,nn+i,None,formTop)
  sheet.write(52,nn+i+1,'Time of illness and training / internal meetings:',formTop2)
  sheet.write(53,nn+i+1,'Hours:'                                           ,form2)
  sheet.write(54,nn+i+1,'Days (here 8 hours are one day):'                 ,form2)
  sheet.write(52,nn+i+2, None                                              ,formRightTop)
  sheet.write(53,nn+i+2,'=AG41+AG42'                                       ,formRight)
  sheet.write(54,nn+i+2,'=AG54/8'                                          ,formRight)
  sheet.merge_range(55,n-8,57,n+1,'If the days are higher than 15 for a '
     'duration of one year the time should be duly justified (Financial '
     'Guide, Version 02/04/2009, page 45).'                                ,formSpecial)

  # finally images
  sheet.insert_image(0, 0,'sev.bmp', {'x_scale': 1.0, 'y_scale': 0.99})
  sheet.insert_image(0,21,'EUB.bmp', {'x_scale': 1.1, 'y_scale': 1.2})
  sheet.insert_image(0,27,'NKS.bmp', {'x_scale': 1.2, 'y_scale': 1.2, 'y_offset': 6})

book.close()

