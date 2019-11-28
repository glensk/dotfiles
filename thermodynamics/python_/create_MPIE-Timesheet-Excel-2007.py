# _*_ coding: utf-8

from xlsxwriter import *
import glob
import sys
import datetime

name  = sys.argv[1].decode('UTF-8')
ID    = sys.argv[2]
proj  = sys.argv[3]
nrProjects = int(sys.argv[4])

book = Workbook('MPIE-Timesheet.xlsx')

f = open('_tmp_soll_und_fehlzeit_stunden','r')

for line in f:
  year     = int(line.split()[0].split('_')[0])-2000
  month    = int(line.split()[0].split('_')[1])
  monthStr = line.split()[0].split('_')[2]
  soll     = float(line.split()[1])
  fehl     = float(line.split()[2])
  arbeit   = soll - fehl

  sheet = book.add_worksheet(str(2000+year)+" "+monthStr)
  sheet.set_margins(left=1,right=.3,top=0.7)
  sheet.set_print_scale(86)

  nrow = 33
  ncol = 28
  sheet.set_column(     0,ncol-2,1.5)
  sheet.set_column(ncol-1,ncol-1,12)
  sheet.set_column(  ncol,  ncol,13)

  r1 = 20
  r2 = 34
  r3 = 23
  sheet.set_row(0,r1)
  sheet.set_row(1,r2)
  sheet.set_row(2,r1)
  sheet.set_row(3,r1)
  sheet.set_row(4,r2)
  for i in range(5,nrow):
    sheet.set_row(i,r3)

  formTitle   = book.add_format({'border': 2, 'font_name': 'Arial',        'font_size': 15, 'align': 'center', 'valign': 'vcenter'})
  formSpecial = book.add_format({'border': 2, 'font_name': 'Arial',        'font_size': 26, 'align': 'center', 'valign': 'vcenter', 'bold': True})
  formCenter  = book.add_format({'border': 2, 'font_name': 'Arial',        'font_size': 10, 'align': 'center', 'valign': 'vcenter', 'text_wrap': True})
  formCenter2 = book.add_format({'border': 2, 'font_name': 'Arial',        'font_size': 12, 'align': 'center', 'valign': 'vcenter', 'text_wrap': True})
  formCenterN = book.add_format({'border': 2, 'font_name': 'Arial Narrow', 'font_size': 10, 'align': 'center', 'valign': 'vcenter', 'text_wrap': True})

  sheet.merge_range( 0,   0, 0,ncol-1,'STUNDENZETTEL'   ,formTitle)
  sheet.merge_range( 1,   0, 1,     1,'KA'              ,formCenter)
  sheet.merge_range( 2,   0, 2,     1,'65'              ,formCenter2)
  sheet.merge_range( 1,   2, 1,     3,'von Firma'       ,formCenterN)
  sheet.merge_range( 2,   2, 2,     3,'90'              ,formCenter2)
  sheet.merge_range( 1,   4, 1,     8,'Personal- Nummer',formCenter)
  sheet.merge_range( 2,   4, 2,     8, ID               ,formCenter2)
  sheet.merge_range( 1,   9, 1,    10,'Mon.'            ,formCenter)
  sheet.merge_range( 2,   9, 2,    10, month            ,formCenter2)
  sheet.merge_range( 1,  11, 1,    12,'Jahr'            ,formCenter)
  sheet.merge_range( 2,  11, 2,    12, year             ,formCenter2)
  sheet.merge_range( 1,  13, 1,    15,'von KST'         ,formCenter)
  sheet.merge_range( 2,  13, 2,    15,''                ,formCenter2)
  sheet.merge_range( 1,  16, 1,    19,'Soll Stunden'    ,formCenter)
  sheet.merge_range( 2,  16, 2,    19, soll             ,formCenter2)
  sheet.merge_range( 1,  20, 1,ncol-1,'Name'            ,formCenter)
  sheet.merge_range( 2,  20, 2,ncol-1, name             ,formCenter2)
  sheet.merge_range( 0,ncol, 2,  ncol,'MPI'             ,formSpecial)
                                                        
  formLeft   = book.add_format({'left': 2, 'right': 1, 'top': 1, 'bottom': 1, 'font_name': 'Arial', 'font_size': 12, 'valign': 'vcenter', 'align': 'center'})
  formRight  = book.add_format({'left': 1, 'right': 2, 'top': 1, 'bottom': 1, 'font_name': 'Arial', 'font_size': 12, 'valign': 'vcenter', 'align': 'center'})
  formDash1  = book.add_format({'left': 1, 'right': 3, 'top': 1, 'bottom': 1, 'font_name': 'Arial', 'font_size': 12, 'valign': 'vcenter', 'align': 'center'})
  formDash2  = book.add_format({'left': 3, 'right': 2, 'top': 1, 'bottom': 1, 'font_name': 'Arial', 'font_size': 12, 'valign': 'vcenter', 'align': 'center'})
  formTop2   = book.add_format({'top': 2,  'font_name': 'Arial', 'font_size': 12, 'valign': 'top'})
  formLeft2  = book.add_format({'left': 2, 'font_name': 'Arial', 'font_size': 12, 'valign': 'top'})
  formRight2 = book.add_format({'right': 2})
  form       = book.add_format({'border': 1, 'font_name': 'Arial', 'font_size': 12, 'valign': 'vcenter', 'align': 'center'})
  form2      = book.add_format({'border': 2, 'font_name': 'Arial', 'font_size': 12, 'valign': 'vcenter'})

  off = ncol - 12
  for i in range(2):
    sheet.merge_range( 4, 0+i*off, 4, 1+i*off,'f'+unichr(0x00FC)+'r Firma'             ,formCenterN)
    sheet.merge_range( 4, 2+i*off, 4, 6+i*off,'f'+unichr(0x00FC)+'r KST    oder KSTTR' ,formCenter)
    sheet.merge_range( 4, 7+i*off, 4,10+i*off,'Std.'                  ,formCenter)
    for j in range(5,nrow):
      sheet.write(j,0+i*off,None,formLeft)
      sheet.write(j,1+i*off,None,formRight)
      sheet.write(j,2+i*off,None,formLeft)
      for k in range(3,6): sheet.write(j,k+i*off,None,form)
      sheet.write(j,6+i*off,None,formRight)
      sheet.write(j,7+i*off,None,formLeft)
      for k in range(8,9): sheet.write(j,k+i*off,None,form)
      sheet.write(j, 9+i*off,None,formDash1)
      sheet.write(j,10+i*off,None,formDash2)

  sheet.merge_range(nrow-1,  0,nrow-1,    6,unichr(0x00DC)+'bertrag Std.'  , form2)
  sheet.merge_range(     5,off,     5,off+6,unichr(0x00DC)+'bertrag Std.'  , form2)
  sheet.merge_range(nrow-3,off,nrow-3,off+6,'Summe Arb. Std.', form2)
  sheet.merge_range(nrow-2,off,nrow-2,off+6,'Fehlzeitstd.'   , form2)
  sheet.merge_range(nrow-1,off,nrow-1,off+6,'Summe Sollstd.' , form2)

  sheet.write(3,0,None,formLeft2)
  for i in range(0,ncol+1): sheet.write(nrow,     i,None,formTop2)
  for i in range(0,nrow):   sheet.write(   i,ncol+1,None,formLeft2)
  dat = datetime.date.today().strftime("%B %d, %Y")
  sheet.merge_range(nrow-10,ncol-1,nrow-10,ncol,'Datum:'       ,formTop2)
  sheet.merge_range(nrow- 9,ncol-1,nrow-8 ,ncol, '        '+dat,formLeft2)
  sheet.merge_range(nrow- 7,ncol-1,nrow-5 ,ncol,'Unterschrift:',formLeft2)
  sheet.merge_range(nrow- 4,ncol-1,nrow-3 ,ncol,'Abt.-Leiter:' ,formTop2)

  sheet.write(5,0,"9",formLeft)
  sheet.write(5,1,"0",formRight)

  if nrProjects == 2:
    sheet.write(6,0,"9",formLeft)
    sheet.write(6,1,"0",formRight)

  for i in range(len(proj)-1,-1,-1):
    if i==0:        sheet.write(5,2+i,proj[i],formLeft)
    if i>0 and i<4: sheet.write(5,2+i,proj[i],form)
    if i==4:        sheet.write(5,2+i,proj[i],formRight)

  if nrProjects == 2:
    proj2 = "77010"  # TIME-BRIDGE Nr hard coded here
    for i in range(len(proj2)-1,-1,-1):
      if i==0:        sheet.write(6,2+i,proj2[i],formLeft)
      if i>0 and i<4: sheet.write(6,2+i,proj2[i],form)
      if i==4:        sheet.write(6,2+i,proj2[i],formRight)

  for j in range(nrProjects):
    strArbeitHalf = '%.1f' % (arbeit/nrProjects)
    l = len(strArbeitHalf)
    for i in range(l):
      if i==l-1: sheet.write(     5+j,10    ,strArbeitHalf[i],formDash2)
      if i==l-3: sheet.write(     5+j, 9    ,strArbeitHalf[i],formDash1)
      if i==l-4: sheet.write(     5+j, 8    ,strArbeitHalf[i],form)
      if i==l-5: sheet.write(     5+j, 7    ,strArbeitHalf[i],formLeft)

  strArbeit = '%.1f' % arbeit
  l = len(strArbeit)
  for i in range(l):
    if i==l-1:
      sheet.write(     5  ,10+off,strArbeit[i],formDash2)
      sheet.write(nrow-1  ,10    ,strArbeit[i],formDash2)
      sheet.write(nrow-3  ,10+off,strArbeit[i],formDash2)
    if i==l-3:
      sheet.write(     5,   9+off,strArbeit[i],formDash1)
      sheet.write(nrow-1,   9    ,strArbeit[i],formDash1)
      sheet.write(nrow-3,   9+off,strArbeit[i],formDash1)
    if i==l-4:
      sheet.write(     5,   8+off,strArbeit[i],form)
      sheet.write(nrow-1,   8    ,strArbeit[i],form)
      sheet.write(nrow-3,   8+off,strArbeit[i],form)
    if i==l-5:
      sheet.write(     5,   7+off,strArbeit[i],formLeft)
      sheet.write(nrow-1,   7    ,strArbeit[i],formLeft)
      sheet.write(nrow-3,   7+off,strArbeit[i],formLeft)

  strFehl = '%.1f' % fehl
  l = len(strFehl)
  for i in range(l):
    if i==l-1: sheet.write(nrow-2,10+off,strFehl[i],formDash2)
    if i==l-3: sheet.write(nrow-2, 9+off,strFehl[i],formDash1)
    if i==l-4: sheet.write(nrow-2, 8+off,strFehl[i],form)
    if i==l-5: sheet.write(nrow-2, 7+off,strFehl[i],formLeft)

  strSoll = '%.1f' % soll
  l = len(strSoll)
  for i in range(l):
    if i==l-1: sheet.write(nrow-1,10+off,strSoll[i],formDash2)
    if i==l-3: sheet.write(nrow-1, 9+off,strSoll[i],formDash1)
    if i==l-4: sheet.write(nrow-1, 8+off,strSoll[i],form)
    if i==l-5: sheet.write(nrow-1, 7+off,strSoll[i],formLeft)

f.close()
book.close()

