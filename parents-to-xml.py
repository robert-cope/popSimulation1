import re,sys,os

#
#    Copyright 2013 Robert Cope cope.robert.c@gmail.com
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#

filedir = './'
infile = sys.argv[1]+'.txt'
print infile
inFile = open(infile,'r')
print infile
outFile = open(filedir + 'TEMP'+sys.argv[1]+'.txt','w')

matchingRegex = re.compile("\s{0,20}(\d{1,6})\ {1,4}\(\d{0,6}\s{0,4}\d{0,6}\)\ {1,4}(\d{1,2})(?:\ {1,4}(\d{1,4})\ {1,4}(\d{1,4}))?\ {1,4}([\?0-9MF]{1,4})\ {1,4}([\?0-9]{1,4})(?:\ {1,4}([\?0-9]{1,4}))?((?:\ {1,8}[\?0-9]{1,8}\.[\?0-9]{1,8}){1,100})")
print >> outFile, "1 25 . temp"
for i in range(25):
  print >> outFile, "abc%d" %(i)
print >> outFile, (len(inFile.readlines())-2)
inFile = open(filedir+'/'+infile,'r')
print infile
for line in inFile:
#match them all
  match = matchingRegex.match(line)
  if match:
  #start the first document and the loci list
    fifth=match.group(7)
    if fifth == None:
      fifth = '?'
    print >> outFile, "    %s %s %s %s%s - %s %s" %(match.group(1), match.group(5), match.group(6),fifth,match.group(8),match.group(3), match.group(4))
  #close everything
outFile.close()

infile = 'TEMP'+sys.argv[1]+'.txt'
print infile
inFile = open(filedir+'/'+infile,'r')
print infile
outFile = open(filedir + sys.argv[1]+'.xml','w')
matchingRegex0 = re.compile("([0-9]{1,4}) ([0-9]{1,4}) \. ([a-zA-Z]{1,12})")
matchingRegex1 = re.compile("([a-zA-Z]{1,6}[a-zA-Z0-9]{1,8})\s{0,80}$")
matchingRegex2 = re.compile("^\s{0,8}([0-9]{1,5})\s{0,80}$")
matchingRegex3 = re.compile("\s{0,20}(\d{1,6})\ {1,8}(\d{1,2})\ {1,4}([\?0-9MF]{1,4})\ {1,4}([\?0-9]{1,4})(?:\ {1,4}([\?0-9]{1,4}))?((?:\ {1,8}[\?0-9]{1,4}\.[\?0-9]{1,4}){1,100})(?:(?:\s{1,4}[0-9.]{3,15})?\ {1,4}\-\ {1,4}(\d{1,4})\ {1,4}(\d{1,4}))?\s{0,80}$")
for line in inFile:
  match0 = matchingRegex0.match(line)
  match1 = matchingRegex1.match(line)
  match2 = matchingRegex2.match(line)
  match3 = matchingRegex3.match(line)
  if match0:
    #start the first document and the loci list
    print >> outFile, "<population>"
    print >> outFile, "    <loci n=\"%s\">" %(match0.group(2))
  if match1:
    #elements of the loci list
    print >> outFile, "        <locus name=\"%s\" />" %(match1.group(1))
  
  if match2:
    print >> outFile, "    </loci>"
    print >> outFile, "    <individuals n=\"%s\">" %(match2.group(1))
  #close the loci list, start the individual list
  
  if match3:
    print 'match3'
    print line
    print match3.group(0)
  #put in all the elements from the individual list
    print >> outFile, "        <ind id=\"%s\" sex=\"%s\" birth=\"%s\" observed=\"%s\" sizeClass=\"%s\" genotype=\"%s\" />" %(match3.group(1), match3.group(2),match3.group(4),match3.group(7), match3.group(8),match3.group(6))
 
print >> outFile, "    </individuals>"
print >> outFile, "</population>"
  #close everything
outFile.close() 
os.remove('TEMP'+sys.argv[1]+'.txt') 
