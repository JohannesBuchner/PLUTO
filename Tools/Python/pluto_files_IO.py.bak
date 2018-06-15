import os
import sys

class PlutoFiles(object):
    def __init__(self, fname):
        self.fname = fname

    def OpenFile(self):
        try:
            self.f = open(self.fname, 'r')
        except IOError:
            print "File : %s not found"%self.fname

    def List2File(self, contents):
        self.fc = open(self.fname, 'w')
        # contents is a List.
        for item in contents:
            self.fc.write("%s"%item)

        self.fc.close()

    def File2List(self):
        NewList = []
        self.OpenFile()
        for line in self.f.readlines():
            NewList.append(line)
        return NewList

    def LocateString(self,string):
        strlist = []
        scrh = self.File2List()
        for item in scrh:
            if len(string.split()) == 1 and string.split()[0] == string: # Look for 1 word
                if string in item.split():
                    strlist.append([scrh.index(item), item])
            
            # Look for the whole string... Careful with \n matches thats why 1 word match is better
            elif len(string.split()) > 1: 
                if string == item:
                    strlist.append([scrh.index(item), item])

            else:
                print "Invalid String in LocateString"
                sys.exit()
                
        return strlist

    def DeleteLines(self, lbeg, lend):
        scrh = self.File2List()
        del scrh[lbeg:lend+1]
        self.List2File(scrh)

    def InsertLine(self, string, lbeg):
        scrh = self.File2List()
        scrh.insert(lbeg, string) #Inserts String at position lbeg.
        self.List2File(scrh)

    def ReplaceLine(self, NewLine, OldLineIndex):
        scrh = self.File2List()
        scrh[OldLineIndex] = NewLine
        self.List2File(scrh)

    def ReplaceWord(self, oldword, newword, DelOld = False):
        scrh = self.LocateString(oldword)
        try:
            scrh[0]
        except IndexError:
            #print 'Word %s not found in file %s'%(oldword, self.fname)
            pass
        else:
            l = scrh[0][1].split(' ')
            l[l.index(oldword)] = newword
            newline = ' '.join(l)
            if DelOld:
                self.DeleteLines(scrh[0][0], scrh[0][0])
                self.InsertLine(newword, scrh[0][0])
            else:
                self.ReplaceLine(newline,scrh[0][0])

    def ReadLines(self, beg, end):
        scrh = self.File2List()
        return scrh[beg:end+1]

    def UpdateListFromFile(self, scrhlist, updatelist):
        self.OpenFile()
        for line in self.f.readlines():
            for e in scrhlist:
                if e in line.split():
                    updatelist[scrhlist.index(e)] = line.split()[-1]
                else:
                    continue


    def ReplaceObsoletes(self):
        oldwords = ['INCLUDE_ROTATION', 'INCLUDE_BODY_FORCE', 'INCLUDE_COOLING', 'INCLUDE_BACKGROUND_FIELD',
                    'TIME_EVOLUTION', 'RAYMOND', 'NEQ', 'CT_VEC_POT_INIT', 'USE_VECTOR_POTENTIAL',
                    'SAVE_VEC_POT']
        newwords = ['ROTATING_FRAME', 'BODY_FORCE', 'COOLING','BACKGROUND_FIELD', 'TIME_STEPPING', 'SNEq',
                    'MINEq', 'USE_VECTOR_POTENTIAL', 'ASSIGN_VECTOR_POTENTIAL', 'UPDATE_VECTOR_POTENTIAL']

        for i in range(len(oldwords)):
            self.ReplaceWord(oldwords[i],newwords[i])
    
